
#' Read CosMx Spatial Transcriptomics Data
#'
#' Loads CosMx flat files from one or more TMA/slide subdirectories, builds
#' globally unique cell IDs across slides, separates negative control and
#' SystemControl probes, and converts pixel coordinates to millimetres.
#'
#' @param experiment_repo Character. Path to the root experiment directory.
#'   Each immediate subdirectory is treated as one TMA/slide and must contain
#'   \code{*_exprMat_file.csv.gz}, \code{*_metadata_file.csv.gz},
#'   and optionally \code{*-polygons.csv.gz}.
#'
#' @return A Seurat object with assays \code{RNA} (gene counts),
#'   \code{negprobes}, and \code{falsecode}, cell metadata, one spatial
#'   FOV per slide containing centroids, polygon coordinates in
#'   code{obj@misc$polygons},
#'   \code{orig.ident} set to \code{slidename}, active idents set to
#'   \code{orig.ident}, and \code{obj@misc$xy_condensed} pre-computed by
#'   \code{condenseTissues()} for immediate use with \code{TmaPlot()}.
#'   
readCosMx <- function(experiment_repo, MM_PER_PX = 0.120280945 / 1e3) {
  
  # 1. Discover TMA directories
  # ===========================
  message("Discovering TMA directories...")
  slide_paths <- list.dirs(experiment_repo, recursive = FALSE)
  if (length(slide_paths) == 0) stop("No TMAs found in the experiment repository.")
  
  slidenames <- basename(slide_paths)
  if (anyDuplicated(slidenames)) {
    stop("Duplicate TMA folder names: ",
         paste(unique(slidenames[duplicated(slidenames)]), collapse = ", "))
  }
  
  # 2. Load each slide
  # ==================
  slide_data <- lapply(seq_along(slide_paths), function(i) {
    current_path <- slide_paths[[i]]
    slidename    <- slidenames[[i]]
    message("Loading slide ", slidename, " (", i, "/", length(slide_paths), ")")
    
    find_file <- function(pattern, required = TRUE) {
      matches <- list.files(current_path, pattern, full.names = TRUE)
      if (length(matches) == 1L) return(matches)
      if (!required && length(matches) == 0L) return(NULL)
      stop("Expected one file matching '", pattern, "' in ", current_path,
           "; found ", length(matches), ".")
    }
    
    # a. Metadata
    message("Metadata...")
    meta_file <- find_file("metadata_file.*\\.csv\\.gz$")
    meta <- data.table::fread(meta_file, showProgress = FALSE)
    meta[, `:=`(
      slidename        = slidename,
      slide_ID_numeric = i,
      global_cell_ID   = paste0("c_", i, "_", fov, "_", cell_ID),
      FOV              = paste0("s", i, "f", fov)
    )]
    
    # b. Polygons
    message("Polygons...")
    poly_file <- find_file("polygons.*\\.csv\\.gz$", required = FALSE)
    polygons  <- NULL
    if (!is.null(poly_file)) {
      polygons <- data.table::fread(poly_file, showProgress = FALSE)
      polygons[, `:=`(
        cell_ID         = paste0("c_", i, "_", fov, "_", cellID),
        cell            = paste0("c_", i, "_", fov, "_", cellID),
        Run_Tissue_name = slidename,
        slidename       = slidename,
        FOV             = paste0("s", i, "f", fov),
        x_slide_mm      = x_global_px * MM_PER_PX,
        y_slide_mm      = y_global_px * MM_PER_PX
      )]
    } else {
      message("No polygon file found for ", slidename)
    }
    
    # c. Counts
    # Read the whole exprMat once (fast, multithreaded, single decompress), then
    # convert the dense gene block to sparse in in-memory row-blocks. Reading
    # once matters because the file is gzip-compressed and not seekable — a
    # skip-based chunked read would re-decompress from the top for every chunk.
    # Blocking the dense->sparse conversion keeps peak memory bounded without
    # ever building a zero-inclusive long-format table for the whole slide.
    message("Counts...")
    counts_path <- find_file("exprMat_file.*\\.csv\\.gz$")

    dt <- data.table::fread(counts_path, showProgress = FALSE)
    stopifnot(
      "exprMat file must contain 'fov' and 'cell_ID' columns" =
        all(c("fov", "cell_ID") %in% colnames(dt))
    )
    gene_cols <- setdiff(colnames(dt), c("fov", "cell_ID"))
    cell_ids  <- paste0("c_", i, "_", dt$fov, "_", dt$cell_ID)

    n_cells    <- nrow(dt)
    block_rows <- max(1L, floor(5e7 / length(gene_cols)))  # cap dense block size
    starts     <- seq.int(1L, n_cells, by = block_rows)
    i_acc <- j_acc <- x_acc <- vector("list", length(starts))
    for (b in seq_along(starts)) {
      rows        <- starts[b]:min(starts[b] + block_rows - 1L, n_cells)
      dense       <- as.matrix(dt[rows, gene_cols, with = FALSE])
      mode(dense) <- "numeric"
      tri         <- methods::as(Matrix::Matrix(dense, sparse = TRUE), "TsparseMatrix")
      i_acc[[b]]  <- tri@i + rows[1]      # tri@i is 0-based; rows[1] shifts to global 1-based
      j_acc[[b]]  <- tri@j + 1L
      x_acc[[b]]  <- tri@x
    }
    counts_mat <- Matrix::sparseMatrix(
      i        = unlist(i_acc),
      j        = unlist(j_acc),
      x        = unlist(x_acc),
      dims     = c(n_cells, length(gene_cols)),
      dimnames = list(cell_ids, gene_cols)
    )
    cell_order <- match(meta$global_cell_ID, rownames(counts_mat))
    if (anyNA(cell_order)) stop("Cell IDs in exprMat do not match metadata for ", slidename, ".")
    counts_mat <- counts_mat[cell_order, , drop = FALSE]
    list(meta = meta, polygons = polygons, counts = counts_mat)
  })
  
  # 3. Combine across slides
  # ========================
  message("Combining across slides...")
  all_genes <- Reduce(union, lapply(slide_data, \(s) colnames(s$counts)))
  if (length(unique(lapply(slide_data, \(s) sort(colnames(s$counts))))) > 1)
    warning("Gene panels differ across TMA slides — missing genes filled with 0.")
  
  per_slide <- lapply(slide_data, function(s) {
    missing_genes <- setdiff(all_genes, colnames(s$counts))
    if (length(missing_genes) > 0) {
      zero_block <- Matrix::sparseMatrix(
        i = integer(0), j = integer(0),
        dims     = c(nrow(s$counts), length(missing_genes)),
        dimnames = list(rownames(s$counts), missing_genes)
      )
      s$counts <- cbind(s$counts, zero_block)
    }
    s$counts[, all_genes, drop = FALSE]
  })
  # Single-slide runs must not go through do.call(rbind, .) — a one-element
  # list can fail to dispatch to the sparse-matrix rbind method.
  counts <- if (length(per_slide) == 1) per_slide[[1]] else do.call(rbind, per_slide)

  metadata <- data.table::rbindlist(lapply(slide_data, \(s) s$meta),  fill = TRUE)
  polygon_list <- Filter(Negate(is.null), lapply(slide_data, \(s) s$polygons))
  polygons <- if (length(polygon_list)) {
    data.table::rbindlist(polygon_list, fill = TRUE)
  } else {
    data.table::data.table()
  }
  
  # 4. Finalize metadata
  # ====================
  metadata[, `:=`(
    cell_ID        = global_cell_ID,
    cell_id        = global_cell_ID,
    x_slide_mm     = CenterX_global_px * MM_PER_PX,
    y_slide_mm     = CenterY_global_px * MM_PER_PX,
    cell           = NULL,
    global_cell_ID = NULL
  )]
  
  # 5. Separate control probes
  # ==========================
  is_neg      <- grepl("^Negative",      colnames(counts))
  is_sys      <- grepl("^SystemControl", colnames(counts))
  negcounts   <- counts[,  is_neg,           drop = FALSE]
  falsecounts <- counts[,  is_sys,           drop = FALSE]
  counts      <- counts[, !is_neg & !is_sys, drop = FALSE]

  # 6. Build Seurat object
  # ======================
  message("Building Seurat object...")
  meta_df           <- as.data.frame(metadata)
  rownames(meta_df) <- meta_df$cell_ID
  meta_df$orig.ident <- meta_df$slidename  # one TMA = one identity

  # counts is cells x genes — Seurat expects genes x cells.
  # Use Matrix::t() explicitly: base t() dispatches to t.default on a sparse
  # matrix unless the Matrix package happens to be attached, which fails with
  # "argument is not a matrix".
  obj <- Seurat::CreateSeuratObject(
    counts    = Matrix::t(counts),
    assay     = "RNA",
    meta.data = meta_df
  )

  # Set active idents to orig.ident (TMA / slide identity)
  Seurat::Idents(obj) <- obj$orig.ident

  # Negative probes and system controls as dedicated assays
  obj[["negprobes"]] <- Seurat::CreateAssayObject(counts = Matrix::t(negcounts))
  obj[["falsecode"]] <- Seurat::CreateAssayObject(counts = Matrix::t(falsecounts))
  
  # Add centroid coordinates. Full polygons remain in misc for QC without
  # creating a very large Seurat segmentation object.
  for (s in slidenames) {
    slide_meta <- metadata[slidename == s]
    cents <- SeuratObject::CreateCentroids(data.frame(
      x    = slide_meta$x_slide_mm,
      y    = slide_meta$y_slide_mm,
      cell = slide_meta$cell_ID
    ))
    obj[[s]] <- SeuratObject::CreateFOV(
      coords = list(centroids = cents),
      type   = "centroids",
      assay  = "RNA"
    )
  }

  # Store polygon local-px coordinates for downstream QC (e.g. computeSplitRatio)
  if (nrow(polygons) > 0)
    obj@misc$polygons <- polygons[, .(cell, slidename, FOV, x_local_px, y_local_px)]

  # Pre-compute condensed tissue layout for use with TmaPlot()
  obj <- condenseTissues(obj)

  obj
}


#' Condense Multiple TMA Tissues into a Compact Layout
#'
#' Rearranges cell coordinates so that spatially distant TMA cores are packed
#' side-by-side into a compact grid. The condensed coordinates are stored in
#' \code{obj@misc$xy_condensed} for downstream use by \code{TmaPlot()}.
#'
#' @param obj A Seurat object with coordinate and tissue columns in
#'   \code{@meta.data}.
#' @param tissue_col Metadata column for tissue/slide grouping.
#'   Default \code{"slidename"}.
#' @param x_col Metadata column for X coordinate in mm. Default \code{"x_slide_mm"}.
#' @param y_col Metadata column for Y coordinate in mm. Default \code{"y_slide_mm"}.
#' @param tissueorder Character vector controlling tissue placement order.
#'   Default \code{NULL} orders by decreasing tissue height.
#' @param buffer Gap in mm between tissue slots. Default \code{0.2}.
#' @param widthheightratio Target width/height ratio of the layout. Default \code{8/3}.
#' @param seed Integer seed for reproducibility. Default \code{1}.
#'
#' @return The Seurat object with \code{obj@misc$xy_condensed} added — a
#'   data frame with columns \code{x_mm}, \code{y_mm}, and \code{tissue},
#'   rownames matching cell IDs.
condenseTissues <- function(obj,
                             tissue_col       = "slidename",
                             x_col            = "x_slide_mm",
                             y_col            = "y_slide_mm",
                             tissueorder      = NULL,
                             buffer           = 0.2,
                             widthheightratio = 8/3,
                             seed             = 1) {

  md <- obj@meta.data
  missing_cols <- setdiff(c(tissue_col, x_col, y_col), colnames(md))
  if (length(missing_cols) > 0)
    stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))

  xy     <- as.matrix(md[, c(x_col, y_col)])
  tissue <- md[[tissue_col]]

  # Compute each tissue's bounding-box dimensions
  tissdf        <- data.frame(tissue = unique(tissue), stringsAsFactors = FALSE)
  tissdf$width  <- sapply(tissdf$tissue, function(t) diff(range(xy[tissue == t, 1])))
  tissdf$height <- sapply(tissdf$tissue, function(t) diff(range(xy[tissue == t, 2])))

  # Tissue order: explicit or by decreasing height
  if (!is.null(tissueorder)) {
    if (length(setdiff(tissdf$tissue, tissueorder)) > 0)
      stop("values in tissue missing from tissueorder")
    if (length(setdiff(tissueorder, tissdf$tissue)) > 0)
      stop("values in tissueorder missing from tissue")
    tissdf$order <- match(tissdf$tissue, tissueorder)
  } else {
    tissdf$order <- order(tissdf$height, decreasing = TRUE)
  }
  tissdf <- tissdf[tissdf$order, ]

  # Number of tissues per row to approximate widthheightratio
  tissuesperrow <- round(sqrt(nrow(tissdf)) * widthheightratio *
                           mean(tissdf$height) / mean(tissdf$width))
  targetwidth   <- sum(tissdf$width[seq_len(tissuesperrow)], na.rm = TRUE) +
                   buffer * (tissuesperrow - 1)

  # Place tissues on shelves
  tissdf$x        <- NA_real_
  tissdf$y        <- NA_real_
  tempx           <- 0
  tempy           <- 0
  tempshelfheight <- 0
  tempshelfwidth  <- 0

  for (i in seq_len(nrow(tissdf))) {
    tissdf$x[i]    <- tempx
    tissdf$y[i]    <- tempy
    tempshelfheight <- max(tempshelfheight, tissdf$height[i])
    tempshelfwidth  <- tempx + tissdf$width[i]
    tempx           <- tempx + tissdf$width[i] + buffer

    if (i < nrow(tissdf)) {
      if (abs(tempshelfwidth - targetwidth) <
          abs(tempshelfwidth + buffer + tissdf$width[i + 1] - targetwidth)) {
        tempy           <- tempy + tempshelfheight + buffer
        tempx           <- 0
        tempshelfheight <- 0
        tempshelfwidth  <- 0
      }
    }
  }

  # Shift each tissue's cells to their assigned slot origin
  set.seed(seed)
  for (t in unique(tissue)) {
    idx        <- tissue == t
    xy[idx, 1] <- xy[idx, 1] - min(xy[idx, 1]) + tissdf$x[tissdf$tissue == t]
    xy[idx, 2] <- xy[idx, 2] - min(xy[idx, 2]) + tissdf$y[tissdf$tissue == t]
  }

  obj@misc$xy_condensed <- data.frame(
    x_mm    = xy[, 1],
    y_mm    = xy[, 2],
    tissue  = tissue,
    row.names = rownames(md)
  )

  message("Condensed layout stored in obj@misc$xy_condensed. ",
          "Use TmaPlot(obj) to visualise.")
  obj
}


#' QC Plots for CosMx Cells
#'
#' Produces eight panels covering the main cell-level QC metrics:
#' histogram and scatter for each of \code{nCount_RNA}/\code{nFeature_RNA},
#' \code{Area}, \code{log2SBR}, and \code{SplitRatioToLocal} (plotted on log2 scale).
#' Optional red threshold lines can be added to any panel via the \code{*_min}/\code{*_max} arguments.
#'
#' @param obj A Seurat object with \code{nCount_RNA}, \code{nFeature_RNA},
#'   \code{Area}, \code{log2SBR}, and \code{SplitRatioToLocal} in \code{@meta.data}.
#' @param bins Integer. Number of histogram bins. Default \code{40}.
#' @param point_size Numeric. Point size for scatter plots. Default \code{1}.
#' @param base_size Numeric. Base font size passed to \code{ggthemes::theme_clean}. Default \code{20}.
#' @param col Character. Colour for all geoms. Default \code{"gray"}.
#' @param count_min Numeric or \code{NULL}. Red vertical line at lower \code{nCount_RNA} bound. Default \code{NULL}.
#' @param count_max Numeric or \code{NULL}. Red vertical line at upper \code{nCount_RNA} bound. Default \code{NULL}.
#' @param area_min Numeric or \code{NULL}. Red vertical line at lower \code{Area} bound. Default \code{NULL}.
#' @param area_max Numeric or \code{NULL}. Red vertical line at upper \code{Area} bound. Default \code{NULL}.
#' @param sbr_min Numeric or \code{NULL}. Red vertical line at lower \code{log2SBR} bound. Default \code{NULL}.
#' @param sbr_max Numeric or \code{NULL}. Red vertical line at upper \code{log2SBR} bound. Default \code{NULL}.
#' @param split_ratio_min Numeric or \code{NULL}. Red vertical line at lower \code{log2(SplitRatioToLocal)} bound. Default \code{NULL}.
#' @param split_ratio_max Numeric or \code{NULL}. Red vertical line at upper \code{log2(SplitRatioToLocal)} bound. Default \code{NULL}.
#'
#' @return A named list of eight ggplot objects:
#'   \code{nCounts_hist}, \code{nCounts_nGenes},
#'   \code{area_hist}, \code{nCounts_area},
#'   \code{log2SBR_hist}, \code{nCounts_log2SBR},
#'   \code{splitRatio_hist}, \code{nCounts_splitRatio}.
plotQCs <- function(obj,
                    bins            = 40,
                    point_size      = 1,
                    base_size       = 20,
                    col             = "gray",
                    count_min       = NULL,
                    count_max       = NULL,
                    area_min        = NULL,
                    area_max        = NULL,
                    sbr_min         = NULL,
                    sbr_max         = NULL,
                    split_ratio_min = NULL,
                    split_ratio_max = NULL) {

  md  <- obj@meta.data
  thr <- function(v) ggplot2::geom_vline(xintercept = v, colour = "red", linewidth = 0.5)
  thr_h <- function(v) ggplot2::geom_hline(yintercept = v, colour = "red", linewidth = 0.5)

  p1 <- data.table::data.table(nCounts = md$nCount_RNA) |>
    ggplot2::ggplot(ggplot2::aes(nCounts)) +
    ggplot2::geom_histogram(col = col, bins = bins) +
    ggplot2::ggtitle("CosMx: nCounts x nCells") +
    ggplot2::ylab("nCells") +
    ggthemes::theme_clean(base_size = base_size) +
    { if (!is.null(count_min)) thr(count_min) } +
    { if (!is.null(count_max)) thr(count_max) }

  p2 <- data.table::data.table(nCounts = md$nCount_RNA, nGenes = md$nFeature_RNA) |>
    ggplot2::ggplot(ggplot2::aes(nCounts, nGenes)) +
    ggplot2::geom_point(col = col, size = point_size) +
    ggplot2::ggtitle("CosMx: nCounts x nGenes") +
    ggthemes::theme_clean(base_size = base_size) +
    { if (!is.null(count_min)) thr(count_min) } +
    { if (!is.null(count_max)) thr(count_max) }

  p3 <- data.table::data.table(cell_area = md$Area) |>
    ggplot2::ggplot(ggplot2::aes(cell_area)) +
    ggplot2::geom_histogram(col = col, bins = bins) +
    ggplot2::ggtitle("CosMx: cellArea x nCells") +
    ggthemes::theme_clean(base_size = base_size) +
    { if (!is.null(area_min)) thr(area_min) } +
    { if (!is.null(area_max)) thr(area_max) }

  p4 <- data.table::data.table(nCounts = md$nCount_RNA, cell_area = md$Area) |>
    ggplot2::ggplot(ggplot2::aes(cell_area, nCounts)) +
    ggplot2::geom_point(col = col, size = point_size) +
    ggplot2::ggtitle("CosMx: cellArea x nCells") +
    ggthemes::theme_clean(base_size = base_size) +
    { if (!is.null(count_min)) thr_h(count_min) } +
    { if (!is.null(count_max)) thr_h(count_max) } +
    { if (!is.null(area_min))  thr(area_min) } +
    { if (!is.null(area_max))  thr(area_max) }

  has_sbr   <- "log2SBR"           %in% colnames(md)
  has_split <- "SplitRatioToLocal" %in% colnames(md)

  if (!has_sbr)   message("log2SBR not found — run computeSBR() first.")
  if (!has_split) message("SplitRatioToLocal not found — run computeSplitRatio() first.")

  p5 <- if (has_sbr) {
    data.table::data.table(log2SBR = md$log2SBR) |>
      ggplot2::ggplot(ggplot2::aes(log2SBR)) +
      ggplot2::geom_histogram(col = col, bins = bins, linewidth = 0.5) +
      ggplot2::ggtitle("CosMx: log2SBR x nCells") +
      ggthemes::theme_clean(base_size = base_size) +
      { if (!is.null(sbr_min)) thr(sbr_min) } +
      { if (!is.null(sbr_max)) thr(sbr_max) }
  } else {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "log2SBR unavailable\nrun computeSBR()") +
      ggplot2::theme_void()
  }

  p6 <- if (has_sbr) {
    data.table::data.table(nCounts = md$nCount_RNA, log2SBR = md$log2SBR) |>
      ggplot2::ggplot(ggplot2::aes(log2SBR, nCounts)) +
      ggplot2::geom_point(col = col, size = point_size) +
      ggplot2::ggtitle("CosMx: log2SBR x nCells") +
      ggthemes::theme_clean(base_size = base_size) +
      { if (!is.null(sbr_min))   thr(sbr_min) } +
      { if (!is.null(sbr_max))   thr(sbr_max) } +
      { if (!is.null(count_min)) thr_h(count_min) } +
      { if (!is.null(count_max)) thr_h(count_max) }
  } else {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "log2SBR unavailable\nrun computeSBR()") +
      ggplot2::theme_void()
  }

  p7 <- if (has_split) {
    data.table::data.table(spRatioToLocal = md$SplitRatioToLocal) |>
      ggplot2::ggplot(ggplot2::aes(log2(spRatioToLocal))) +
      ggplot2::geom_histogram(col = col, bins = bins) +
      ggplot2::ggtitle("CosMx: log2(SplitRatio) x nCells") +
      ggplot2::ylab("nCells") +
      ggthemes::theme_clean(base_size = base_size) +
      { if (!is.null(split_ratio_min)) thr(split_ratio_min) } +
      { if (!is.null(split_ratio_max)) thr(split_ratio_max) }
  } else {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "SplitRatioToLocal unavailable\nrun computeSplitRatio()") +
      ggplot2::theme_void()
  }

  p8 <- if (has_split) {
    data.table::data.table(nCounts = md$nCount_RNA, spRatioToLocal = md$SplitRatioToLocal) |>
      ggplot2::ggplot(ggplot2::aes(log2(spRatioToLocal), nCounts)) +
      ggplot2::geom_point(col = col, size = point_size) +
      ggplot2::ggtitle("CosMx: log2(SplitRatio) x nCells") +
      ggthemes::theme_clean(base_size = base_size) +
      { if (!is.null(split_ratio_min)) thr(split_ratio_min) } +
      { if (!is.null(split_ratio_max)) thr(split_ratio_max) } +
      { if (!is.null(count_min))       thr_h(count_min) } +
      { if (!is.null(count_max))       thr_h(count_max) }
  } else {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "SplitRatioToLocal unavailable\nrun computeSplitRatio()") +
      ggplot2::theme_void()
  }

  list(
    nCounts_hist       = p1,
    nCounts_nGenes     = p2,
    area_hist          = p3,
    nCounts_area       = p4,
    log2SBR_hist       = p5,
    nCounts_log2SBR    = p6,
    splitRatio_hist    = p7,
    nCounts_splitRatio = p8
  )
}



#'
#' Produces four \code{FeatureScatter} panels covering the main cell-level QC
#' metrics, with optional red threshold lines. Requires
#' \code{computeSplitRatio()} and \code{computeSBR()} to have been run first.
#'
#' @param obj A Seurat object built by \code{readCosMx()}.
#' @param count_min Numeric or \code{NULL}. Adds a vertical/horizontal red line
#'   at this \code{nCount_RNA} lower bound. \code{NULL} = no line. Default \code{NULL}.
#' @param count_max Numeric or \code{NULL}. Upper \code{nCount_RNA} line.
#'   Default \code{NULL}.
#' @param area_max Numeric or \code{NULL}. Upper \code{Area} line. Default
#'   \code{NULL}.
#' @param area_min Numeric or \code{NULL}. Lower \code{Area} line. Default
#'   \code{NULL}.
#' @param sbr_min Numeric or \code{NULL}. Lower \code{log2SBR} threshold line.
#'   Default \code{NULL}.
#' @param sbr_max Numeric or \code{NULL}. Upper \code{log2SBR} threshold line.
#'   Default \code{NULL}.
#' @param split_ratio_min Numeric or \code{NULL}. Lower bound of the flagged
#'   \code{SplitRatioToLocal} range (exclusive). \code{NULL} disables. Default \code{0}.
#' @param split_ratio_max Numeric or \code{NULL}. Upper bound of the flagged
#'   \code{SplitRatioToLocal} range (exclusive). \code{NULL} disables. Default \code{NULL}.
#'   Cells with \code{split_ratio_min < SplitRatioToLocal < split_ratio_max} are flagged.
#' @param title Character. Overall plot title. Default \code{"CosMx QC metrics"}.
#'
#' @return Invisibly returns a named list of ggplot objects:
#'   \code{nCount_nFeature}, \code{nCount_Area}, \code{SBR_nCount},
#'   \code{SplitRatio_Area}. The combined patchwork is printed as a side effect.
QCplots <- function(obj,
                    count_min        = NULL,
                    count_max        = NULL,
                    area_max         = NULL,
                    area_min         = NULL,
                    sbr_min          = NULL,
                    sbr_max          = NULL,
                    split_ratio_min  = NULL,
                    split_ratio_max  = NULL,
                    title            = "CosMx QC metrics") {

  thr_line <- function(...)
    ggplot2::geom_vline(..., colour = "red", linewidth = 0.5)
  thr_hline <- function(...)
    ggplot2::geom_hline(..., colour = "red", linewidth = 0.5)

  p1 <- Seurat::FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
      raster = FALSE) +
    ggplot2::ggtitle("nCount ~ nFeature") +
    { if (!is.null(count_min)) thr_line(xintercept = count_min) } +
    { if (!is.null(count_max)) thr_line(xintercept = count_max) }

  p2 <- Seurat::FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "Area",
      raster = FALSE) +
    ggplot2::ggtitle("nCount ~ Area") +
    { if (!is.null(count_min)) thr_line(xintercept  = count_min) } +
    { if (!is.null(count_max)) thr_line(xintercept  = count_max) } +
    { if (!is.null(area_max))  thr_hline(yintercept = area_max)  } +
    { if (!is.null(area_min))  thr_hline(yintercept = area_min)  }

  has_sbr   <- "log2SBR"           %in% colnames(obj@meta.data)
  has_split <- "SplitRatioToLocal" %in% colnames(obj@meta.data)

  if (!has_sbr)   message("log2SBR not found — run computeSBR() to add threshold line.")
  if (!has_split) message("SplitRatioToLocal not found — run computeSplitRatio() to add threshold line.")

  p3 <- Seurat::FeatureScatter(obj, raster = FALSE,
      feature1 = if (has_sbr) "log2SBR"    else "nCount_RNA",
      feature2 = if (has_sbr) "nCount_RNA" else "nFeature_RNA") +
    ggplot2::ggtitle(if (has_sbr) "SBR ~ nCount" else "nCount ~ nFeature (SBR unavailable)") +
    { if (has_sbr && !is.null(sbr_min)) thr_line(xintercept  = sbr_min)  } +
    { if (has_sbr && !is.null(sbr_max)) thr_line(xintercept  = sbr_max)  } +
    { if (has_sbr && !is.null(count_min)) thr_hline(yintercept = count_min) } +
    { if (has_sbr && !is.null(count_max)) thr_hline(yintercept = count_max) }

  p4 <- Seurat::FeatureScatter(obj, raster = FALSE,
      feature1 = if (has_split) "Area"              else "nCount_RNA",
      feature2 = if (has_split) "SplitRatioToLocal" else "Area") +
    ggplot2::ggtitle(if (has_split) "SplitRatio ~ Area" else "nCount ~ Area (SplitRatio unavailable)") +
    { if (has_split && !is.null(area_max))             thr_line(xintercept  = area_max)             } +
    { if (has_split && !is.null(area_min))             thr_line(xintercept  = area_min)             } +
    { if (has_split && !is.null(split_ratio_min))      thr_hline(yintercept = split_ratio_min)      } +
    { if (has_split && !is.null(split_ratio_max))      thr_hline(yintercept = split_ratio_max)      }

  plots <- list(
    nCount_nFeature = p1,
    nCount_Area     = p2,
    SBR_nCount      = p3,
    SplitRatio_Area = p4
  )

  combined <- Reduce(`+`, lapply(plots, function(p) p + Seurat::NoLegend())) +
    patchwork::plot_layout(ncol = 2) +
    patchwork::plot_annotation(title = title)

  print(combined)
  invisible(plots)
}


#' Per-TMA QC Overview (distribution + scatter panels)
#'
#' Produces cell-level QC histograms and scatter plots for \code{nCount_RNA},
#' \code{nFeature_RNA}, \code{Area}, \code{log2SBR}, and
#' \code{log2(SplitRatioToLocal)} — but \strong{faceted by TMA / slide}
#' (\code{group.by}) so slides can be compared side by side. No spatial maps.
#'
#' The regional and boundary metrics are treated as \strong{prerequisites}: run
#' \code{computeSBR()} and \code{computeSplitRatio()} \emph{before} calling this
#' function. If a required metric column is missing, the function stops with the
#' exact call to run — it never silently recomputes heavy metrics.
#'
#' @param obj A Seurat object built by \code{readCosMx()}, with QC metrics
#'   already computed.
#' @param group.by Metadata column defining the TMA/slide facet. Default
#'   \code{"slidename"}.
#' @param bins Integer. Histogram bins. Default \code{40}.
#' @param point_size Numeric. Scatter point size. Default \code{0.4}.
#' @param base_size Numeric. Base font size. Default \code{14}.
#' @param col Character. Colour for geoms. Default \code{"gray40"}.
#' @param count_min,count_max Optional red threshold lines for \code{nCount_RNA}.
#' @param feature_min,feature_max Optional red threshold lines for
#'   \code{nFeature_RNA}.
#' @param area_min,area_max Optional red threshold lines for \code{Area}.
#' @param sbr_min,sbr_max Optional red threshold lines for \code{log2SBR}.
#' @param split_ratio_min,split_ratio_max Optional red threshold lines for
#'   \code{log2(SplitRatioToLocal)} (border cells only).
#'
#' @return A named list of eight ggplot objects (one set of panels, each faceted
#'   by \code{group.by}): \code{nCounts_hist}, \code{nCounts_nGenes},
#'   \code{area_hist}, \code{nCounts_area}, \code{log2SBR_hist},
#'   \code{nCounts_log2SBR}, \code{splitRatio_hist}, \code{nCounts_splitRatio}.
QCoverview <- function(obj,
                       group.by        = "slidename",
                       bins            = 40,
                       point_size      = 0.4,
                       base_size       = 14,
                       col             = "gray40",
                       count_min       = NULL, count_max       = NULL,
                       feature_min     = NULL, feature_max     = NULL,
                       area_min        = NULL, area_max        = NULL,
                       sbr_min         = NULL, sbr_max         = NULL,
                       split_ratio_min = NULL, split_ratio_max = NULL,
                       neg_max         = NULL,
                       density_min     = NULL, density_max     = NULL,
                       extranuclear_max = NULL) {

  md <- obj@meta.data

  # --- Prerequisites -------------------------------------------------------
  if (!group.by %in% colnames(md))
    stop("group.by column '", group.by, "' not found in @meta.data.")
  miss_basic <- setdiff(c("nCount_RNA", "nFeature_RNA", "Area"), colnames(md))
  if (length(miss_basic))
    stop("Missing basic metric columns: ", paste(miss_basic, collapse = ", "))
  if (!"log2SBR" %in% colnames(md))
    stop("log2SBR not found — run `obj <- computeSBR(obj)` before QCoverview().")
  if (!"SplitRatioToLocal" %in% colnames(md))
    stop("SplitRatioToLocal not found — run `obj <- computeSplitRatio(obj)` before QCoverview().")

  # --- Helpers -------------------------------------------------------------
  thr    <- function(v) if (!is.null(v)) ggplot2::geom_vline(xintercept = v, colour = "red", linewidth = 0.5)
  thr_h  <- function(v) if (!is.null(v)) ggplot2::geom_hline(yintercept = v, colour = "red", linewidth = 0.5)
  facetL <- ggplot2::facet_wrap(ggplot2::vars(.data[[".grp"]]), scales = "free")
  theme0 <- ggthemes::theme_clean(base_size = base_size)

  # All cells (for count/area/SBR panels)
  df <- data.table::data.table(
    .grp       = as.factor(md[[group.by]]),
    nCounts    = md$nCount_RNA,
    nGenes     = md$nFeature_RNA,
    area       = md$Area,
    log2SBR    = md$log2SBR,
    splitRatio = md$SplitRatioToLocal
  )
  # Optional extra metrics (only plotted when the compute* helper has been run)
  if ("pct_false"        %in% colnames(md)) df[, pct_false        := md$pct_false]
  if ("counts_per_um2"   %in% colnames(md)) df[, counts_per_um2   := md$counts_per_um2]
  if ("pct_extranuclear" %in% colnames(md)) df[, pct_extranuclear := md$pct_extranuclear]

  # Border cells only (SplitRatio > 0), on log2 scale
  bd <- df[splitRatio > 0]
  if (nrow(bd) > 0) bd[, log2SplitRatio := log2(splitRatio)]

  # Faceted-histogram helper for the optional single-metric panels
  histf <- function(xcol, title, xlab, lo = NULL, hi = NULL)
    ggplot2::ggplot(df, ggplot2::aes(.data[[xcol]])) +
      ggplot2::geom_histogram(col = col, bins = bins) + facetL +
      ggplot2::labs(title = title, x = xlab, y = "nCells") + theme0 +
      thr(lo) + thr(hi)

  # --- Panels (mirror plotQCs, faceted by TMA) ----------------------------
  p1 <- ggplot2::ggplot(df, ggplot2::aes(nCounts)) +
    ggplot2::geom_histogram(col = col, bins = bins) + facetL +
    ggplot2::labs(title = "nCounts x nCells", y = "nCells") + theme0 +
    thr(count_min) + thr(count_max)

  p2 <- ggplot2::ggplot(df, ggplot2::aes(nCounts, nGenes)) +
    ggplot2::geom_point(col = col, size = point_size) + facetL +
    ggplot2::labs(title = "nCount_RNA x nFeature_RNA") + theme0 +
    thr(count_min) + thr(count_max) + thr_h(feature_min) + thr_h(feature_max)

  p3 <- ggplot2::ggplot(df, ggplot2::aes(area)) +
    ggplot2::geom_histogram(col = col, bins = bins) + facetL +
    ggplot2::labs(title = "cellArea x nCells", x = "cell_area", y = "nCells") + theme0 +
    thr(area_min) + thr(area_max)

  p4 <- ggplot2::ggplot(df, ggplot2::aes(area, nCounts)) +
    ggplot2::geom_point(col = col, size = point_size) + facetL +
    ggplot2::labs(title = "cellArea x nCounts", x = "cell_area") + theme0 +
    thr(area_min) + thr(area_max) + thr_h(count_min) + thr_h(count_max)

  p5 <- ggplot2::ggplot(df, ggplot2::aes(log2SBR)) +
    ggplot2::geom_histogram(col = col, bins = bins) + facetL +
    ggplot2::labs(title = "log2SBR x nCells", y = "nCells") + theme0 +
    thr(sbr_min) + thr(sbr_max)

  p6 <- ggplot2::ggplot(df, ggplot2::aes(log2SBR, nCounts)) +
    ggplot2::geom_point(col = col, size = point_size) + facetL +
    ggplot2::labs(title = "log2SBR x nCounts") + theme0 +
    thr(sbr_min) + thr(sbr_max) + thr_h(count_min) + thr_h(count_max)

  if (nrow(bd) > 0) {
    p7 <- ggplot2::ggplot(bd, ggplot2::aes(log2SplitRatio)) +
      ggplot2::geom_histogram(col = col, bins = bins) + facetL +
      ggplot2::labs(title = "log2(SplitRatio) x nCells (border cells)", y = "nCells") + theme0 +
      thr(split_ratio_min) + thr(split_ratio_max)

    p8 <- ggplot2::ggplot(bd, ggplot2::aes(log2SplitRatio, nCounts)) +
      ggplot2::geom_point(col = col, size = point_size) + facetL +
      ggplot2::labs(title = "log2(SplitRatio) x nCounts (border cells)") + theme0 +
      thr(split_ratio_min) + thr(split_ratio_max) + thr_h(count_min) + thr_h(count_max)
  } else {
    message("No border cells (SplitRatioToLocal > 0) — SplitRatio panels are empty placeholders.")
    p7 <- p8 <- ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No border cells")
  }

  panels <- list(
    nCounts_hist       = p1,
    nCounts_nGenes     = p2,
    area_hist          = p3,
    nCounts_area       = p4,
    log2SBR_hist       = p5,
    nCounts_log2SBR    = p6,
    splitRatio_hist    = p7,
    nCounts_splitRatio = p8
  )

  # Optional extra-metric panels (appended only if the metric was computed)
  if ("pct_false" %in% colnames(df))
    panels$negFraction_hist <- histf("pct_false",
      "neg-probe % of RNA counts x nCells", "pct_false", lo = neg_max)
  if ("counts_per_um2" %in% colnames(df))
    panels$density_hist <- histf("counts_per_um2",
      "counts per um2 x nCells", "counts_per_um2", lo = density_min, hi = density_max)
  if ("pct_extranuclear" %in% colnames(df))
    panels$extranuclear_hist <- histf("pct_extranuclear",
      "extranuclear transcript % x nCells", "pct_extranuclear", hi = extranuclear_max)

  panels
}


#' Plot TMA Cores with Flexible Metadata Colouring
#'
#' The primary spatial visualisation function for CosMx TMA data. Plots cells
#' at their condensed mm coordinates, faceted by tissue/TMA core, and coloured
#' by any metadata column. Works on the full object or any \code{subset()}.
#' Automatically runs \code{condenseTissues()} if coordinates are not yet
#' computed.
#'
#' @param obj A Seurat object. Must contain \code{x_slide_mm} and
#'   \code{y_slide_mm} in \code{@meta.data} (set by \code{readCosMx()}).
#' @param col.by Character. Name of a \code{@meta.data} column to colour cells
#'   by. Continuous columns use a viridis gradient; discrete columns use
#'   ggplot2's default categorical palette. Default \code{NULL} plots all cells
#'   in a single colour.
#' @param pt.size Numeric. Point size. Default \code{0.5}.
#' @param main Character. Plot title. Default \code{NULL} uses \code{col.by}
#'   as title, or \code{"TMA Layout"} when \code{col.by} is \code{NULL}.
#' @param dark Logical. Use a dark background theme. Default \code{FALSE}.
#' @param cols Character vector. Custom colours. For discrete \code{col.by},
#'   one colour per level; for continuous, passed to
#'   \code{scale_colour_gradientn()}. Default \code{NULL}.
#' @param subsample_frac Numeric in (0, 1]. Fraction of cells to plot.
#'   Default \code{1} (all cells). Reduce for faster interactive previews.
#' @param legend.max.levels Integer. Maximum number of discrete levels before
#'   the legend is automatically hidden. Default \code{50}.
#' @param facet Logical. Facet by tissue. Default \code{TRUE}. Set \code{FALSE}
#'   to plot all tissues on a single panel (uses condensed coordinates).
#' @param label Logical. Overlay a centred label for each discrete \code{col.by}
#'   group (analogous to \code{DimPlot(label = TRUE)} in Seurat). Ignored for
#'   continuous \code{col.by}. Default \code{FALSE}.
#' @param label.size Numeric. Font size for group labels (passed to
#'   \code{geom_label} / \code{geom_label_repel}). Default \code{3}.
#' @param label.repel Logical. When \code{TRUE} (the default) and
#'   \code{ggrepel} is installed, uses \code{ggrepel::geom_label_repel()} to
#'   avoid overlapping labels. Falls back to \code{geom_label()} when
#'   \code{ggrepel} is not available.
#' @param seed Integer. Random seed for reproducible subsampling. Default \code{1}.
#'
#' @return A \code{ggplot} object.
TmaPlot <- function(obj,
                    col.by            = "orig.ident",
                    pt.size           = 0.01,
                    main              = NULL,
                    dark              = TRUE,
                    cols              = NULL,
                    facet             = TRUE,
                    subsample_frac    = 1,
                    legend.max.levels = 50,
                    label             = FALSE,
                    label.size        = 3,
                    label.repel       = TRUE,
                    seed              = 1) {

  # 1. Get / compute condensed coordinates
  if (is.null(obj@misc$xy_condensed)) {
    message("obj@misc$xy_condensed not found — running condenseTissues() automatically.")
    obj <- condenseTissues(obj)
  }
  xy_cond <- obj@misc$xy_condensed[
    intersect(rownames(obj@misc$xy_condensed), colnames(obj)), , drop = FALSE
  ]

  # 2. Attach colour metadata
  if (!col.by %in% c(colnames(obj@meta.data), colnames(xy_cond)))
    stop("col.by '", col.by, "' not found in obj@meta.data.")
  if (!col.by %in% colnames(xy_cond))
    xy_cond[[col.by]] <- obj@meta.data[rownames(xy_cond), col.by]

  # 3. Subsample
  set.seed(seed)
  if (subsample_frac < 1)
    xy_cond <- xy_cond[sample(nrow(xy_cond), round(nrow(xy_cond) * subsample_frac)), ]

  # 4. Title
  if (is.null(main)) main <- col.by

  # 5. Build plot
  p <- ggplot2::ggplot(xy_cond,
                       ggplot2::aes(x = x_mm, y = y_mm, colour = .data[[col.by]])) +
    ggplot2::geom_point(size = pt.size) +
    ggplot2::labs(title = main, y = "y (mm)", x = "x (mm)", colour = col.by)

  if (facet) p <- p + ggplot2::facet_wrap(~ tissue, scales = "free")

  # 6. Colour scale
  is_continuous <- is.numeric(xy_cond[[col.by]])
  fov_signal_cols <- c(
    "<-2"   = "#2166AC",
    "-2:-1" = "#92C5DE",
    "-1:0"  = "#D1E5F0",
    "0:1"   = "#FDDBC7",
    "1:2"   = "#F4A582",
    ">2"    = "#D6604D"
  )
  if (is_continuous) {
    p <- p + if (!is.null(cols))
      ggplot2::scale_colour_gradientn(colours = cols)
    else
      ggplot2::scale_colour_viridis_c()
  } else {
    n_levels  <- length(unique(xy_cond[[col.by]]))
    col_vals  <- if (!is.null(cols)) cols else
      if (col.by == "fov_signal_loss_cat") fov_signal_cols else
        setNames(scales::hue_pal()(n_levels), levels(factor(xy_cond[[col.by]])))
    p <- p + ggplot2::scale_colour_manual(values = col_vals)
    # scale_fill_manual only matters for the label layer (the only fill user);
    # adding it unconditionally warns "no shared levels" on point-only plots.
    if (label) p <- p + ggplot2::scale_fill_manual(values = col_vals)
    if (n_levels > legend.max.levels) {
      message("col.by '", col.by, "' has ", n_levels, " levels — legend hidden ",
              "(increase legend.max.levels to show).")
      p <- p + ggplot2::guides(colour = "none")
    }
  }

  # 7. Group labels (discrete col.by only)
  if (label && !is_continuous) {
    group_cols <- if (facet) c("tissue", col.by) else col.by
    label_df <- do.call(
      data.frame,
      c(
        lapply(
          setNames(group_cols, group_cols),
          function(col) xy_cond[[col]]
        ),
        list(
          y_mm = xy_cond$y_mm,
          x_mm = xy_cond$x_mm
        )
      )
    )
    label_df <- aggregate(
      cbind(y_mm, x_mm) ~ .,
      data  = label_df,
      FUN   = median
    )

    use_repel <- label.repel && requireNamespace("ggrepel", quietly = TRUE)
    label_aes <- ggplot2::aes(
      x     = x_mm,
      y     = y_mm,
      label = .data[[col.by]],
      fill  = .data[[col.by]]
    )
    label_base <- list(
      data        = label_df,
      mapping     = label_aes,
      size        = label.size,
      colour      = "white",
      alpha       = 0.7,
      show.legend = FALSE
    )
    if (use_repel) {
      p <- p + do.call(ggrepel::geom_label_repel,
                       c(label_base, list(label.size        = NA,
                                          min.segment.length = 0,
                                          box.padding        = 0.25)))
    } else {
      p <- p + do.call(ggplot2::geom_label,
                       c(label_base, list(linewidth = 0)))
    }
  }

  # 8. Theme
  dark_theme <- ggplot2::theme_dark() +
    ggplot2::theme(
      aspect.ratio          = NULL,
      plot.background       = ggplot2::element_rect(fill = "#1a1a1a", colour = NA),
      panel.background      = ggplot2::element_rect(fill = "#1a1a1a"),
      panel.grid.major      = ggplot2::element_line(colour = "#333333"),
      panel.grid.minor      = ggplot2::element_blank(),
      strip.background      = ggplot2::element_rect(fill = "#333333"),
      strip.text            = ggplot2::element_text(colour = "white"),
      axis.text             = ggplot2::element_text(colour = "grey70"),
      axis.title            = ggplot2::element_text(colour = "grey70"),
      plot.title            = ggplot2::element_text(colour = "white"),
      legend.background     = ggplot2::element_rect(fill = "#1a1a1a", colour = NA),
      legend.box.background = ggplot2::element_rect(fill = "#1a1a1a", colour = NA),
      legend.key            = ggplot2::element_rect(fill = "#1a1a1a"),
      legend.text           = ggplot2::element_text(colour = "grey70"),
      legend.title          = ggplot2::element_text(colour = "grey70")
    )

  p <- p + if (dark) dark_theme else ggplot2::theme_bw()

  # Return the plot as a VISIBLE value so it auto-prints exactly once:
  #   TmaPlot(obj)                 -> draws the plot
  #   TmaPlot(obj) + theme_clean() -> draws the plot with the extra theme, still once
  # Do NOT print(p) here as well: an internal print + a theme added outside would
  # draw the plot twice (once without the theme, once with it).
  p
}



#' Compute SplitRatioToLocal for CosMx Cells
#'
#' Identifies cells whose polygon vertices touch the edge of their slide's local
#' FOV coordinate space and computes a ratio of their area relative to the
#' mean cell area within the same FOV. Non-boundary cells receive a value of 0.
#'
#' Follows the same approach as Giuseppe's \code{dataprep_cosmx()}: boundary
#' detection uses \code{x_local_px} / \code{y_local_px} stored in
#' \code{obj@misc$polygons}, computed independently per slide so that each
#' slide's own coordinate extremes are used as the FOV frame boundary.
#'
#' @param obj A Seurat object built by \code{readCosMx()}, which stores polygon
#'   vertex coordinates in \code{obj@misc$polygons} (columns \code{cell},
#'   \code{slidename}, \code{FOV}, \code{x_local_px}, \code{y_local_px}). The
#'   \code{@meta.data} slot must contain \code{cell_id}, \code{fov}, and
#'   \code{Area}.
#'
#' @return The Seurat object with \code{SplitRatioToLocal} added to
#'   \code{@meta.data}. Value is \code{0} for non-boundary cells and
#'   \code{round(Area / mean_fov_area, 2)} for boundary cells.
#'   Values \code{> 1} are strong filter candidates.
computeSplitRatio <- function(obj) {

  md <- obj@meta.data

  # Warn if already present and non-trivial — mirrors Giuseppe's guard
  if ("SplitRatioToLocal" %in% colnames(md) && !all(is.na(md$SplitRatioToLocal)))
    warning("SplitRatioToLocal already present in metadata and will be overwritten.")

  if (is.null(obj@misc$polygons))
    stop("obj@misc$polygons not found. Ensure the Seurat object was built with ",
         "readCosMx(), which stores x_local_px and y_local_px in obj@misc$polygons.")

  polygons <- obj@misc$polygons

  # Boundary detection and area normalization must be independent for each FOV.
  # ----------------------------------------------------------------------------
  boundary_cells <- unique(unlist(lapply(
    split(polygons, interaction(polygons$slidename, polygons$FOV, drop = TRUE)),
    function(s) {
      s$cell[
        s$x_local_px %in% c(min(s$x_local_px), max(s$x_local_px)) |
        s$y_local_px %in% c(min(s$y_local_px), max(s$y_local_px))
      ]
    }
  )))
  
  # 2. Per-FOV mean area for normalisation
  has_boundary  <- md$cell_id %in% boundary_cells
  mean_area     <- ave(md$Area, interaction(md$slidename, md$fov), FUN = mean)
  
  ratio         <- ifelse(has_boundary,
                          round(md$Area / mean_area, 2),
                          0)
  names(ratio)  <- rownames(md)
  
  message("Boundary cells detected: ", sum(has_boundary),
          " / ", nrow(md), " total cells.")
  
  obj@meta.data$SplitRatioToLocal <- ratio[rownames(obj@meta.data)]
  obj
}


#' Compute Smoothed Signal-to-Background Ratio (SBR) for CosMx Cells
#'
#' For each cell, computes a spatially smoothed signal-to-background ratio using
#' negative probe counts as the background estimate. A Gaussian kernel weights
#' each cell's neighbourhood; the SBR is the ratio of smoothed mean gene counts
#' to smoothed mean negative probe counts. Low \code{log2(SBR)} values indicate
#' cells in tissue regions where background dominates signal.
#'
#' Mirrors the regional QC step in \code{gbspatial::run_spatial_qc()}.
#' Requires \code{dbscan} and \code{Matrix}.
#'
#' @param obj A Seurat object built by \code{readCosMx()}. Must contain
#'   \code{FOV} in \code{@meta.data} and assays \code{RNA} and \code{negprobes}.
#'   \code{obj@misc$xy_condensed} is used for spatial coordinates (matching
#'   Giuseppe's pipeline); if absent, \code{condenseTissues()} is called
#'   automatically.
#' @param bandwidth Numeric. Gaussian kernel bandwidth in mm. Default \code{0.01}.
#' @param weight_cutoff Numeric. Minimum kernel weight below which neighbours
#'   are excluded (controls search radius). Default \code{0.08}.
#' @param bg_quantile Numeric in [0, 1). The smoothed background is floored at
#'   this quantile of its positive values before forming the ratio. This stops
#'   cells in near-zero-background neighbourhoods from producing an exploded
#'   \code{SBR} (a \code{smoothed_neg} of ~0 otherwise gave \code{SBR} up to
#'   ~1e9, an artifactual high-\code{log2SBR} second mode). Only the lowest
#'   \code{bg_quantile} of background estimates are affected; the rest of the
#'   distribution is unchanged. Set to \code{0} to restore the old behaviour.
#'   Default \code{0.01}.
#'
#' @return The Seurat object with \code{SBR} and \code{log2SBR} added to
#'   \code{@meta.data}.
computeSBR <- function(obj, bandwidth = 0.01, weight_cutoff = 0.08, bg_quantile = 0.01) {

  if (!requireNamespace("dbscan", quietly = TRUE))
    stop("Package 'dbscan' is required. Install with install.packages('dbscan').")

  # Use condensed tissue coordinates — mirrors Giuseppe's run_spatial_qc(),
  # which always receives condensed xy from dataprep_cosmx().
  if (is.null(obj@misc$xy_condensed)) {
    message("obj@misc$xy_condensed not found — running condenseTissues() first.")
    obj <- condenseTissues(obj)
  }

  md      <- obj@meta.data
  xy      <- as.matrix(obj@misc$xy_condensed[rownames(md), c("x_mm", "y_mm")])
  counts    <- Matrix::t(GetAssayData(obj, assay = "RNA",       layer = "counts"))
  negcounts <- Matrix::t(GetAssayData(obj, assay = "negprobes", layer = "counts"))

  message("Building spatial neighbourhood graph (bandwidth = ", bandwidth, " mm)...")
  max_dist   <- sqrt(-2 * bandwidth^2 * log(weight_cutoff))
  nn         <- dbscan::frNN(xy, eps = max_dist)
  n_neighbors <- sapply(nn$id, length)
  i_idx      <- rep(seq_len(nrow(xy)), times = n_neighbors)
  j_idx      <- unlist(nn$id)
  weights    <- exp(-(unlist(nn$dist)^2) / (2 * bandwidth^2))

  conn <- Matrix::sparseMatrix(
    i    = i_idx, j = j_idx, x = weights,
    dims = c(nrow(xy), nrow(xy))
  )
  Matrix::diag(conn) <- 1
  conn <- Matrix::Diagonal(x = 1 / Matrix::rowSums(conn)) %*% conn

  # Avoid materialising conn %*% counts (n_cells × n_genes dense matrix, ~6 GB).
  # rowMeans(conn %*% counts) == conn %*% rowMeans(counts) — cheap vector op.
  message("Smoothing counts and negprobes...")
  smoothed_signal <- as.vector(conn %*% Matrix::rowMeans(counts))
  smoothed_neg    <- as.vector(conn %*% Matrix::rowMeans(negcounts))

  # Floor the smoothed background so near-zero-background neighbourhoods don't
  # blow the ratio up to ~1e9. The floor is a low quantile of the positive
  # background values, so only the pathological low tail is affected.
  pos_neg  <- smoothed_neg[smoothed_neg > 0]
  bg_floor <- if (bg_quantile > 0 && length(pos_neg) > 0)
    as.numeric(stats::quantile(pos_neg, bg_quantile, na.rm = TRUE)) else 1e-9
  n_floored <- sum(smoothed_neg < bg_floor)
  message("Background floor (q", bg_quantile, " of positive smoothed_neg) = ",
          signif(bg_floor, 3), "; cells floored: ", n_floored,
          " (", round(100 * n_floored / length(smoothed_neg), 2), "%)")

  sbr             <- smoothed_signal / pmax(smoothed_neg, bg_floor)
  names(sbr)      <- rownames(md)

  obj@meta.data$SBR     <- sbr[rownames(obj@meta.data)]
  obj@meta.data$log2SBR <- log2(obj@meta.data$SBR)

  n_low <- sum(obj@meta.data$log2SBR < 0, na.rm = TRUE)
  message("Cells with log2(SBR) < 0: ", n_low, " / ", nrow(md),
          " (", round(100 * n_low / nrow(md), 1), "%)")
  obj
}



#' Compute FOV Imaging Integrity QC for CosMx Data
#'
#' Reproduces the FOV barcode QC step from \code{gbspatial::run_spatial_qc()},
#' adapted to take a Seurat object as input and write all metrics back to
#' \code{@meta.data}. Internally calls \code{gbspatial:::runFOVQC()}.
#'
#' Matches Giuseppe's implementation exactly: uses \strong{condensed tissue
#' coordinates} (from \code{obj@misc$xy_condensed}) so that FOV neighbourhood
#' lookup in \code{runFOVQC} uses the same coordinate space as
#' \code{gbspatial::run_spatial_qc()}. If \code{obj@misc$xy_condensed} is
#' absent, \code{condenseTissues()} is called automatically.
#'
#' Two failure modes are checked per FOV (via \code{runFOVQC}):
#' \itemize{
#'   \item \strong{Barcode bias}: >50\% of grids within the FOV show >
#'     \code{max_prop_loss} dropout for one or more reporter cycles.
#'   \item \strong{Total count loss}: >75\% of grids show total counts
#'     \code{max_totalcounts_loss} below neighbouring FOVs.
#' }
#'
#' Two metrics are added to \code{@meta.data}:
#' \itemize{
#'   \item \code{flag_fov_integrity}: logical. \code{TRUE} = cell belongs to a
#'     failed FOV. Mirrors Giuseppe's \code{flag_fovqc}.
#'   \item \code{fov_signal_loss}: numeric. Per-cell log2 fold-change of total
#'     counts relative to neighbouring grid squares. Mirrors the colour axis of
#'     Giuseppe's \code{FOVSignalLossSpatialPlot}. Cells in sparse grid squares
#'     (<10 cells) are \code{NA}.
#' }
#'
#' @param obj A Seurat object built by \code{readCosMx()}. Must contain
#'   \code{FOV} in \code{@meta.data} and assay \code{RNA}.
#'   \code{obj@misc$xy_condensed} is used for spatial coordinates; if absent,
#'   \code{condenseTissues()} is called automatically.
#' @param panel_name Character. CosMx panel name. One of \code{"Hs_6k"},
#'   \code{"Hs_IO"}, \code{"Hs_UCC"}, \code{"Hs_WTX"}, \code{"Mm_Neuro"},
#'   \code{"Mm_UCC"}. Default \code{"Hs_6k"}.
#' @param max_prop_loss Numeric in (0,1). Maximum fraction of barcode positions
#'   allowed to drop out before an FOV is flagged. Default \code{0.6}.
#' @param max_totalcounts_loss Numeric in (0,1). Maximum fractional total count
#'   loss relative to neighbouring FOVs before an FOV is flagged. Default \code{0.6}.
#'
#' @return The Seurat object with \code{flag_fov_integrity} and
#'   \code{fov_signal_loss} added to \code{@meta.data}, and the full
#'   \code{runFOVQC} result stored in \code{obj@misc$fov_integrity}.
computeFOVintegrity <- function(obj,
                                panel_name           = NULL,
                                max_prop_loss        = 0.6,
                                max_totalcounts_loss = 0.6) {

  if (!requireNamespace("gbspatial", quietly = TRUE))
    stop("Package 'gbspatial' is required.")

  # Use condensed tissue coordinates — mirrors Giuseppe's run_spatial_qc(),
  # which always receives condensed xy from dataprep_cosmx().
  if (is.null(obj@misc$xy_condensed)) {
    message("obj@misc$xy_condensed not found — running condenseTissues() first.")
    obj <- condenseTissues(obj)
  }

  md         <- obj@meta.data
  counts_mat <- Matrix::t(GetAssayData(obj, assay = "RNA", layer = "counts"))

  # xy column names must match runFOVQC expectations (x_slide_mm / y_slide_mm)
  xy_mat           <- as.matrix(obj@misc$xy_condensed[rownames(md), c("x_mm", "y_mm")])
  colnames(xy_mat) <- c("x_slide_mm", "y_slide_mm")

  # Slide identity keeps FOV IDs unique in multi-slide experiments.
  tissue_vec <- md$slidename

  # panel_name: inferred from unique Panel values unless explicitly supplied
  if (is.null(panel_name))
    stop("panel_name must be supplied. Choose from: ",
         paste(names(gbspatial:::barcodes_by_panel), collapse = ", "))

  barcodemap <- gbspatial:::barcodes_by_panel[[panel_name]]
  if (is.null(barcodemap))
    stop("panel_name '", panel_name, "' not found. Choose from: ",
         paste(names(gbspatial:::barcodes_by_panel), collapse = ", "))

  message("Running FOV integrity QC (panel: ", panel_name, ")...")
  res <- gbspatial:::runFOVQC(
    counts               = counts_mat,
    xy                   = xy_mat,
    fov                  = md$FOV,
    tissue               = tissue_vec,
    barcodemap           = barcodemap,
    max_prop_loss        = max_prop_loss,
    max_totalcounts_loss = max_totalcounts_loss
  )

  # Per-cell tissue+fov ID — matches the format used in res$flaggedfovs
  cell_fovid <- paste0(tissue_vec, md$FOV)

  # Three flag columns mirroring the three flaggedfovs outputs from runFOVQC:
  #   flag_fov_integrity   — combined flag (union of the two below)
  #   flag_fov_totalcounts — FOV total counts too low vs spatial neighbours
  #   flag_fov_bias        — barcode-channel bias detected in FOV
  obj@meta.data$flag_fov_integrity   <- cell_fovid %in% res$flaggedfovs
  obj@meta.data$flag_fov_totalcounts <- cell_fovid %in% res$flaggedfovs_fortotalcounts
  obj@meta.data$flag_fov_bias        <- cell_fovid %in% res$flaggedfovs_forbias

  # Per-cell log2 FC of total counts vs neighbouring grid squares.
  # res$gridinfo$gridid is positionally aligned with rows of counts_mat;
  # assign cell barcodes as names so metadata lookup works correctly.
  cell_scores <- res$totalcountsresids[res$gridinfo$gridid]
  names(cell_scores) <- rownames(md)
  obj@meta.data$fov_signal_loss <- cell_scores[rownames(md)]
  obj@meta.data$fov_signal_loss_cat <- cut(
    cell_scores[rownames(md)],
    breaks = c(-Inf, -2, -1, 0, 1, 2, Inf),
    labels = c("<-2", "-2:-1", "-1:0", "0:1", "1:2", ">2"),
    right  = TRUE
  )

  obj@misc$fov_integrity <- res

  n_fov     <- length(unique(cell_fovid))
  n_flagged <- length(res$flaggedfovs)
  n_cells   <- sum(obj@meta.data$flag_fov_integrity)
  message("FOVs flagged: ", n_flagged, " / ", n_fov,
          " — affects ", n_cells, " cells (",
          round(100 * n_cells / nrow(md), 1), "%)")
  if (length(res$flaggedfovs_fortotalcounts) > 0)
    message("  Total-counts failures: ",
            paste(res$flaggedfovs_fortotalcounts, collapse = ", "))
  if (length(res$flaggedfovs_forbias) > 0)
    message("  Barcode-bias failures: ",
            paste(res$flaggedfovs_forbias, collapse = ", "))

  # ── Diagnostic plots ────────────────────────────────────────────────────────
  message("Generating diagnostic plots...")
  .colramp <- colorRampPalette(c("darkblue", "blue", "grey80", "red", "darkred"))(101)

  # Helper: open a null device, draw, record, close
  .rec <- function(draw_fn) {
    pdf(NULL)
    dev.control(displaylist = "enable")
    draw_fn()
    p <- recordPlot()
    dev.off()
    p
  }

  plots <- list()

  # 1. Map of flagged FOVs (all FOVs in blue, flagged in red)
  plots$flagged_fovs <- .rec(function() {
    plot(res$xy, cex = 0.1, asp = 1, pch = 16, col = "grey80", main = "Flagged FOVs")
    for (f in unique(res$fov)) {
      inds <- res$fov == f
      rect(min(res$xy[inds, 1]), min(res$xy[inds, 2]),
           max(res$xy[inds, 1]), max(res$xy[inds, 2]),
           col = adjustcolor("dodgerblue2", alpha.f = 0.5))
    }
    for (f in res$flaggedfovs) {
      inds <- res$fov == f
      rect(min(res$xy[inds, 1]), min(res$xy[inds, 2]),
           max(res$xy[inds, 1]), max(res$xy[inds, 2]),
           col = adjustcolor("red", alpha.f = 0.5))
      text(median(range(res$xy[inds, 1])), median(range(res$xy[inds, 2])), f, col = "green")
    }
  })

  # 2. Spatial log2 fold-change in total counts vs comparable regions
  plots$signal_loss <- .rec(function() {
    plot(res$xy, cex = 0.2, asp = 1, pch = 16,
         col = .colramp[pmax(pmin(
           51 + res$totalcountsresids[
             match(res$gridinfo$gridid, names(res$totalcountsresids))] * 25,
           101), 1)],
         main = "Log2 fold-change in total counts vs comparable regions")
    for (f in unique(res$fov)) {
      inds <- res$fov == f
      rect(min(res$xy[inds, 1]), min(res$xy[inds, 2]),
           max(res$xy[inds, 1]), max(res$xy[inds, 2]), border = "black")
    }
    for (f in res$flaggedfovs_fortotalcounts) {
      inds <- res$fov == f
      rect(min(res$xy[inds, 1]), min(res$xy[inds, 2]),
           max(res$xy[inds, 1]), max(res$xy[inds, 2]), border = "yellow", lwd = 2)
      text(median(range(res$xy[inds, 1])), median(range(res$xy[inds, 2])), f, col = "green")
    }
    legend("right", pch = 16,
           col    = rev(c("darkblue", "blue", "grey80", "red", "darkred")),
           legend = rev(c("< -2", -1, 0, 1, "> 2")))
  })

  # 3. Heatmap of per-FOV barcode-bit bias (pheatmap returns a grob directly)
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    plots$bias_heatmap <- pheatmap::pheatmap(
      res$fovstats$bias * res$fovstats$flag,
      col    = colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(100),
      breaks = seq(-2, 2, length.out = 101),
      main   = "FOV bias: log2(fold-change) from comparable regions",
      silent = TRUE
    )
  } else {
    message("  pheatmap not available — skipping bias heatmap.")
  }

  # 4. Spatial plots for each flagged reporter-cycle × channel combination
  bitnames   <- colnames(res$fovstats$p)
  colorvals  <- unique(substr(bitnames, nchar(bitnames), nchar(bitnames)))
  flagged_rc <- colnames(res$flags_per_fov_x_reportercycle)[
    colSums(res$flags_per_fov_x_reportercycle >= 0.5) > 0]

  if (length(flagged_rc) > 0) {
    bits_to_plot <- match(
      paste0(rep(flagged_rc, each = length(colorvals)),
             rep(colorvals, length(flagged_rc))),
      colnames(res$resid))
    bits_to_plot <- bits_to_plot[!is.na(bits_to_plot)]

    bit_plots <- lapply(bits_to_plot, function(i) {
      .rec(function() {
        par(mar = c(0, 0, 2, 0))
        plot(res$xy, cex = 0.2, asp = 1, pch = 16,
             col = .colramp[pmax(pmin(
               51 + res$resid[match(res$gridinfo$gridid, rownames(res$resid)), i] * 50,
               101), 1)],
             main = paste0(colnames(res$resid)[i],
                           ": log2(fold-change)\nfrom comparable regions"))
        for (f in unique(res$fov)) {
          inds <- res$fov == f
          rect(min(res$xy[inds, 1]), min(res$xy[inds, 2]),
               max(res$xy[inds, 1]), max(res$xy[inds, 2]), border = "black")
        }
        for (f in rownames(res$fovstats$flag)[res$fovstats$flag[, i] > 0]) {
          inds <- res$fov == f
          rect(min(res$xy[inds, 1]), min(res$xy[inds, 2]),
               max(res$xy[inds, 1]), max(res$xy[inds, 2]), lwd = 2, border = "yellow")
        }
        legend("right", pch = 16,
               col    = rev(c("darkblue", "blue", "grey80", "red", "darkred")),
               legend = rev(c("< -1", -0.5, 0, 0.5, "> 1")))
      })
    })
    names(bit_plots) <- colnames(res$resid)[bits_to_plot]
    plots$bit_effects <- bit_plots
  } else {
    plots$bit_effects <- list()
  }

  obj@misc$fov_integrity_plots <- plots

  obj
}


#' Run Full CosMx QC Pipeline on a Seurat Object
#'
#' Sequentially computes all five QC metrics used by Giuseppe's
#' \code{gbspatial::run_spatial_qc()} pipeline, adds per-cell flags to
#' \code{@meta.data}, prints a summary table, and returns the filtered object.
#'
#' Steps performed (can be toggled individually):
#' \enumerate{
#'   \item \strong{nCount_RNA} — flag cells outside [\code{count_min}, \code{count_max}].
#'   \item \strong{Cell area} — flag cells outside [\code{area_min}, \code{area_max}].
#'   \item \strong{FOV integrity} — flag cells in FOVs with barcode signal dropout
#'     via \code{computeFOVintegrity()}.
#'   \item \strong{FOV boundary} — flag partially cropped cells via
#'     \code{computeSplitRatio()} (skipped if already computed).
#'   \item \strong{Regional SBR} — flag cells in high-background regions via
#'     \code{computeSBR()} (skipped if already computed).
#'   \item \strong{Density} (optional, \code{do_density}) — flag cells with
#'     \code{counts_per_um2} outside [\code{density_min}, \code{density_max}];
#'     computed via \code{computeDensity()} if missing.
#'   \item \strong{Extranuclear} (optional, \code{do_extranuclear}) — flag cells
#'     with \code{pct_extranuclear > extranuclear_max} (segmentation leakage /
#'     no-nucleus cells). Requires \code{computeTranscriptComposition()} first;
#'     off by default (review in \code{QCoverview()} before enabling).
#' }
#'
#' All steps only \strong{flag}; a single removal happens at the end via
#' \code{flag_overall} when \code{filter = TRUE}. Set \code{filter = FALSE} to
#' add flags without removing any cells.
#'
#' @param obj A Seurat object built by \code{readCosMx()}.
#' @param count_min Integer or \code{NULL}. Flag cells with \code{nCount_RNA <
#'   count_min}. \code{NULL} disables the lower bound. Default \code{20}.
#' @param count_max Integer or \code{NULL}. Flag cells with \code{nCount_RNA >
#'   count_max}. \code{NULL} disables the upper bound. Default \code{NULL}.
#' @param area_min Numeric or \code{NULL}. Flag cells with \code{Area <
#'   area_min}. \code{NULL} disables the lower bound. Default \code{NULL}.
#' @param area_max Numeric or \code{NULL}. Flag cells with \code{Area >
#'   area_max}. \code{NULL} disables the upper bound. Default \code{30000}.
#' @param split_ratio_min Numeric or \code{NULL}. Lower bound of the flagged
#'   \code{SplitRatioToLocal} range (exclusive). Default \code{0}.
#' @param split_ratio_max Numeric or \code{NULL}. Upper bound of the flagged
#'   \code{SplitRatioToLocal} range (exclusive). Default \code{0.5}.
#' @param sbr_min Numeric or \code{NULL}. Flag cells with \code{log2SBR < sbr_min}. \code{NULL} disables the lower bound. Default \code{0}.
#' @param sbr_max Numeric or \code{NULL}. Flag cells with \code{log2SBR > sbr_max}. \code{NULL} disables the upper bound. Default \code{NULL}.
#' @param panel_name Character. CosMx panel for FOV integrity QC. Default \code{"Hs_6k"}.
#' @param fov_integrity_threshold Numeric in (0,1). Fraction of barcode/count
#'   signal loss above which an FOV is flagged. Default \code{0.6} (60\%).
#' @param do_nCount Logical. Run nCount filter. Default \code{TRUE}.
#' @param do_area Logical. Run cell area filter. Default \code{TRUE}.
#' @param do_fov_integrity Logical. Run FOV integrity filter. Default \code{TRUE}.
#' @param do_boundary Logical. Run FOV boundary filter. Default \code{TRUE}.
#' @param do_sbr Logical. Run regional SBR filter. Default \code{TRUE}.
#' @param density_min,density_max Numeric or \code{NULL}. Flag cells with
#'   \code{counts_per_um2} below/above these bounds. Default \code{NULL}.
#' @param extranuclear_max Numeric or \code{NULL}. Flag cells with
#'   \code{pct_extranuclear} above this (segmentation leakage). Default \code{95}.
#' @param do_density Logical. Run transcript-density filter. Default \code{FALSE}.
#' @param do_extranuclear Logical. Run extranuclear-fraction filter (needs
#'   \code{computeTranscriptComposition()} first). Default \code{FALSE}.
#' @param upset_nintersects Number of intersection columns to show, largest
#'   first. Default \code{20} keeps the plot readable by trimming the long tail
#'   of tiny (often 1-cell) combinations; the per-filter totals remain complete
#'   in the set-size bars and \code{summary} table. Use \code{NA} to show every
#'   combination.
#' @param upset_text_scale Numeric scalar or length-6 vector passed to
#'   \code{UpSetR::upset(text.scale=)} (intersection title, intersection ticks,
#'   set-size title, set-size ticks, set names, bar numbers). Default \code{1.3}.
#' @param upset_point_size,upset_line_size Numeric. Dot/line size of the UpSet
#'   intersection matrix. Defaults \code{2.2} and \code{0.8}.
#' @param upset_mb_ratio Length-2 numeric \code{c(bars, matrix)} giving the
#'   height split between the intersection barplot and the dot matrix. Default
#'   \code{c(0.6, 0.4)} (taller bars than the UpSetR default).
#' @param filter Logical. Return filtered object. If \code{FALSE}, returns the
#'   object with flags added but no cells removed. Default \code{TRUE}.
#'
#' @return The Seurat object with per-cell flag columns added to
#'   \code{@meta.data} (\code{flag_nCount}, \code{flag_area},
#'   \code{flag_fov_integrity}, \code{flag_boundary}, \code{flag_sbr},
#'   \code{flag_overall}) and, if \code{filter = TRUE}, flagged cells removed.
QCfiltering <- function(obj,
                        count_min             = 30,
                        count_max             = NULL,
                        feature_min           = NULL,
                        feature_max           = NULL,
                        area_min              = NULL,
                        area_max              = 80e3,
                        split_ratio_min       = log2(0.5),
                        split_ratio_max       = log2(4),
                        sbr_min               = 0,
                        sbr_max               = NULL,
                        density_min           = NULL,
                        density_max           = NULL,
                        extranuclear_max      = 95,
                        panel_name            = "Hs_6k",
                        fov_integrity_threshold = 0.6,
                        do_nCount             = TRUE,
                        do_nFeature           = TRUE,
                        do_area               = TRUE,
                        do_fov_integrity      = TRUE,
                        do_boundary           = TRUE,
                        do_sbr                = TRUE,
                        do_density            = FALSE,
                        do_extranuclear       = FALSE,
                        upset_nintersects     = 20,
                        upset_text_scale      = 1.3,
                        upset_point_size      = 2.2,
                        upset_line_size       = 0.8,
                        upset_mb_ratio        = c(0.6, 0.4),
                        filter                = TRUE) {

  n_start <- ncol(obj)

  # 1. nCount_RNA — flag below min and/or above max (NULL = no bound)
  obj@meta.data$flag_nCount <- if (do_nCount) {
    flag <- rep(FALSE, ncol(obj))
    if (!is.null(count_min)) flag <- flag | obj@meta.data$nCount_RNA < count_min
    if (!is.null(count_max)) flag <- flag | obj@meta.data$nCount_RNA > count_max
    flag
  } else FALSE

  # 1b. nFeature_RNA (gene complexity) — WTx guide recommends complexity filtering
  #     alongside counts (examine the joint distribution). NULL bound = no cutoff.
  obj@meta.data$flag_nFeature <- if (do_nFeature) {
    flag <- rep(FALSE, ncol(obj))
    if (!is.null(feature_min)) flag <- flag | obj@meta.data$nFeature_RNA < feature_min
    if (!is.null(feature_max)) flag <- flag | obj@meta.data$nFeature_RNA > feature_max
    flag
  } else FALSE

  # 2. Cell area — flag below min and/or above max (NULL = no bound)
  obj@meta.data$flag_area <- if (do_area) {
    flag <- rep(FALSE, ncol(obj))
    if (!is.null(area_min)) flag <- flag | obj@meta.data$Area < area_min
    if (!is.null(area_max)) flag <- flag | obj@meta.data$Area > area_max
    flag
  } else FALSE

  # 3. FOV integrity (barcode QC)
  if (do_fov_integrity) {
    if (is.null(obj@misc$fov_integrity)) {
      obj <- computeFOVintegrity(obj, panel_name = panel_name,
                                   max_prop_loss        = fov_integrity_threshold,
                                   max_totalcounts_loss = fov_integrity_threshold)
    } else {
      cell_fovid <- paste0(obj@meta.data$slidename, obj@meta.data$FOV)
      obj@meta.data$flag_fov_integrity <- cell_fovid %in% obj@misc$fov_integrity$flaggedfovs
    }
  } else {
    obj@meta.data$flag_fov_integrity <- FALSE
  }

  # 4. FOV boundary (SplitRatioToLocal = area / local FOV mean) — only border
  #    cells (SplitRatio > 0) are eligible; flag those whose log2(ratio) falls
  #    outside [min, max]. Default split_ratio_min = log2(0.5) matches the WTx
  #    guide (remove border cells < 50% of local mean area, i.e. truncated
  #    fragments); split_ratio_max = log2(4) additionally drops oversized edge
  #    cells. Non-border cells (ratio == 0) are never flagged. Flagging only;
  #    removal happens once, via flag_overall, at the end.
  if (do_boundary) {
    if (!"SplitRatioToLocal" %in% colnames(obj@meta.data)) {
      obj <- computeSplitRatio(obj)
    }
    sr     <- obj@meta.data$SplitRatioToLocal
    border <- sr > 0
    l2     <- suppressWarnings(log2(sr))
    flag   <- rep(FALSE, ncol(obj))
    if (!is.null(split_ratio_min)) flag <- flag | (border & l2 < split_ratio_min)
    if (!is.null(split_ratio_max)) flag <- flag | (border & l2 > split_ratio_max)
    obj@meta.data$flag_boundary <- flag
  } else {
    obj@meta.data$flag_boundary <- FALSE
  }


  # 5. Regional SBR
  if (do_sbr) {
    if (!"log2SBR" %in% colnames(obj@meta.data)) {
      obj <- computeSBR(obj)
    }
    obj@meta.data$flag_sbr <- {
      flag <- rep(FALSE, ncol(obj))
      if (!is.null(sbr_min)) flag <- flag | obj@meta.data$log2SBR < sbr_min
      if (!is.null(sbr_max)) flag <- flag | obj@meta.data$log2SBR > sbr_max
      flag
    }
  } else {
    obj@meta.data$flag_sbr <- FALSE
  }

  # 6. Transcript density (counts per um^2) — computed if missing (cheap)
  if (do_density) {
    if (!"counts_per_um2" %in% colnames(obj@meta.data)) obj <- computeDensity(obj)
    flag <- rep(FALSE, ncol(obj)); d <- obj@meta.data$counts_per_um2
    if (!is.null(density_min)) flag <- flag | d < density_min
    if (!is.null(density_max)) flag <- flag | d > density_max
    obj@meta.data$flag_density <- flag
  } else {
    obj@meta.data$flag_density <- FALSE
  }

  # 7. Extranuclear transcript fraction (segmentation leakage). NOT auto-computed
  #    (needs the transcript file) — run computeTranscriptComposition() first.
  if (do_extranuclear) {
    if (!"pct_extranuclear" %in% colnames(obj@meta.data))
      stop("pct_extranuclear not found — run computeTranscriptComposition(obj, tx_file) ",
           "before QCfiltering(do_extranuclear = TRUE).")
    obj@meta.data$flag_extranuclear <-
      if (!is.null(extranuclear_max)) obj@meta.data$pct_extranuclear > extranuclear_max
      else rep(FALSE, ncol(obj))
  } else {
    obj@meta.data$flag_extranuclear <- FALSE
  }

  # Overall flag
  obj@meta.data$flag_overall <-
    obj@meta.data$flag_nCount        |
    obj@meta.data$flag_nFeature      |
    obj@meta.data$flag_area          |
    obj@meta.data$flag_fov_integrity |
    obj@meta.data$flag_boundary      |
    obj@meta.data$flag_sbr           |
    obj@meta.data$flag_density       |
    obj@meta.data$flag_extranuclear

  # Summary table
  md      <- obj@meta.data
  n_total <- nrow(md)

  flag_rule <- function(mn, mx) {
    lo <- if (!is.null(mn)) paste0("< ", signif(mn, 4)) else NULL
    hi <- if (!is.null(mx)) paste0("> ", signif(mx, 4)) else NULL
    if (is.null(lo) && is.null(hi)) "—" else paste(c(lo, hi), collapse = " or ")
  }

  flags <- list(
    md$flag_nCount, md$flag_nFeature, md$flag_area, md$flag_fov_integrity,
    md$flag_boundary, md$flag_sbr, md$flag_density, md$flag_extranuclear
  )
  active <- c(do_nCount, do_nFeature, do_area, do_fov_integrity,
              do_boundary, do_sbr, do_density, do_extranuclear)
  filter_names <- c("RNA counts", "Detected genes", "Cell area", "FOV integrity",
                    "FOV boundary", "Regional SBR", "Transcript density",
                    "Extranuclear fraction")
  rules <- c(
    paste("nCount_RNA", flag_rule(count_min, count_max)),
    paste("nFeature_RNA", flag_rule(feature_min, feature_max)),
    paste("Area", flag_rule(area_min, area_max)),
    paste0("signal loss > ", fov_integrity_threshold),
    paste("log2 ratio", flag_rule(split_ratio_min, split_ratio_max)),
    paste("log2 SBR", flag_rule(sbr_min, sbr_max)),
    paste("density", flag_rule(density_min, density_max)),
    paste("extranuclear %", flag_rule(NULL, extranuclear_max))
  )

  summary_df <- data.frame(
    Filter    = filter_names[active],
    Flag_rule = rules[active],
    Flagged   = vapply(flags[active], sum, numeric(1)),
    Pct       = round(100 * vapply(flags[active], mean, numeric(1)), 2)
  )
  summary_df <- rbind(summary_df, data.frame(
    Filter = "Overall", Flag_rule = "any filter",
    Flagged = sum(md$flag_overall), Pct = round(100 * mean(md$flag_overall), 2)
  ))
  message("\n--- QC summary (", n_total, " cells) ---")
  print(summary_df, row.names = FALSE)

  # UpSet plot — mirrors Giuseppe's run_spatial_qc() approach:
  # build a named list of cell IDs per filter, use UpSetR::fromList(),
  # and add an overall pass/fail pie chart inset via patchwork.
  if (requireNamespace("UpSetR", quietly = TRUE) &&
      requireNamespace("patchwork", quietly = TRUE)) {

    cell_ids <- rownames(md)
    filter_list <- list(
      `RNA counts (nCount_RNA)`       = cell_ids[md$flag_nCount],
      `Detected genes (nFeature_RNA)` = cell_ids[md$flag_nFeature],
      `Cell area`                      = cell_ids[md$flag_area],
      `FOV integrity`                  = cell_ids[md$flag_fov_integrity],
      `FOV boundary`                   = cell_ids[md$flag_boundary],
      `Regional SBR`                   = cell_ids[md$flag_sbr],
      `Transcript density`             = cell_ids[md$flag_density],
      `Extranuclear fraction`          = cell_ids[md$flag_extranuclear]
    )
    # keep only filters that flagged at least one cell (mirrors Giuseppe's purrr::keep)
    filter_list <- Filter(function(x) length(x) > 0, filter_list)

    if (length(filter_list) > 1) {
      # Per-filter totals ("cells lost") are the set-size bars + summary table and
      # are always complete; nintersects only caps the combination columns to keep
      # the plot readable (order.by="freq" keeps the largest). Use NA to show all.
      u_plot <- UpSetR::upset(
        UpSetR::fromList(filter_list),
        nintersects = upset_nintersects,
        order.by    = "freq",
        nsets       = length(filter_list),
        text.scale  = upset_text_scale,
        point.size  = upset_point_size,
        line.size   = upset_line_size,
        mb.ratio    = upset_mb_ratio
      )
      p_pie_overall <- {
        df <- data.frame(Category = ifelse(md$flag_overall, "Flagged", "Kept")) |>
          (\(d) { d$n <- ave(d$Category, d$Category, FUN = length); d })() |>
          unique()
        df$n         <- as.integer(df$n)
        df$Pct       <- round(100 * df$n / sum(df$n), 1)
        df$Label     <- paste0(df$Category, "\nn=", df$n, "\n(", df$Pct, "%)")
        ggplot2::ggplot(df, ggplot2::aes(x = 2, y = n, fill = Category)) +
          ggplot2::geom_bar(stat = "identity", width = 1, color = "white") +
          ggplot2::coord_polar("y", start = 0) +
          ggplot2::scale_fill_manual(
            values = c(Flagged = "#D73027", Kept = "gray90")) +
          ggplot2::geom_text(ggplot2::aes(label = Label),
            position = ggplot2::position_stack(vjust = 0.5),
            size = 2.5, fontface = "bold") +
          ggplot2::theme_void() +
          ggplot2::xlim(0.5, 2.5) +
          ggplot2::theme(legend.position = "none")
      }
      p_upset <- suppressWarnings(
        patchwork::wrap_elements(grid::grid.grabExpr(print(u_plot))) +
          patchwork::inset_element(p_pie_overall,
            left = 0.65, bottom = 0.65, right = 1, top = 1)
      )
    } else {
      message("Fewer than 2 active filters — skipping UpSet plot.")
      p_upset <- NULL
    }
  } else {
    message("UpSetR or patchwork not available — skipping UpSet plot.")
    p_upset <- NULL
  }

  if (is.null(obj@misc$qc_plots)) obj@misc$qc_plots <- list()
  obj@misc$qc_plots$upset   <- p_upset
  obj@misc$qc_plots$summary <- summary_df

  if (filter) {
    keep  <- rownames(md)[!md$flag_overall]
    obj   <- suppressWarnings(subset(obj, cells = keep))
    message("\nRetained: ", ncol(obj), " / ", n_start,
            " cells (", round(100 * ncol(obj) / n_start, 1), "%)")
  }

  obj
}


#' Data-driven QC cutoffs (WTx guide: 2.5 SD in log10 space)
#'
#' Computes adaptive lower/upper cutoffs for count/complexity metrics as
#' \code{10^(mean(log10 x) +/- nsd * sd(log10 x))}, matching the NanoString WTx
#' quality-processing guide's thresholding for \code{nCount_RNA} /
#' \code{nFeature_RNA}. Feed the results into \code{QCfiltering()} (e.g.
#' \code{count_min = cutoffs$nCount_RNA["min"]}) to replace fixed thresholds.
#'
#' @param obj A Seurat object.
#' @param metrics Character vector of \code{@meta.data} columns. Default
#'   \code{c("nCount_RNA", "nFeature_RNA")}.
#' @param nsd Number of standard deviations. Default \code{2.5} (guide value).
#' @return Named list; each element a numeric \code{c(min, max)} on the raw scale.
qc_cutoffs_2p5sd <- function(obj, metrics = c("nCount_RNA", "nFeature_RNA"), nsd = 2.5) {
  md  <- obj@meta.data
  out <- lapply(metrics, function(m) {
    if (!m %in% colnames(md)) stop("metric '", m, "' not in meta.data")
    x  <- md[[m]]; x <- x[is.finite(x) & x > 0]
    lx <- log10(x); mu <- mean(lx); s <- stats::sd(lx)
    c(min = 10^(mu - nsd * s), max = 10^(mu + nsd * s))
  })
  names(out) <- metrics
  out
}


#' Parse a TMA Map Matrix from an Excel File
#'
#' Reproduces Giuseppe's approach from \code{gb_model.Rmd}: reads the entire
#' sheet as a raw character matrix, strips everything after the first comma in
#' each cell (e.g. \code{"37166,B1"} → \code{"37166"}), then extracts the
#' patient-ID rows.
#'
#' For the TA649 \strong{constructionKey} layout (6 annotation rows per TMA
#' row: patient-ID, year, block, core-position, tissue, disease), the
#' patient-ID rows fall at Excel rows 2, 8, 14, … — i.e. every 6th row
#' starting from row 2.  The function detects this automatically from
#' \code{rows_per_core} (default \code{6}) and \code{first_pid_row} (default
#' \code{2}).
#'
#' For the \strong{PrintingKey} layout (3 rows per TMA row) pass
#' \code{rows_per_core = 3}.
#'
#' @param path Path to the \code{.xlsx} file.
#' @param sheet Sheet name or index. Default \code{"TA649 constructionKey"}.
#' @param first_col,last_col Integer column indices (1-based) of the first and
#'   last TMA column in the sheet. Default \code{2} and \code{11} (columns
#'   2–11 for a 10-column TMA map).
#' @param rows_per_core Integer. Number of Excel rows per TMA row. Default
#'   \code{6} (constructionKey). Use \code{3} for PrintingKey.
#' @param first_pid_row Integer. Excel row index of the first patient-ID row.
#'   Default \code{2}.
#' @param flip_rows Logical. Reverse the row order of the extracted matrix
#'   (as Giuseppe does for some slides). Default \code{FALSE}.
#' @param flip_cols Logical. Reverse the column order. Default \code{FALSE}.
#'
#' @return A character matrix (\code{n_tma_rows × n_tma_cols}) of patient IDs.
#'   Empty / NA cells remain \code{NA}.
parseTMAmap <- function(path,
                        sheet        = "TA649 constructionKey",
                        first_col    = 2L,
                        last_col     = 11L,
                        rows_per_core = 6L,
                        first_pid_row = 2L,
                        flip_rows    = FALSE,
                        flip_cols    = FALSE) {

  if (!requireNamespace("openxlsx", quietly = TRUE))
    stop("Package 'openxlsx' is required.")

  # Read full sheet as raw matrix (Giuseppe's exact call)
  raw <- as.matrix(openxlsx::read.xlsx(path, sheet = sheet, colNames = FALSE))

  # Strip everything after the first comma in each cell (Giuseppe's sub(",.*","",x))
  raw <- matrix(sub(",.*", "", raw), nrow = nrow(raw), ncol = ncol(raw))

  # Identify patient-ID rows: every rows_per_core-th row starting at first_pid_row
  pid_rows <- seq(first_pid_row, nrow(raw), by = rows_per_core)

  tma_mat <- raw[pid_rows, first_col:last_col, drop = FALSE]

  # Replace empty strings / "NA" strings with NA
  tma_mat[tma_mat == "" | tma_mat == "NA" | tma_mat == "NULL"] <- NA

  if (flip_rows) tma_mat <- tma_mat[nrow(tma_mat):1, , drop = FALSE]
  if (flip_cols) tma_mat <- tma_mat[, ncol(tma_mat):1, drop = FALSE]

  rownames(tma_mat) <- NULL
  colnames(tma_mat) <- NULL
  tma_mat
}


#' Resize the in-plot text labels of a ggplot
#'
#' The core-label size in the \code{fov_core_plots} produced by
#' \code{assignCoresToCells()} (via \code{gbspatial}) is baked into the
#' \code{geom_text}/\code{geom_label} layer at creation, so \code{theme()} cannot
#' change it. This helper overrides the \code{size} of every text/label layer in
#' a ggplot, returning the modified plot. Point/line layers are untouched.
#'
#' @param p A ggplot object (e.g. \code{obj@misc$fov_core_plots[[1]]$Whole_Slide_by_Row}).
#' @param size Numeric label text size (geom units). Default \code{2}.
#' @return The ggplot with text/label layers resized.
setPlotLabelSize <- function(p, size = 2) {
  text_geoms <- c("GeomText", "GeomLabel", "GeomTextRepel", "GeomLabelRepel")
  for (i in seq_along(p$layers))
    if (inherits(p$layers[[i]]$geom, text_geoms))
      p$layers[[i]]$aes_params$size <- size
  p
}


#' Warn if a TMA map's orientation does not match the slide's FOV layout
#'
#' Compares the shape of the supplied \code{tma_map} (rows x cols) with the
#' physical FOV grid inferred from cell centroids. TMA cores are laid out on a
#' roughly regular grid, so the map's column:row aspect should match the slide's
#' width:height aspect. When it does not — or when the map has more columns than
#' the slide has FOV-columns (or more rows than FOV-rows) — the map most likely
#' needs transposing (\code{t()}) and/or flipping. Emits a \code{warning()} only;
#' never stops. Flips (row/column reversal) cannot be detected automatically and
#' must be confirmed against \code{obj@misc$fov_core_plots}.
#'
#' @param tma_map Matrix/data.frame as passed to \code{assignCoresToCells()}.
#' @param cell_df Data.frame with \code{CenterX_global_px}, \code{CenterY_global_px}.
#' @param fov_size FOV width/height in pixels.
#' @param slide Character label for messages.
#' @return Invisibly \code{NULL}; called for its warning side effect.
checkTMAorientation <- function(tma_map, cell_df, fov_size, slide = "") {
  x <- cell_df$CenterX_global_px; y <- cell_df$CenterY_global_px
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 2) return(invisible(NULL))
  x <- x[ok]; y <- y[ok]

  # Physical FOV grid: number of distinct FOV columns / rows occupied on the slide
  n_fov_cols <- length(unique(round((x - min(x)) / fov_size)))
  n_fov_rows <- length(unique(round((y - min(y)) / fov_size)))
  R <- nrow(tma_map); C <- ncol(tma_map)
  tag <- if (nzchar(slide)) paste0("[", slide, "] ") else ""

  # (1) Hard check: a core grid can't be finer than the FOV grid on either axis.
  if (C > n_fov_cols || R > n_fov_rows) {
    swap_fits <- (R <= n_fov_cols && C <= n_fov_rows)
    warning(tag, "TMA map is ", R, " rows x ", C, " cols, but the slide spans only ",
            "~", n_fov_cols, " FOV-columns x ", n_fov_rows, " FOV-rows. ",
            "Orientation/shape looks wrong",
            if (swap_fits) " — the transposed map fits, try t(tma_map)." else
              " — check the map extraction (rows_per_core / first_col / last_col).",
            call. = FALSE)
    return(invisible(NULL))
  }

  # (2) Soft check: does the transposed aspect match the slide much better?
  asp_phys  <- n_fov_cols / n_fov_rows
  d_given   <- abs(log((C / R) / asp_phys))
  d_swapped <- abs(log((R / C) / asp_phys))
  if (R != C && d_swapped + log(1.3) < d_given)
    warning(tag, "TMA map orientation (", R, " x ", C, ") fits the slide layout ",
            "(~", n_fov_cols, " cols x ", n_fov_rows, " rows) poorly; the transposed ",
            "map matches much better — consider t(tma_map). Verify against ",
            "obj@misc$fov_core_plots.", call. = FALSE)
  invisible(NULL)
}


#' Assign TMA Core Identities to Cells via FOV-to-Core Mapping
#'
#' Wraps \code{gbspatial::assign_fovs_to_cores()} to map each FOV to a TMA
#' core grid position, then propagates the assignment to every cell in the
#' Seurat object. Strictly reproduces Giuseppe's approach but takes a Seurat
#' object as input and returns a Seurat object.
#'
#' \strong{Algorithm (mirrors Giuseppe exactly):}
#' \enumerate{
#'   \item Build a regular grid of \code{n_rows × n_cols} anchor points from
#'     the extent of the FOV positions and the TMA map dimensions.
#'   \item Assign each cell to the nearest anchor (initial assignment).
#'   \item Refine anchors using cell-centroid means; discard anchors with
#'     \code{< min_cells} cells or excessive drift.
#'   \item Re-assign cells to the nearest \emph{valid} (refined) anchor.
#'   \item Assign each FOV to a core by majority vote of its cells; smooth
#'     boundary FOVs by neighbourhood majority.
#'   \item Look up the sample ID from the TMA map matrix using the core
#'     grid position (row, col).
#'   \item Join FOV-level core assignments back to all cells.
#' }
#'
#' @param obj A Seurat object, typically after \code{QCfiltering()}. Must
#'   contain \code{slidename}, \code{fov} (integer), \code{CenterX_global_px},
#'   and \code{CenterY_global_px} in \code{@meta.data}.
#' @param tma_map A matrix or data.frame where rows = TMA rows, columns = TMA
#'   columns, and each cell contains the sample identifier (e.g. patient ID).
#'   For multi-slide objects, pass a named list of such matrices — one per
#'   slide in the same order as \code{unique(obj$slidename)} — or a single
#'   matrix recycled for all slides.
#' @param fov_positions A data.frame with columns \code{FOV} (integer),
#'   \code{x_global_px}, \code{y_global_px} (top-left corner of each FOV box),
#'   or a file path to \code{*_fov_positions_file.csv.gz}. For multi-slide
#'   objects, pass a list of data.frames or paths in slide order.
#' @param fov_size Integer. FOV width/height in pixels. Default \code{4256}.
#' @param core_drift_tolerance Numeric in (0,1). Maximum allowed fraction of
#'   a core-grid cell width/height that a refined anchor may drift from the
#'   initial grid anchor before the core is discarded. Default \code{0.4}.
#' @param min_cells Integer. Minimum number of cells required for a core to be
#'   considered valid. Default \code{50}.
#' @param sample_id_col Character. Name of the new \code{@meta.data} column
#'   that will receive the sample identifier from \code{tma_map}. Default
#'   \code{"sample_id"}.
#'
#' @return The Seurat object with four new \code{@meta.data} columns:
#'   \code{core_str} (e.g. \code{"C3R2"}), \code{core_col} (integer),
#'   \code{core_row} (integer), and \code{sample_id_col} (the value from
#'   \code{tma_map}). Cells whose FOV could not be assigned receive \code{NA}.
#'   The full FOV-level mapping is stored in \code{obj@misc$fov_core_mapping}
#'   and diagnostic plots in \code{obj@misc$fov_core_plots}.
assignCoresToCells <- function(obj,
                               tma_map,
                               fov_positions,
                               fov_size             = 4256,
                               core_drift_tolerance = 0.4,
                               min_cells            = 50,
                               sample_id_col        = "sample_id") {

  if (!requireNamespace("gbspatial", quietly = TRUE))
    stop("Package 'gbspatial' is required.")

  md       <- obj@meta.data
  slides   <- unique(md$slidename)
  n_slides <- length(slides)

  # Normalise tma_map → list of length n_slides
  if (is.matrix(tma_map) || is.data.frame(tma_map))
    tma_map_list <- rep(list(tma_map), n_slides)
  else if (is.list(tma_map))
    tma_map_list <- tma_map
  else
    stop("tma_map must be a matrix, data.frame, or list thereof.")
  if (length(tma_map_list) != n_slides)
    stop("tma_map must have 1 entry or one entry per slide (", n_slides, " slides).")

  # Normalise fov_positions → list of length n_slides
  if (is.data.frame(fov_positions) || (is.character(fov_positions) && length(fov_positions) == 1))
    fov_list <- rep(list(fov_positions), n_slides)
  else if (is.list(fov_positions) || (is.character(fov_positions) && length(fov_positions) > 1))
    fov_list <- as.list(fov_positions)
  else
    stop("fov_positions must be a data.frame, file path, or list thereof.")
  if (length(fov_list) != n_slides)
    stop("fov_positions must have 1 entry or one entry per slide (", n_slides, " slides).")

  # Build cell_input list — one data.frame per slide
  # assign_fovs_to_cores requires columns: fov, CenterX_global_px, CenterY_global_px
  cell_list <- lapply(slides, function(s) {
    md[md$slidename == s, c("fov", "CenterX_global_px", "CenterY_global_px"), drop = FALSE]
  })

  # Orientation sanity check: warn (per slide) if the TMA map likely needs
  # transposing/reorienting to match the physical FOV layout on the slide.
  for (i in seq_along(slides))
    checkTMAorientation(tma_map_list[[i]], cell_list[[i]], fov_size, slide = slides[i])

  message("Assigning FOVs to TMA cores...")
  result <- gbspatial::assign_fovs_to_cores(
    fov_input            = fov_list,
    cell_input           = cell_list,
    tma_map_input        = tma_map_list,
    fov_size             = fov_size,
    core_drift_tolerance = core_drift_tolerance,
    min_cells            = min_cells,
    slidelabels          = slides
  )

  # Join FOV-level assignments back to cells via slidename + fov integer.
  # assign_fovs_to_cores returns mapped_data with original_FOV (character of
  # the integer FOV number) and slidename — use both to form a unique join key.
  fov_map          <- result$mapped_data
  fov_map$join_key <- paste0(fov_map$slidename, "_", as.integer(fov_map$original_FOV))
  cell_keys        <- paste0(md$slidename, "_", md$fov)
  idx              <- match(cell_keys, fov_map$join_key)

  obj@meta.data$core_str         <- fov_map$core_str[idx]
  obj@meta.data$core_col         <- fov_map$core_col[idx]
  obj@meta.data$core_row         <- fov_map$core_row[idx]
  obj@meta.data[[sample_id_col]] <- fov_map$id[idx]

  obj@misc$fov_core_mapping <- result$mapped_data
  obj@misc$fov_core_plots   <- result$plots

  n_assigned   <- sum(!is.na(obj@meta.data$core_str))
  n_unassigned <- ncol(obj) - n_assigned
  # round to 2 dp so a handful of unassigned cells never displays as a clean 100%
  message("Core assignment complete: ", n_assigned, " / ", ncol(obj),
          " cells assigned (", round(100 * n_assigned / ncol(obj), 2), "%); ",
          n_unassigned, " unassigned.")
  message("Unique cores assigned: ",
          length(unique(na.omit(obj@meta.data$core_str))))

  # Surface unassigned cells explicitly — small FOVs (< min_cells) or dropped
  # anchors are otherwise hidden by rounding (e.g. FOV 182 in TA649).
  if (n_unassigned > 0) {
    na_fovs <- sort(unique(
      paste0(obj@meta.data$slidename, ":FOV", obj@meta.data$fov)[is.na(obj@meta.data$core_str)]))
    warning(n_unassigned, " cell(s) have no core assignment (NA) — usually FOVs with ",
            "< min_cells (", min_cells, ") cells or anchors dropped by ",
            "core_drift_tolerance. Affected: ", paste(head(na_fovs, 20), collapse = ", "),
            if (length(na_fovs) > 20) ", ..." else "",
            ". Lower min_cells, raise core_drift_tolerance, or assign these FOVs manually.",
            call. = FALSE)
  }
  obj
}


#' Manually correct a core / sample assignment for one or more FOVs
#'
#' After \code{assignCoresToCells()}, a small or boundary FOV can be mis-assigned
#' (e.g. FOV 182 in TA649). This overwrites \code{core_str} and/or the sample-ID
#' column for every cell in the given FOV(s), and — when the new \code{core_str}
#' looks like \code{"C<col>R<row>"} — keeps \code{core_col}/\code{core_row}
#' consistent (which a bare \code{ifelse} overwrite does not).
#'
#' @param obj A Seurat object after \code{assignCoresToCells()}.
#' @param fov Integer vector of FOV id(s) to correct (matched against
#'   \code{@meta.data$fov}).
#' @param core_str New core label (e.g. \code{"C2R3"}), or \code{NULL} to leave
#'   \code{core_str}/\code{core_col}/\code{core_row} unchanged.
#' @param sample_id New sample identifier, or \code{NULL} to leave unchanged.
#' @param slide Optional slide name to disambiguate FOV ids in multi-slide
#'   objects (matched against \code{@meta.data$slidename}). Default \code{NULL}.
#' @param sample_id_col Name of the sample-ID column. Default \code{"sample_id"}.
#' @return The Seurat object with the corrected cells.
reassignCore <- function(obj, fov, core_str = NULL, sample_id = NULL,
                         slide = NULL, sample_id_col = "sample_id") {
  md  <- obj@meta.data
  sel <- md$fov %in% fov
  if (!is.null(slide)) sel <- sel & md$slidename %in% slide
  if (!any(sel)) {
    warning("No cells matched fov = ", paste(fov, collapse = ", "),
            if (!is.null(slide)) paste0(" on slide ", slide) else "", call. = FALSE)
    return(obj)
  }

  if (!is.null(core_str)) {
    obj@meta.data$core_str[sel] <- core_str
    # keep core_col / core_row consistent when core_str is "C<col>R<row>"
    m <- regmatches(core_str, regexec("^C(\\d+)R(\\d+)$", core_str))[[1]]
    if (length(m) == 3L) {
      obj@meta.data$core_col[sel] <- as.integer(m[2])
      obj@meta.data$core_row[sel] <- as.integer(m[3])
    } else {
      warning("core_str '", core_str, "' is not 'C<col>R<row>' — core_col/core_row left unchanged.",
              call. = FALSE)
    }
  }
  if (!is.null(sample_id)) obj@meta.data[[sample_id_col]][sel] <- sample_id

  message("Reassigned ", sum(sel), " cells (fov ", paste(fov, collapse = ", "), ") -> core_str=",
          if (is.null(core_str)) "(unchanged)" else core_str, ", ", sample_id_col, "=",
          if (is.null(sample_id)) "(unchanged)" else sample_id)
  obj
}


#' Batch-corrected Dimensionality Reduction with scPearsonPCA
#'
#' Runs the batch-aware PCA -> UMAP -> graph -> clustering pipeline used
#' throughout the analysis: builds a batch variable from one or more
#' metadata columns, computes per-batch gene frequencies
#' (\code{scPearsonPCA::gene_frequency()}), selects variable features
#' (\code{Seurat::FindVariableFeatures()}), computes a batch-corrected
#' quasi-Poisson PCA (\code{sparse_quasipoisson_pca_seurat_batch()}), derives
#' a UMAP + SNN graph (\code{scPearsonPCA::make_umap()}), and clusters with
#' \code{Seurat::FindClusters()}.
#'
#' @param obj A Seurat object.
#' @param label Character. Short label used to name the reductions stored on
#'   the object, e.g. \code{label = "despotx"} stores \code{despotx_pca},
#'   \code{despotx_umap}, and \code{despotx_graph}.
#' @param assay Character. Assay to run the reduction on. Default is the
#'   object's current \code{DefaultAssay(obj)}.
#' @param batch_vars Character vector of \code{obj@meta.data} columns
#'   combined (dash-separated) into the batch variable used for gene
#'   frequency and PCA batch correction. Default \code{c("sample_id", "fov")}.
#' @param nfeatures Integer. Number of variable features passed to
#'   \code{FindVariableFeatures()}. Default \code{2000}.
#' @param resolution Numeric (optionally vector). Resolution(s) passed to
#'   \code{Seurat::FindClusters()}. Default \code{1}. When a vector is given,
#'   \code{clusters} is derived from the last resolution in the vector.
#' @param algorithm Integer. Clustering algorithm passed to
#'   \code{Seurat::FindClusters()}. Default \code{4} (SLM).
#' @param cell_id_col Character. Column in \code{obj@meta.data} holding the
#'   unique cell ID used by the \code{scPearsonPCA} functions. Default
#'   \code{"cell_id"}.
#' @param scale.max,do.scale,do.center,ncores Passed through to
#'   \code{sparse_quasipoisson_pca_seurat_batch()}.
#'
#' @return \code{obj} with reductions \code{<label>_pca}, \code{<label>_umap},
#'   \code{<label>_graph}; the clustering columns added by
#'   \code{FindClusters()} (one \code{<label>_graph_res.<r>} per resolution);
#'   and a \code{clusters} metadata column (\code{"G<n>"} factor) from the
#'   last resolution.
runScPearsonDimRed <- function(obj,
                                label,
                                assay       = Seurat::DefaultAssay(obj),
                                batch_vars  = c("sample_id", "fov"),
                                nfeatures   = 2000,
                                resolution  = 1,
                                algorithm   = 4,
                                cell_id_col = "cell_id",
                                scale.max   = 10,
                                do.scale    = TRUE,
                                do.center   = TRUE,
                                ncores      = 8) {

  Seurat::DefaultAssay(obj) <- assay

  # 1. Batch variable
  # ==================
  batch_col <- paste0(label, "_batch")
  obj@meta.data[[batch_col]] <- do.call(paste, c(obj@meta.data[batch_vars], sep = "-"))
  obs_dt <- data.table::as.data.table(obj@meta.data[, c(cell_id_col, batch_col)])

  # 2. Gene frequency & variable features
  # ======================================
  message("Computing gene frequency for assay '", assay, "'...")
  genefreq_batch <- scPearsonPCA::gene_frequency(
    obj[[assay]]$counts,
    obs            = obs_dt,
    cellid_colname = cell_id_col,
    batch_variable = batch_col
  )

  obj  <- Seurat::FindVariableFeatures(obj, nfeatures = nfeatures)
  hvgs <- Seurat::VariableFeatures(obj, assay = assay)

  # 3. Batch-corrected PCA
  # =======================
  message("Computing batch-corrected PCA ('", label, "')...")
  tc <- Matrix::colSums(obj[[assay]]$counts)
  pcaobj_batch <- sparse_quasipoisson_pca_seurat_batch(
    obj[[assay]]$counts[hvgs, ],
    totalcounts    = tc,
    grate          = genefreq_batch[hvgs, ],
    obs            = obs_dt,
    batch_variable = batch_col,
    cellid_colname = cell_id_col,
    scale.max      = scale.max,
    do.scale       = do.scale,
    do.center      = do.center,
    ncores         = ncores
  )

  # 4. UMAP + graph
  # ================
  message("Computing UMAP ('", label, "')...")
  umapobj <- scPearsonPCA::make_umap(pcaobj_batch)

  pca_name   <- paste0(label, "_pca")
  umap_name  <- paste0(label, "_umap")
  graph_name <- paste0(label, "_graph")

  obj[[pca_name]]   <- pcaobj_batch$reduction.data
  obj[[umap_name]]  <- umapobj$ump
  obj[[graph_name]] <- Seurat::as.Graph(umapobj$grph)

  # 5. Clustering
  # ==============
  message("Clustering ('", label, "', resolution = ", paste(resolution, collapse = ", "), ")...")
  obj <- Seurat::FindClusters(
    obj,
    resolution = resolution,
    graph      = graph_name,
    algorithm  = algorithm
  )

  res_col <- paste0(graph_name, "_res.", tail(resolution, 1))
  obj@meta.data$clusters <- factor(paste0("G", obj@meta.data[[res_col]]))

  obj
}
