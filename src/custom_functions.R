
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
#'   FOV per slide containing centroids and (if available) segmentation polygons,
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
    
    files <- dir(current_path)
    
    # a. Metadata
    message("Metadata...")
    meta_file <- files[grepl("metadata_file", files)]
    if (length(meta_file) == 0) stop("No metadata file found for: ", slidename)
    meta <- data.table::fread(
      file.path(current_path, meta_file), 
      showProgress = FALSE)
    meta[, `:=`(
      slidename        = slidename,
      slide_ID_numeric = i,
      global_cell_ID   = paste0("c_", i, "_", fov, "_", cell_ID),
      FOV              = paste0("s", i, "f", fov)
    )]
    
    # b. Polygons
    message("Polygons...")
    poly_file <- files[grepl("polygons", files)]
    polygons  <- NULL
    if (length(poly_file) > 0) {
      polygons <- data.table::fread(
        file.path(current_path, poly_file), 
        showProgress = FALSE)
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
    message("Counts...")
    counts_file <- files[grepl("exprMat_file", files)]
    if (length(counts_file) == 0) stop("No exprMat file found for: ", slidename)
    dt <- data.table::fread(
      file.path(current_path, counts_file),
      showProgress = FALSE)
    gene_cols <- setdiff(colnames(dt), c("fov", "cell_ID"))
    cell_ids  <- paste0("c_", i, "_", dt$fov, "_", dt$cell_ID)
    
    # Melt to triplet (non-zero only), build sparse matrix without dense intermediate
    long <- data.table::melt(dt[, c("fov", "cell_ID", gene_cols), with = FALSE],
                             id.vars         = c("fov", "cell_ID"),
                             variable.name   = "gene",
                             value.name      = "count",
                             variable.factor = FALSE)[count > 0]
    counts_mat <- Matrix::sparseMatrix(
      i        = match(paste0("c_", i, "_", long$fov, "_", long$cell_ID), cell_ids),
      j        = match(long$gene, gene_cols),
      x        = long$count,
      dims     = c(length(cell_ids), length(gene_cols)),
      dimnames = list(cell_ids, gene_cols)
    )
    counts_mat <- counts_mat[match(meta$global_cell_ID, rownames(counts_mat)), ]
    list(meta = meta, polygons = polygons, counts = counts_mat)
  })
  
  # 3. Combine across slides
  # ========================
  message("Combining across slides...")
  all_genes <- Reduce(union, lapply(slide_data, \(s) colnames(s$counts)))
  if (length(unique(lapply(slide_data, \(s) sort(colnames(s$counts))))) > 1)
    warning("Gene panels differ across TMA slides — missing genes filled with 0.")
  
  counts <- do.call(rbind, lapply(slide_data, function(s) {
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
  }))
  
  metadata <- data.table::rbindlist(lapply(slide_data, \(s) s$meta),  fill = TRUE)
  polygons <- data.table::rbindlist(lapply(slide_data, \(s) s$polygons), fill = TRUE)
  
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
  is_neg      <- grepl("Negative",      colnames(counts))
  is_sys      <- grepl("SystemControl", colnames(counts))
  negcounts   <- counts[,  is_neg,           drop = FALSE]
  falsecounts <- counts[,  is_sys,           drop = FALSE]
  counts      <- counts[, !is_neg & !is_sys, drop = FALSE]
  
  # 6. Build Seurat object
  # ======================
  message("Building Seurat object...")
  meta_df           <- as.data.frame(metadata)
  rownames(meta_df) <- meta_df$cell_ID
  meta_df$orig.ident <- meta_df$slidename  # one TMA = one identity

  # counts is cells x genes — Seurat expects genes x cells
  obj <- Seurat::CreateSeuratObject(
    counts    = t(counts),
    assay     = "RNA",
    meta.data = meta_df
  )

  # Set active idents to orig.ident (TMA / slide identity)
  Seurat::Idents(obj) <- obj$orig.ident

  # Negative probes and system controls as dedicated assays
  obj[["negprobes"]] <- Seurat::CreateAssayObject(counts = t(negcounts))
  obj[["falsecode"]] <- Seurat::CreateAssayObject(counts = t(falsecounts))
  
  # Add spatial FOV per slide
  for (s in slidenames) {
    slide_meta <- metadata[slidename == s]
    slide_poly <- polygons[Run_Tissue_name == s]
    
    cents <- SeuratObject::CreateCentroids(data.frame(
      x    = slide_meta$x_slide_mm,
      y    = slide_meta$y_slide_mm,
      cell = slide_meta$cell_ID
    ))
    
    if (!is.null(slide_poly) && nrow(slide_poly) > 0) {
      segs <- SeuratObject::CreateSegmentation(data.frame(
        x    = slide_poly$x_slide_mm,
        y    = slide_poly$y_slide_mm,
        cell = slide_poly$cell_ID
      ))
      fov_coords <- list(centroids = cents, segmentation = segs)
      fov_types  <- c("segmentation", "centroids")
    } else {
      fov_coords <- list(centroids = cents)
      fov_types  <- "centroids"
    }
    
    obj[[s]] <- SeuratObject::CreateFOV(
      coords = fov_coords,
      type   = fov_types,
      assay  = "RNA"
    )
  }

  # Store polygon local-px coordinates for downstream QC (e.g. computeSplitRatio)
  if (nrow(polygons) > 0)
    obj@misc$polygons <- polygons[, .(cell, slidename, x_local_px, y_local_px)]

  # Pre-compute condensed tissue layout for use with TmaPlot()
  obj <- condenseTissues(obj)

  return(obj)
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
    p <- p + ggplot2::scale_colour_manual(values = col_vals) +
             ggplot2::scale_fill_manual(values = col_vals)
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

  print(p)
  invisible(p)
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
#'   \code{slidename}, \code{x_local_px}, \code{y_local_px}). The
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

  # 1. Boundary detection — per slide, using each slide's own local-px extremes
  # ----------------------------------------------------------------------------
  boundary_cells <- unique(unlist(lapply(
    split(polygons, polygons$slidename),
    function(s) {
      s$cell[
        s$x_local_px %in% c(min(s$x_local_px), max(s$x_local_px)) |
        s$y_local_px %in% c(min(s$y_local_px), max(s$y_local_px))
      ]
    }
  )))
  
  # 2. Per-FOV mean area for normalisation
  has_boundary  <- md$cell_id %in% boundary_cells
  mean_area     <- ave(md$Area, md$fov, FUN = mean)
  
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
#'
#' @return The Seurat object with \code{SBR} and \code{log2SBR} added to
#'   \code{@meta.data}.
computeSBR <- function(obj, bandwidth = 0.01, weight_cutoff = 0.08) {

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

  sbr             <- smoothed_signal / (smoothed_neg + 1e-9)
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

  # tissue: per-cell Panel label from metadata. runFOVQC uses this to build
  # unique FOV IDs as paste0(tissue, fov) and to separate multi-slide runs.
  # The unique value also serves as the barcode map key (panel_name).
  if (!"Panel" %in% colnames(md))
    stop("No 'Panel' column found in @meta.data. ",
         "Ensure the Seurat object was built with readCosMx().")
  tissue_vec <- md$Panel

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
  cell_scores[is.na(cell_scores)] <- 0
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
#' }
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
#' @param filter Logical. Return filtered object. If \code{FALSE}, returns the
#'   object with flags added but no cells removed. Default \code{TRUE}.
#'
#' @return The Seurat object with per-cell flag columns added to
#'   \code{@meta.data} (\code{flag_nCount}, \code{flag_area},
#'   \code{flag_fov_integrity}, \code{flag_boundary}, \code{flag_sbr},
#'   \code{flag_overall}) and, if \code{filter = TRUE}, flagged cells removed.
QCfiltering <- function(obj,
                        count_min             = 20,
                        count_max             = NULL,
                        area_min              = NULL,
                        area_max              = 30000,
                        split_ratio_min       = NULL,
                        split_ratio_max       = abs(log2(0.5)),
                        sbr_min               = 0,
                        sbr_max               = NULL,
                        panel_name            = "Hs_6k",
                        fov_integrity_threshold = 0.6,
                        do_nCount             = TRUE,
                        do_area               = TRUE,
                        do_fov_integrity      = TRUE,
                        do_boundary           = TRUE,
                        do_sbr                = TRUE,
                        filter                = TRUE) {

  n_start <- ncol(obj)

  # 1. nCount_RNA — flag below min and/or above max (NULL = no bound)
  obj@meta.data$flag_nCount <- if (do_nCount) {
    flag <- rep(FALSE, ncol(obj))
    if (!is.null(count_min)) flag <- flag | obj@meta.data$nCount_RNA < count_min
    if (!is.null(count_max)) flag <- flag | obj@meta.data$nCount_RNA > count_max
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
      message("FOV integrity already computed — re-applying threshold.")
      cell_fovid <- paste0(obj@meta.data$Panel, obj@meta.data$FOV)
      obj@meta.data$flag_fov_integrity <- cell_fovid %in% obj@misc$fov_integrity$flaggedfovs
    }
  } else {
    obj@meta.data$flag_fov_integrity <- FALSE
  }

  # 4. FOV boundary (SplitRatioToLocal)
  if (do_boundary) {
    if (!"SplitRatioToLocal" %in% colnames(obj@meta.data)) {
      obj <- computeSplitRatio(obj)
    } else {
      message("SplitRatioToLocal already computed — re-applying threshold.")
    }
    obj@meta.data$flag_boundary <- {
      flag <- rep(FALSE, ncol(obj))
      sr   <- log2(obj@meta.data$SplitRatioToLocal)
      if (!is.null(split_ratio_min)) flag <- flag | sr > split_ratio_min
      if (!is.null(split_ratio_max)) flag <- flag & sr < split_ratio_max
      flag
    }
  } else {
    obj@meta.data$flag_boundary <- FALSE
  }
  if (!is.null(split_ratio_max)) obj = subset(obj, log2(SplitRatioToLocal)<split_ratio_max) %>% suppressWarnings()
  

  # 5. Regional SBR
  if (do_sbr) {
    if (!"log2SBR" %in% colnames(obj@meta.data)) {
      obj <- computeSBR(obj)
    } else {
      message("log2SBR already computed — re-applying threshold.")
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

  # Overall flag
  obj@meta.data$flag_overall <-
    obj@meta.data$flag_nCount       |
    obj@meta.data$flag_area         |
    obj@meta.data$flag_fov_integrity |
    obj@meta.data$flag_boundary     |
    obj@meta.data$flag_sbr

  # Summary table
  md      <- obj@meta.data
  n_total <- nrow(md)

  fmt_thr <- function(mn, mx, fmt = "%g") {
    lo <- if (!is.null(mn)) paste0(">", sprintf(fmt, mn)) else NULL
    hi <- if (!is.null(mx)) paste0("<", sprintf(fmt, mx)) else NULL
    if (is.null(lo) && is.null(hi)) "—" else paste(c(lo, hi), collapse = " & ")
  }

  flags <- list(md$flag_nCount, md$flag_area, md$flag_fov_integrity,
                md$flag_boundary, md$flag_sbr)

  # unique: cells flagged by this filter only (not already caught by prior filters)
  already <- rep(FALSE, n_total)
  unique_counts <- sapply(flags, function(f) {
    u <- sum(f & !already)
    already <<- already | f
    u
  })

  summary_df <- data.frame(
    Filter    = c("nCount_RNA", "Cell area", "FOV integrity",
                  "FOV boundary", "Regional SBR", "Overall"),
    Threshold = c(
      fmt_thr(count_min, count_max),
      fmt_thr(area_min, area_max),
      paste0("prop_loss>", fov_integrity_threshold),
      fmt_thr(split_ratio_min, split_ratio_max),
      fmt_thr(sbr_min, sbr_max),
      "—"
    ),
    Flagged      = c(sapply(flags, sum), sum(md$flag_overall)),
    Flagged_only = c(unique_counts, sum(md$flag_overall)),
    Pct          = round(100 * c(sapply(flags, mean), mean(md$flag_overall)), 2)
  )
  message("\n--- QC summary (", n_total, " cells) ---")
  print(summary_df, row.names = FALSE)

  # UpSet plot — mirrors Giuseppe's run_spatial_qc() approach:
  # build a named list of cell IDs per filter, use UpSetR::fromList(),
  # and add an overall pass/fail pie chart inset via patchwork.
  if (requireNamespace("UpSetR", quietly = TRUE) &&
      requireNamespace("patchwork", quietly = TRUE)) {

    cell_ids <- rownames(md)
    filter_list <- list(
      `Low Counts`      = cell_ids[md$flag_nCount],
      `High Cell Area`  = cell_ids[md$flag_area],
      `FOV QC`          = cell_ids[md$flag_fov_integrity],
      `Low SplitRatio`  = cell_ids[md$flag_boundary],
      `Low SBR`         = cell_ids[md$flag_sbr]
    )
    # keep only filters that flagged at least one cell (mirrors Giuseppe's purrr::keep)
    filter_list <- Filter(function(x) length(x) > 0, filter_list)

    if (length(filter_list) > 1) {
      u_plot <- UpSetR::upset(
        UpSetR::fromList(filter_list),
        nintersects = 10,
        order.by    = "freq",
        nsets       = length(filter_list),
        text.scale  = 1.5
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

  # Re-orient plots to match TmaPlot convention: aes(x = y_global, y = x_global).
  # Giuseppe's plots use aes(x = CenterX_global_px, y = CenterY_global_px);
  obj@misc$fov_core_plots <- result$plots

  n_assigned <- sum(!is.na(obj@meta.data$core_str))
  message("Core assignment complete: ", n_assigned, " / ", ncol(obj),
          " cells assigned (", round(100 * n_assigned / ncol(obj), 1), "%)")
  message("Unique cores assigned: ",
          length(unique(na.omit(obj@meta.data$core_str))))
  obj
}
