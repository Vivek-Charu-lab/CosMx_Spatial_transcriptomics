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

