

dataprep_cosmx <- function(myflatfiledir, plot_tissues = FALSE) {
  
  # Automatically get slide paths from the provided directories
  slide_paths <- character()
  for (dir_path in myflatfiledir) {
    s_names <- dir(dir_path)
    # Ensure we only keep directories (ignoring loose files in the parent dir)
    valid_dirs <- s_names[dir.exists(file.path(dir_path, s_names))]
    slide_paths <- c(slide_paths, file.path(dir_path, valid_dirs))
  }
  
  if (length(slide_paths) == 0) {
    stop("No slide directories found in the provided path(s).")
  }
  
  slidenames <- basename(slide_paths)
  #make unique slidenames if duplicates
  slidenames=ave(slidenames, slidenames, FUN = \(v) ifelse(seq_along(v) == 1,v,paste0(v, "-", seq_along(v))))
  
  # Lists to collect the counts matrices, metadata, and polygons, one per slide
  countlist <- vector(mode = 'list', length = length(slide_paths)) 
  metadatalist <- vector(mode = 'list', length = length(slide_paths)) 
  polygonlist <- vector(mode = 'list', length = length(slide_paths))
  _
  for(i in seq_along(slide_paths)) {
    
    current_path <- slide_paths[i]
    slidename <- slidenames[i] 
    
    msg <- paste0("Loading slide ", slidename, ", ", i, "/", length(slide_paths), ".")
    message(msg)    
    
    # slide-specific files:
    thisslidesfiles <- dir(current_path)
    
    # load in metadata:
    thisslidesmetadata <- thisslidesfiles[grepl("metadata\\_file", thisslidesfiles)]
    if (length(thisslidesmetadata) == 0) stop(paste("No metadata file found for", slidename))
    
    tempdatatable <- data.table::fread(file.path(current_path, thisslidesmetadata))
    
    # Use base R assignment instead of data.table := to avoid namespace errors
    tempdatatable$slidename <- slidename
    
    # numeric slide ID 
    slide_ID_numeric <- i
    tempdatatable$slide_ID_numeric <- i
    
    # global cell ID 
    tempdatatable$global_cell_ID <- paste0("c_", slide_ID_numeric, "_", tempdatatable$fov, "_", tempdatatable$cell_ID)
    
    # ALWAYS load polygons if available
    thisslidespolygon <- thisslidesfiles[grepl("polygons", thisslidesfiles)]
    if (length(thisslidespolygon) != 0) {
      polygons <- data.table::fread(file.path(current_path, thisslidespolygon))
      
      # Modify cell_ID to match metadata global_cell_ID format
      polygons$cell_ID <- paste0("c_", slide_ID_numeric, "_", polygons$fov, "_", polygons$cellID)
      
      # Add Run_Tissue_name
      polygons$Run_Tissue_name <- slidename
      
      # Add global FOV
      polygons$FOV <- paste0("s", slide_ID_numeric, "f", polygons$fov)
      
      # Save to list
      polygonlist[[i]] <- polygons
      
      # Create SplitRatioToLocal if it does not exist 
      if (!"SplitRatioToLocal" %in% names(tempdatatable) || all(is.na(tempdatatable$SplitRatioToLocal))) {
        boundarycells=unique(polygons$cell[polygons$x_local_px %in% c(min(polygons$x_local_px), max(polygons$x_local_px)) | polygons$y_local_px %in% c(min(polygons$y_local_px), max(polygons$y_local_px))])
        has_boundary <- tempdatatable$cell_id %in% boundarycells 
        mean_area <- ave(tempdatatable$Area, tempdatatable$fov, FUN = mean) 
        tempdatatable$SplitRatioToLocal <- ifelse(has_boundary, round(tempdatatable$Area / mean_area, 2), 0)
      }
      
    } else {
      # Fallback if no polygon file exists
      message(paste("No polygon file found for", slidename))
      if (!"SplitRatioToLocal" %in% names(tempdatatable) || all(is.na(tempdatatable$SplitRatioToLocal))) {
        message("Cannot generate SplitRatioToLocal because the polygon file is missing.")
      }
    }
    
    # load in counts as a data table:
    thisslidescounts <- thisslidesfiles[grepl("exprMat\\_file", thisslidesfiles)]
    if (length(thisslidescounts) == 0) stop(paste("No exprMat file found for", slidename))
    
    countsfile <- file.path(current_path, thisslidescounts)
    nonzero_elements_perchunk <- 5 * 10^7
    
    ### Safely read in the dense (0-filled) counts matrices in chunks.
    lastchunk <- FALSE 
    skiprows <- 0
    chunkid <- 1
    
    required_cols <- data.table::fread(countsfile, select=c("fov", "cell_ID"))
    stopifnot("columns 'fov' and 'cell_ID' are required, but not found in the counts file" = 
                all(c("cell_ID", "fov") %in% colnames(required_cols)))
    number_of_cells <- nrow(required_cols)
    
    number_of_cols <- ncol(data.table::fread(countsfile, nrows = 2))
    number_of_chunks <- ceiling(number_of_cols * number_of_cells / nonzero_elements_perchunk)
    chunk_size <- floor(number_of_cells / number_of_chunks)
    sub_counts_matrix <- vector(mode = 'list', length = number_of_chunks)
    
    pb <- txtProgressBar(min = 0, max = number_of_chunks, initial = 0, char = "=", style = 3)
    cellcount <- 0
    
    while(lastchunk == FALSE) {
      read_header <- ifelse(chunkid == 1, TRUE, FALSE)
      
      countsdatatable <- data.table::fread(countsfile,
                                           nrows = chunk_size,
                                           skip = skiprows + (chunkid > 1),
                                           header = read_header)
      if(chunkid == 1) {
        header <- colnames(countsdatatable)
      } else {
        colnames(countsdatatable) <- header
      }
      
      cellcount <- nrow(countsdatatable) + cellcount     
      if(cellcount == number_of_cells) lastchunk <- TRUE
      
      skiprows <- skiprows + chunk_size
      slide_fov_cell_counts <- paste0("c_", slide_ID_numeric, "_", countsdatatable$fov, "_", countsdatatable$cell_ID)
      
      # Define columns to keep by subtracting fov and cell_ID safely
      cols_to_keep <- setdiff(colnames(countsdatatable), c("fov", "cell_ID"))
      
      # Create base matrix, enforce numeric mode, and use Matrix() constructor
      dense_mat <- as.matrix(countsdatatable[, cols_to_keep, with = FALSE])
      mode(dense_mat) <- "numeric"
      sub_counts_matrix[[chunkid]] <- Matrix::Matrix(dense_mat, sparse = TRUE)
      rownames(sub_counts_matrix[[chunkid]]) <- slide_fov_cell_counts 
      
      setTxtProgressBar(pb, chunkid)
      chunkid <- chunkid + 1
    }
    
    close(pb)   
    
    countlist[[i]] <- do.call(rbind, sub_counts_matrix) 
    
    # ensure that cell-order in counts matches cell-order in metadata   
    countlist[[i]] <- countlist[[i]][match(tempdatatable$global_cell_ID, rownames(countlist[[i]])), ] 
    metadatalist[[i]] <- tempdatatable 
    
    # track common genes and common metadata columns across slides
    if(i == 1) {
      sharedgenes <- colnames(countlist[[i]]) 
      sharedcolumns <- colnames(tempdatatable)
    } else {
      sharedgenes <- intersect(sharedgenes, colnames(countlist[[i]]))
      sharedcolumns <- intersect(sharedcolumns, colnames(tempdatatable))
    }
  }
  
  # reduce to shared metadata columns and shared genes
  for(i in seq_along(slide_paths)) {
    metadatalist[[i]] <- metadatalist[[i]][, sharedcolumns, with = FALSE]
    countlist[[i]] <- countlist[[i]][, sharedgenes, drop = FALSE]
  }
  
  counts <- do.call(rbind, countlist)
  metadata <- data.table::rbindlist(metadatalist)
  polygons_all <- data.table::rbindlist(polygonlist, fill = TRUE)
  
  # add to metadata: add a global non-slide-specific FOV ID:
  metadata$FOV <- paste0("s", metadata$slide_ID_numeric, "f", metadata$fov)
  
  # remove cell_ID metadata column, which only identifies cell within slides, not across slides:
  metadata$cell_ID <- metadata$global_cell_ID
  metadata$cell <- NULL
  metadata$cell_id <- metadata$cell_ID
  metadata$global_cell_ID <- NULL
  
  # add coordinates in mm
  um_per_px <- 0.120280945 # dimension of each pixel in micro-meter
  metadata$x_slide_mm <- um_per_px * metadata$CenterX_global_px / 1e3 
  metadata$y_slide_mm <- um_per_px * metadata$CenterY_global_px / 1e3
  polygons_all$x_slide_mm <- um_per_px * polygons_all$x_global_px / 1e3 
  polygons_all$y_slide_mm <- um_per_px * polygons_all$y_global_px / 1e3
  
  # isolate negative control matrices:
  negcounts <- counts[, grepl("Negative", colnames(counts)), drop = FALSE]
  falsecounts <- counts[, grepl("SystemControl", colnames(counts)), drop = FALSE]
  
  # reduce counts matrix to only genes:
  counts <- counts[, !grepl("Negative", colnames(counts)) & !grepl("SystemControl", colnames(counts)), drop = FALSE]
  
  # Then break out cells' xy positions in a distinct data object:
  xy <- as.matrix(metadata[, c("CenterX_global_px", "CenterY_global_px"), with = FALSE])
  
  # Use the generated global_cell_ID for rownames to ensure accuracy 
  rownames(xy) <- metadata$cell_ID
  
  # rescale to mm:
  thisinstrument_nanometers_per_pixel = 120.280945   
  xy <- xy * thisinstrument_nanometers_per_pixel / 1000000
  colnames(xy) <- paste0(c("x", "y"), "_mm")
  
  # Condense tissues
  set.seed(1)
  xy <- condenseTissues(xy = xy, 
                        tissue = metadata$Run_Tissue_name, 
                        tissueorder = NULL,  
                        buffer = 1, 
                        widthheightratio = 8/3) 
  
  # Optional: Plot tissues if requested
  if (plot_tissues) {
    sub <- sample(1:nrow(xy), round(nrow(xy) / 20), replace = FALSE)
    plot(xy[sub, ], pch = 16, cex = 0.2, 
         asp = 1, 
         col = as.numeric(as.factor(metadata$Run_Tissue_name[sub])),
         main = "Condensed Tissues Layout")
    
    for (s_name in unique(metadata$Run_Tissue_name)) {
      text(median(xy[metadata$Run_Tissue_name == s_name, 1]), 
           max(xy[metadata$Run_Tissue_name == s_name, 2]), 
           s_name)
    }
    condensed_tissue_plot <- recordPlot()
  }
  
  # Return final list object
  result_obj <- list(
    counts = counts,
    negcounts = negcounts,
    falsecounts = falsecounts,
    metadata = metadata,
    xy = xy,
    polygons = polygons_all,
    condensed_tissue_plot = if (exists("condensed_tissue_plot")) condensed_tissue_plot else NULL
  )
  
  return(result_obj)
}