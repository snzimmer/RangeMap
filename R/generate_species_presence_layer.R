#' Generate raster layers of species presence.
#' 
#' -- Warning! --
#' This feature should be used to analyze areas of presence within a species' known range, not to determine a species' range.
#' RangeMap is likely to indicate presence extending beyond the known range of a species; these results should be viewed with caution!
#' 
#' 
#' @param raster_path Path to RangeMap raster file for a single year
#' @param attributes_path Path to attributes table file. It is preferable to use the RangeMap_Attributes.csv so full field names are preserved, but a tif.vat.dbf file associated with one year's raster may be used
#' @param species_codes Exact species code(s) for which to generate a species presence raster layer. A single raster layer will be generated showing presence of any of the species code(s) provided
#' @param AOI Area of interest. Spatial area for generating raster layers. This can be a polygon or raster loaded into environment (as terra of sf object), or a file path to a .shp or .tif
#' @param output_directory Full directory path for the output rasters (not a file path). Does not need to already exist
#' @param n_cores Optional. Sets the number of cores to use (defaults to 40% of total cores)
#' @param tile_size_adjustment Optional. Adjust sizing of tiles run in parallel. Set this to less than 1 if raster generation fails
#' @return Raster file showing presence of any of the species codes provided
#' @export
#'
#'
generate_species_presence_layer<- function(raster_path,
                                     attributes_path,
                                     species_codes,
                                     AOI,
                                     output_directory,
                                     n_cores = NULL,
                                     tile_size_adjustment = 1){
  options(warn = 1)
  #
  message("Prepping data")
  # Check output directory path
  if (!dir.exists(output_directory)) dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  
  # Set up persistent temp directory for raster tiles
  tile_temp_dir<- file.path(tempdir(), paste0("tiles_", format(Sys.time(), "%Y%m%d%H%M%S")))
  dir.create(tile_temp_dir, showWarnings = FALSE)
  
  # if (file_test("-f", output_directory)) {
  #   stop("output_directory must be a directory, not a file")
  # }
  # #
  # if(!dir.exists(output_directory)){
  #   dir.create(output_directory, recursive = T)
  # }
  
  
  # Check if the attributes path is csv or dbf, then load attributes
  if(endsWith(attributes_path, ".dbf")){
    attributes<- foreign::read.dbf(attributes_path)
    warning("Attributes are from dbf file. Recommend using Attributes csv for full attribute names")
  } else if(endsWith(attributes_path, ".csv")){
    attributes<- utils::read.csv(attributes_path, check.names = F)
  } else{
    message("Expected .csv or .dbf file, something else provided")
  }
  
  
  # Generate new column indicating presence or absence of any of the provided species code(s)
  attributes_species<- attributes[,startsWith(names(attributes), "Species")]
  #
  species_code_pattern <- pattern <- paste0("\\b(", paste(species_codes, collapse = "|"), ")\\b")
  #
  species_presence<- attributes_species |>
    dplyr::rowwise() |>
    dplyr::mutate(presence = any(stringr::str_detect(dplyr::c_across(dplyr::everything()), species_code_pattern), na.rm = TRUE)) |>
    dplyr::ungroup()
  #
  species_presence$presence_numeric<- NA
  species_presence$presence_numeric[species_presence$presence=="TRUE"]<- 1
  
  # Filter to RM_ID and species presence
  fields<- data.frame("RM_ID" = attributes$RM_ID, 
                      "Presence" = species_presence$presence_numeric)
  #
  
  # Write functions to process the files in parallel, with automatic retry of failed tiles
  process_tile <- function(tile_path){
    out_path <- file.path(tile_temp_dir, paste0("processed_", basename(tile_path)))
    
    tryCatch({
      r <- terra::rast(tile_path)
      # Extract raw pixel values as a single vector
      pixel_ids <- terra::values(r, mat = FALSE)
      # Generate attribute lookup table
      id_col_name<- names(final_fields)[1]
      attr_df<- final_fields[, -1, drop = FALSE]
      
      # Run vectorized lookup (maps every pixel ID to corresponding row index in final_fields)
      row_idx <- match(pixel_ids, final_fields[[id_col_name]])
      
      # Create the mapped matrix (rows = pixels, cols = attribute layers)
      mapped_matrix <- as.matrix(attr_df[row_idx, , drop = FALSE])
      
      # Build output raster
      out_raster <- terra::rast(r, nlyrs = ncol(mapped_matrix))
      names(out_raster) <- colnames(mapped_matrix)
      #
      terra::values(out_raster) <- mapped_matrix
      #
      terra::writeRaster(out_raster, out_path,
                         overwrite = TRUE, datatype = "INT4S", NAflag = -9999)
      
      rm(r, out_raster, pixel_ids, row_idx, mapped_matrix); gc()
      
      list(tile = tile_path, status = "success", error = NA_character_)
      
    }, error = function(e){
      # Remove any partial output left behind by a failed tile being written
      if (file.exists(out_path)) file.remove(out_path)
      list(tile = tile_path, status = "failed", error = conditionMessage(e))
    })
  }
  #
  run_tiles <- function(tiles_to_process, max_retries = 5){
    # Track failed tiles so they can be re-run
    all_failed_log <- list()
    
    for (attempt in seq_len(max_retries + 1)){
      if (length(tiles_to_process) == 0) break
      
      cores <- min(n_cores, length(tiles_to_process))
      
      if (cores > 1){
        cl <- parallel::makeCluster(cores)
        parallel::clusterExport(cl, varlist = c("final_fields", "tile_temp_dir"), envir = environment())
        parallel::clusterEvalQ(cl, terra::setGDALconfig("GDAL_PAM_ENABLED", "NO"))
        results <- tryCatch(
          parallel::clusterApplyLB(cl, tiles_to_process, process_tile),
          finally = try(parallel::stopCluster(cl), silent = TRUE)
        )
      } else {
        results <- lapply(tiles_to_process, process_tile)
      }
      
      results_df <- do.call(rbind, lapply(results, as.data.frame, stringsAsFactors = FALSE))
      failed <- results_df$tile[results_df$status == "failed"]
      all_failed_log[[attempt]] <- results_df[results_df$status == "failed", ]
      
      if (length(failed) == 0){
        return(list(failed_tiles = character(0)))
      }
      failed_errors <- results_df$error[results_df$status == "failed"]
      tiles_to_process <- failed
    }
    
    list(failed_tiles = tiles_to_process, log = do.call(rbind, all_failed_log))
  }
  
  # Start running raster data here
  message("Generating species presence raster")
  
  # Load raster
  ras<- terra::rast(raster_path)
  
  # Reproject AOI
  # If AOI is a file path to .shp or .tif, read these in
  if(inherits(AOI, "character") && endsWith(AOI, ".shp")){
    AOI_proj<- terra::project(terra::vect(AOI), terra::crs(ras))
  } else if (inherits(AOI, "character") && endsWith(AOI, ".tif")){
    AOI_proj<- terra::ext(terra::project(terra::rast(AOI), terra::crs(ras)))
  } else if(class(AOI)[1] == "SpatVector"){
    AOI_proj<- terra::project(AOI, terra::crs(ras))
  } else if (class(AOI)[1] %in% c("sf", "SpatialPolygons")){
    AOI_proj<- terra::project(terra::vect(AOI), terra::crs(ras))
  } else if(class(AOI)[1] == "SpatRaster"){
    AOI_proj<- terra::ext(terra::project(AOI, terra::crs(ras)))
  } else{
    message("Expected AOI to be either a path to a .shp or .tif,
            or a terra or sf object. Please provide one of these.")
  }
  
  # Mask / crop the raster to the AOI
  if(class(AOI_proj)[1] == "SpatExtent"){
    ras<- terra::crop(ras, AOI_proj)
  } else if(class(AOI_proj)[1] == "SpatVector"){
    ras<- terra::crop(ras, AOI_proj, mask=T)
  }
  
  # Begin to work on raster generation
  if(is.null(n_cores)){
    n_cores <- max(1, ceiling(parallel::detectCores() * 0.4))
  }
  #
  freeRAM_mb <- terra::free_RAM() / 1024
  budget_bytes <- ((freeRAM_mb * 0.5) / n_cores) * 1024^2
  max_cells <- budget_bytes / (32 * 20)
  tile_dim <- max(floor(sqrt(max_cells)), 500) * tile_size_adjustment
  
  # Generate tiles, only needs to be done once regardless of the number of attributes being generated
  terra::setGDALconfig("GDAL_PAM_ENABLED", "NO")
  tile_files <- terra::makeTiles(
    ras,
    y = c(tile_dim, tile_dim),
    filename = file.path(tile_temp_dir, "tile_.tif"),
    na.rm = TRUE
  )

    #
    # Define the fields (attributes) for this batch
    final_fields<- fields
    #
    tile_result <- run_tiles(tile_files, max_retries = 5)
    
    if (length(tile_result$failed_tiles) > 0){
      stop("Tiles are failing during processing. Try reducing tile_size_adjustment to less than 1!")
    }
    
    # Mosaic all processed tiles into one raster (multiband if mulutple attributes selected)
    # Direct to processed tiles
    processed_tiles <- list.files(
      tile_temp_dir,
      pattern = "processed_tile",
      full.names = TRUE
    )
    
    message("Saving species presence raster")
    terra::setGDALconfig("GDAL_MAX_DATASET_POOL_SIZE", "1000")
    terra::setGDALconfig("GDAL_CACHEMAX","4000")
    
    # Export single mosaicked raster of processed attributes
    vrt_file <- file.path(tile_temp_dir, "vrt.vrt")
    terra::vrt(processed_tiles,vrt_file, set_names = TRUE, overwrite=T)
    r <- terra::rast(vrt_file)
    names(r)<- paste0("Presence_", paste0(paste(species_codes, collapse = "_")))
    #
    out_tif<- file.path(output_directory, paste0("Presence_", paste0(paste(species_codes, collapse = "_"), ".tif")))
    #
    dt<- "INT1U"
    nodata<- 255
      
    if(file.exists(out_tif)){
        warning(paste0(out_tif), " already exists and will be overwritten")
    }
    terra::writeRaster(
        r,
        out_tif,
        datatype = dt,
        NAflag = nodata,
        overwrite = T,
        gdal = c("COMPRESS=DEFLATE", "ZLEVEL=8", "PREDICTOR=2",
                 "TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512",
                 "NUM_THREADS=ALL_CPUS", "SPARSE_OK=YES", "BIGTIFF=YES"))
    
    # Clean up processed tiles
    file.remove(processed_tiles)
    file.remove(vrt_file)
    unlink(tile_temp_dir, recursive = TRUE)
}
    
   

