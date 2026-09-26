#' Internal helper to tile rasters
#' @noRd
manual_make_tiles<- function(
    ras,
    tile_dim,
    tile_temp_dir,
    datatype = "INT4S",
    gdal_opts = c(
      "COMPRESS=DEFLATE", "ZLEVEL=8", "PREDICTOR=2",
      "TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512",
      "NUM_THREADS=ALL_CPUS", "SPARSE_OK=YES", "BIGTIFF=YES"),
    na_rm = TRUE) {

  terra::setGDALconfig("GDAL_PAM_ENABLED", "NO")
  #
  nr<- terra::nrow(ras)
  nc<- terra::ncol(ras)
  #
  row_starts<- seq(1, nr, by = tile_dim)
  col_starts<- seq(1, nc, by = tile_dim)
  #
  # Grid expand to vectorize grid creation instead of double nested for-loop
  grid<- expand.grid(r0 = row_starts, c0 = col_starts)
  #
  tile_files<- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    r0<- grid$r0[i]
    c0<- grid$c0[i]

    r1<- min(r0 + tile_dim - 1, nr)
    c1<- min(c0 + tile_dim - 1, nc)

    # Crop raster slice directly by row/col ranges
    tile_ras<- ras[r0:r1, c0:c1, drop = FALSE]

    # Skip empty/all-NA tiles
    if (na_rm && terra::hasValues(tile_ras) && all(is.na(terra::values(tile_ras)))) {
      next
    }
    out_file<- file.path(tile_temp_dir, sprintf("tile_%03d.tif", i))
    #
    terra::writeRaster(
      tile_ras,
      filename  = out_file,
      datatype  = datatype,
      gdal      = gdal_opts,
      overwrite = TRUE
    )
    tile_files[[i]]<- out_file
  }

  # Return non-null created file paths
  unlist(tile_files, use.names = FALSE)
}
