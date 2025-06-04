# Doc
# https://r-lidar.github.io/lidRbook/index.html
# https://r-lidar.github.io/lasR/articles/tutorial.html

# Install lasR package if needed
installed_packages <- as.data.frame(installed.packages())$Package
if(!("lasR" %in% installed_packages)) {
  install.packages("lasR", repos="https://r-lidar.r-universe.dev")
}

# Load libraries
library(lidR)
library(lasR)
#library(lidRviewer)
library(here)
library(glue)
library(sf)
library(terra)
#library(ggplot2)
library(readr)
library(dplyr)
library(tmap)

# The Agathis folder
# Agathis folder in project directory is a symlink to /media/ghislain/Agathis

# Output directory
out_dir_chm <- here("Agathis/outputs/outputs-lidar-chm")
dir.create(out_dir_chm, showWarnings=FALSE)

# Import lidar database
ifile <- here("Agathis", "lidar-flights-newcal.csv")
lidar_df <- readr::read_csv(ifile, show_col_types=FALSE)
nsites <- nrow(lidar_df)

# Compute chm at 1m resolution for each site
for (i in 1:nsites) {
  # Site name and source
  site <- lidar_df$site[i]
  source <- lidar_df$source[i]
  # CHM file
  chm_file <- here(out_dir_chm, glue("chm_{site}.tif"))
  # Run only if chm file does not exist
  if (!file.exists(chm_file)) {
    # LAS files
    las_files <- lidar_df$path_las_files[i]
    if (grepl(";", las_files)) {
      las_files <- unlist(strsplit(las_files, split=";"))
    }
    # Get catalogue of files
    # Here lidR is only used to get files
    # Could be replaced with list files and pattern (cf. help of lasR::exec).
    ctg <- lidR::readLAScatalog(here("Agathis", las_files))
    las_files <- ctg$filename
    
    # DTM and DSM
    # Attention to terminology here: chm() funtion
    # provides a digital surface model (DSM).
    read <- reader_las()
    clean <- delete_noise()
    classify <- classify_with_csf(
      slope_smooth=TRUE,
      class_threshold=0.5,
      cloth_resolution=0.5)
    chm <- chm(1, ofile=here(out_dir_chm, glue("dsm_{site}.tif")))
    dtm <- dtm(1, ofile=here(out_dir_chm, glue("dtm_{site}.tif")))
    ## if (source == "amap") {
      pipeline <- read + clean + classify + chm + dtm
    ## } else {
    ##   pipeline <- read + chm + dtm
    ## }
    set_exec_options(
      progress=TRUE,
      chunk=250,
      buffer=30,
      ncores=concurrent_points(half_cores()))
    exec(pipeline, on=las_files)

    # CHM
    dsm <- terra::rast(here(out_dir_chm, glue("dsm_{site}.tif")))
    dtm <- terra::rast(here(out_dir_chm, glue("dtm_{site}.tif")))
    chm <- dsm - dtm
    chm[chm < 0] <- 0
    varnames(chm) <- "chm"
    names(chm) <- "height"
    # Save CHM raster file
    terra::writeRaster(
      chm,
      here(out_dir_chm, glue("chm_{site}.tif")),
      datatype="INT1U",
      filetype="GTiff",
      overwrite=TRUE,
      NAflag=255)
  }
}
