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
library(lidRviewer)
library(here)
library(glue)
library(sf)
library(terra)
library(ggplot2)

# las_file
las_file <- paste0(
  "/home/ghislain/ExtDrives/nas-amap1.ird.nc/",
  "USBDisk1/LIDAR_NC/lidar_pouembout/IRD_CIRAD/",
  "lidar/Mission LIDAR ADMIR M1 M2 M3 M4/lidars/terra_las/",
  "cloud202bea4d2fe2f57b.las")

# Catalogue
ctg <- readLAScatalog(las_file)
ctg
plot(ctg)
sf_pb <- st_sf(name="Pouembout", geometry=ctg$geometry)
dir_out <- "pb_outputs"
dir.create(here(dir_out))
write_sf(sf_pb, here(dir_out, "pouembout.gpkg"))

las_check(ctg)
# Inconsistent offsets
# Some tiles seem to overlap each other

# Chunk size
opt_chunk_size(ctg) <- 150
plot(ctg, chunk=TRUE)

# DTM and DSM
pipeline <- reader_las() +
  chm(1, ofile=file.path(out_dir_ssd, "dsm_2P_full.tif")) +
  dtm(1, ofile=file.path(out_dir_ssd, "dtm_2P_full.tif"))
set_exec_options(progress=TRUE, ncores=sequential())
exec(pipeline, on=ctg)

# CHM
dsm <- terra::rast(file.path(out_dir_ssd, "dsm_2P_aoi.tif"))
dtm <- terra::rast(file.path(out_dir_ssd, "dtm_2P_aoi.tif"))
chm <- dsm - dtm
chm[chm < 0] <- 0
varnames(chm) <- "chm"
names(chm) <- "height"
# Save CHM raster file
terra::writeRaster(
  chm,
  file.path(out_dir_ssd, "chm_2P_aoi.tif"),
  datatype="INT1U",
  filetype="GTiff",
  overwrite=TRUE,
  NAflag=255)

# End
