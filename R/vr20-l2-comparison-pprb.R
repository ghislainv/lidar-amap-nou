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

# Variables
las_file_l2 <- "/media/ghislain/Agathis/20241209_IRD_CIRAD_LIVRABLES_PPRB_LDR_L2_ECOTONE/IRD_CIRAD_PPRB_LDR_L2_M1/terra_las/cloud_merged.las"
las_file_vx20 <- "/media/ghislain/Agathis/20241212-YS-PPRB-P1-test/export/20241212-YS-PPRB-P1-test-20241217-053955-F001.laz"

# Output directory
out_dir_ssd <- "/media/ghislain/Agathis/outputs/l2-vx20-pprb"
dir.create(out_dir_ssd, showWarnings=FALSE)

# ========================================
# CHM with DJI-L2
# ========================================

# Catalogue
ctg <- readLAScatalog(las_file_l2)
ctg
plot(ctg)

# Chunk size
opt_chunk_size(ctg) <- 150
plot(ctg, chunk=TRUE)

# DTM and DSM
pipeline <- reader_las() +
  chm(1, ofile=file.path(out_dir_ssd, "dsm_l2.tif")) +
  dtm(1, ofile=file.path(out_dir_ssd, "dtm_l2.tif"))
set_exec_options(progress=TRUE, ncores=sequential())
exec(pipeline, on=ctg)

# CHM
dsm <- terra::rast(file.path(out_dir_ssd, "dsm_l2.tif"))
dtm <- terra::rast(file.path(out_dir_ssd, "dtm_l2.tif"))
chm <- dsm - dtm
chm[chm < 0] <- 0
varnames(chm) <- "chm"
names(chm) <- "height"
# Save CHM raster file
terra::writeRaster(
  chm,
  file.path(out_dir_ssd, "chm_l2.tif"),
  datatype="INT1U",
  filetype="GTiff",
  overwrite=TRUE,
  NAflag=255)

# ========================================
# CHM with YS-Vx20
# ========================================

# Catalogue
ctg <- readLAScatalog(las_file_vx20)
ctg
plot(ctg)
las_check(ctg)

# Plot transect
las_vx20_pprb <- readLAS(las_file_vx20, select="xyzc")
print(las_vx20_pprb)
p1 <- c(470977, 232806)
p2 <- c(470999, 232667)
las_tr <- clip_transect(las_vx20_pprb, p1, p2, width=2, xz=TRUE)
p <- ggplot(payload(las_tr), aes(X, Z, color = Z)) + 
  geom_point(size = 0.5) + 
  coord_equal() + 
  theme_minimal() +
  scale_color_gradientn(colours = height.colors(50))
ggsave(file.path(out_dir_ssd, "transect.png"), p, width=10, height=2)

# Chunk
opt_chunk_size(ctg) <- 160
plot(ctg, chunk=TRUE)

# Classifying points
ofile <- file.path(out_dir_ssd, "ys-vx20-classified.las")
pipeline <- reader_las() +
  classify_with_csf(
    slope_smooth=TRUE,
    class_threshold=0.5,
    cloth_resolution=0.5) +
  write_las(ofile)
exec(pipeline, on=ctg)

# Plot ground points
las_vx20_pprb <- readLAS(ofile, select="xyzc")
print(las_vx20_pprb)
p1 <- c(470977, 232806)
p2 <- c(470999, 232667)
las_tr <- clip_transect(las_vx20_pprb, p1, p2, width=2, xz=TRUE)
p <- ggplot(payload(las_tr), aes(X, Z, color = Classification)) + 
  geom_point(size = 0.5) + 
  coord_equal() + 
  theme_minimal()
ggsave(file.path(out_dir_ssd, "tr-ground-points.png"), p, width=10, height=2)

# DTM and DSM
pipeline <- reader_las() +
  chm(1, ofile=file.path(out_dir_ssd, "dsm_vx20.tif")) +
  dtm(1, ofile=file.path(out_dir_ssd, "dtm_vx20.tif"))
set_exec_options(progress=TRUE, ncores=sequential())
exec(pipeline, on=ofile)

# CHM
dsm <- terra::rast(file.path(out_dir_ssd, "dsm_vx20.tif"))
dtm <- terra::rast(file.path(out_dir_ssd, "dtm_vx20.tif"))
chm <- dsm - dtm
chm[chm < 0] <- 0
varnames(chm) <- "chm"
names(chm) <- "height"
# Save CHM raster file
terra::writeRaster(
  chm,
  file.path(out_dir_ssd, "chm_vx20.tif"),
  datatype="INT1U",
  filetype="GTiff",
  overwrite=TRUE,
  NAflag=255)

# ========================================
# CHM comparison at 1m resolution
# ========================================




# ========================================
# CHM comparison at 10m resolution
# ========================================

# PPRB P1 plot
source("~/Code/forest-plot-grid/forest-plot-grid.R")
get_grid(coords=c(166.688126, -22.106773),
         azimut=306,
         proj="epsg:3163",
         distance=100,
         bin=10,
         output_dir=out_dir_ssd)
ifile <- file.path(out_dir_ssd, "grid.gpkg")
p1_pprb <- sf::st_read(ifile, quiet=TRUE) |>
  sf::st_transform(crs="epsg:3163")

# Zonal statistics in subplot
chm_l2 <- terra::rast(file.path(out_dir_ssd, "chm_l2.tif"))
chm_vx20 <- terra::rast(file.path(out_dir_ssd, "chm_vx20.tif"))
p1_pprb_l2 <- terra::zonal(
  x=chm_l2,
  z=terra::vect(p1_pprb),
  fun="mean", na.rm=TRUE,
  as.polygons=TRUE)
p1_pprb_l2_vx20 <- terra::zonal(
  x=chm_vx20,
  z=p1_pprb_l2,
  fun="mean", na.rm=TRUE,
  as.polygons=TRUE)
df <- as.data.frame(p1_pprb_l2_vx20) |>
  dplyr::rename(subplot=name, height_l2=2, height_vx20=3)

# Plotting
p <- ggplot(df, aes(x=height_l2, y=height_vx20)) +
  geom_point() +
  xlab("Mean canopy height (m), DJI-L2, 10 x 10 m") +
  ylab("Mean canopy height (m), YS-Vx20, 10 x 10 m") +
  geom_abline(slope=1, intercept=0, colour="red") +
  theme_bw()
ggsave(file.path(out_dir_ssd, "chm_comp.png"), p)

# End
