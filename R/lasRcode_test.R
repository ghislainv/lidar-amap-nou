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

# Variables
las_P1_dir <- "/media/ghislain/Agathis/IRD_CIRAD_LIVRABLES_PPGF_LDR_L2_ECOTONE/IRD_CIRAD_PPGF_LDR_L2_P1"
las_P2_dir <- "/media/ghislain/Agathis/IRD_CIRAD_LIVRABLES_PPGF_LDR_L2_ECOTONE/IRD_CIRAD_PPGF_LDR_L2_P2"

# Files
las_files_P1 <- list.files(las_P1_dir, pattern="*.las", full.names=TRUE)
las_files_P2 <- list.files(las_P2_dir, pattern="*.las", full.names=TRUE)

# Output directory
out_dir_ssd <- "/media/ghislain/Agathis/outputs"
dir.create(out_dir_ssd, showWarnings=FALSE)

# PPGF bounding box
ifile <- "/home/ghislain/Documents/Research_current/PROJECTS/MELANOBS/Data/PPGF/parcelle_ppgf.gpkg"
ppgf <- sf::read_sf(ifile)
bbox_ppgf <- sf::st_bbox(ppgf)

# ========================================
# lidR
# ========================================

# Catalogue
ctg <- readLAScatalog(las_P1_dir)
ctg
plot(ctg)
las_check(ctg)
# Inconsistent offsets
# Some tiles seem to overlap each other

# Clip catalog to ppgf
opt_output_files(ctg) <- paste0("/media/ghislain/Agathis/outputs/ppgf_", "{ID}")
ppgf_las <- clip_roi(ctg, geometry=bbox_ppgf)

# Load las ppgf
ifile <- "/media/ghislain/Agathis/outputs/ppgf_1.las"
las_ppgf <- readLAS(ifile)

# View
plot(las_ppgf, backend="lidRviewer")

# ========================================
# CHM with one pass
# ========================================

# DTM
ifile <- "/media/ghislain/Agathis/outputs/ppgf_1.las"
pipeline <- reader_las(filter=keep_ground()) +
  triangulate() +
  rasterize(1, ofile=file.path(out_dir_ssd, "dtm.tif"))
set_exec_options(progress=TRUE, ncores=sequential())
exec(pipeline, on=ifile)

# DTM and DSM
ifile <- "/media/ghislain/Agathis/outputs/ppgf_1.las"
pipeline <- reader_las() +
  chm(1, ofile=file.path(out_dir_ssd, "dsm.tif")) +
  dtm(1, ofile=file.path(out_dir_ssd, "dtm.tif"))
set_exec_options(progress=TRUE, ncores=sequential())
exec(pipeline, on=ifile)

# CHM
dsm <- terra::rast(file.path(out_dir_ssd, "dsm.tif"))
dtm <- terra::rast(file.path(out_dir_ssd, "dtm.tif"))
chm <- dsm - dtm
varnames(chm) <- "chm"
names(chm) <- "height"
# Save CHM raster file
terra::writeRaster(
  chm,
  file.path(out_dir_ssd, "chm.tif"),
  datatype="INT1U",
  filetype="GTiff",
  overwrite=TRUE,
  NAflag=255)

# Plot CHM
df_chm <- as.data.frame(chm, xy=TRUE)
p_chm <-  ggplot() +
  geom_raster(data=df_chm, 
              aes(x=x, y=y, fill=height)) + 
  scale_fill_gradientn(name="Canopy height (m)", colors=terrain.colors(10, rev=TRUE)) +
  coord_equal()
ggsave(file.path(out_dir_ssd, "chm.pdf"), p_chm)

# ========================================
# CHM with two passes
# ========================================

# Catalogue
ctg <- readLAScatalog(c(las_files_P1, las_files_P2))
ctg
plot(ctg)
las_check(ctg)

# Clip catalogue to ppgf
opt_output_files(ctg) <- paste0("/media/ghislain/Agathis/outputs/ppgf_2P_", "{ID}")
ppgf_las <- clip_roi(ctg, geometry=bbox_ppgf)

# Load las ppgf
ifile <- "/media/ghislain/Agathis/outputs/ppgf_2P_1.las"
las_ppgf <- readLAS(ifile)

# DTM and DSM
ifile <- "/media/ghislain/Agathis/outputs/ppgf_2P_1.las"
pipeline <- reader_las() +
  chm(1, ofile=file.path(out_dir_ssd, "dsm_2P.tif")) +
  dtm(1, ofile=file.path(out_dir_ssd, "dtm_2P.tif"))
set_exec_options(progress=TRUE, ncores=sequential())
exec(pipeline, on=ifile)

# CHM
dsm <- terra::rast(file.path(out_dir_ssd, "dsm_2P.tif"))
dtm <- terra::rast(file.path(out_dir_ssd, "dtm_2P.tif"))
chm <- dsm - dtm
varnames(chm) <- "chm"
names(chm) <- "height"
# Save CHM raster file
terra::writeRaster(
  chm,
  file.path(out_dir_ssd, "chm_2P.tif"),
  datatype="INT1U",
  filetype="GTiff",
  overwrite=TRUE,
  NAflag=255)

# ========================================
# Plot CHM relationship
# ========================================

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# Plot relationship
chm_1P <- terra::rast(file.path(out_dir_ssd, "chm.tif"))
chm_2P <- terra::rast(file.path(out_dir_ssd, "chm_2P.tif"))

# Data-frames
df_1P <- as.data.frame(chm_1P, xy=TRUE)
df_1P$loc <- paste0(df_1P$x, df_1P$y)
df_2P <- as.data.frame(chm_2P, xy=TRUE)
df_2P$loc <- paste0(df_2P$x, df_2P$y)
df_comp <- merge(df_1P, df_2P, by="loc", all.x=TRUE, suffixes=c("_P1", "_P2"))
df_comp$density <- nrow(df_comp) * get_density(df_comp$height_P1, df_comp$height_P2, n=100)

# Plotting
p_comp <- ggplot(aes(height_P1, height_P2, color=density), data=df_comp) +
  geom_point() +
  scale_fill_continuous(type="viridis") +
  xlab("Canopy height 1 pass") +
  ylab("Canopy height 2 passes") +
  coord_equal() +
  theme_bw()
ggsave(file.path(out_dir_ssd, "comparison_1P_2P.pdf"), p_comp)

# ========================================
# After resampling at 10m
# ========================================

# Resample
chm_1P_resamp <- terra::aggregate(x=chm_1P, fact=10, fun="mean")
chm_2P_resamp <- terra::aggregate(x=chm_2P, fact=10, fun="mean")
# Data frame
df_comp <- data.frame(height_P1=as.numeric(values(chm_1P_resamp)),
                      height_P2=as.numeric(values(chm_2P_resamp)))

# Plotting
p_comp <- ggplot(aes(height_P1, height_P2), data=df_comp) +
  geom_point() +
  ggtitle("10m mean canopy height") +
  xlab("Canopy height 1 pass") +
  ylab("Canopy height 2 passes") +
  coord_equal() +
  theme_bw()
ggsave(file.path(out_dir_ssd, "comparison_1P_2P_10m.pdf"), p_comp)

# ========================================
# Full analysis
# ========================================

# Catalogue
ctg <- readLAScatalog(c(las_files_P1, las_files_P2))
ctg
plot(ctg)

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
