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
out_dir_lidar <- here("Agathis/outputs/outputs-lidar-sites")
dir.create(out_dir_lidar, showWarnings=FALSE)

# Import lidar database
ifile <- here("Agathis", "lidar-flights-newcal.csv")
lidar_df <- readr::read_csv(ifile, show_col_types=FALSE)
nsites <- nrow(lidar_df)

# Data frame to hold results
df <- data.frame(
  name=lidar_df$site, X=NA, Y=NA,
  area=NA, geometry=NA)

# Get boundaries from las files and export as gpkg
for (i in 1:nsites) {
  site <- lidar_df$site[i]
  ofile <- here(out_dir_lidar, glue("{site}.gpkg"))
  # Only execute if GPKG file does not exist
  if (!file.exists(ofile)) {
    # LAS files
    las_files <- lidar_df$path_las_files[i]
    if (grepl(";", las_files)) {
      las_files <- unlist(strsplit(las_files, split=";"))
    }
    # Get catalogue of files
    ctg <- lidR::readLAScatalog(here("Agathis", las_files))
    ## # Get boundaries using convex hull (not efficient, use lasboundary)
    ## ctg_boundaries <- lidR::catalog_boundaries(
    ##   ctg,
    ##   concavity=1,
    ##   length_threshold=15)
    # Get crs of las files
    crs <- ctg$CRS[1]
    if (site == "ouemo") {crs <- "epsg:3163"}
    if (site == "grand-lac") {crs <- "epsg:32758"}
    # Number of files
    nfiles <- length(ctg$filename)
    # Convert catalogue to sf
    sf_files <- sf::st_sf(file=paste0("file-", 1:nfiles),
                          geometry=ctg$geometry)
    sf_files <- sf::st_set_crs(sf_files, crs)
    # Union of features
    sf_files_union <- sf::st_union(sf_files)
    if (sf::st_crs(sf_files_union) != sf::st_crs("epsg:32758")) {
      sf_files_union <- sf::st_transform(
        sf_files_union, "epsg:32758")
    }
    # Compute centroid coordinates and area (in ha)
    cntrd <- sf::st_centroid(sf_files_union)
    cntrd_coord <- sf::st_coordinates(cntrd)
    area <- sf::st_area(sf_files_union)
    # Filling in the data frame
    area <- round(as.numeric(area) / 10000, 2)
    X <- cntrd_coord[1]
    Y <- cntrd_coord[2]
    geom <- sf_files_union
    df[i, 2:4] <- c(X, Y, area)
    df$geometry[i] <- geom
  }
}

# Create sf and export as gpkg
ofile <- here(out_dir_lidar, "lidar-sites.gpkg")
sf_lidar_sites <- sf::st_as_sf(df) |>
  sf::st_set_crs("epsg:32758") |>
  write_sf(ofile)

# Create point data and export as gpkg
ofile <- here(out_dir_lidar, "lidar-sites-pts.gpkg")
df_pts <- df |>
  select(-geometry) |>
  sf::st_as_sf(
    coords=c("X", "Y"),
    remove=FALSE,
    crs=sf::st_crs("epsg:32758")) |>
  write_sf(ofile)

# =====================================
# Map of lidar sites
# =====================================

# Lidar sites
ifile <- here(out_dir_lidar, "lidar-sites-pts.gpkg")
lidar_sites <- sf::read_sf(ifile)

# Borders from GADM
ifile <- here("gisdata", "borders-newcal.gpkg")
borders <- sf::read_sf(ifile) |>
  sf::st_transform("epsg:32758")

# Forest from TMFv2023
ifile <- here("gisdata", "forest-tmf-2000-2010-2020.tif")
ofile <- here("gisdata", "forest-newcal-tmf-2020.tif")
gdal_utils(
  util="warp",
  source=ifile,
  destination=ofile,
  options=c(
    "-b", "3",
    "-tr", "500", "500",
    "-t_srs", "EPSG:32758",
    "-r", "mode",
    "-srcnodata", "None",
    "-dstnodata", "0",
    "-ot", "Byte",
    "-tap", "-overwrite",
    "-co", "COMPRESS=DEFLATE"))
forest <- terra::rast(ofile)

# Options
tmap_opt <- function(npix=1e5, ...) {
  tmap_options_reset()
  tmap_options(raster.max_cells = c(plot=npix, view=npix), ...)
  cat("raster.max_cells set to: ", npix, "\n")
}
tmap_opt(npix=1e7)
green <- rgb(34, 139, 34, 255, maxColorValue=255)
ofile <- here(out_dir_lidar, "lidar-sites.png")
textwidth <- 16.6  # textwidth (in cm) for figure width

# Map
tm_lidar <- 
  tm_shape(forest) +
  tm_raster(
    col.scale=tm_scale_categorical(values=c(green)),
    col_alpha=0.5,
    col.legend=tm_legend(
      # frame=FALSE, # not working
      frame.lwd=0,
      text.color="transparent",
      title="Forest cover in 2020\n(source TMFv2023)", 
      #orientation="landscape",
      position=tm_pos_in(0.02, 0.24, just.h="left", just.v="bottom"))
  ) +
  tm_shape(borders, unit="km") +
  tm_borders(col="black", lwd=0.75) +
  tm_scalebar(
    breaks=c(0, 50, 100), text.size=0.7,
    position=tm_pos_in(0.5, 1, just.h="center", just.v="top")
  ) +
  tm_shape(lidar_sites) +
  tm_squares(
    col="blue", fill="white", fill_alpha=0, size="area",
    size.legend=tm_legend(
      frame=FALSE,
      title="Area (ha)", 
      orientation="landscape",
      position=tm_pos_in(0.02, 0.08, just.h="left", just.v="bottom"))
  ) +
  tm_layout(inner.margins=rep(0.02, 4),
            outer.margins=rep(0.02, 4)) +
  tm_title(text="AMAP lidar sites in New Caledonia",
           position=tm_pos_in(0.02, 0.02, just.h="left", just.v="bottom"),
           size=1.2)

# Save map
tmap_save(tm_lidar, filename=ofile, width=textwidth, units="cm", dpi=300, scale=1)

# End
