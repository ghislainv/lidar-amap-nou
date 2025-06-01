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
library(readr)
library(dplyr)

# The Agathis folder
# Agathis folder in project directory is a symlink to /media/ghislain/Agathis

# Output directory
out_dir_ssd <- here("Agathis/outputs/outputs-lidar-sites")
dir.create(out_dir_ssd, showWarnings=FALSE)

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
  ofile <- here(out_dir_ssd, glue("{site}.gpkg"))
  # Only execute if GPKG file does not exist
  if (!file.exists(ofile) & site != "dzumac-2") {
    # LAS files
    las_files <- lidar_df$path_las_files[i]
    if (grepl(";", las_files)) {
      las_files <- unlist(strsplit(las_files, split=";"))
    }
    # Get catalogue of file
    ctg <- lidR::readLAScatalog(here("Agathis", las_files))
    ## # Get boundaries using convex hull (not efficient)
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
    # Add columns to data frame
    area <- round(as.numeric(area) / 10000, 2)
    X <- cntrd_coord[1]
    Y <- cntrd_coord[2]
    geom <- sf_files_union
    df[i, 2:4] <- c(X, Y, area)
    df$geometry[i] <- geom
    ## sf_lidar_area <- sf::st_as_sf(
    ##   lidar_data,
    ##   geometry=sf_files_union)
    ## # Export as GPKG file
    ## write_sf(sf_lidar_area, ofile)
  }
}

# Create sf and export
sf_lidar_sites <- sf::st_as_sf(df) |>
  sf::st_set_crs("epsg:32758")

# Export as GPKG file
ofile <- here(out_dir_ssd, "lidar-sites.gpkg")
write_sf(sf_lidar_sites, ofile)

# End
