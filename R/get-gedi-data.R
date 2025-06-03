# Load libaries
library(rGEDI)
library(here)
library(glue)
library(dplyr)

# Borders from GADM
ifile <- here("gisdata", "borders-newcal.gpkg")
borders <- sf::read_sf(ifile)

# =======================================
# Find GEDI data within your study area (GEDI finder tool)
# =======================================

# Study area boundary box coordinates
ul_lat <- -19.5249919999999015
ul_lon <- 163.5694430000000068
lr_lat <- -22.8480570010000008
lr_lon <- 168.1347199999999873

# Specifying the date range
daterange <- c("2024-01-01","2025-12-31")

# Get path to GEDI data
# L2A, relative height metrics
gLevel2A <- gedifinder(
  product="GEDI02_A", ul_lat, ul_lon, lr_lat, lr_lon,
  version="002", daterange=daterange)

# Output directory
out_dir_ssd <- here("Agathis/outputs/outputs-gedi/downloads")
dir.create(out_dir_ssd, recursive=TRUE, showWarnings=FALSE)

# =======================================
# Downloading GEDI data
# =======================================

# Downloading GEDI data
rGEDI::gediDownload(filepath=gLevel2A,outdir=out_dir_ssd)

# =======================================
# Reading GEDI data
# =======================================

file_list <- list.files(out_dir_ssd, pattern="^.*\\.h5$")
nfiles <- length(file_list)

for (i in 1:nfiles) {
  # Message
  cat(glue("Processing file {file_list[i]}: {i}/{nfiles}\n\n"))
  ifile <- here(out_dir_ssd, file_list[i])
  gedilevel2a <- rGEDI::readLevel2A(level2Apath=ifile)
  level2a_metrics <- rGEDI::getLevel2AM(level2a=gedilevel2a)
  # Converting level2a data.table to sf
  level2a_sf <- level2a_metrics |>
    dplyr::select(1:11, rh0, rh95, rh98, rh100) |>
    dplyr::filter(!is.na(lon_lowestmode) & !is.na(lat_lowestmode)) |>
    sf::st_as_sf(
      coords=c("lon_lowestmode", "lat_lowestmode"),
      crs="epsg:4326") |>
    sf::st_intersection(borders) |>
    sf::st_transform("epsg:32758") |>
    sf::st_write(here(out_dir_ssd, "level2a.gpkg"), append=TRUE)
}

# End
