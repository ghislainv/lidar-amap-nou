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
start <- "2022-01-01"
end <- "2025-12-31"
daterange <- c("2022-01-01", "2025-12-31")

# Get path to GEDI data
# L2A, relative height metrics
gLevel2A <- gedifinder(
  product="GEDI02_A", ul_lat, ul_lon, lr_lat, lr_lon,
  version="002", daterange=daterange)

# Output directory
out_dir_download <- here("Agathis/outputs/outputs-gedi/downloads")
dir.create(out_dir_download, recursive=TRUE, showWarnings=FALSE)
out_dir_gedi <- here("Agathis/outputs/outputs-gedi")
out_dir_lidar <- here("Agathis/outputs/outputs-lidar-sites")

# =======================================
# Downloading GEDI data
# =======================================

## # Downloading GEDI data
## rGEDI::gediDownload(filepath=gLevel2A, outdir=out_dir_download)

# =======================================
# Downloading and reading GEDI data
# =======================================

# file_list <- list.files(out_dir_download, pattern="^.*\\.h5$")
nfiles <- length(gLevel2A)

for (i in 1:nfiles) {
  # Message
  file_name <- basename(gLevel2A[i])
  ifile <- here(out_dir_download, file_name)
  ofile <- here(out_dir_download, glue("level2a_{start}_{end}_{i}.gpkg"))
  if (!file.exists(ofile)) {
    cat(glue("Processing file {file_name}: {i}/{nfiles}\n\n"))
    # Downloading GEDI data
    rGEDI::gediDownload(filepath=gLevel2A[i], outdir=out_dir_download)
    # Read metrics
    gedilevel2a <- rGEDI::readLevel2A(level2Apath=ifile)
    level2a_metrics <- rGEDI::getLevel2AM(level2a=gedilevel2a)
    # Converting level2a data.table to small sf
    level2a_sf <- level2a_metrics |>
      dplyr::select(1:11, rh0, rh25, rh50,
                    rh75, rh90, rh95, rh98, rh100) |>
      dplyr::filter(
        !is.na(lon_lowestmode) & !is.na(lat_lowestmode)) |>
      sf::st_as_sf(
        coords=c("lon_lowestmode", "lat_lowestmode"),
        crs="epsg:4326") |>
      sf::st_intersection(borders) |>
      dplyr::select(-GID_0, -COUNTRY) |>
      sf::st_transform("epsg:32758") |>
      sf::st_write(ofile)
    # Delete file on disk
    file.remove(ifile)
  }
}

# =======================================
# Combine gpkg files in one big gpkg file
# =======================================

file_list <- list.files(out_dir_download, pattern="^.*\\.gpkg$")
nfiles <- length(file_list)
gedi <- sf::st_read(here(out_dir_download, file_list[1]))
for (i in 2:nfiles) {
  f_bind <- sf::st_read(here(out_dir_download, file_list[i]))
  gedi <- gedi |> 
    dplyr::bind_rows(f_bind)
}
sf::st_write(gedi, here(out_dir_gedi, "gedi-rh-nc.gpkg"))
remove(gedi)

# =======================================
# Select gedi shots based on date
# =======================================

library(lubridate)

# Load data
gedi <- sf::st_read(here(out_dir_gedi, "gedi-rh-nc.gpkg"))

# Select data for year 2022
gedi_2022 <- gedi |>
  mutate(date=date(ymd_hms("2018-01-01 00:00:00") + dseconds(delta_time))) |>
  mutate(year=as.numeric(year(date))) |>
  filter(year==2022 & quality_flag==0 & degrade_flag==0) |>
  sf::st_write(here(out_dir_gedi, "gedi-rh-nc-2022.gpkg"), append=FALSE)

# =============================================
# Keep gedi shots intersecting with LiDAR areas
# =============================================

# Load lidar areas
lidar <- sf::st_read(here(out_dir_lidar, "lidar-sites.gpkg"))
if (!("gedi" %in% ls())) gedi <- sf::st_read(here(out_dir_gedi, "gedi-rh-nc.gpkg"))

gedi_lidar <- gedi |>
  mutate(date=date(ymd_hms("2018-01-01 00:00:00") + dseconds(delta_time))) |>
  mutate(year=as.numeric(year(date))) |>
  filter(quality_flag==0 & degrade_flag==0) |>
  sf::st_intersection(lidar) |>
  sf::st_write(here(out_dir_gedi, "gedi-lidar.gpkg"), append=FALSE)

# Buffer of 12.5m around center of gedi shot
gedi_lidar_25m <- gedi_lidar |>
  sf::st_buffer(12.5) |>
  sf::st_write(here(out_dir_gedi, "gedi-lidar-25m.gpkg"), append=FALSE)

# End
