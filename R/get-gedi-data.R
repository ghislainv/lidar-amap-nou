# Load libaries
library(terra)
library(rGEDI)
library(here)
library(glue)
library(dplyr)
library(sf)
library(ggplot2)

# Borders from GADM
ifile <- here("gisdata", "borders-newcal.gpkg")
borders <- sf::read_sf(ifile)

# Output directory
out_dir_download <- here("Agathis/outputs/outputs-gedi/downloads")
dir.create(out_dir_download, recursive=TRUE, showWarnings=FALSE)
out_dir_gedi <- here("Agathis/outputs/outputs-gedi")
out_dir_lidar <- here("Agathis/outputs/outputs-lidar-sites")
out_dir_chm <- here("Agathis/outputs/outputs-lidar-chm")

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

# Compute canopy height ALS derived for each gedi shot
ifile <- here(out_dir_gedi, "gedi-lidar-25m.gpkg")
gedi_lidar_25m_SpatVect <- terra::vect(ifile)
gedi_lidar_25m <- sf::st_read(ifile) |>
  mutate(als_height_mean=NA,
         als_height_max=NA,
         als_site=NA)
als_chm_files <- list.files(path=out_dir_chm, pattern="^chm.*\\.tif$")
n_als_sites <- length(als_chm_files)
for (i in 1:n_als_sites) {
  site_name <- substr(als_chm_files[i], 5, nchar(als_chm_files[i]) - 4)
  als_chm <- terra::rast(here(out_dir_chm, als_chm_files[i]))
  # Zonal statistics, mean
  # Use na.rm=FALSE to remove shots where
  # als include NAs (eg. borders of als lidar sites) 
  als_mean_canopy_height <- terra::zonal(
    x=als_chm,
    z=gedi_lidar_25m_SpatVect,
    fun="mean", na.rm=FALSE) |>
    mutate(height=ifelse(is.nan(height), NA, height))
  als_mean_canopy_height <- als_mean_canopy_height$height
  # Zonal statistics, max
  # Use na.rm=FALSE to remove shots where
  # als include NAs (eg. borders of als lidar sites) 
  als_max_canopy_height <- terra::zonal(
    x=als_chm,
    z=gedi_lidar_25m_SpatVect,
    fun="max", na.rm=FALSE) |>
    mutate(height=ifelse(is.nan(height), NA, height))
  als_max_canopy_height <- als_max_canopy_height$height
  # Add values to attribute table of gedi shots
  gedi_lidar_25m <- gedi_lidar_25m |>
    mutate(als_height_mean=ifelse(
      !is.na(als_mean_canopy_height),
      als_mean_canopy_height,
      als_height_mean)) |>
    mutate(als_height_max=ifelse(
      !is.na(als_max_canopy_height),
      als_max_canopy_height,
      als_height_max)) |>
    mutate(als_site=ifelse(
      !is.na(als_mean_canopy_height),
      site_name,
      als_site))
}

# Remove gedi shots with unrealistic values
gedi_lidar_25m_filter <- gedi_lidar_25m |>
  dplyr::filter(rh98 != 0) |>
  dplyr::filter(!is.na(als_site)) |>
  sf::st_write(here(out_dir_gedi, "gedi-lidar-filter.gpkg"),
               append=FALSE)

# Plot relationship between GEDI rh98 and ALS canopy height
# rh98-mean
p <- gedi_lidar_25m_filter |>
  ggplot(aes(als_height_mean, rh98, col=als_site)) +
  geom_point() +
  geom_smooth(formula=y~x, method="loess", col="red") +
  theme_bw()
ggsave(here(out_dir_gedi, "gedi-rh98-als-mean.png"))

# rh98-max
p <- gedi_lidar_25m_filter |>
  ggplot(aes(als_height_max, rh98, col=als_site)) +
  geom_point() +
  geom_smooth(formula=y~x, method="loess", col="red") +
  theme_bw()
ggsave(here(out_dir_gedi, "gedi-rh98-als-max.png"))

# mean-max
h_max <- max(gedi_lidar_25m_filter$als_height_mean)
p <- gedi_lidar_25m_filter |>
  ggplot(aes(als_height_mean, als_height_max, col=als_site)) +
  geom_point() +
  geom_smooth(formula=y~x, method="loess", col="red") +
  geom_segment(x=0, y=0, xend=h_max, yend=h_max,
               col="black") +
  theme_bw()
ggsave(here(out_dir_gedi, "als-max-mean.png"))

# End
