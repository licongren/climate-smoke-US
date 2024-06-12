library(raster)
library(gstat)
library(sf)
library(dplyr)   
library(terra)

rm(list=ls())
gc()

input_path <- "~/oak_space/mhqiu/smokePM_expsoure/input"
output_path <- "~/oak_space/mhqiu/smokePM_expsoure/output"
setwd(output_path)

#states_map <- map_data("state")
us_shp <- vect(paste0(input_path, "/grid_10km_wgs84/grid_10km_wgs84.shp"))
population = rast(paste0(input_path, "/gpw_v4_population_count_rev11_2020_2pt5_min.tif"))
# Crop population raster to CONUS
population = raster::crop(population, extent(-124.848974, -66.885444, 24.396308, 49.384358))

smoke_pm_raw <- readRDS(paste0(input_path, "/epa_station_day_totalPM_smokePM_panel_20000101-20230705.rds"))

date_range <- c(seq(as.Date("2006-01-01"), as.Date("2023-07-05"), by="days"))

############ interpolating smoke PM2.5
data_summ_combined <- NULL
for (pick_date in date_range){
  print(.Date(pick_date))
  smoke_pm <-   smoke_pm_raw %>%
    filter(date == pick_date, !is.na(smokePM))

  smoke_data_1 <- smoke_pm %>%
    rename(x=lon, y=lat, z=smokePM) %>%
    dplyr::select(x, y, z)

  idw <- interpIDW(population, as.matrix(smoke_data_1),
                          radius=2.5, power=5, smooth=0.5)
  names(idw) <- "interp_smokePM"

  ### mask to US range
  idw <- terra::mask(idw, us_shp)
  population <- terra::mask(population, us_shp)
  writeRaster(idw, paste0(output_path, "/interpolation/interpolated_smokePM_",
                          .Date(pick_date),".tif"), overwrite=T)

  df_smoke_pop <- data.frame(as.matrix(idw), as.matrix(population)) %>%
    rename(population=gpw_v4_population_count_rev11_2020_2pt5_min) %>%
    filter(!is.na(population), !is.na(interp_smokePM))

  data_summ <- df_smoke_pop %>%
    summarise(
      ppl_times_smoke=sum(population*interp_smokePM),
      pop_weight_smoke=weighted.mean(interp_smokePM, w=population),
      ppl_above_50=sum(population*as.numeric(interp_smokePM>50)),
      ppl_above_75=sum(population*as.numeric(interp_smokePM>75)),
      ppl_above_100=sum(population*as.numeric(interp_smokePM>100))) %>%
    mutate(date=.Date(pick_date))

  data_summ_combined <- bind_rows(data_summ_combined, data_summ)
}

data_summ_combined$date <- .Date(data_summ_combined$date)
saveRDS(data_summ_combined,  paste0(output_path, "/df_ppl_smoke_2006_20230705.rds"))

############ interpolating total PM2.5
data_summ_combined <- NULL
for (pick_date in date_range){
  print(.Date(pick_date))
  smoke_pm <-   smoke_pm_raw %>% 
    filter(date == pick_date, !is.na(pm25))
  
 pm25_data_1 <- smoke_pm %>%
    rename(x=lon, y=lat, z=pm25) %>%
    dplyr::select(x, y, z)
  
  idw <- interpIDW(population, as.matrix(pm25_data_1), 
                   radius=2.5, power=5, smooth=0.5) 
  names(idw) <- "interp_totalPM"
  
  ### mask to US range
  idw <- terra::mask(idw, us_shp)
  population <- terra::mask(population, us_shp)
  writeRaster(idw, paste0(output_path, "/interpolation/totalPM/interpolated_totalPM_",
                          .Date(pick_date),".tif"), overwrite=T)
  
  df_total_pop <- data.frame(as.matrix(idw), as.matrix(population)) %>%
    rename(population=gpw_v4_population_count_rev11_2020_2pt5_min) %>%
    filter(!is.na(population), !is.na(interp_totalPM))
  
  data_summ <- df_total_pop %>%
    summarise( 
      ppl_times_total=sum(population*interp_totalPM),
      pop_weight_total=weighted.mean(interp_totalPM, w=population),
      ppl_above_50=sum(population*as.numeric(interp_totalPM>50)),
      ppl_above_75=sum(population*as.numeric(interp_totalPM>75)),
      ppl_above_100=sum(population*as.numeric(interp_totalPM>100))) %>%
    mutate(date=.Date(pick_date))
  
  data_summ_combined <- bind_rows(data_summ_combined, data_summ)
}

data_summ_combined$date <- .Date(data_summ_combined$date)
saveRDS(data_summ_combined,  paste0(output_path, "/df_ppl_totalPM_2006_20230705.rds"))

# 
