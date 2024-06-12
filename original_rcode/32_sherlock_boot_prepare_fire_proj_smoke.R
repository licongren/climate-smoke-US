library(terra);library(ncdf4)
library(sf)
library(dplyr)   
library(fixest)
library(tidyverse)
library(foreach)
library(doParallel)
library(doSNOW)

rm(list=ls())
gc()

input_path <- "~/oak_space/mhqiu/smoke_climate/input_data"
output_path <- "~/oak_space/mhqiu/smoke_climate/output"

#-------------------------------------------------------------------------------
# Summarise the projected emissions to the smoke cell level and by distance and wind
# Written by Minghao
# Last edited Nov 2023
#-------------------------------------------------------------------------------
smoke_pm <- readRDS(paste0(input_path, "/smokePM_predictions_20060101_20201231.rds"))

grid_10km_list <- sort(unique(smoke_pm$grid_id_10km))#[1:25000] #50001:100156

era5_month <-  readRDS(paste0(input_path, "/era5_10km_grid_monthly_2006_2020.rds")) %>% 
  filter(grid_id_10km %in% grid_10km_list) %>%
  #select(grid_id_10km, u10m, v10m, year, month) %>% 
  group_by(grid_id_10km, month) %>%
  summarise(u10m=mean(u10m, na.rm=T), 
            v10m=mean(v10m, na.rm=T)) %>% ungroup()

grid_10km <- st_read(paste0(input_path, "/10km_grid_shp/grid_10km_wgs84.shp")) %>%
  st_centroid()
coord <- st_coordinates(grid_10km)
grid_10km <- grid_10km %>% mutate(lon=coord[,1], lat=coord[,2])
grid_10km_df <- as.data.frame(grid_10km) %>% mutate(grid_id_10km=ID) %>% select(-ID, -COORDX, -COORDY, -geometry)

getDistanceFromLatLonInKm  <- function(lat1=lat1,lon1=lon1,lat2=lat2,lon2=lon2) {
  R <-  6371; #Radius of the earth in km
  dLat <- deg2rad(lat2-lat1) # deg2rad below
  dLon <- deg2rad(lon2-lon1);
  a <- sin(dLat/2) * sin(dLat/2) + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(dLon/2) * sin(dLon/2)
  c <- 2 * atan2(sqrt(a), sqrt(1-a));
  d <- R * c; # Distance in km
  return(d)
}

deg2rad <- function(deg=deg) {
  return (deg*pi/180)
}

getangle_latlon <- function(lat1=lat1,lon1=lon1,lat2=lat2,lon2=lon2) {
  dLat <- deg2rad(lat2-lat1) # deg2rad below
  dLon <- deg2rad(lon2-lon1)
  angle <- atan2(sin(dLon)*cos(deg2rad(lat2)), (cos(deg2rad(lat1))*sin(deg2rad(lat2)) - sin(deg2rad(lat1))*cos(deg2rad(lat2))*cos(dLon)))
  angle1 <- (angle/pi*180 + 360) %% 360
  return(angle1)
}

########  load the projected fire emissions at the gfed cell level #############
gfed_data <- readRDS(paste0(input_path,"/climate-fire_training_data_2001-2021_clean_varname.rds"))
gfed_coord  <- gfed_data %>%
  distinct(gfed_cell_id, lat, lon) %>%
  mutate(lat_gfed=lat, lon_gfed=lon) %>% select(-lat, -lon)

scen_name <- "ssp370"
#fire_proj <-  readRDS(paste0(input_path,"/downscaled_emis_final_threshold_05.rds"))  %>% filter(scenario == scen_name, year==2055) 
#fire_proj <-  readRDS(paste0(input_path,"/downscaled_emis_2046_2055_years.rds"))  %>% filter(scenario == scen_name) 
##downscaled_emis_2046_2055_years.rds
#fire_proj <-  readRDS(paste0(input_path,"/downscaled_emis_uncertainty_boot.rds")) %>% filter(scenario == scen_name) 
fire_proj <-  readRDS(paste0(input_path,"/downscaled_emis_MC_100_fire_gcm.rds")) %>% 
  filter(scenario == scen_name, mc_id %in% 1:50) 

#### downscale the fire proj to the monthly-level
gfed_month <- gfed_data %>% group_by(gfed_cell_id, month) %>%
  summarise(DM_kg=sum(DM_kg, na.rm=T)) %>%
  ungroup() %>%
  group_by(gfed_cell_id) %>%
  mutate(DM_month_ratio=DM_kg/sum(DM_kg)) %>%
  ungroup() %>% select(-DM_kg)

fire_proj_month <- fire_proj %>% slice(rep(1:n(), each = 12)) %>%
  mutate(month=rep(1:12, nrow(fire_proj))) %>%
  left_join(gfed_month) %>%
  mutate(pred_DM_grid= pred_DM_grid * DM_month_ratio)

##### the function to calculate emission for each smoke grid cell by wind and distance
get_gfed_emis_monthwind <- function(grid_id=grid_id, fire_proj_data=fire_proj_data){
  tmp_10km <-  filter(grid_10km_df, grid_id_10km==grid_id)
  distance_tmp <- gfed_coord %>% mutate(lat_smoke=tmp_10km$lat, lon_smoke=tmp_10km$lon)
  distance_tmp <- distance_tmp %>% mutate(distance=getDistanceFromLatLonInKm(lat1=lat_gfed, lon1=lon_gfed, lat2=lat_smoke, lon2=lon_smoke),
                                          angle=getangle_latlon(lat1=lat_gfed, lon1=lon_gfed, lat2=lat_smoke, lon2=lon_smoke))
  
  fire_proj_tmp <- fire_proj_data %>%
    left_join(distance_tmp[,c("gfed_cell_id","distance","angle")], by="gfed_cell_id") %>%
    select(year, month, mc_id, scenario, gfed_cell_id, pred_DM_grid, distance, angle)
  
  era5_tmp <- era5_month %>%  filter(grid_id_10km==grid_id) %>%
    mutate(wind_direction=(atan2(u10m,v10m)*180/pi+360)%%360,
           wind_direct1=(wind_direction-45)%%360,
           wind_direct2=(wind_direction+45)%%360) %>%
    mutate(wind_direct3=(wind_direct1+180)%%360,
           wind_direct4=(wind_direct2+180)%%360)
  
  fire_era5_tmp <- left_join(fire_proj_tmp, era5_tmp, by="month") %>%
    mutate(upwind=2) %>%     ##### 0: downwind-90, 1:upwind-90; 2: the other 180
    mutate(upwind=replace(upwind, (wind_direct1<wind_direct2 &  angle>wind_direct1 & angle<wind_direct2), 1)) %>%
    mutate(upwind=replace(upwind, (wind_direct1>wind_direct2 & (angle>wind_direct1 | angle<wind_direct2)), 1)) %>%
    mutate(upwind=replace(upwind, (wind_direct3<wind_direct4 &  angle>wind_direct3 & angle<wind_direct4), 0)) %>%
    mutate(upwind=replace(upwind, (wind_direct3>wind_direct4 & (angle>wind_direct3 | angle<wind_direct4)), 0))
  
  fire_era5_tmp <- fire_era5_tmp %>% pivot_longer(cols=pred_DM_grid) %>% 
    group_by(grid_id_10km, year, month, mc_id, scenario, name, upwind) %>%
    summarise(within50km=sum(value*as.numeric(distance<50), na.rm=T),
              dist50_100km=sum(value*as.numeric(distance>50 & distance<100), na.rm=T),
              dist100_200km=sum(value*as.numeric(distance>100 & distance<200), na.rm=T),
              dist200_350km=sum(value*as.numeric(distance>200 & distance<350), na.rm=T),
              dist350_500km=sum(value*as.numeric(distance>350 & distance<500), na.rm=T),
              dist500_750km=sum(value*as.numeric(distance>500 & distance<750), na.rm=T),
              dist750_1000km=sum(value*as.numeric(distance>750 & distance<1000), na.rm=T),
              dist1000_1500km=sum(value*as.numeric(distance>1000 & distance<1500), na.rm=T),
              dist1500_2000km=sum(value*as.numeric(distance>1500 & distance<2000), na.rm=T),
              above_2000km=sum(value*as.numeric(distance>2000), na.rm=T)) %>% ungroup() %>%
    select(-name)
  
  return(list(fire_era5_tmp))
}

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}


  fire_proj_combine  <- NULL
  cl <- makeCluster(15)  #### change the number of clusters based on the configuration of the computing environment
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(min = 1, max =  length(grid_10km_list), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  fire_proj_combine  <- foreach(ii=1:length(grid_10km_list),.options.snow = opts,.packages=c("dplyr","tidyr"),
                                .combine='comb', .multicombine=TRUE,
                                .init=list(list())) %dopar% {
                                  get_gfed_emis_monthwind(grid_id=grid_10km_list[ii], fire_proj_data=fire_proj_month)
                                }
  stopCluster(cl)
  
  fire_combine <- dplyr::bind_rows(fire_proj_combine[[1]])
  #saveRDS(fire_combine, paste0(output_path,"/downscaled_emission/emis_proj_2050_",scen_name,"_boot.rds"))
  saveRDS(fire_combine, paste0(output_path,"/downscaled_emission/emis_proj_2050_",scen_name,"_MC_1_50.rds"))
  


