library(terra);library(ncdf4)
library(openxlsx)
library(ggplot2)
library(sf)
library(dplyr)   
library(fixest)
library(tidyverse)
library(MetBrewer)
library(dplyr)
library(purrr)
library(data.table)
library(mlr3verse)
library(xgboost)
library(future)
library(broom)
library(zoo)

rm(list=ls())
gc()

code_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/rcode"
data_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data"
figure_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/figures"
result_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/results/fire_model"

setwd(figure_path)
states_map <- map_data("state")

db_path <- "/Users/mhqiu/BurkeLab Dropbox/Projects/smoke-climate"

#-------------------------------------------------------------------------------
# Use the different models to predict levels of emissions
# Written by Minghao
# Last edited July 2022
#-------------------------------------------------------------------------------
gfed_grid <- readRDS(paste0(data_path, "/GFED4S_NorthAmerica_025deg_grid_EPSG4326.rds"))

gfed_ecoregion_link <- readRDS(paste0(db_path, "/crosswalk_GFED_NA_CEC_Eco_Level2_3.rds")) %>%
  mutate(intersection_area=as.numeric(intersection_area)) %>%   group_by(gfed_cell_id)  %>%
  arrange(desc(intersection_area)) %>%
  mutate(area_weight=intersection_area/sum(intersection_area)) %>% ungroup() # %>% distinct(gfed_cell_id, .keep_all=T)

gfed_wang_region <- readRDS(paste0(data_path, "/crosswalk_GFED_wang_regions.rds")) 

gfed_data <- readRDS(file.path(db_path, "climate-fire_training_data_2001-2021_clean_varname.rds"))
gfed_data <- left_join(gfed_data, gfed_wang_region, by="gfed_cell_id")
gfed_data <- left_join(gfed_data, gfed_ecoregion_link, by="gfed_cell_id")
gfed_data$area <- as.numeric(gfed_data$area) 

vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
                   "narr_vpd", "narr_runoff",
                   "nldas_soilm")

vlist_crossec <- c("barren", "cropland", "forest", "grassland", 
                   "lichen", "shrubland", "snowice", "urban", "water","wetland",
                   "elevation", "slope")  

### create crosswalk between ecoregions and 1deg cells of CMIP6 
crosswalk_gfed_1deg <- readRDS(paste0(data_path, "/crosswalk_gfed_1deg.rds"))
crosswalk_ecoregion <- left_join(crosswalk_gfed_1deg, gfed_ecoregion_link) %>%
  left_join(gfed_wang_region) %>%
  filter(gfed_cell_id %in% unique(gfed_data$gfed_cell_id)) ### drop all grid cells that do not have wildfire for CMIP6 projections

gfed_veg <- gfed_data %>% distinct(gfed_cell_id, .keep_all = T) %>%
  mutate(veg_frac = forest * as.numeric(region %in% c("Western US", "Canada-Alaska")) +
        (forest + grassland + cropland) * as.numeric(!region %in% c("Western US", "Canada-Alaska"))) %>%
  select(gfed_cell_id, veg_frac)

crosswalk_ecoregion <- left_join(crosswalk_ecoregion, gfed_veg)

#### generate the cross-sectional data for the cross-section features
gfed_2021 <- filter(gfed_data, year==2021)

crossec_region <- gfed_2021 %>% group_by(region) %>%
  summarise(across(vlist_crossec,  weighted.mean, w=intersection_area)) %>%  ungroup()

crossec_eco2 <- gfed_2021 %>% group_by(region, NA_L2KEY) %>%
  summarise(across(vlist_crossec, weighted.mean, w=intersection_area)) %>%  ungroup() 

crossec_eco3 <- gfed_2021 %>% group_by(region, NA_L3KEY) %>%
  summarise(across(vlist_crossec, weighted.mean, w=intersection_area)) %>%  ungroup() 

crossec_1deg <- gfed_2021 %>% 
  left_join(crosswalk_gfed_1deg) %>%
  group_by(region, grid_1deg_id) %>%
  summarise(across(vlist_crossec, weighted.mean, w=intersection_area)) %>%  ungroup() 

#---------------------
## Prepare the training data from CMIP to Eco-3 level
#---------------------
var_long_list <- c("temperature_surface", "precipitation", "RH_surface", "vpd",
                   "soil_moisture",  "soil_moisture", "runoff", "wind_speed")
var_labels <- c("narr_temp", "narr_precip", "narr_rhum", "narr_vpd", 
                "narr_soilm", "nldas_soilm", "narr_runoff", "narr_wspd")

cmip_eco3 <- NULL
for (vv in 1: length(var_long_list)){
  var_long <- var_long_list[vv]  
  var_name <- var_labels[vv]  
  print(var_name)
  
  cmip_debias <- readRDS(paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/",
                                var_long, "/CMIP6_annual_debias_1deg.rds")) 
  
  if (var_name == "nldas_soilm"){
    cmip_debias <- readRDS(paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/",
                                  var_long, "/CMIP6_annual_debias_1deg_nldas.rds")) 
  }
  
  model_list <- unique(cmip_debias$model)
  
  #   print(length(model_list))
  # }
  
  cmip_data_combined <- NULL
  for (mm in model_list){
    print(mm)
    cmip_data <-  filter(cmip_debias, model == mm) %>% 
      left_join(crosswalk_ecoregion) %>%
      filter(!is.na(NA_L3KEY)) %>%
      group_by(model, scenario, year, region, NA_L3KEY) %>%
      summarise(cmip_debias = weighted.mean(cmip_debias, w= intersection_area, na.rm=T)) %>%
      rename(!!as.name(var_name):= cmip_debias)
    
    cmip_data_combined <- bind_rows(cmip_data_combined, cmip_data)
  }
  
  if (vv ==1){
    cmip_eco3 <- cmip_data_combined
  } else{
    cmip_eco3 <- left_join(cmip_eco3, cmip_data_combined)
  }
}

cmip_eco3_1 <- left_join(cmip_eco3, crossec_eco3)
saveRDS(cmip_eco3_1, paste0(result_path,"/CMIP6_fire_proj/final_cmip_eco3_data.rds"))

##### generate extra features for eco3 data and de-mean based on mean of subregions ##### 
cmip_eco3 <- readRDS(paste0(result_path,"/CMIP6_fire_proj/final_cmip_eco3_data.rds")) %>%
  ungroup() %>% rename(region_sub = NA_L3KEY)

vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
                   "narr_vpd", "narr_runoff","nldas_soilm")

cmip_eco3_full <- cmip_eco3 %>%
  mutate_at(vlist_climate, list(forest = function(x){x*cmip_eco3$forest},
                                crop = function(x){x*cmip_eco3$cropland},
                                grass = function(x){x*cmip_eco3$grassland})) %>%
  mutate_at(vlist_climate, list(temp = function(x){x*cmip_eco3$narr_temp})) %>%
  mutate_at(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff", "nldas_soilm"), 
            list(precip = function(x){x*cmip_eco3$narr_precip})) %>%
  mutate_at(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"),
            list(rhum = function(x){x*cmip_eco3$narr_rhum})) %>%  
  mutate_at(c("narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"),
            list(wspd = function(x){x*cmip_eco3$narr_wspd})) %>% 
  mutate_at(c("narr_vpd", "narr_runoff", "nldas_soilm"),
            list(vpd = function(x){x*cmip_eco3$narr_vpd})) %>% 
  mutate_at(c("narr_runoff", "nldas_soilm"),
            list(runoff = function(x){x*cmip_eco3$narr_runoff})) %>% 
  mutate(nldas_soilm_soilm = nldas_soilm^2) %>%                              
  select(-barren, -cropland, -forest, -grassland, -shrubland,
         -snowice, -urban, -water, -wetland, -lichen, -slope, -elevation)

var_list <- c(vlist_climate, paste(rep(vlist_climate, each = 3),
                                   c("forest","crop", "grass"), sep="_"),
              paste0(vlist_climate, "_temp"),
              paste0(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff", "nldas_soilm"), "_precip"),
              paste0(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_rhum"),
              paste0(c("narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_wspd"),
              paste0(c("narr_vpd", "narr_runoff", "nldas_soilm"), "_vpd"),
              paste0(c("narr_runoff", "nldas_soilm"), "_runoff"),
              paste0(c("nldas_soilm"), "_soilm"))

##### load the mean of log-based dataset to de-mean
region_mean_log <- readRDS(paste0(result_path, "/final_models/training_data/obs_subregion_mean_log.rds"))$data_eco3

cmip_eco3_log <- left_join(cmip_eco3_full, region_mean_log[,c("region", "region_sub", "area", "outcome")]) %>%
      rename(outcome_mean_log=outcome) %>%
  filter(!is.na(outcome_mean_log))

for (vv in var_list){
region_mean_tmp <- region_mean_log[,c("region", "region_sub",vv)] %>%
  rename(regionsub_mean = !!as.name(vv))

cmip_eco3_log <- left_join(cmip_eco3_log, region_mean_tmp) %>%
  mutate(!!as.name(vv):= !!as.name(vv) - regionsub_mean) %>% select(-regionsub_mean)
}
saveRDS(cmip_eco3_log, paste0(result_path,"/CMIP6_fire_proj/cmip_eco3_full_demean_log.rds"))

##### load the mean of level-based dataset to de-mean
region_mean_level <- readRDS(paste0(result_path, "/final_models/training_data/obs_subregion_mean_level.rds"))$data_eco3

cmip_eco3_level <- left_join(cmip_eco3_full, region_mean_level[,c("region", "region_sub", "area", "outcome")]) %>%
  rename(outcome_mean_level=outcome) %>%
  filter(!is.na(outcome_mean_level))

for (vv in var_list){
  region_mean_tmp <- region_mean_level[,c("region", "region_sub",vv)] %>%
    rename(regionsub_mean = !!as.name(vv))
  
  cmip_eco3_level <- left_join(cmip_eco3_level, region_mean_tmp) %>%
    mutate(!!as.name(vv):= !!as.name(vv) - regionsub_mean) %>% select(-regionsub_mean)
}
saveRDS(cmip_eco3_level, paste0(result_path,"/CMIP6_fire_proj/cmip_eco3_full_demean_level.rds"))

#---------------------
## Prepare the training data from CMIP to Eco-2 level
#---------------------
var_long_list <- c("temperature_surface", "precipitation", "RH_surface", "vpd",
                   "soil_moisture",  "soil_moisture", "runoff", "wind_speed")
var_labels <- c("narr_temp", "narr_precip", "narr_rhum", "narr_vpd", 
                "narr_soilm", "nldas_soilm", "narr_runoff", "narr_wspd")

cmip_eco2 <- NULL
for (vv in 1: length(var_long_list)){
  var_long <- var_long_list[vv]  
  var_name <- var_labels[vv]  
  print(var_name)
  
  cmip_debias <- readRDS(paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/",
                                var_long, "/CMIP6_annual_debias_1deg.rds")) 
  
  if (var_name == "nldas_soilm"){
    cmip_debias <- readRDS(paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/",
                                  var_long, "/CMIP6_annual_debias_1deg_nldas.rds")) 
  }
  
  model_list <- unique(cmip_debias$model)
  
  #   print(length(model_list))
  # }
  
  cmip_data_combined <- NULL
  for (mm in model_list){
    print(mm)
    cmip_data <-  filter(cmip_debias, model == mm) %>% 
      left_join(crosswalk_ecoregion) %>%
      filter(!is.na(NA_L2KEY)) %>%
      group_by(model, scenario, year, region, NA_L2KEY) %>%
      summarise(cmip_debias = weighted.mean(cmip_debias, w= intersection_area, na.rm=T)) %>%
      rename(!!as.name(var_name):= cmip_debias)
    
    cmip_data_combined <- bind_rows(cmip_data_combined, cmip_data)
  }
  
  if (vv ==1){
    cmip_eco2 <- cmip_data_combined
  } else{
    cmip_eco2 <- left_join(cmip_eco2, cmip_data_combined)
  }
}

cmip_eco2_1 <- left_join(cmip_eco2, crossec_eco2)
saveRDS(cmip_eco2_1, paste0(result_path,"/CMIP6_fire_proj/final_cmip_eco2_data.rds"))

##### generate extra features for eco2 data and de-mean based on mean of subregions ##### 
cmip_eco2 <- readRDS(paste0(result_path,"/CMIP6_fire_proj/final_cmip_eco2_data.rds")) %>%
  ungroup() %>% rename(region_sub = NA_L2KEY)

vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
                   "narr_vpd", "narr_runoff","nldas_soilm")

cmip_eco2_full <- cmip_eco2 %>%
  mutate_at(vlist_climate, list(forest = function(x){x*cmip_eco2$forest},
                                crop = function(x){x*cmip_eco2$cropland},
                                grass = function(x){x*cmip_eco2$grassland})) %>%
  mutate_at(vlist_climate, list(temp = function(x){x*cmip_eco2$narr_temp})) %>%
  mutate_at(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff", "nldas_soilm"), 
            list(precip = function(x){x*cmip_eco2$narr_precip})) %>%
  mutate_at(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"),
            list(rhum = function(x){x*cmip_eco2$narr_rhum})) %>%  
  mutate_at(c("narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"),
            list(wspd = function(x){x*cmip_eco2$narr_wspd})) %>% 
  mutate_at(c("narr_vpd", "narr_runoff", "nldas_soilm"),
            list(vpd = function(x){x*cmip_eco2$narr_vpd})) %>% 
  mutate_at(c("narr_runoff", "nldas_soilm"),
            list(runoff = function(x){x*cmip_eco2$narr_runoff})) %>% 
  mutate(nldas_soilm_soilm = nldas_soilm^2) %>%                              
  select(-barren, -cropland, -forest, -grassland, -shrubland,
         -snowice, -urban, -water, -wetland, -lichen, -slope, -elevation)

var_list <- c(vlist_climate, paste(rep(vlist_climate, each = 3),
                                   c("forest","crop", "grass"), sep="_"),
              paste0(vlist_climate, "_temp"),
              paste0(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff", "nldas_soilm"), "_precip"),
              paste0(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_rhum"),
              paste0(c("narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_wspd"),
              paste0(c("narr_vpd", "narr_runoff", "nldas_soilm"), "_vpd"),
              paste0(c("narr_runoff", "nldas_soilm"), "_runoff"),
              paste0(c("nldas_soilm"), "_soilm"))

##### load the mean of log-based dataset to de-mean
region_mean_log <- readRDS(paste0(result_path, "/final_models/training_data/obs_subregion_mean_log.rds"))$data_eco2

cmip_eco2_log <- left_join(cmip_eco2_full, region_mean_log[,c("region", "region_sub", "area", "outcome")]) %>%
  rename(outcome_mean_log=outcome) %>%
  filter(!is.na(outcome_mean_log))

for (vv in var_list){
  region_mean_tmp <- region_mean_log[,c("region", "region_sub",vv)] %>%
    rename(regionsub_mean = !!as.name(vv))
  
  cmip_eco2_log <- left_join(cmip_eco2_log, region_mean_tmp) %>%
    mutate(!!as.name(vv):= !!as.name(vv) - regionsub_mean) %>% select(-regionsub_mean)
}
saveRDS(cmip_eco2_log, paste0(result_path,"/CMIP6_fire_proj/cmip_eco2_full_demean_log.rds"))

##### load the mean of level-based dataset to de-mean
region_mean_level <- readRDS(paste0(result_path, "/final_models/training_data/obs_subregion_mean_level.rds"))$data_eco2

cmip_eco2_level <- left_join(cmip_eco2_full, region_mean_level[,c("region", "region_sub", "area", "outcome")]) %>%
  rename(outcome_mean_level=outcome) %>%
  filter(!is.na(outcome_mean_level))

for (vv in var_list){
  region_mean_tmp <- region_mean_level[,c("region", "region_sub",vv)] %>%
    rename(regionsub_mean = !!as.name(vv))
  
  cmip_eco2_level <- left_join(cmip_eco2_level, region_mean_tmp) %>%
    mutate(!!as.name(vv):= !!as.name(vv) - regionsub_mean) %>% select(-regionsub_mean)
}
saveRDS(cmip_eco2_level, paste0(result_path,"/CMIP6_fire_proj/cmip_eco2_full_demean_level.rds"))


#---------------------
## Prepare the training data from CMIP to regional level
#---------------------
cmip_region <- NULL
for (vv in 1: length(var_long_list)){
  var_long <- var_long_list[vv]  
  var_name <- var_labels[vv]  
  print(var_name)
  
  cmip_debias <- readRDS(paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/",
                                var_long, "/CMIP6_annual_debias_1deg.rds")) 
  
  if (var_name == "nldas_soilm"){
    cmip_debias <- readRDS(paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/",
                                  var_long, "/CMIP6_annual_debias_1deg_nldas.rds")) 
  }
  
  model_list <- unique(cmip_debias$model)
  
  cmip_data_combined <- NULL
  for (mm in model_list){
    print(mm)
    cmip_data <-  filter(cmip_debias, model == mm) %>% 
      left_join(crosswalk_ecoregion) %>%
      filter(!is.na(region)) %>%
      group_by(model, scenario, year, region) %>%
      summarise(cmip_debias = weighted.mean(cmip_debias, w= intersection_area*veg_frac, na.rm=T)) %>%
      rename(!!as.name(var_name):= cmip_debias)
    
    cmip_data_combined <- bind_rows(cmip_data_combined, cmip_data)
  }
  
  if (vv == 1){
    cmip_region <- cmip_data_combined
  } else{
    cmip_region <- left_join(cmip_region, cmip_data_combined)
  }
}
saveRDS(cmip_region, paste0(result_path,"/CMIP6_fire_proj/final_cmip_region_data.rds"))

cmip_eco2_1 <- left_join(cmip_eco2, crossec_eco2)
saveRDS(cmip_eco2_1, paste0(result_path,"/CMIP6_fire_proj/final_cmip_eco2_data.rds"))

##### generate extra features for regional data and de-mean based on mean of subregions ##### 
cmip_region <- readRDS(paste0(result_path,"/CMIP6_fire_proj/final_cmip_region_data_veg.rds")) %>% ungroup() 

vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
                   "narr_vpd", "narr_runoff","nldas_soilm")

cmip_region_full <- cmip_region %>%
  mutate_at(vlist_climate, list(temp = function(x){x*cmip_region$narr_temp})) %>%
  mutate_at(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff", "nldas_soilm"), 
            list(precip = function(x){x*cmip_region$narr_precip})) %>%
  mutate_at(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"),
            list(rhum = function(x){x*cmip_region$narr_rhum})) %>%  
  mutate_at(c("narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"),
            list(wspd = function(x){x*cmip_region$narr_wspd})) %>% 
  mutate_at(c("narr_vpd", "narr_runoff", "nldas_soilm"),
            list(vpd = function(x){x*cmip_region$narr_vpd})) %>% 
  mutate_at(c("narr_runoff", "nldas_soilm"),
            list(runoff = function(x){x*cmip_region$narr_runoff})) %>% 
  mutate(nldas_soilm_soilm = nldas_soilm^2)

var_list <- c(vlist_climate,
              paste0(vlist_climate, "_temp"),
              paste0(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff", "nldas_soilm"), "_precip"),
              paste0(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_rhum"),
              paste0(c("narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_wspd"),
              paste0(c("narr_vpd", "narr_runoff", "nldas_soilm"), "_vpd"),
              paste0(c("narr_runoff", "nldas_soilm"), "_runoff"),
              paste0(c("nldas_soilm"), "_soilm"))

##### load the mean of log-based dataset to de-mean
region_mean_log <- readRDS(paste0(result_path, "/final_models/training_data/obs_subregion_mean_log.rds"))$data_region

cmip_region_log <- left_join(cmip_region_full, region_mean_log[,c("region", "area", "outcome")]) %>%
  rename(outcome_mean_log=outcome) %>%
  filter(!is.na(outcome_mean_log))

for (vv in var_list){
  region_mean_tmp <- region_mean_log[,c("region",vv)] %>%
    rename(regionsub_mean = !!as.name(vv))
  
  cmip_region_log <- left_join(cmip_region_log, region_mean_tmp) %>%
    mutate(!!as.name(vv):= !!as.name(vv) - regionsub_mean) %>% select(-regionsub_mean)
}
saveRDS(cmip_region_log, paste0(result_path,"/CMIP6_fire_proj/cmip_region_full_demean_log.rds"))

##### load the mean of level-based dataset to de-mean
region_mean_level <- readRDS(paste0(result_path, "/final_models/training_data/obs_subregion_mean_level.rds"))$data_region

cmip_region_level <- left_join(cmip_region_full, region_mean_level[,c("region", "area", "outcome")]) %>%
  rename(outcome_mean_level=outcome) %>%
  filter(!is.na(outcome_mean_level))

for (vv in var_list){
  region_mean_tmp <- region_mean_level[,c("region", vv)] %>%
    rename(regionsub_mean = !!as.name(vv))
  
  cmip_region_level <- left_join(cmip_region_level, region_mean_tmp) %>%
    mutate(!!as.name(vv):= !!as.name(vv) - regionsub_mean) %>% select(-regionsub_mean)
}
saveRDS(cmip_region_level, paste0(result_path,"/CMIP6_fire_proj/cmip_region_full_demean_level.rds"))


#---------------------
## Prepare the training data at the original CMIP6 grid level (just differ by regions)
#---------------------
var_long_list <- c("temperature_surface", "precipitation", "RH_surface", "vpd",
                   "soil_moisture",  "soil_moisture", "runoff", "wind_speed")
var_labels <- c("narr_temp", "narr_precip", "narr_rhum", "narr_vpd", 
                "narr_soilm", "nldas_soilm", "narr_runoff", "narr_wspd")

cmip_region <- NULL
for (vv in 1: length(var_long_list)){
  var_long <- var_long_list[vv]  
  var_name <- var_labels[vv]  
  print(var_name)
  
  cmip_debias <- readRDS(paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/",
                                var_long, "/CMIP6_annual_debias_1deg.rds")) 
  
  if (var_name == "nldas_soilm"){
    cmip_debias <- readRDS(paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/",
                                  var_long, "/CMIP6_annual_debias_1deg_nldas.rds")) 
  }
  
  model_list <- unique(cmip_debias$model)
  
  cmip_data_combined <- NULL
  for (mm in model_list){
    print(mm)
    cmip_data <-  filter(cmip_debias, model == mm) %>% 
      left_join(crosswalk_ecoregion) %>%
      filter(!is.na(region)) %>%
      filter(year > 2070)    %>%
      group_by(model, scenario, region, year, grid_1deg_id) %>%
      mutate(area = as.numeric(area)) %>%
      summarise(cmip_debias = weighted.mean(cmip_debias, w= area, na.rm=T)) %>%
      rename(!!as.name(var_name):= cmip_debias) 
    
    cmip_data_combined <- bind_rows(cmip_data_combined, cmip_data)
  }
  
  if (vv == 1){
    cmip_region <- cmip_data_combined
  } else{
    cmip_region <- left_join(cmip_region, cmip_data_combined)
  }
}
saveRDS(cmip_region, paste0(result_path,"/CMIP6_fire_proj/final_cmip_region_1deg_data_2071_2099.rds"))

# gfed_1deg <- readRDS(paste0(result_path,"/final_models/training_data/final_train_region_FE_level_1deg_grid.rds"))
# cmip_eco2_1 <- left_join(cmip_eco2, crossec_eco2)
# saveRDS(cmip_eco2_1, paste0(result_path,"/CMIP6_fire_proj/final_cmip_eco2_data.rds"))
##### generate extra features for regional data and de-mean based on mean of subregions ##### 
cmip_region <- readRDS(paste0(result_path,"/CMIP6_fire_proj/final_cmip_region_1deg_data.rds")) %>% 
  ungroup() 
cmip_region1 <- readRDS(paste0(result_path,"/CMIP6_fire_proj/final_cmip_region_1deg_data_2071_2099.rds")) %>% 
  ungroup() 
cmip_region <- bind_rows(cmip_region, cmip_region1)

cmip_region <- left_join(cmip_region, crossec_1deg)

vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
                   "narr_vpd", "narr_runoff","nldas_soilm")

cmip_region_full <- cmip_region %>%
   mutate_at(vlist_climate, list(forest = function(x){x*cmip_region$forest},
                                crop = function(x){x*cmip_region$cropland},
                                grass = function(x){x*cmip_region$grassland})) %>%
  mutate_at(vlist_climate, list(temp = function(x){x*cmip_region$narr_temp})) %>%
  mutate_at(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff", "nldas_soilm"), 
            list(precip = function(x){x*cmip_region$narr_precip})) %>%
  mutate_at(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"),
            list(rhum = function(x){x*cmip_region$narr_rhum})) %>%  
  mutate_at(c("narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"),
            list(wspd = function(x){x*cmip_region$narr_wspd})) %>% 
  mutate_at(c("narr_vpd", "narr_runoff", "nldas_soilm"),
            list(vpd = function(x){x*cmip_region$narr_vpd})) %>% 
  mutate_at(c("narr_runoff", "nldas_soilm"),
            list(runoff = function(x){x*cmip_region$narr_runoff})) %>% 
  mutate(nldas_soilm_soilm = nldas_soilm^2)

var_list <- c(vlist_climate)
              , paste(rep(vlist_climate, each = 3),
                                   c("forest","crop", "grass"), sep="_"),
              paste0(vlist_climate, "_temp"),
              paste0(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff", "nldas_soilm"), "_precip"),
              paste0(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_rhum"),
              paste0(c("narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_wspd"),
              paste0(c("narr_vpd", "narr_runoff", "nldas_soilm"), "_vpd"),
              paste0(c("narr_runoff", "nldas_soilm"), "_runoff"),
              paste0(c("nldas_soilm"), "_soilm"))

##### load the mean of log-based dataset to de-mean
region_mean_log <- readRDS(paste0(result_path, "/final_models/training_data/obs_grid_1deg_mean_log.rds"))
region_mean_log$grid_1deg_id <- as.numeric(as.character(region_mean_log$grid_1deg_id))

cmip_region_log <- left_join(cmip_region_full[,c("model", "scenario","region", "year","grid_1deg_id", var_list)], 
                             region_mean_log[,c("region", "grid_1deg_id","area", "outcome")]) %>%
  rename(outcome_mean_log=outcome) %>%
  filter(!is.na(outcome_mean_log))

for (vv in var_list){
  region_mean_tmp <- region_mean_log[,c("region","grid_1deg_id",vv)] %>%
    rename(regionsub_mean = !!as.name(vv))
  
  cmip_region_log <- left_join(cmip_region_log, region_mean_tmp) %>%
    mutate(!!as.name(vv):= !!as.name(vv) - regionsub_mean) %>% select(-regionsub_mean)
}
saveRDS(cmip_region_log, paste0(result_path,"/CMIP6_fire_proj/cmip_region_1deg_full_demean_log_2099.rds"))

##### load the mean of level-based dataset to de-mean
region_mean_level <- readRDS(paste0(result_path, "/final_models/training_data/obs_grid_1deg_mean_level.rds"))
region_mean_level$grid_1deg_id <- as.numeric(as.character(region_mean_level$grid_1deg_id))

cmip_region_level <- left_join(cmip_region_full[,c("model", "scenario","region", "year","grid_1deg_id", var_list)], 
                               region_mean_level[,c("region", "grid_1deg_id","area", "outcome")]) %>%
  rename(outcome_mean_level=outcome) %>%
  filter(!is.na(outcome_mean_level))

for (vv in var_list){
  region_mean_tmp <- region_mean_level[,c("region","grid_1deg_id",vv)] %>%
    rename(regionsub_mean = !!as.name(vv))
  
  cmip_region_level <- left_join(cmip_region_level, region_mean_tmp) %>%
    mutate(!!as.name(vv):= !!as.name(vv) - regionsub_mean) %>% select(-regionsub_mean)
}
saveRDS(cmip_region_level, paste0(result_path,"/CMIP6_fire_proj/cmip_region_1deg_full_demean_level_2099.rds"))


# #---------------------
# ## Prepare the training data from CMIP to GFED grid cell level
# #---------------------
# var_long_list <- c("temperature_surface", "precipitation", "RH_surface", "vpd",
#                    "soil_moisture",  "soil_moisture", "runoff", "wind_speed")
# var_labels <- c("narr_temp", "narr_precip", "narr_rhum", "narr_vpd", 
#                 "narr_soilm", "nldas_soilm", "narr_runoff", "narr_wspd")
# 
# #### Only doing this for Mexico now
# cmip_grid <- NULL
# for (vv in 1: length(var_long_list)){
#   var_long <- var_long_list[vv]  
#   var_name <- var_labels[vv]  
#   print(var_name)
#   
#   cmip_debias <- readRDS(paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/",
#                                 var_long, "/CMIP6_annual_debias_1deg.rds")) 
#   
#   if (var_name == "nldas_soilm"){
#     cmip_debias <- readRDS(paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/",
#                                   var_long, "/CMIP6_annual_debias_1deg_nldas.rds")) 
#   }
#   
#   model_list <- unique(cmip_debias$model)
#   
#   cmip_data_combined <- NULL
#   for (mm in model_list){
#     print(mm)
#     cmip_data <-  filter(cmip_debias, model == mm) %>% 
#       left_join(crosswalk_ecoregion) %>%
#       filter(region=="Mexico") %>%
#       distinct(model, scenario, year, cmip_debias, gfed_cell_id, region) %>%
#       rename(!!as.name(var_name):= cmip_debias)
#     
#     cmip_data_combined <- bind_rows(cmip_data_combined, cmip_data)
#   }
#   
#   if (vv == 1){
#     cmip_grid <- cmip_data_combined
#   } else{
#     cmip_grid <- left_join(cmip_grid, cmip_data_combined)
#   }
# }
# saveRDS(cmip_grid, paste0(result_path,"/CMIP6_fire_proj/meixco_cmip_gfed_grid_data.rds"))
# 
# 
# ##### generate extra features for regional data and de-mean based on mean of subregions ##### 
# cmip_grid <- readRDS(paste0(result_path,"/CMIP6_fire_proj/meixco_cmip_gfed_grid_data.rds")) %>% ungroup() 
# 
# vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
#                    "narr_vpd", "narr_runoff")
# 
# ##### load the mean of level-based dataset to de-mean
# grid_mean_level <- readRDS(paste0(result_path, "/final_models/training_data/obs_subregion_mean_level.rds"))$data_grid %>%
#   filter(region == "Mexico") %>%
#   mutate(gfed_cell_id = as.numeric(as.character(gfed_cell_id)))
# 
# cmip_grid_level <- left_join(cmip_grid, grid_mean_level[,c("gfed_cell_id", "area", "outcome")]) %>%
#   rename(outcome_mean_level=outcome) %>%
#   filter(!is.na(outcome_mean_level))
# 
# var_list <- vlist_climate
# 
# for (vv in var_list){
#   grid_mean_tmp <- grid_mean_level[,c("gfed_cell_id", vv)] %>%
#     rename(gridsub_mean = !!as.name(vv))
#   
#   cmip_grid_level <- left_join(cmip_grid_level, grid_mean_tmp) %>%
#     mutate(!!as.name(vv):= !!as.name(vv) - gridsub_mean) %>% select(-gridsub_mean)
# }
# #cmip_grid_level <- cmip_grid_level %>% select(-nldas_soilm, -narr_soilm)
# saveRDS(cmip_grid_level, paste0(result_path,"/CMIP6_fire_proj/mexico_cmip_gfed_grid_demean_level.rds"))






#------------------------------------------------
####### Plot the CMIP6 climate after aggregations
#------------------------------------------------
color_map <- c("grey65",
               rev(met.brewer(name="Homer1", n=8, type="discrete"))[c(2,6,8)])

cmip_region <- readRDS(paste0(result_path,"/CMIP6_fire_proj/final_cmip_region_data_veg.rds")) %>%
  pivot_longer(cols = narr_temp:narr_wspd) %>%
  rename(gcm = model) %>% ungroup()

cmip_2001_2021 <- filter(cmip_region, scenario%in%c("historical", "ssp370"), year %in% seq(2001,2021)) 
cmip_2001_2021_mean <- cmip_2001_2021 %>%
  group_by(gcm, region, name) %>%
  summarise(mean_value= mean(value, na.rm=T)) %>% ungroup()

cmip_region <-  bind_rows(filter(cmip_region, year>2021),
                          cmip_2001_2021 %>% mutate(scenario = "historical"),
                          cmip_2001_2021 %>% mutate(scenario = "ssp126"),
                          cmip_2001_2021 %>% mutate(scenario = "ssp245"),
                          cmip_2001_2021 %>% mutate(scenario = "ssp370"),
                          cmip_2001_2021 %>% mutate(scenario = "ssp585")) %>%
  left_join(cmip_2001_2021_mean) %>%
  mutate(anomaly = value - mean_value)

cmip_10yr_smooth <- cmip_region %>%
  group_by(scenario, region, gcm, name) %>%
  arrange(year) %>%
  mutate(anomaly_10yr=rollapplyr(anomaly, width = 10, FUN=mean, partial=T),
         value_10yr=rollapplyr(value, width = 10, FUN=mean, partial=T)) %>%
  ungroup()

cmip_10yr_smooth <- cmip_10yr_smooth %>%
  group_by(scenario, year, region, name) %>%
  summarise(anomaly_median=median(anomaly_10yr, na.rm=T), 
            anomaly_10=quantile(anomaly_10yr, 0.1, na.rm=T),
            anomaly_90=quantile(anomaly_10yr, 0.9, na.rm=T),
            anomaly_25=quantile(anomaly_10yr, 0.25, na.rm=T),
            anomaly_75=quantile(anomaly_10yr, 0.75, na.rm=T)) %>% ungroup() %>%
  filter(scenario=="historical" | year>2021)

## just to connect 2014 and 2015 (for visualization purpose)
cmip_10yr_smooth <- bind_rows(cmip_10yr_smooth,
                              filter(cmip_10yr_smooth, year==2021) %>% mutate(scenario = "ssp126"),
                              filter(cmip_10yr_smooth, year==2021) %>% mutate(scenario = "ssp245"),
                              filter(cmip_10yr_smooth, year==2021) %>% mutate(scenario = "ssp370"),
                              filter(cmip_10yr_smooth, year==2021) %>% mutate(scenario = "ssp585"))

# cmip_10yr_smooth$region <- factor(cmip_10yr_smooth$region,
#                                   levels=c("Canada-Alaska", "Northeastern US", "Southeastern US",
#                                            "Western US", "Mexico"))

rrr <- "Western US"
ggplot(filter(cmip_10yr_smooth, year<2065, scenario!="ssp585", region==rrr), 
       aes(x=year, y= anomaly_median, colour=scenario, fill=scenario)) +
  geom_line(size=1) +
  geom_vline(xintercept = c(2047, 2049), linetype="dashed") +
  geom_ribbon(aes(ymin=anomaly_25, ymax=anomaly_75), alpha=0.3, colour=NA) +
  scale_fill_manual(values = color_map) +
  scale_colour_manual(values = color_map) +
  theme_classic() + theme(text = element_text(size=16)) +
  labs(x="", fill="", colour="") +
  facet_wrap(~name, scales = "free")
ggsave("west_CMIP6_debias_10yr_smooth.png", width=9, height=6.5)

ggplot(filter(cmip_10yr_smooth, year<2065, scenario=="ssp370", region==rrr, name=="nldas_soilm"), 
       aes(x=year, y= anomaly_median, colour=scenario, fill=scenario)) +
  geom_line(size=1) +
  geom_vline(xintercept = c(2046, 2049), linetype="dashed") +
  geom_ribbon(aes(ymin=anomaly_25, ymax=anomaly_75), alpha=0.3, colour=NA) +
  scale_fill_manual(values = color_map) +
  scale_colour_manual(values = color_map) +
  theme_classic() + theme(text = element_text(size=16)) +
  labs(x="", fill="", colour="") +
  facet_wrap(~name, scales = "free")
ggsave("west_CMIP6_debias_10yr_smooth.png", width=9, height=6.5)

#individual models

model_list <- c("CanESM5", "EC-Earth3","MRI-ESM2-0","MIROC-ES2L")
ggplot(filter(cmip_10yr_smooth, year<2065, region==rrr, gcm%in%model_list, scenario!="ssp585", name!="narr_soilm"), 
       aes(x=year, y= value_10yr, colour=scenario)) +
  geom_line(size=1) +
  #geom_vline(xintercept = c(2046, 2049), linetype="dashed") +
  #geom_ribbon(aes(ymin=anomaly_25, ymax=anomaly_75), alpha=0.3, colour=NA) +
  # scale_fill_manual(values = color_map) +
  scale_colour_manual(values = color_map) +
  theme_classic() + theme(text = element_text(size=16)) +
  labs(x="", fill="", colour="") +
  facet_grid(name~gcm, scales = "free")
ggsave("west_gcm_diag.png", width=10, height=7.5)

# cmip_region <- readRDS(paste0(result_path,"/CMIP6_fire_proj/final_cmip_region_data_veg.rds"))
# normalize_baseline <- readRDS(paste0(result_path,"/final_models/west_regional_data_norm_veg.rds"))
# 
# vlist <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
#            "narr_vpd", "narr_runoff","nldas_soilm") 
# 
# for (vv in vlist){
#   mean_tmp <- filter(normalize_baseline, variable == vv, stat=="mean")$value
#   sd_tmp <- filter(normalize_baseline, variable == vv, stat=="sd")$value
#   cmip_region[,vv] <- (cmip_region[,vv] - mean_tmp)/sd_tmp
# }
# saveRDS(cmip_region, paste0(result_path,"/CMIP6_fire_proj/final_cmip_region_data_normalize_veg.rds"))


