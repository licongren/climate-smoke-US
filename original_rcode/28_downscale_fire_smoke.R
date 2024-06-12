library(terra);library(ncdf4)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(plotly)
library(sp);library(sf)
library(dplyr)   
library(fixest)
library(tidyverse)
library(MetBrewer)
library(stringr)
library(splines2)
library(exactextractr)
library(zoo)

rm(list=ls())
gc()

code_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/rcode"
data_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data"
figure_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/figures"
result_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/results/fire_model"

db_path <- "/Users/mhqiu/BurkeLab Dropbox/Projects/smoke-climate"

setwd(figure_path)
states_map <- map_data("state")

color_map <- c("grey65",
               rev(met.brewer(name="Homer1", n=8, type="discrete"))[c(2,6,8)])

#-------------------------------------------------------------------------------
# Downscale the projected fire DM to grid cell level
# Written by Minghao
# Last edited October 2022
#-------------------------------------------------------------------------------
## load historical gfed emissions
gfed_grid <- readRDS(paste0(data_path, "/GFED4S_NorthAmerica_025deg_grid_EPSG4326.rds"))

gfed_ecoregion_link <- readRDS(paste0(db_path, "/crosswalk_GFED_NA_CEC_Eco_Level2_3.rds")) %>%
  mutate(intersection_area=as.numeric(intersection_area)) %>%   group_by(gfed_cell_id)  %>%
  arrange(desc(intersection_area)) %>%
  mutate(area_weight=intersection_area/sum(intersection_area)) %>% ungroup() # %>% distinct(gfed_cell_id, .keep_all=T)

gfed_wang_region <- readRDS(paste0(data_path, "/crosswalk_GFED_wang_regions.rds")) 
crosswalk_1deg <-  readRDS(paste0(data_path, "/crosswalk_gfed_1deg.rds"))

gfed_data <- readRDS(file.path(db_path, "climate-fire_training_data_2001-2021_clean_varname.rds"))
gfed_data <- left_join(gfed_data, gfed_wang_region, by="gfed_cell_id")
gfed_data <- left_join(gfed_data, gfed_ecoregion_link, by="gfed_cell_id")
gfed_data <- left_join(gfed_data, crosswalk_1deg[,c("gfed_cell_id", "coverage_area", "grid_1deg_id")], 
                       by="gfed_cell_id")

######## scale using aggregated emissions over the entire period
#downscale <- "average" ### by year
gfed_grid_annual <- gfed_data %>% 
  group_by(gfed_cell_id, region, NA_L2KEY, NA_L3KEY, grid_1deg_id) %>%
  summarise(DM_kg=sum(DM_kg*area_weight, na.rm=T)) %>%
  ungroup() 
  
gfed_grid_annual$grid_1deg_id <- as.character(gfed_grid_annual$grid_1deg_id)
#------------------------------
# load future CMIP6 fire projections
#------------------------------
# fire_pred <- readRDS(paste0(result_path,"/CMIP6_fire_proj/final_pred_best_models_threshold_05.rds"))

# sample_proj1 <- fire_pred %>% filter(algorithm_outcome=="Linear, level",
#                                     region=="Western US", year==2030)
# 
# sample_proj2 <- fire_pred %>% filter(algorithm_outcome=="Linear, level",
#                                      region=="Northeastern US", year==2030)
# 
# fire_proj = sample_proj1
# fire_var = "pred_final_10yr_positive"
# scale_data= gfed_grid_annual

#------------------------------
# function for downscaling
#------------------------------
downscale_fire_emis <- function(fire_proj=fire_proj, fire_var=fire_var, scale_data=scale_data){
  spatial_res = unique(fire_proj$spatial_res)
  # if (spatial_res != agg_level){
  #   print("spatial scales not matched")
  #   break
  # }
  if (spatial_res == "regional"){
    agg_level = "region"
  }  else if (spatial_res == "eco2"){
    agg_level = "NA_L2KEY"
  } else if (spatial_res == "eco3"){
    agg_level = "NA_L3KEY"
  } else if (spatial_res == "grid"){
    agg_level = "grid_1deg_id"
  } 
  
  if (spatial_res!="regional"){
      scale_data = scale_data %>% group_by(gfed_cell_id, !!as.name(agg_level), region) %>%
      summarise(grid_emis=sum(DM_kg, na.rm=T)) %>% ungroup() %>%
      group_by(!!as.name(agg_level), region) %>%
      mutate(total_emis=sum(grid_emis, na.rm=T),
             grid_ratio=grid_emis/total_emis) %>% ungroup() %>%
      select(-grid_emis, -total_emis)
    
      fire_proj = fire_proj %>%
        rename(!!as.name(agg_level) := region_sub) 
      
    fire_scale = full_join(fire_proj, scale_data, by=c("region", agg_level)) %>%
      mutate(pred_DM_grid = !!as.name(fire_var)  * grid_ratio) %>%
      filter(!is.na(pred_DM_grid))
  }
  else{
    scale_data = scale_data %>% group_by(gfed_cell_id, region) %>%
      summarise(grid_emis=sum(DM_kg, na.rm=T)) %>% ungroup() %>%
      group_by(region) %>%
      mutate(total_emis=sum(grid_emis, na.rm=T),
             grid_ratio=grid_emis/total_emis) %>% ungroup() %>%
      select(-grid_emis, -total_emis)
    
    fire_scale = full_join(fire_proj, scale_data, by="region") %>%
      mutate(pred_DM_grid = !!as.name(fire_var) * grid_ratio) %>%
      filter(!is.na(pred_DM_grid))
  }
  fire_scale = fire_scale %>% select(gcm, scenario, region, year, algorithm_outcome,
                                       gfed_cell_id, pred_DM_grid)
  return(fire_scale)
}

#------------------------------
# downscaling emissions 
#------------------------------
#fire_pred <- readRDS(paste0(result_path,"/CMIP6_fire_proj/final_pred_best_models_subregion.rds"))
fire_pred <- readRDS(paste0(result_path,"/CMIP6_fire_proj/final_pred_fire_2001_2099_annualmax20.rds"))

region_algorithm <- fire_pred %>% distinct(region, algorithm_outcome)
year_list <-  seq(2015,2095,by=10) # 2055 #c(2015, 2025, 2035, 2045, 2055)
fire_grid_combined <- NULL
for (ii in 1:nrow(region_algorithm)){
  print(ii)
  region_tmp <- region_algorithm[ii, "region"] %>% unlist() %>% as.vector()
  algorithm_tmp <- region_algorithm[ii, "algorithm_outcome"] %>% unlist() %>% as.vector() 

  for (yy in year_list){  
  fire_scale_tmp <- downscale_fire_emis(fire_proj = filter(fire_pred, 
                                                           region == region_tmp,
                                                           algorithm_outcome == algorithm_tmp,
                                                           year==yy),
                                        fire_var = "pred_final_10yr_positive",
                                        scale_data= gfed_grid_annual)
  
  fire_grid_combined <- bind_rows(fire_grid_combined, fire_scale_tmp)
  }
}
# saveRDS(fire_grid_combined, paste0(result_path,"/CMIP6_fire_proj/downscaled_emis_final_best_models_threshold_05.rds"))
# 
# fire_grid_combined <- readRDS(paste0(result_path,"/CMIP6_fire_proj/downscaled_emis_final_best_models_threshold_05.rds"))

fire_grid_combined <- fire_grid_combined %>% 
  group_by(gcm, scenario, gfed_cell_id, year, algorithm_outcome) %>% ## sum across cells that cover more than 1 region
  summarise(pred_DM_grid=sum(pred_DM_grid)) %>% ungroup() %>%
  group_by(gcm, scenario, gfed_cell_id, year) %>% ## ##mean across algorithms
  summarise(pred_DM_grid=mean(pred_DM_grid)) %>% ungroup()
#saveRDS(fire_grid_combined, paste0(result_path,"/CMIP6_fire_proj/downscaled_emis_final_threshold_05.rds"))
#saveRDS(fire_grid_combined, paste0(result_path,"/CMIP6_fire_proj/downscaled_emis_final_subregion.rds"))
saveRDS(fire_grid_combined, paste0(result_path,"/CMIP6_fire_proj/downscaled_emis_10yr_2010_2090.rds"))

######################## save all 2025-2055 years
region_algorithm <- fire_pred %>% distinct(region, algorithm_outcome)
fire_pred[fire_pred$pred_DM_final<0, "pred_DM_final"] <- 0

year_list <- c(2025:2045) #c(2046:2055)
fire_grid_combined <- NULL
for (ii in 1:nrow(region_algorithm)){
  print(ii)
  region_tmp <- region_algorithm[ii, "region"] %>% unlist() %>% as.vector()
  algorithm_tmp <- region_algorithm[ii, "algorithm_outcome"] %>% unlist() %>% as.vector() 
  
  for (yy in year_list){  
    fire_scale_tmp <- downscale_fire_emis(fire_proj = filter(fire_pred, 
                                                             region == region_tmp,
                                                             algorithm_outcome == algorithm_tmp,
                                                             year==yy),
                                          fire_var = "pred_DM_final",
                                          scale_data= gfed_grid_annual)
    
    fire_grid_combined <- bind_rows(fire_grid_combined, fire_scale_tmp)
  }
}
# saveRDS(fire_grid_combined, paste0(result_path,"/CMIP6_fire_proj/downscaled_emis_2046_2055_years.rds"))
# 
# fire_grid_combined <- readRDS(paste0(result_path,"/CMIP6_fire_proj/downscaled_emis_2046_2055_years.rds"))

fire_grid_combined <- fire_grid_combined %>% 
  group_by(gcm, scenario, gfed_cell_id, year, algorithm_outcome) %>%
  summarise(pred_DM_grid=sum(pred_DM_grid)) %>% ungroup() %>%
  group_by(gcm, scenario, gfed_cell_id, year) %>%
  summarise(pred_DM_grid=mean(pred_DM_grid)) %>% ungroup()
saveRDS(fire_grid_combined, paste0(result_path,"/CMIP6_fire_proj/downscaled_emis_2025_2045_years.rds"))


######################## down scale the 2055 emissions geenrated by the boostraps runs 
fire_pred <- readRDS(paste0(result_path,"/CMIP6_fire_proj/pred_fire_model_uncertainty_boot_2050.rds")) 
region_algorithm <- fire_pred %>% distinct(region, algorithm_outcome)
nboot <- max(fire_pred$bootid)

fire_grid_combined <- NULL
for (bbb in 1:nboot){  
  print(bbb)

fire_grid_tmp <- NULL
for (ii in 1:nrow(region_algorithm)){
  print(ii)
  region_tmp <- region_algorithm[ii, "region"] %>% unlist() %>% as.vector()
  algorithm_tmp <- region_algorithm[ii, "algorithm_outcome"] %>% unlist() %>% as.vector() 

  fire_proj_tmp <- filter(fire_pred,  bootid == bbb,
                          region == region_tmp, algorithm_outcome == algorithm_tmp)
  if (nrow(fire_proj_tmp) == 0){next}
  fire_scale_tmp <- downscale_fire_emis(fire_proj = fire_proj_tmp ,
                                        fire_var = "pred_final_10yr_positive",
                                        scale_data= gfed_grid_annual)
  fire_grid_tmp <- bind_rows(fire_grid_tmp, fire_scale_tmp)
}
    
fire_grid_tmp <-  fire_grid_tmp %>% 
  group_by(gcm, scenario, region, year, gfed_cell_id) %>% ## mean across algorithms 
              summarise(pred_DM_grid = mean(pred_DM_grid, na.rm=T)) %>% ungroup() %>%
  group_by(scenario, region, year, gfed_cell_id) %>% ## mean across gcms
  summarise(pred_DM_grid = mean(pred_DM_grid, na.rm=T)) %>% ungroup() %>%
  group_by(scenario, year, gfed_cell_id) %>% ## summ across grid that could span more than one region
  summarise(pred_DM_grid = mean(pred_DM_grid, na.rm=T)) %>% ungroup() %>%
  mutate(bootid=bbb)
  
fire_grid_combined <- bind_rows(fire_grid_combined, fire_grid_tmp)
}

## only keep the bootruns with all method available
# boot_list <- fire_grid_combined %>% distinct(bootid, region, algorithm_outcome) %>%
#   group_by(bootid) %>%
#   summarise(n=length(bootid)) %>% ungroup()

saveRDS(fire_grid_combined, paste0(result_path,"/CMIP6_fire_proj/downscaled_emis_uncertainty_boot.rds"))

######################## down scale the 2055, 
#------   but with a MC draw of bootid and GCMs
fire_pred <- readRDS(paste0(result_path,"/CMIP6_fire_proj/pred_fire_model_uncertainty_boot_2050.rds")) 
region_algorithm <- fire_pred %>% distinct(region, algorithm_outcome)
nboot <- max(fire_pred$bootid)

boot_gcm <- fire_pred %>% distinct(bootid, gcm)

set.seed(100)
boot_gcm_100 <- sample_n(boot_gcm, size=100, replace = F)
boot_gcm_100$mc_id = 1:100
saveRDS(boot_gcm_100, paste0(result_path,"/CMIP6_fire_proj/MC_100_fire_gcm_list.rds"))


fire_grid_combined <- NULL
for (bbb in 1:nrow(boot_gcm_100)){  
  print(bbb)
  boot_tmp <- as.numeric(boot_gcm_100[bbb,"bootid"])
  gcm_tmp  <- as.character(boot_gcm_100[bbb,"gcm"])
  
  fire_grid_tmp <- NULL
  for (ii in 1:nrow(region_algorithm)){
    print(ii)
    region_tmp <- region_algorithm[ii, "region"] %>% unlist() %>% as.vector()
    algorithm_tmp <- region_algorithm[ii, "algorithm_outcome"] %>% unlist() %>% as.vector() 
    
    fire_proj_tmp <- filter(fire_pred,  bootid == boot_tmp, gcm== gcm_tmp,
                            region == region_tmp, algorithm_outcome == algorithm_tmp)
    if (nrow(fire_proj_tmp) == 0){next}
    fire_scale_tmp <- downscale_fire_emis(fire_proj = fire_proj_tmp ,
                                          fire_var = "pred_final_10yr_positive",
                                          scale_data= gfed_grid_annual)
    fire_grid_tmp <- bind_rows(fire_grid_tmp, fire_scale_tmp)
  }
  
  fire_grid_tmp <-  fire_grid_tmp %>% 
    group_by(gcm, scenario, region, year, gfed_cell_id) %>% ## mean across algorithms 
    summarise(pred_DM_grid = mean(pred_DM_grid, na.rm=T)) %>% ungroup() %>%
    group_by(scenario, region, year, gfed_cell_id) %>% ## mean across gcms
    summarise(pred_DM_grid = mean(pred_DM_grid, na.rm=T)) %>% ungroup() %>%
    group_by(scenario, year, gfed_cell_id) %>% ## summ across grid that could span more than one region
    summarise(pred_DM_grid = mean(pred_DM_grid, na.rm=T)) %>% ungroup() %>%
    mutate(mc_id=bbb)
  
  fire_grid_combined <- bind_rows(fire_grid_combined, fire_grid_tmp)
}

## only keep the bootruns with all method available
# boot_list <- fire_grid_combined %>% distinct(bootid, region, algorithm_outcome) %>%
#   group_by(bootid) %>%
#   summarise(n=length(bootid)) %>% ungroup()

saveRDS(fire_grid_combined, paste0(result_path,"/CMIP6_fire_proj/downscaled_emis_MC_100_fire_gcm.rds"))




# aa <- fire_grid_combined %>% 
#   filter(year==2055) %>%
#   group_by(region, algorithm_outcome) %>%
#   summarise(tot_emis=sum(pred_DM_grid)) %>% ungroup()
# 
# bb <- fire_pred %>% 
#   filter(year==2055) %>%
#   group_by(region, algorithm_outcome) %>%
#   summarise(tot_emis=sum(pred_final_10yr_positive)) %>% ungroup()

