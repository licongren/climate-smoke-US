library(terra);library(ncdf4)
library(sf)
library(dplyr)   
library(fixest)
library(tidyverse)
library(purrr)
library(data.table)
library(mlr3verse)
library(xgboost)
library(nnet)
library(ranger)
library(future)

rm(list=ls())
gc()

input_path <- "~/oak_space/mhqiu/smoke_climate/input_data"
output_path <- "~/oak_space/mhqiu/smoke_climate/output"

gfed_grid <- readRDS(paste0(input_path, "/GFED4S_NorthAmerica_025deg_grid_EPSG4326.rds"))

gfed_ecoregion_link <- readRDS(paste0(input_path, "/crosswalk_GFED_NA_CEC_Eco_Level2_3.rds")) %>%
  mutate(intersection_area=as.numeric(intersection_area)) %>%   group_by(gfed_cell_id)  %>%
  arrange(desc(intersection_area)) %>%
  mutate(area_weight=intersection_area/sum(intersection_area)) %>% ungroup() # %>% distinct(gfed_cell_id, .keep_all=T)

gfed_region <- readRDS(paste0(input_path, "/crosswalk_GFED_wang_regions.rds"))

gfed_data <- readRDS(file.path(input_path, "climate-fire_training_data_2001-2021_clean_varname.rds")) 

gfed_data <- left_join(gfed_data, gfed_region, by="gfed_cell_id")
gfed_data <- left_join(gfed_data, gfed_ecoregion_link, by="gfed_cell_id")

# gfed_agg_data <- readRDS(paste0(output_path, "/gfed_agg_data.rds"))
# data_region <- gfed_agg_data$data_region 
# data_eco3 <- gfed_agg_data$data_eco3 
# data_eco2 <- gfed_agg_data$data_eco2 
# data_grid <- gfed_agg_data$data_grid 
# 
# ##### drop the eco-3 and eco-2 regions with very few grid cells
# eco3_summ <- gfed_data %>% 
#   group_by(region, NA_L3KEY) %>% summarise(n=length(region)/252) %>%
#   arrange(n) %>% ungroup() %>% mutate(region_sub=factor(NA_L3KEY)) %>% select(-NA_L3KEY)
# 
# eco2_summ <- gfed_data %>% 
#   group_by(region, NA_L2KEY) %>% summarise(n=length(region)/252) %>%
#   arrange(n) %>% ungroup() %>% mutate(region_sub=factor(NA_L2KEY)) %>% select(-NA_L2KEY)
# 
# ### only select regions with more than 10 grid cells in the region
# data_eco3 <- left_join(data_eco3, eco3_summ) %>% filter(n>10) %>% select(-n)
# 
# data_eco2 <- left_join(data_eco2, eco2_summ) %>% filter(n>10) %>% select(-n)

### CV linear model
#region_list <- unique(gfed_data$region) #c("Western US", "Northeastern US", "Southeastern US")  #
LOG <- F
area_normalize <- T
region_FE <- T

# for (experiment in c("eco3", "eco2","grid","regional")){ 
#   print(paste("exp:", experiment))
#   if (experiment == "regional"){
#     data_tmp <- data_region
#   }else if (experiment == "eco2"){
#     data_tmp <- data_eco2
#   }else if (experiment == "eco3"){
#     data_tmp <- data_eco3
#   }else if (experiment == "grid"){
#     data_tmp <- data_grid 
#   }

#### de-mean the data of 1deg cell  
data_tmp <- readRDS(paste0(output_path, "/gfed_1deg_data.rds"))

experiment <- "grid"  
  if (area_normalize){
    data_tmp <- data_tmp %>% mutate(outcome = DM_kg / area)
  } else{
    data_tmp <- data_tmp %>% mutate(outcome = DM_kg)
  }
    
    if (LOG){
      data_tmp <- data_tmp %>% filter(outcome>0) %>% mutate(outcome=log(outcome))
    }
  
  vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
                     "narr_vpd", "narr_runoff","nldas_soilm")

    if (experiment != "regional"){
      data_tmp <- data_tmp %>%
        mutate_at(vlist_climate, list(forest = function(x){x*data_tmp$forest},
                                      crop = function(x){x*data_tmp$cropland},
                                      grass = function(x){x*data_tmp$grassland})) %>%
        mutate_at(vlist_climate, list(temp = function(x){x*data_tmp$narr_temp})) %>%
        mutate_at(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff", "nldas_soilm"), 
                  list(precip = function(x){x*data_tmp$narr_precip})) %>%
        mutate_at(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"),
                  list(rhum = function(x){x*data_tmp$narr_rhum})) %>%  
        mutate_at(c("narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"),
                  list(wspd = function(x){x*data_tmp$narr_wspd})) %>% 
        mutate_at(c("narr_vpd", "narr_runoff", "nldas_soilm"),
                  list(vpd = function(x){x*data_tmp$narr_vpd})) %>% 
        mutate_at(c("narr_runoff", "nldas_soilm"),
                  list(runoff = function(x){x*data_tmp$narr_runoff})) %>% 
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
    } else{
      data_tmp <- data_tmp %>%
        mutate_at(vlist_climate, list(temp = function(x){x*data_tmp$narr_temp})) %>%
        mutate_at(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff", "nldas_soilm"), 
                  list(precip = function(x){x*data_tmp$narr_precip})) %>%
        mutate_at(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"),
                  list(rhum = function(x){x*data_tmp$narr_rhum})) %>%  
        mutate_at(c("narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"),
                  list(wspd = function(x){x*data_tmp$narr_wspd})) %>% 
        mutate_at(c("narr_vpd", "narr_runoff", "nldas_soilm"),
                  list(vpd = function(x){x*data_tmp$narr_vpd})) %>% 
        mutate_at(c("narr_runoff", "nldas_soilm"),
                  list(runoff = function(x){x*data_tmp$narr_runoff})) %>% 
        mutate(nldas_soilm_soilm = nldas_soilm^2)
      
      var_list <- c(vlist_climate,
                    paste0(vlist_climate, "_temp"),
                    paste0(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff", "nldas_soilm"), "_precip"),
                    paste0(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_rhum"),
                    paste0(c("narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_wspd"),
                    paste0(c("narr_vpd", "narr_runoff", "nldas_soilm"), "_vpd"),
                    paste0(c("narr_runoff", "nldas_soilm"), "_runoff"),
                    paste0(c("nldas_soilm"), "_soilm"))
    }
    
    ##### to de-mean features
    # if (experiment =="grid"){
    #     data_tmp <- data_tmp %>%
    #       mutate(grid_1deg_id = factor(grid_1deg_id)) %>%
    #       group_by(region, grid_1deg_id) %>%
    #       mutate(sd_DM=sd(DM_kg, na.rm=T)) %>%
    #       filter(!is.na(sd_DM)) %>%
    #       mutate(outcome_mean = mean(outcome, na.rm=T),
    #              outcome_demean = outcome - outcome_mean) %>%
    #       mutate_at(var_list, .funs = function(x){x-mean(x, na.rm=T)}) %>% ungroup()
    # }
    #
    #  if (experiment !="grid"){
    #     data_tmp <- data_tmp %>%
    #       mutate(region_sub = factor(region_sub)) %>%
    #       group_by(region, region_sub) %>%
    #       mutate(sd_DM=sd(DM_kg, na.rm=T)) %>%
    #       filter(!is.na(sd_DM)) %>%
    #       mutate(outcome_mean = mean(outcome, na.rm=T),
    #              outcome_demean = outcome - outcome_mean) %>%
    #       mutate_at(var_list, .funs = function(x){x-mean(x, na.rm=T)}) %>% ungroup()
    # }
    # 
    ##### to get the mean of features
    if (experiment =="grid"){
        data_tmp <- data_tmp %>%
          mutate(grid_1deg_id = factor(grid_1deg_id)) %>%
          group_by(region, grid_1deg_id) %>%
          mutate(sd_DM=sd(DM_kg, na.rm=T)) %>%
          filter(!is.na(sd_DM)) %>%
          summarise_at(c(var_list,"outcome","area"), .funs = mean, na.rm=T) %>% ungroup()
    }

    # if (experiment !="grid"){
    #     data_tmp <- data_tmp %>%
    #       mutate(region_sub = factor(region_sub)) %>%
    #       group_by(region, region_sub) %>%
    #       mutate(sd_DM=sd(DM_kg, na.rm=T)) %>%
    #       filter(!is.na(sd_DM)) %>%
    #       summarise_at(c(var_list,"outcome","area"), .funs = mean, na.rm=T) %>% ungroup()
    # }

#   if (experiment == "regional"){
#     data_region <- data_tmp
#   }else if (experiment == "eco2"){
#     data_eco2 <- data_tmp
#   }else if (experiment == "eco3"){
#     data_eco3 <- data_tmp
#   }else if (experiment == "grid"){
#     data_grid <- data_tmp
#   }
# }
#   
# gfed_agg_data <- list(data_region = data_region, data_eco2 = data_eco2, data_eco3 = data_eco3, data_grid = data_grid)
#saveRDS(gfed_agg_data, paste0(output_path, "/final_train_region_FE_level.rds"))
  #saveRDS(data_tmp, paste0(output_path, "/final_train_region_FE_log_1deg_grid.rds"))
  saveRDS(data_tmp, paste0(output_path, "/obs_grid_1deg_mean_level.rds"))
  