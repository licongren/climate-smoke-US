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

nthread <- 3

#-------------------------------------------------------------------------------
# Use MLR3 linear models to predict GFED DM emissions at the regional levels
# Written by Minghao
# Last edited Mar 2022
#-------------------------------------------------------------------------------
gfed_grid <- readRDS(paste0(input_path, "/GFED4S_NorthAmerica_025deg_grid_EPSG4326.rds"))

gfed_ecoregion_link <- readRDS(paste0(input_path, "/crosswalk_GFED_NA_CEC_Eco_Level2_3.rds")) %>%
  mutate(intersection_area=as.numeric(intersection_area)) %>%   group_by(gfed_cell_id)  %>%
  arrange(desc(intersection_area)) %>%
  mutate(area_weight=intersection_area/sum(intersection_area)) %>% ungroup() # %>% distinct(gfed_cell_id, .keep_all=T)

gfed_region <- readRDS(paste0(input_path, "/crosswalk_GFED_wang_regions.rds"))

gfed_data <- readRDS(file.path(input_path, "climate-fire_training_data_2001-2021_clean_varname.rds")) 

gfed_data <- left_join(gfed_data, gfed_region, by="gfed_cell_id")
gfed_data <- left_join(gfed_data, gfed_ecoregion_link, by="gfed_cell_id")

### CV linear model
num_temporal_folds = length(unique(gfed_data$temporal_fold))
outer_resampling = rsmp("cv", folds = num_temporal_folds)

region_list <- unique(gfed_data$region)
test <- NULL 
LOG <- T
LOOCV <- T

### CV linear model
if (LOOCV){
num_temporal_folds = length(unique(gfed_data$year))
outer_resampling = rsmp("cv", folds = num_temporal_folds)
}

##### Load training data #####
if (LOG){
  data_full <- readRDS(paste0(output_path, "/final_train_region_FE_log.rds"))
  data_grid <- readRDS(paste0(output_path, "/final_train_region_FE_log_1deg_grid.rds"))
} else{
  data_full <- readRDS(paste0(output_path, "/final_train_region_FE_level.rds"))
  data_grid <- readRDS(paste0(output_path, "/final_train_region_FE_level_1deg_grid.rds"))
}

for (experiment in c("eco3", "eco2","grid","regional")){ 
  print(paste("exp:", experiment))
  if (experiment == "regional"){
    data <- data_full$data_region
  }else if (experiment == "eco2"){
    data <- data_full$data_eco2
  }else if (experiment == "eco3"){
    data <- data_full$data_eco3
  }else if (experiment == "grid"){
    #data <- data_full$data_grid 
    data <- data_grid 
  }
  
  for (rrr in region_list){
    print(rrr)
    data_train <- data %>% filter(region == rrr)
    
    vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
                       "narr_vpd", "narr_runoff",
                       "nldas_soilm")
    
    if (rrr %in% c("Canada-Alaska","Mexico")){
      data_train <- data_train %>% select(-nldas_soilm) 
      vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
                         "narr_vpd", "narr_runoff")
    }
    
        if (experiment =="grid"){
          data_folds <- data_train %>%
            mutate(row_ids = row_number()) %>%
            select(row_ids, year, temporal_fold, region, grid_1deg_id, area, outcome_mean)
          
           data_train <- data_train[,c("year", "temporal_fold", "region", "grid_1deg_id", 
                                       "outcome_demean",
                                       vlist_climate)]

          task = as_task_regr(data_train, target = "outcome_demean")

          task$set_col_roles(
            "year", roles = "group")$set_col_roles(
              c("temporal_fold","region", "grid_1deg_id"),remove_from = "feature")

         #task$set_col_roles("sd_DM",roles = "weight")
        }
    
     if (experiment !="grid"){
       data_folds <- data_train %>%
         mutate(row_ids = row_number()) %>%
         select(row_ids, year, temporal_fold, region, region_sub, area, outcome_mean)
       
       data_train <- data_train[,c("year", "temporal_fold", "region", "region_sub",
                                   "outcome_demean",
                                   vlist_climate)]

      task = as_task_regr(data_train, target = "outcome_demean")
      
      task$set_col_roles(
        "year", roles = "group")$set_col_roles(
          c("temporal_fold","region", "region_sub"), remove_from = "feature")
      
      #task$set_col_roles("sd_DM",roles = "weight")
    }
    
    rr = mlr3::resample(task = task, 
                        learner = lrn("regr.lm"), 
                        resampling = outer_resampling, 
                        store_models = T)
    
    # Extract inner tuning results
    predictions <- rr$prediction() %>% as.data.table() %>%
      left_join(data_folds, by="row_ids")
    
      predictions <- predictions %>% 
        mutate(truth = truth + outcome_mean,
               response = response + outcome_mean)  
    
    if (LOG){
      predictions <- predictions %>% mutate(truth=exp(truth),
                                            response=exp(response)) 
    }
    
    predictions <- predictions %>% 
      mutate(truth = truth * area,
             response = response * area,
             region=rrr, experiment = experiment)
    test <- bind_rows(test, predictions)
  }
}

saveRDS(test, paste0(output_path, "/linear_model_log_loocv.rds"))

# for (experiment in c("eco3", "eco2","grid","regional")){ 
#   print(paste("exp:", experiment))
#   if (experiment == "regional"){
#     data <- data_region
#   }else if (experiment == "eco2"){
#     data <- data_eco2
#   }else if (experiment == "eco3"){
#     data <- data_eco3
#   }else if (experiment == "grid"){
#     data <- data_grid 
#   }
#   
#   if (experiment != "regional"){
#   data <- data %>% select(-barren, -cropland, -forest, -grassland, -shrubland,
#                           -snowice, -urban, -water, -wetland, -lichen, -slope, -elevation)
#   }
#   
#   if (area_normalize){
#     data <- data %>% mutate(outcome = DM_kg / area)
#   }
#   else{
#     data <- data %>% mutate(outcome = DM_kg)
#   }
#   
#   for (rrr in region_list){
#     print(rrr)
#     data_train <- data %>% filter(region == rrr)
#     
#     vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
#                        "narr_vpd", "narr_runoff",
#                        "nldas_soilm")
#     
#     if (rrr %in% c("Canada-Alaska","Mexico")){
#       data_train <- data_train %>% select(-nldas_soilm) 
#       vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
#                          "narr_vpd", "narr_runoff")
#     }
#     
#     if (LOG){
#       data_train <- data_train %>% filter(outcome>0) %>% mutate(outcome=log(outcome))
#     }else{
#       data_train <- data_train 
#     }
# 
#     if (experiment =="grid"){
#       if (region_FE){
#         data_train <- data_train %>% 
#           mutate(gfed_cell_id = factor(gfed_cell_id)) %>%
#           group_by(gfed_cell_id) %>%
#           mutate(sd_DM=sd(DM_kg, na.rm=T)) %>%
#           filter(!is.na(sd_DM)) %>%
#           mutate(outcome_mean = mean(outcome, na.rm=T),
#                  outcome = outcome - outcome_mean) %>% 
#           mutate_at(vlist_climate, .funs = function(x){x-mean(x, na.rm=T)}) %>% ungroup()
#         
#         data_folds <- data_train %>% 
#           mutate(row_ids = row_number()) %>% 
#           select(row_ids, year, temporal_fold, region, gfed_cell_id, area, outcome_mean)
#       } else{
#         data_folds <- data_train %>% 
#           mutate(row_ids = row_number()) %>% 
#           select(row_ids, year, temporal_fold, region, gfed_cell_id, area)
#       }
#       task = as_task_regr(data_train, target = "outcome")
#       
#       task$set_col_roles(
#         "temporal_fold",
#         roles = "group")$set_col_roles(
#           c("year","region", "gfed_cell_id","area","DM_kg"),
#           remove_from = "feature")
#       
#       task$set_col_roles("sd_DM",roles = "weight")
#     }
#     if (experiment !="grid"){
#       if (region_FE){
#         data_train <- data_train %>%
#           mutate(region_sub = factor(region_sub)) %>%
#           group_by(region_sub) %>%
#           mutate(sd_DM=sd(DM_kg, na.rm=T)) %>%
#           filter(!is.na(sd_DM)) %>%
#           mutate(outcome_mean = mean(outcome, na.rm=T),
#                  outcome = outcome - outcome_mean) %>% 
#           mutate_at(vlist_climate, .funs = function(x){x-mean(x, na.rm=T)}) %>%
#           ungroup()
# 
#         data_folds <- data_train %>%
#           mutate(row_ids = row_number()) %>%
#           select(row_ids, year, temporal_fold, region, region_sub, area, outcome_mean)
#       } else{
#         data_folds <- data_train %>%
#           mutate(row_ids = row_number()) %>%
#           select(row_ids, year, temporal_fold, region, region_sub, area)
#       }
#       task = as_task_regr(data_train, target = "outcome")
# 
#       task$set_col_roles(
#         "temporal_fold", roles = "group")$set_col_roles(
#         c("year","region", "region_sub","area", "DM_kg"), remove_from = "feature")
#       
#       task$set_col_roles("sd_DM",roles = "weight")
#     }
#   
#     rr = mlr3::resample(task = task, 
#                         learner = lrn("regr.lm"), 
#                         resampling = outer_resampling, 
#                         store_models = T)
#     
#     # Extract inner tuning results
#     predictions <- rr$prediction() %>% as.data.table() %>%
#       left_join(data_folds, by="row_ids")
#     
#     if (region_FE){
#       predictions <- predictions %>% 
#         mutate(truth = truth + outcome_mean,
#                response = response + outcome_mean)  
#     }
#     
#     if (LOG){
#       predictions <- predictions %>% mutate(truth=exp(truth),
#                                             response=exp(response)) 
#     }
#     
#     predictions <- predictions %>% 
#       mutate(truth = truth * area,
#              response = response * area,
#              region=rrr, experiment = experiment)
#     test <- bind_rows(test, predictions)
#   }
# }
# 
# saveRDS(test, paste0(output_path, "/linear_model_log_FE_weighted.rds"))


### archived for FE
# if (experiment !="grid"){
#   if (region_FE){
#     data_train <- data_train %>% 
#       mutate(region_sub = factor(region_sub)) %>%
#       group_by(region_sub) %>%
#       mutate(outcome_mean = mean(outcome, na.rm=T),
#              outcome = outcome - outcome_mean) %>% ungroup()
#     
#     data_folds <- data_train %>% 
#       mutate(row_ids = row_number()) %>% 
#       select(row_ids, year, temporal_fold, region, region_sub, area, outcome_mean)
#   } else{
#     
#     data_folds <- data_train %>% 
#       mutate(row_ids = row_number()) %>% 
#       select(row_ids, year, temporal_fold, region, region_sub, area)
#   }
#   task = as_task_regr(data_train, target = "outcome")
#   
#   task$set_col_roles(
#     "temporal_fold", 
#     roles = "group")$set_col_roles(
#       c("year","region", "region_sub","area", "outcome_mean"), 
#       remove_from = "feature")  
# }

# if (experiment !="grid"){
#   if (region_FE){
#     data_train <- data_train %>% 
#       mutate(region_sub = factor(region_sub)) %>%
#       group_by(region_sub) %>%
#       mutate(outcome_mean = mean(outcome, na.rm=T),
#              outcome = outcome - outcome_mean) %>% ungroup()
#     
#     data_folds <- data_train %>% 
#       mutate(row_ids = row_number()) %>% 
#       select(row_ids, year, temporal_fold, region, region_sub, area, outcome_mean)
#   } else{
#     
#     data_folds <- data_train %>% 
#       mutate(row_ids = row_number()) %>% 
#       select(row_ids, year, temporal_fold, region, region_sub, area)
#   }
#   task = as_task_regr(data_train, target = "outcome")
#   
#   task$set_col_roles(
#     "temporal_fold", 
#     roles = "group")$set_col_roles(
#       c("year","region", "region_sub","area", "outcome_mean"), 
#       remove_from = "feature")  
# }

