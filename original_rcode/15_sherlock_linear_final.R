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
library(glmnet)

rm(list=ls())
gc()

input_path <- "~/oak_space/mhqiu/smoke_climate/input_data"
output_path <- "~/oak_space/mhqiu/smoke_climate/output"

nthread <- 1

#-------------------------------------------------------------------------------
# Get the final linear models
# Written by Minghao
# Last edited Mar 2022
#-------------------------------------------------------------------------------
region_list <- c("Western US", "Canada-Alaska", "Mexico", "Northeastern US", "Southeastern US") ##unique(gfed_data$region)
region_label_list <- c("west", "canada", "mexico", "northeast", "southeast")
LOG <- T
final_train <- T

##### Load training data #####
if (LOG){
  data_full <- readRDS(paste0(output_path, "/final_train_region_FE_log.rds"))
  data_grid <- readRDS(paste0(output_path, "/final_train_region_FE_log_1deg_grid.rds"))
} else{
  data_full <- readRDS(paste0(output_path, "/final_train_region_FE_level.rds"))
  data_grid <- readRDS(paste0(output_path, "/final_train_region_FE_level_1deg_grid.rds"))
}

###### if train the model ###### 
if (final_train){
  for (rrr in 1:length(region_list)){
    print(region_list[rrr])
    region_tmp <-   region_list[rrr] 
    region_label <- region_label_list[rrr]
    
    for (experiment in c("eco3", "eco2","regional", "grid")){ 
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
      
      data_train <- filter(data, region == region_tmp)
      
      if (region_tmp %in% c("Canada-Alaska","Mexico")){
        vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
                           "narr_vpd", "narr_runoff")
              } 
      
      if (!region_tmp %in% c("Canada-Alaska","Mexico")){        
        vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
                           "narr_vpd", "narr_runoff", "nldas_soilm")
      }
      
     lm_fml  <-  paste("outcome_demean ~ ", paste(vlist_climate, collapse = " + "))
     model_tmp <- lm(as.formula(lm_fml), data= data_train,
                     weights = data_train$sd_DM)
        
      if (experiment == "regional"){
        linear_region <- model_tmp
      }else if (experiment == "eco2"){
        linear_eco2 <- model_tmp
      }else if (experiment == "eco3"){
        linear_eco3 <- model_tmp
      }else if (experiment == "grid"){
        linear_grid <- model_tmp 
      }
    }
    final_models <- list(linear_region = linear_region, 
                         linear_eco2 = linear_eco2, 
                         linear_eco3 = linear_eco3,
                         linear_grid = linear_grid)
    
    saveRDS(final_models, file = paste0(output_path, "/",region_label,"_linear_log_weighted_final_models.rds"))
  }
}
