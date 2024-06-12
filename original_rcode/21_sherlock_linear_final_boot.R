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
# Get the final linear models and use bootstrap to get uncertainty
# Written by Minghao
# Last edited Jan 2024
#-------------------------------------------------------------------------------
region_list <- c("Western US", "Canada-Alaska", "Mexico", "Northeastern US", "Southeastern US") ##unique(gfed_data$region)
region_label_list <- c("west", "canada", "mexico", "northeast", "southeast")
LOG <- F
final_train <- T
nboot <- 40

##### Load training data #####
if (LOG){
  data_full <- readRDS(paste0(output_path, "/final_train_region_FE_log.rds"))
  data_grid <- readRDS(paste0(output_path, "/final_train_region_FE_log_1deg_grid.rds"))
} else{
  data_full <- readRDS(paste0(output_path, "/final_train_region_FE_level.rds"))
  data_grid <- readRDS(paste0(output_path, "/final_train_region_FE_level_1deg_grid.rds"))
}

###### if train the model ###### 
final_models <- vector(mode='list', length=4*nboot*length(region_list))

if (final_train){
  index <- 1
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
      
      if (region_tmp %in% c("Canada-Alaska","Mexico")){
        vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
                           "narr_vpd", "narr_runoff")
      } 
      
      if (!region_tmp %in% c("Canada-Alaska","Mexico")){        
        vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
                           "narr_vpd", "narr_runoff", "nldas_soilm")
      }
      
  for (bbb in 1:nboot){
        #print(bbb)
        print(index)
        set.seed(bbb)
        
      data_train <- filter(data, region == region_tmp)
      data_train <- data_train[sample(1:nrow(data_train), replace = T), ]
      
      lm_fml  <-  paste("outcome_demean ~ ", paste(vlist_climate, collapse = " + "))
      model_tmp <- lm(as.formula(lm_fml), data= data_train)
      
    final_models[[index]] <- model_tmp
    names(final_models)[index] <- paste(region_label, experiment, bbb, sep="_")
    index <- index + 1
  }
  }
  }
    saveRDS(final_models, file = paste0(output_path, "/linear_level_final_models_boot.rds"))
}
