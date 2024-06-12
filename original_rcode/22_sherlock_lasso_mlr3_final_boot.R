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

nthread <- 7

#-------------------------------------------------------------------------------
# Use MLR3 linear models to predict GFED DM emissions at the regional levels
# Written by Minghao
# Last edited Mar 2022
#-------------------------------------------------------------------------------
# Set up measure
measure = msr("regr.mse")
# Set up terminator
terminator = trm("none")
# Set up tuner
tuner = tnr("grid_search")

#rrr <- "Western US"
#region_label <- "west"
region_list <- c("Western US", "Canada-Alaska", "Mexico", "Northeastern US", "Southeastern US") ##unique(gfed_data$region)
region_label_list <- c("west", "canada", "mexico", "northeast", "southeast")
MODEL <- "LASSO"
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

nboot <- 40
final_models <- vector(mode='list', length=4*nboot*length(region_list))

###### If train the model ###### 
if (final_train){
  index <- 1
  for (rrr in 1:length(region_list)){
    print(region_list[rrr])
    region_tmp <-   region_list[rrr] 
    region_label <- region_label_list[rrr]
    ### If MODEL=="LASSO", the learners are set in the training (as need to get the range of lambdas ex-ante)
    for (experiment in c("regional","eco2","eco3","grid")){ 
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
        if (experiment != "regional"){
          var_list <- c(vlist_climate, paste(rep(vlist_climate, each = 3),
                                             c("forest","crop", "grass"), sep="_"),
                        paste0(vlist_climate, "_temp"),
                        paste0(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff"), "_precip"),
                        paste0(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff"), "_rhum"),
                        paste0(c("narr_wspd",  "narr_vpd", "narr_runoff"), "_wspd"),
                        paste0(c("narr_vpd", "narr_runoff"), "_vpd"),
                        paste0(c("narr_runoff"), "_runoff"))
        } else {
          var_list <- c(vlist_climate,
                        paste0(vlist_climate, "_temp"),
                        paste0(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff"), "_precip"),
                        paste0(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff"), "_rhum"),
                        paste0(c("narr_wspd",  "narr_vpd", "narr_runoff"), "_wspd"),
                        paste0(c("narr_vpd", "narr_runoff"), "_vpd"),
                        paste0(c("narr_runoff"), "_runoff"))
        }
      } 
      
      if (!region_tmp %in% c("Canada-Alaska","Mexico")){        
        vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
                           "narr_vpd", "narr_runoff", "nldas_soilm")
        if (experiment != "regional"){
          var_list <- c(vlist_climate, paste(rep(vlist_climate, each = 3),
                                             c("forest","crop", "grass"), sep="_"),
                        paste0(vlist_climate, "_temp"),
                        paste0(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff", "nldas_soilm"), "_precip"),
                        paste0(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_rhum"),
                        paste0(c("narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_wspd"),
                        paste0(c("narr_vpd", "narr_runoff", "nldas_soilm"), "_vpd"),
                        paste0(c("narr_runoff", "nldas_soilm"), "_runoff"),
                        paste0(c("nldas_soilm"), "_soilm"))
        } else {
          var_list <- c(vlist_climate,
                        paste0(vlist_climate, "_temp"),
                        paste0(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff", "nldas_soilm"), "_precip"),
                        paste0(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_rhum"),
                        paste0(c("narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_wspd"),
                        paste0(c("narr_vpd", "narr_runoff", "nldas_soilm"), "_vpd"),
                        paste0(c("narr_runoff", "nldas_soilm"), "_runoff"),
                        paste0(c("nldas_soilm"), "_soilm"))
        }
      }
      
  for (bbb in 1:nboot){
    print(index)
    data_train <- filter(data, region == region_tmp)
    data_train <- data_train[sample(1:nrow(data_train), replace = T), ]
    
    n_fold <- length(unique(data_train$year))
    
      if (experiment =="grid"){
        data_folds <- data_train %>% 
          mutate(row_ids = row_number()) %>% 
          select(row_ids, year, temporal_fold, region, grid_1deg_id, area, outcome_mean)
        
        data_train <- data_train[,c("year", "temporal_fold", "region", "grid_1deg_id", 
                                    "outcome_demean", "sd_DM",
                                    var_list)]
        
        if (MODEL =="LASSO"){
          x <- data.matrix(data_train[, var_list])
          cv_model <- cv.glmnet(x, data_train$outcome_demean, alpha = 1, nfolds =n_fold)
          
          full_learner = lrn("regr.glmnet",
                             lambda = to_tune(c(
                               c(1000,100,50,10,5,2)*max(cv_model$lambda),
                               cv_model$lambda,
                               c(0.5,0.2,0.1,0.02,0.01,0.001)*min(cv_model$lambda))),
                             predict_sets=c("train","test")) #, 200, 400
        }
        
        task = as_task_regr(data_train, target = "outcome_demean")
        
        task$set_col_roles(
          "year", 
          roles = "group")$set_col_roles(
            c("temporal_fold","region", "grid_1deg_id"), 
            remove_from = "feature") 
        
        task$set_col_roles("sd_DM",roles = "weight")
      }
      
      if (experiment !="grid"){
        data_folds <- data_train %>%
          mutate(row_ids = row_number()) %>%
          select(row_ids, year, temporal_fold, region, region_sub, area, outcome_mean)
        
        data_train <- data_train[,c("year", "temporal_fold", "region", "region_sub", 
                                    "outcome_demean", "sd_DM",
                                    var_list)]
        
        if (MODEL =="LASSO"){
          x <- data.matrix(data_train[, var_list])
          cv_model <- cv.glmnet(x, data_train$outcome_demean, alpha = 1, nfolds =n_fold)
          
          full_learner = lrn("regr.glmnet",
                             lambda = to_tune(c(
                               c(1000,100,50,10,5,2)*max(cv_model$lambda),
                               cv_model$lambda,
                               c(0.5,0.2,0.1,0.02,0.01,0.001)*min(cv_model$lambda))),
                             predict_sets=c("train","test")) #, 200, 400
        }
        task = as_task_regr(data_train, target = "outcome_demean")
        
        task$set_col_roles(
          "year", 
          roles = "group")$set_col_roles(
            c("temporal_fold","region", "region_sub"), 
            remove_from = "feature")  
        
        task$set_col_roles("sd_DM",roles = "weight")
      }
      
      #### the following step is to get a potential range for the lambda values for tuning
      #perform k-fold cross-validation to find optimal lambda value
      set_threads(full_learner, nthread)
      
      outer_resampling = rsmp("cv", folds = n_fold)
      instance <- ti(task = task,
                     learner = full_learner,
                     resampling = outer_resampling,
                     measure = measure,
                     terminator = terminator,
                     store_models = T)
      
      tr <- tuner$optimize(instance)
      
      # Train final model
      tuned <- lrn("regr.glmnet",
                   lambda = 1,
                   predict_sets=c("test"))## these are just initial values, will be replaced by optimal hyperparameters
      
      tuned$param_set$values <- instance$result_learner_param_vals
      tuned$train(task)
      
    final_models[[index]] <- tuned
    names(final_models)[index] <- paste(region_label, experiment, bbb, sep="_")
    index <- index + 1
     }
    }
  }
    saveRDS(final_models, file = paste0(output_path, "/LASSO_log_weighted_final_models_boot.rds"))
  }

