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

nthread <- 8

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

## Set up resampling folds
# num_temporal_folds = length(unique(gfed_data$temporal_fold))
# outer_resampling = rsmp("cv", folds = num_temporal_folds)
# inner_resampling = rsmp("cv", folds = num_temporal_folds-1)

LOOCV <- T
if (LOOCV){
  #num_temporal_folds = length(unique(gfed_data$temporal_fold))
  outer_resampling = rsmp("cv", folds = 21)
  inner_resampling = rsmp("cv", folds = 5)
}
# Set up measure
measure = msr("regr.mse")
# Set up terminator
terminator = trm("none")
# Set up tuner
tuner = tnr("grid_search")

region_list <- unique(gfed_data$region)
test <- NULL 
train <- NULL 
tuning_combined <- NULL
LOG <- T
nestedCV <- T
final_train <- F

##### Load training data #####
if (LOG){
  data_full <- readRDS(paste0(output_path, "/final_train_region_FE_log.rds"))
  data_grid <- readRDS(paste0(output_path, "/final_train_region_FE_log_1deg_grid.rds"))
} else{
  data_full <- readRDS(paste0(output_path, "/final_train_region_FE_level.rds"))
  data_grid <- readRDS(paste0(output_path, "/final_train_region_FE_level_1deg_grid.rds"))
}

MODEL <- "LASSO"# "LASSO"
## set up learners 
if (MODEL =="XGBOOST"){
  full_learner = lrn("regr.xgboost",
                     eta = to_tune(c(0.01, 0.1, 0.4)),  #0.01, 0.1, 0.4
                     max_depth = to_tune(c(2, 4, 8)), #2, 4, 8
                     nrounds = to_tune(c(50, 100, 200, 400)), #50, 100, 200, 400
                     subsample = to_tune(c(0.2, 0.5, 0.8)), #0.2, 0.5, 0.8
                     lambda   = to_tune(c(3, 5, 10)), #3, 5, 10
                     predict_sets=c("train","test")) #, 200, 400             
} else if (MODEL =="NN"){
  full_learner = lrn("regr.nnet",
                     size = to_tune(c(5, 10, 20, 50, 100)),  #5, 10, 20, 50, 100
                     maxit = to_tune(c(10, 20, 100, 200, 400)), #10, 20, 50, 100, 200, 400
                     MaxNWts = 50000,
                     predict_sets=c("train","test")) #, 200, 400             
}

###### If using nested CV to evaluate the model ###### 
if (nestedCV){
### If MODEL=="LASSO", the learners are set in the training (as need to get the range of lambdas ex-ante)
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
        
if (rrr %in% c("Canada-Alaska","Mexico")){
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
        
if (!rrr %in% c("Canada-Alaska","Mexico")){        
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

    if (experiment =="grid"){
        data_folds <- data_train %>% 
          mutate(row_ids = row_number()) %>% 
          select(row_ids, year, temporal_fold, region, grid_1deg_id, area, outcome_mean)

        data_train <- data_train[,c("year", "temporal_fold", "region", "grid_1deg_id", 
                                    "outcome_demean", "sd_DM",
                                    var_list)]
        
      if (MODEL =="LASSO"){
      x <- data.matrix(data_train[, var_list])
      cv_model <- cv.glmnet(x, data_train$outcome_demean, alpha = 1, nfolds =21)
      
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
      cv_model <- cv.glmnet(x, data_train$outcome_demean, alpha = 1, nfolds =21)
      
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
      
    at = auto_tuner(method = tuner,
                      learner = full_learner,
                      resampling = inner_resampling,
                      measure = measure,
                      terminator = terminator,
                      store_models = T)
    
    # Run nested resampling
    rr = mlr3::resample(task = task,
                        learner = at,
                        resampling = outer_resampling,
                        store_models = T)
    
    archives = extract_inner_tuning_archives(rr)
    hyperparam <- data.frame(archives[,1:3]) %>% #nn:4 xgboost:7
      mutate(region=rrr, experiment = experiment)
    tuning_combined <- bind_rows(tuning_combined, hyperparam)
    
    # Extract inner tuning results
    predictions <- rr$prediction() %>% as.data.table() %>%
      left_join(data_folds, by="row_ids")

      predictions <- predictions %>% 
        mutate(truth = truth + outcome_mean,
               response = response + outcome_mean)  
    
    if (LOG){
      predictions <- predictions %>% mutate(truth=exp(truth),
                                            response=exp(response),
                                            model="log") 
    }
    
    predictions <- predictions %>% 
      mutate(truth = truth * area,
             response = response * area,
             region=rrr, experiment = experiment)
    test <- bind_rows(test, predictions)
    
    train_tmp <- rr$prediction(predict_sets="train") %>% 
      as.data.table() %>%
      left_join(data_folds, by="row_ids")
    
    train_tmp <- train_tmp %>% 
      mutate(truth = truth + outcome_mean,
             response = response + outcome_mean)  
    
    if (LOG){
      train_tmp <- train_tmp %>% mutate(truth=exp(truth),
                                        response=exp(response),
                                        model="log") 
    }
    
    train_tmp <- train_tmp %>% 
      mutate(truth = truth * area,
             response = response * area,
             region=rrr, experiment = experiment)
    train <- bind_rows(train, train_tmp)
    
    }
 }

saveRDS(test, paste0(output_path, "/LASSO_full_log_loocv.rds"))
saveRDS(train, paste0(output_path, "/LASSO_full_log_loocv_train.rds"))
saveRDS(tuning_combined, paste0(output_path, "/LASSO_full_log_loocv_tuning.rds"))
}

# saveRDS(test, paste0(output_path, "/",region_label,"_",MODEL,"_full_log_FE_weighted.rds"))
# saveRDS(tuning_combined, paste0(output_path, "/",region_label,"_",MODEL,"_full_log_FE_weighted_tuning.rds"))


#          ARCHIVE scripts
# 
# vlist <- c("barren", "cropland", "forest", "grassland", "lichen", "shrubland", "snowice", "urban", "water","wetland",
#            "elevation", "slope", 
#            "narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
#            "narr_vpd", "narr_runoff",
#            "nldas_soilm") 
# 
# vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
#                    "narr_vpd", "narr_runoff",
#                    "nldas_soilm")
# 
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
# 
# ## Set up resampling folds
# num_temporal_folds = length(unique(gfed_data$temporal_fold))
# outer_resampling = rsmp("cv", folds = num_temporal_folds)
# inner_resampling = rsmp("cv", folds = num_temporal_folds-1)
# # Set up measure
# measure = msr("regr.mse")
# # Set up terminator
# terminator = trm("none")
# # Set up tuner
# tuner = tnr("grid_search")
# 
# #region_list <- unique(gfed_data$region)
# rrr <- "Western US"
# region_label <- "west"
# test <- NULL 
# tuning_combined <- NULL
# LOG <- T
# area_normalize <- T
# region_FE <- T
# 
# MODEL <- "NN"# "LASSO"
# ## set up learners 
# if (MODEL =="XGBOOST"){
#   full_learner = lrn("regr.xgboost",
#                      eta = to_tune(c(0.01, 0.1, 0.4)),  #0.01, 0.1, 0.4
#                      max_depth = to_tune(c(2, 4, 8)), #2, 4, 8
#                      nrounds = to_tune(c(50, 100, 200, 400)), #50, 100, 200, 400
#                      subsample = to_tune(c(0.2, 0.5, 0.8)), #0.2, 0.5, 0.8
#                      lambda   = to_tune(c(3, 5, 10)), #3, 5, 10
#                      predict_sets=c("train","test")) #, 200, 400             
# } else if (MODEL =="NN"){
#   full_learner = lrn("regr.nnet",
#                      size = to_tune(c(5, 10, 20, 50, 100)),  #5, 10, 20, 50, 100
#                      maxit = to_tune(c(10, 20, 100, 200, 400)), #10, 20, 50, 100, 200, 400
#                      MaxNWts = 50000,
#                      predict_sets=c("train","test")) #, 200, 400             
# }
# set_threads(full_learner, nthread)
# 
# ### If MODEL=="LASSO", the learners are set in the training (as need to get the range of lambdas ex-ante)
# 
# for (experiment in c("eco3", "eco2","regional")){ 
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
#   if (area_normalize){
#     data <- data %>% mutate(outcome = DM_kg / area)
#   }
#   else{
#     data <- data %>% mutate(outcome = DM_kg)
#   }
#   
#   #for (rrr in region_list){
#   # print(rrr)
#   data_train <- data %>% filter(region == rrr) 
#   
#   vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
#                      "narr_vpd", "narr_runoff",
#                      "nldas_soilm")
#   
#   if (rrr %in% c("Canada-Alaska","Mexico")){
#     data_train <- data_train %>% select(-nldas_soilm) 
#     vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
#                        "narr_vpd", "narr_runoff")
#   }
#   
#   if (LOG){
#     data_train <- data_train %>% filter(outcome>0) %>% mutate(outcome=log(outcome))
#   }else{
#     data_train <- data_train 
#   }
#   
#   ### include more interactive features
#   if (experiment != "regional"){
#     data_train <- data_train %>%
#       mutate_at(vlist_climate, list(forest = function(x){x*data_train$forest},
#                                     crop = function(x){x*data_train$cropland},
#                                     grass = function(x){x*data_train$grassland})) %>%
#       mutate_at(vlist_climate, list(temp = function(x){x*data_train$narr_temp})) %>%
#       mutate_at(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff", "nldas_soilm"), 
#                 list(precip = function(x){x*data_train$narr_precip})) %>%
#       mutate_at(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"),
#                 list(rhum = function(x){x*data_train$narr_rhum})) %>%  
#       mutate_at(c("narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"),
#                 list(wspd = function(x){x*data_train$narr_wspd})) %>% 
#       mutate_at(c("narr_vpd", "narr_runoff", "nldas_soilm"),
#                 list(vpd = function(x){x*data_train$narr_vpd})) %>% 
#       mutate_at(c("narr_runoff", "nldas_soilm"),
#                 list(runoff = function(x){x*data_train$narr_runoff})) %>% 
#       mutate(nldas_soilm_soilm = nldas_soilm^2) %>%                              
#       select(-barren, -cropland, -forest, -grassland, -shrubland,
#              -snowice, -urban, -water, -wetland, -lichen, -slope, -elevation)
#     
#     var_list <- c(vlist_climate, paste(rep(vlist_climate, each = 3),
#                                        c("forest","crop", "grass"), sep="_"),
#                   paste0(vlist_climate, "_temp"),
#                   paste0(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff", "nldas_soilm"), "_precip"),
#                   paste0(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_rhum"),
#                   paste0(c("narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_wspd"),
#                   paste0(c("narr_vpd", "narr_runoff", "nldas_soilm"), "_vpd"),
#                   paste0(c("narr_runoff", "nldas_soilm"), "_runoff"),
#                   paste0(c("nldas_soilm"), "_soilm"))
#     
#   } else{
#     data_train <- data_train %>%
#       mutate_at(vlist_climate, list(temp = function(x){x*data_train$narr_temp})) %>%
#       mutate_at(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff", "nldas_soilm"), 
#                 list(precip = function(x){x*data_train$narr_precip})) %>%
#       mutate_at(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"),
#                 list(rhum = function(x){x*data_train$narr_rhum})) %>%  
#       mutate_at(c("narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"),
#                 list(wspd = function(x){x*data_train$narr_wspd})) %>% 
#       mutate_at(c("narr_vpd", "narr_runoff", "nldas_soilm"),
#                 list(vpd = function(x){x*data_train$narr_vpd})) %>% 
#       mutate_at(c("narr_runoff", "nldas_soilm"),
#                 list(runoff = function(x){x*data_train$narr_runoff})) %>% 
#       mutate(nldas_soilm_soilm = nldas_soilm^2)
#     
#     var_list <- c(vlist_climate,
#                   paste0(vlist_climate, "_temp"),
#                   paste0(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff", "nldas_soilm"), "_precip"),
#                   paste0(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_rhum"),
#                   paste0(c("narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_wspd"),
#                   paste0(c("narr_vpd", "narr_runoff", "nldas_soilm"), "_vpd"),
#                   paste0(c("narr_runoff", "nldas_soilm"), "_runoff"),
#                   paste0(c("nldas_soilm"), "_soilm"))
#   }
#   
#   if (experiment =="grid"){
#     if (region_FE){
#       data_train <- data_train %>% 
#         mutate(gfed_cell_id = factor(gfed_cell_id)) %>%
#         group_by(gfed_cell_id) %>%
#         mutate(sd_DM=sd(DM_kg, na.rm=T)) %>% 
#         filter(!is.na(sd_DM)) %>%
#         mutate(outcome_mean = mean(outcome, na.rm=T),
#                outcome = outcome - outcome_mean) %>% 
#         mutate_at(var_list, .funs = function(x){x-mean(x, na.rm=T)}) %>% ungroup()
#       
#       data_folds <- data_train %>% 
#         mutate(row_ids = row_number()) %>% 
#         select(row_ids, year, temporal_fold, region, gfed_cell_id, area, outcome_mean)
#     } else{
#       data_folds <- data_train %>% 
#         mutate(row_ids = row_number()) %>% 
#         select(row_ids, year, temporal_fold, region, gfed_cell_id, area)
#     }
#     
#     if (MODEL =="LASSO"){
#       x <- data.matrix(data_train[, var_list])
#       cv_model <- cv.glmnet(x, data_train$outcome, alpha = 1, nfolds =7)
#       
#       full_learner = lrn("regr.glmnet",
#                          lambda = to_tune(c(
#                            c(1000,100,50,10,5,2)*max(cv_model$lambda),
#                            cv_model$lambda,
#                            c(0.5,0.2,0.1,0.02,0.01,0.001)*min(cv_model$lambda))),
#                          predict_sets=c("train","test")) #, 200, 400
#       set_threads(full_learner, nthread)
#     }
#     
#     task = as_task_regr(data_train, target = "outcome")
#     
#     task$set_col_roles(
#       "temporal_fold", 
#       roles = "group")$set_col_roles(
#         c("year","region", "gfed_cell_id","area", "outcome_mean", "lat", "lon", "DM_kg"), 
#         remove_from = "feature") 
#     
#     task$set_col_roles("sd_DM",roles = "weight")
#   }
#   
#   if (experiment !="grid"){
#     if (region_FE){
#       data_train <- data_train %>%
#         mutate(region_sub = factor(region_sub)) %>%
#         group_by(region_sub) %>%
#         mutate(sd_DM=sd(DM_kg, na.rm=T)) %>% 
#         filter(!is.na(sd_DM)) %>%
#         mutate(outcome_mean = mean(outcome, na.rm=T),
#                outcome = outcome - outcome_mean) %>% 
#         mutate_at(var_list, .funs = function(x){x-mean(x, na.rm=T)}) %>%
#         ungroup()
#       
#       data_folds <- data_train %>%
#         mutate(row_ids = row_number()) %>%
#         select(row_ids, year, temporal_fold, region, region_sub, area, outcome_mean)
#     } else{
#       data_folds <- data_train %>%
#         mutate(row_ids = row_number()) %>%
#         select(row_ids, year, temporal_fold, region, region_sub, area)
#     }
#     
#     if (MODEL =="LASSO"){
#       x <- data.matrix(data_train[, var_list])
#       cv_model <- cv.glmnet(x, data_train$outcome, alpha = 1, nfolds =7)
#       
#       full_learner = lrn("regr.glmnet",
#                          lambda = to_tune(c(
#                            c(1000,100,50,10,5,2)*max(cv_model$lambda),
#                            cv_model$lambda,
#                            c(0.5,0.2,0.1,0.02,0.01,0.001)*min(cv_model$lambda))),
#                          predict_sets=c("train","test")) #, 200, 400
#       
#       set_threads(full_learner, nthread)
#     }
#     task = as_task_regr(data_train, target = "outcome")
#     
#     task$set_col_roles(
#       "temporal_fold", 
#       roles = "group")$set_col_roles(
#         c("year","region", "region_sub","area", "outcome_mean","DM_kg"), 
#         remove_from = "feature")  
#     
#     task$set_col_roles("sd_DM",roles = "weight")
#   }
#   
#   #### the following step is to get a potential range for the lambda values for tuning
#   #perform k-fold cross-validation to find optimal lambda value
#   at = auto_tuner(method = tuner,
#                   learner = full_learner,
#                   resampling = inner_resampling,
#                   measure = measure,
#                   terminator = terminator,
#                   store_models = T)
#   
#   # Run nested resampling
#   rr = mlr3::resample(task = task,
#                       learner = at,
#                       resampling = outer_resampling,
#                       store_models = T)
#   
#   archives = extract_inner_tuning_archives(rr)
#   hyperparam <- data.frame(archives[,1:4]) %>% #nn:4 xgboost:7
#     mutate(region=rrr, experiment = experiment)
#   tuning_combined <- bind_rows(tuning_combined, hyperparam)
#   
#   # Extract inner tuning results
#   predictions <- rr$prediction() %>% as.data.table() %>%
#     left_join(data_folds, by="row_ids")
#   
#   if (region_FE){
#     predictions <- predictions %>% 
#       mutate(truth = truth + outcome_mean,
#              response = response + outcome_mean)  
#   }
#   
#   if (LOG){
#     predictions <- predictions %>% mutate(truth=exp(truth),
#                                           response=exp(response),
#                                           model="log") 
#   }else{
#     predictions$model <- "level" 
#   }
#   
#   predictions <- predictions %>% 
#     mutate(truth = truth * area,
#            response = response * area,
#            region=rrr, experiment = experiment)
#   test <- bind_rows(test, predictions)
# }
# 
# 
# saveRDS(test, paste0(output_path, "/",region_label,"_",MODEL,"_full_log_FE_weighted.rds"))
# saveRDS(tuning_combined, paste0(output_path, "/",region_label,"_",MODEL,"_full_log_FE_weighted_tuning.rds"))
