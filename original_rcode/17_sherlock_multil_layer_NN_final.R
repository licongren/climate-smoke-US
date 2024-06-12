library(terra);library(ncdf4)
library(sf)
library(dplyr)   
library(tidyverse)
library(purrr)
library(data.table)
library(mlr3verse)
library(future)
library(caret)
library(neuralnet)
library(RSNNS)
library(doSNOW)
library(openxlsx)
library(h2o)
library(RSNNS)
library(doSNOW)

rm(list=ls())
gc()

input_path <- "~/oak_space/mhqiu/smoke_climate/input_data"
output_path <- "~/oak_space/mhqiu/smoke_climate/output"

nthread <- 1

#-------------------------------------------------------------------------------
# Use MLR3 linear models to predict GFED DM emissions at the regional levels
# Written by Minghao
# Last edited Mar 2022
#-------------------------------------------------------------------------------
best_models <- read.xlsx(paste0(output_path, "/best_metrics_rmse_bias_6models_LOOCV.xlsx"))
best_models1 <- read.xlsx(paste0(output_path, "/best_metrics_rmse_bias_6models_subregion_LOOCV.xlsx"))

best_models <- bind_rows(best_models, best_models1) %>% 
  filter(model=="Neural Net", selected==1, metric %in% c("5 years", "10 years")) %>%
  distinct(region, model, experiment, outcome, spec, .keep_all = T)

## Set up resampling folds
# final_train <- T
# projection <- T
nnet_cmip_proj <- NULL
nnet_history <- NULL
feature_normalize <- T
region_list <- c("Western US", "Canada-Alaska", "Mexico", "Southeastern US", "Northeastern US") ##unique(gfed_data$region)
region_label_list <- c("west", "canada", "mexico", "southeast", "northeast")

#--------------------------------------
###       Using H2O
#--------------------------------------
set.seed(111)

h2o.init(max_mem_size = "100g",
         nthreads = nthread)

#if (final_train){
  final_model_list <- vector(mode='list', length=nrow(best_models))
  data_norm_list <- vector(mode='list', length=nrow(best_models))
  
  for (iii in 1:nrow(best_models)){
  
    experiment <- best_models[iii,"experiment"]
    spec <-  best_models[iii,"spec"]
    outcome <- best_models[iii,"outcome"]
    region_tmp <-  best_models[iii,"region"]
    
    print(paste(region_tmp, experiment, outcome, spec, sep="_"))
    
   ## determine the outcome variable
      if (outcome == "log"){
        data_full <- readRDS(paste0(output_path, "/final_train_region_FE_log.rds"))
        data_grid <- readRDS(paste0(output_path, "/final_train_region_FE_log_1deg_grid.rds"))
      } else{
        data_full <- readRDS(paste0(output_path, "/final_train_region_FE_level.rds"))
        data_grid <- readRDS(paste0(output_path, "/final_train_region_FE_level_1deg_grid.rds"))
      }

      ## determine spatial resolution of the model
    if (experiment == "regional"){
      data <- data_full$data_region
    }else if (experiment == "eco2"){
      data <- data_full$data_eco2
    }else if (experiment == "eco3"){
      data <- data_full$data_eco3
    }else if (experiment == "grid"){
      data <- data_grid 
    }
    
      data_train <- data %>% filter(region == region_tmp) %>%
        mutate(row_ids = row_number(),
               year = factor(year))
      
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
      
      if (spec == "simple"){
        var_list <- vlist_climate
      }
  
      if (feature_normalize){
        ###### compute the mean and sd of each variable ###### 
        data_norm_mean <-  data_train %>% 
          summarise_at(var_list, .funs = mean, na.rm=T) %>% t() %>% as.data.frame()
        data_norm_mean <- data_norm_mean %>% mutate(variable= row.names(data_norm_mean),stat ="mean") %>%
          rename(value = V1)
        data_norm_sd <-  data_train %>% 
          summarise_at(var_list, .funs = sd, na.rm=T) %>% t() %>% as.data.frame()
        data_norm_sd <- data_norm_sd %>% mutate(variable= row.names(data_norm_sd),stat ="sd") %>%
          rename(value = V1)
        
        data_norm <- bind_rows(data_norm_mean, data_norm_sd)
        
        data_train <- data_train %>% 
          mutate_at(var_list, .fun=function(x){return((x-mean(x, na.rm=T))/sd(x, na.rm=T))})
      }
    
      ##### set up training, tune the hyper-param with outfolds 
      data_train.hex <- as.h2o(data_train)
      
      grid_id_tmp <- paste(experiment, region_tmp, outcome, spec, sep="_")
      
      ##### set up the number pf hidden neurons for tuning
      if (experiment == "grid"){
      tmp <- expand.grid(c(10, 25, 50),
                         c(10, 25, 50)) %>% as.matrix()
      hidden_neurons <- vector(mode='list', length=nrow(tmp))
      for (kk in  1:nrow(tmp)){
        hidden_neurons[[kk]] <- as.vector(tmp[kk,])
      }
      
      hyper_params <- list(
        activation = c("Rectifier"), 
        hidden = hidden_neurons,
        epochs = c(50, 200),
        rate = c(0.01, 0.001))
      }
      
      if (experiment!="grid"){
        tmp <- expand.grid(c(2,8,16,32),
                           c(2,8,16,32)) %>% as.matrix()
        hidden_neurons <- vector(mode='list', length=nrow(tmp))
        for (kk in  1:nrow(tmp)){
          hidden_neurons[[kk]] <- as.vector(tmp[kk,])
        }
        
        hyper_params <- list(
          activation = c("Rectifier"), 
          hidden = hidden_neurons,
          epochs = c(50, 100, 200),
          rate = c(0.01, 0.005, 0.001))
      }
      
      search_criteria <- list(strategy = "RandomDiscrete",
                              stopping_metric="RMSE",  
                              stopping_tolerance = 0.01,
                              stopping_rounds = 3,
                              seed = 100)
      
      dl_grid <- h2o.grid(algorithm = "deeplearning", 
                          x = var_list,
                          y = "outcome_demean",
                          weights_column = "sd_DM",
                          grid_id = grid_id_tmp,
                          training_frame = data_train.hex,
                          #validation_frame = data_test.hex,
                          fold_column = "year",
                          hyper_params = hyper_params,
                          regression_stop=1e-2,
                          stopping_metric="RMSE",      ## logloss is directly optimized by Deep Learning
                          stopping_tolerance=1e-2,        ## stop when validation logloss does not improve by >=1% for 2 scoring events
                          stopping_rounds=5,
                          search_criteria = search_criteria
      )
      
      dl_grid_performance <- h2o.getGrid(grid_id_tmp, sort_by = "rmse")
      
      final_model <- h2o.getModel(dl_grid_performance@model_ids[[1]])
      
      ##### using the final model to project the fire emission level under CMIP6 scenarios
   #   if (projection){
        proj_path <- "~/oak_space/mhqiu/smoke_climate/fire_ML_model/final_models/proj_data"
        
        if (outcome == "level"){
          history_data_full <- readRDS(paste0(output_path, "/final_train_region_FE_level.rds"))
          history_data_grid <- readRDS(paste0(output_path, "/final_train_region_FE_level_1deg_grid.rds"))
          
          if (experiment == "regional"){
            cmip_data <- readRDS(paste0(proj_path, "/cmip_region_full_demean_level.rds")) %>%
              mutate(region_sub = region)
            history_data <- history_data_full$data_region
          }else if (experiment == "eco2"){
            cmip_data <- readRDS(paste0(proj_path, "/cmip_eco2_full_demean_level.rds"))
            history_data <- history_data_full$data_eco2
          }else if (experiment == "eco3"){
            cmip_data <- readRDS(paste0(proj_path, "/cmip_eco3_full_demean_level.rds"))
            history_data <- history_data_full$data_eco3
          }else if (experiment == "grid"){
            cmip_data <- readRDS(paste0(proj_path, "/cmip_region_1deg_full_demean_level.rds")) %>%
              rename(region_sub = grid_1deg_id) %>%
              mutate(region_sub = as.character(region_sub)) %>%
              filter(scenario != "ssp585" & year<2071) ## save memory
            
            history_data <- history_data_grid %>%
              rename(region_sub = grid_1deg_id) %>%
              mutate(region_sub = as.character(region_sub))
          }
          
          #### normalize the features in CMIP data
          for (vv in unique(data_norm$variable)){
            mean_tmp <- filter(data_norm, variable == vv, stat=="mean")$value
            sd_tmp <- filter(data_norm, variable == vv, stat=="sd")$value
            cmip_data[,vv] <- (cmip_data[,vv] - mean_tmp)/sd_tmp
            history_data[,vv] <- (history_data[,vv] - mean_tmp)/sd_tmp
          }
          
          cmip_proj <- cmip_data %>%
            filter(region == region_tmp)  %>%
            filter_at(vlist_climate, all_vars(!is.na(.))) %>%
            mutate(year=as.factor(year))
          
          cmip_proj.hex <- as.h2o(cmip_proj)
          cmip_proj$pred_DM <-  h2o.predict(object = final_model, newdata = cmip_proj.hex) %>% as.vector()
             
          cmip_proj <- cmip_proj %>% mutate(pred_DM = (pred_DM + outcome_mean_level) * area)  %>%
            select(model, scenario, year, region, region_sub, pred_DM,
                   narr_temp, narr_precip, narr_rhum, narr_vpd, narr_runoff, narr_wspd, nldas_soilm) %>%
            mutate(algorithm = "Neural Net", spatial_res=experiment, outcome=outcome, spec=spec)
          
          ### use the real observations to predict 
          history_data <- history_data %>% filter(region == region_tmp)  %>%
            mutate(year=as.factor(year))
          history_data.hex <- as.h2o(history_data)
          history_data$pred_DM <-  h2o.predict(object = final_model, newdata =  history_data.hex) %>% as.vector()
          
          history_data <- history_data %>% mutate(pred_DM = (pred_DM + outcome_mean) * area)  %>%
            select(year, region, region_sub, pred_DM, 
                   narr_temp, narr_precip, narr_rhum, narr_vpd, narr_runoff, narr_wspd, nldas_soilm) %>%
            mutate(algorithm = "Neural Net", spatial_res=experiment, outcome=outcome, spec=spec, scenario="observation")
        }
        
        #history_data %>% group_by(year) %>% summarise(pred_DM=sum(pred_DM)/1e9) %>% ungroup()
        
        if (outcome == "log"){
          history_data_full <- readRDS(paste0(output_path,"/final_train_region_FE_log.rds"))
          history_data_grid <- readRDS(paste0(output_path,"/final_train_region_FE_log_1deg_grid.rds"))
          
          if (experiment == "regional"){
            cmip_data <- readRDS(paste0(proj_path,"/cmip_region_full_demean_log.rds")) %>%
              mutate(region_sub = region)
            history_data <- history_data_full$data_region
          }else if (experiment == "eco2"){
            cmip_data <- readRDS(paste0(proj_path,"/cmip_eco2_full_demean_log.rds"))
            history_data <- history_data_full$data_eco2
          }else if (experiment == "eco3"){
            cmip_data <- readRDS(paste0(proj_path,"/cmip_eco3_full_demean_log.rds"))
            history_data <- history_data_full$data_eco3
          }else if (experiment == "grid"){
            cmip_data <- readRDS(paste0(proj_path,"/cmip_region_1deg_full_demean_log.rds")) %>%
              rename(region_sub = grid_1deg_id) %>%
              mutate(region_sub = as.character(region_sub)) %>%
              filter(scenario != "ssp585" & year<2071) ## save memory
            
            history_data <- history_data_grid %>%
              rename(region_sub = grid_1deg_id) %>%
              mutate(region_sub = as.character(region_sub))
          }
          
          #### normalize the features in CMIP data
          for (vv in unique(data_norm$variable)){
            mean_tmp <- filter(data_norm, variable == vv, stat=="mean")$value
            sd_tmp <- filter(data_norm, variable == vv, stat=="sd")$value
            cmip_data[,vv] <- (cmip_data[,vv] - mean_tmp)/sd_tmp
            history_data[,vv] <- (history_data[,vv] - mean_tmp)/sd_tmp
          }
          
          cmip_proj <- cmip_data %>%
            filter(region == region_tmp)  %>%
            filter_at(vlist_climate, all_vars(!is.na(.))) %>%
            mutate(year=as.factor(year))
          
          cmip_proj.hex <- as.h2o(cmip_proj)
          cmip_proj$pred_DM <-  h2o.predict(object = final_model, newdata = cmip_proj.hex) %>% as.vector()
          
          cmip_proj <- cmip_proj %>% mutate(pred_DM = exp(pred_DM + outcome_mean_log) * area) %>%
            select(model, scenario, year, region, region_sub, pred_DM,
                   narr_temp, narr_precip, narr_rhum, narr_vpd, narr_runoff, narr_wspd, nldas_soilm) %>%
            mutate(algorithm = "Neural Net", spatial_res=experiment, outcome=outcome, spec=spec)
          
          ### use the real observations to predict 
          history_data <- history_data %>% filter(region == region_tmp) %>%
            mutate(year=as.factor(year))
          history_data.hex <- as.h2o(history_data)
          history_data$pred_DM <-  h2o.predict(object = final_model, newdata =  history_data.hex) %>% as.vector()
          
          history_data <- history_data %>% mutate(pred_DM = exp(pred_DM + outcome_mean) * area)  %>%
            select(year, region, region_sub, pred_DM,
                   narr_temp, narr_precip, narr_rhum, narr_vpd, narr_runoff, narr_wspd, nldas_soilm) %>%
            mutate(algorithm = "Neural Net", spatial_res=experiment, outcome=outcome, spec=spec, scenario="observation")
        }
        
        nnet_cmip_proj <- bind_rows(nnet_cmip_proj,  cmip_proj)
        nnet_history <- bind_rows(nnet_history,  history_data)
        
        final_model_list[[iii]] <- final_model
        names(final_model_list)[iii] <- paste(region_tmp, experiment, outcome, spec, sep="_")
        
        data_norm_list[[iii]] <- data_norm
        names(data_norm_list)[iii] <- paste(region_tmp, experiment, outcome, spec, sep="_")
        
      }
      # h2o.saveModel(final_model,
      #               path="oak_space/mhqiu/smoke_climate/output",
      #               force=T,
      #               filename = paste0(paste("NN",region_label,experiment,outcome,spec, sep="_"),".zip"))
      # h2o.save_mojo(final_model,
      #                path="oak_space/mhqiu/smoke_climate/output",
      #                force=T,
      #                filename = paste0(paste("NN",region_label,experiment,outcome,spec, sep="_"),".zip"))
      
        # }
  #       }
  saveRDS(final_model_list, file = paste0(output_path, "/NN_h2o_final_models.rds"))
  saveRDS(data_norm_list, file = paste0(output_path, "/NN_norm_data.rds"))

  saveRDS(nnet_cmip_proj, paste0(output_path,"/nnet_final_pred.rds"))
  saveRDS(nnet_history, paste0(output_path,"/nnet_final_pred_history.rds"))
#  }

h2o.shutdown(prompt = F)

#--------------------------------------
###       Using caret
#--------------------------------------
# if (final_train){
#   for (rrr in 1:length(region_list)){
#     print(region_list[rrr])
#     region_tmp <-   region_list[rrr] 
#     region_label <- region_label_list[rrr]
#     
#     best_model_tmp <- filter(best_models, region == region_tmp)     
#     
#     final_model_list <- vector(mode='list', length=nrow(best_model_tmp))
#     data_norm_list <- vector(mode='list', length=nrow(best_model_tmp))
#     for (ii in 1:nrow(best_model_tmp)){
#       experiment <- best_model_tmp[ii,"experiment"]
#       threshold_tmp <- best_model_tmp[ii,"threshold"]
#       spec <-  best_model_tmp[ii,"spec"]
#       outcome <- best_model_tmp[ii,"outcome"]
#       
#       ## determine the outcome variable
#       if (outcome =="log"){
#         data_full <- readRDS(paste0(output_path, "/final_train_region_FE_log.rds"))
#       } else{
#         data_full <- readRDS(paste0(output_path, "/final_train_region_FE_level.rds"))
#       }
#       
#       ## determine whether use full or simple features
#       simple_features <- as.numeric(spec == "2 layers simple")
#       
#       ## determine spatial resolution of the model
#       if (experiment == "regional"){
#         data <- data_full$data_region
#       }else if (experiment == "eco2"){
#         data <- data_full$data_eco2
#       }else if (experiment == "eco3"){
#         data <- data_full$data_eco3
#       }else if (experiment == "grid"){
#         data <- data_full$data_grid 
#       }
#       
#       data_train <- data %>% filter(region == region_tmp) %>%
#         mutate(row_ids = row_number())
#       
#       if (region_tmp %in% c("Canada-Alaska","Mexico")){
#         vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
#                            "narr_vpd", "narr_runoff")
#         if (experiment != "regional"){
#           var_list <- c(vlist_climate, paste(rep(vlist_climate, each = 3),
#                                              c("forest","crop", "grass"), sep="_"),
#                         paste0(vlist_climate, "_temp"),
#                         paste0(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff"), "_precip"),
#                         paste0(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff"), "_rhum"),
#                         paste0(c("narr_wspd",  "narr_vpd", "narr_runoff"), "_wspd"),
#                         paste0(c("narr_vpd", "narr_runoff"), "_vpd"),
#                         paste0(c("narr_runoff"), "_runoff"))
#         } else {
#           var_list <- c(vlist_climate,
#                         paste0(vlist_climate, "_temp"),
#                         paste0(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff"), "_precip"),
#                         paste0(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff"), "_rhum"),
#                         paste0(c("narr_wspd",  "narr_vpd", "narr_runoff"), "_wspd"),
#                         paste0(c("narr_vpd", "narr_runoff"), "_vpd"),
#                         paste0(c("narr_runoff"), "_runoff"))
#         }
#       } 
#       
#       if (!region_tmp %in% c("Canada-Alaska","Mexico")){        
#         vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
#                            "narr_vpd", "narr_runoff", "nldas_soilm")
#         if (experiment != "regional"){
#           var_list <- c(vlist_climate, paste(rep(vlist_climate, each = 3),
#                                              c("forest","crop", "grass"), sep="_"),
#                         paste0(vlist_climate, "_temp"),
#                         paste0(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff", "nldas_soilm"), "_precip"),
#                         paste0(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_rhum"),
#                         paste0(c("narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_wspd"),
#                         paste0(c("narr_vpd", "narr_runoff", "nldas_soilm"), "_vpd"),
#                         paste0(c("narr_runoff", "nldas_soilm"), "_runoff"),
#                         paste0(c("nldas_soilm"), "_soilm"))
#         } else {
#           var_list <- c(vlist_climate,
#                         paste0(vlist_climate, "_temp"),
#                         paste0(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff", "nldas_soilm"), "_precip"),
#                         paste0(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_rhum"),
#                         paste0(c("narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_wspd"),
#                         paste0(c("narr_vpd", "narr_runoff", "nldas_soilm"), "_vpd"),
#                         paste0(c("narr_runoff", "nldas_soilm"), "_runoff"),
#                         paste0(c("nldas_soilm"), "_soilm"))
#         }
#       }
#       
#       if (simple_features){
#         var_list <- vlist_climate
#       }
#       
#       if (feature_normalize){
#         ###### compute the mean and sd of each variable ###### 
#         data_norm_mean <-  data_train %>% 
#           summarise_at(var_list, .funs = mean, na.rm=T) %>% t() %>% as.data.frame()
#         data_norm_mean <- data_norm_mean %>% mutate(variable= row.names(data_norm_mean),stat ="mean") %>%
#           rename(value = V1)
#         data_norm_sd <-  data_train %>% 
#           summarise_at(var_list, .funs = sd, na.rm=T) %>% t() %>% as.data.frame()
#         data_norm_sd <- data_norm_sd %>% mutate(variable= row.names(data_norm_sd),stat ="sd") %>%
#           rename(value = V1)
#         
#         data_norm <- bind_rows(data_norm_mean, data_norm_sd)
#         
#         data_train <- data_train %>% 
#           mutate_at(var_list, .fun=function(x){return((x-mean(x, na.rm=T))/sd(x, na.rm=T))})
#       }
#       
#       ##### set up training, tune the hyper-param with outfolds 
#       outer_sampling <-  split(data_train$row_ids, data_train$temporal_fold)
#       
#       grid <-  expand.grid(layer1 = c(1, 2, 4, 8, 16),
#                            layer2 = c(2, 4, 8, 16, 32, 64),
#                            layer3 = 1)
#       
#       final_model <- caret::train(outcome_demean ~ .,
#                                   data = data_train[,c("outcome_demean", var_list)],
#                                   method = "neuralnet",
#                                   tuneGrid = grid,
#                                   metric = "RMSE",
#                                   act.fct = "tanh",
#                                   linear.output= T,
#                                   threshold = threshold_tmp,
#                                   stepmax = 2e+07,
#                                   weights = data_train$sd_DM,
#                                   #preProc = c("center", "scale"), #good idea to do this with neural nets - your error is due to non scaled data
#                                   trControl = trainControl(
#                                     method = "cv",
#                                     index = outer_sampling,
#                                     allowParallel= F,
#                                     verboseIter = T))
#       
#       final_model_list[[ii]] <- final_model
#       names(final_model_list)[ii] <- paste(experiment, outcome, sep="_")
#       
#       data_norm_list[[ii]] <- data_norm
#       names(data_norm_list)[ii] <- paste(experiment, outcome, sep="_")
#     }
#     
#     saveRDS(final_model_list, file = paste0(output_path, "/",region_label,"_NN_final_normalize_models.rds"))
#     saveRDS(data_norm_list, file = paste0(output_path, "/",region_label,"_NN_norm_data.rds"))
#   }
# }

