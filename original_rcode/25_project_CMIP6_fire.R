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
library(h2o)

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
###---- Load the area of different regions to be used for predictions and training purposes
region_list <- c("Western US", "Canada-Alaska", "Mexico", "Southeastern US", "Northeastern US") ##unique(gfed_data$region)
region_label_list <- c("west", "canada", "mexico", "southeast", "northeast")

best_models <- read.xlsx(paste0(result_path,"/model_eval_final/best_metrics_rmse_bias_6models_LOOCV.xlsx"))
best_models1 <- read.xlsx(paste0(result_path,"/model_eval_final/best_metrics_rmse_bias_6models_subregion_LOOCV.xlsx"))

best_models <- bind_rows(best_models, best_models1) %>% 
  filter(selected==1, metric %in% c("5 years", "10 years")) %>%
  distinct(region, model, experiment, outcome, spec, .keep_all = T)

#### to get the complete GCM list 
cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_eco3_full_demean_level.rds")) %>%
  filter(year == 1997) %>%
  filter_at(vars(narr_temp:narr_wspd), all_vars(!is.na(.)))

gcm_list <- unique(cmip_data$model)

# historical_train <-  readRDS(paste0(result_path,"/final_models/training_data/final_train_region_FE_log.rds"))
# cmip_data1 <- historical_train$data_region %>% filter(region == "Western US")

#---------------------------------------------------
# Project future fire emissions using linear models
#---------------------------------------------------
best_linear_models <- filter(best_models, model=="Linear") %>%
  distinct(region, model, experiment, outcome, spec)

linear_cmip_proj <- NULL
linear_history_proj <- NULL
for (ii in 1:nrow(best_linear_models)){
  
  region_tmp <-   best_linear_models[ii, "region"] 
  region_label <- region_label_list[which(region_list == region_tmp)]
  
  experiment <-  best_linear_models[ii, "experiment"] 
  outcome <- best_linear_models[ii, "outcome"] 
  spec <- best_linear_models[ii, "spec"] 
  
  print(paste(region_label, experiment, outcome, spec))
  
  if (experiment == "regional"){
    experiment <- "region"
  }
  
  if (region_tmp %in%c("Canada-Alaska","Mexico")){        
    vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
                       "narr_vpd", "narr_runoff")
  } else{
    vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
                       "narr_vpd", "narr_runoff", "nldas_soilm")
  }
    
 if (outcome == "level"){
   history_data <- readRDS(paste0(result_path,"/final_models/training_data/final_train_region_FE_level.rds"))
   history_data_grid <- readRDS(paste0(result_path,"/final_models/training_data/final_train_region_FE_level_1deg_grid.rds"))
   
   if (spec == "weighted FE"){
     model_tmp <- readRDS(paste0(result_path, "/final_models/linear/",region_label,"_linear_level_weighted_final_models.rds"))[[paste0("linear_",experiment)]]
   } else{
     model_tmp <- readRDS(paste0(result_path, "/final_models/linear/",region_label,"_linear_level_final_models.rds"))[[paste0("linear_",experiment)]]
   }
   
    if (experiment == "region"){
      cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_region_full_demean_level.rds")) %>%
        mutate(region_sub = region)
      history_data <- history_data$data_region
    }else if (experiment == "eco2"){
      cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_eco2_full_demean_level.rds"))
      history_data <- history_data$data_eco2
    }else if (experiment == "eco3"){
      cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_eco3_full_demean_level.rds"))
      history_data <- history_data$data_eco3
    }else if (experiment == "grid"){
      cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_region_1deg_full_demean_level_2099.rds")) %>%
        rename(region_sub = grid_1deg_id) %>%
        mutate(region_sub = as.character(region_sub))
      # %>%
      #   filter(scenario != "ssp585" & year<2071) ## save memory
      
      history_data <- history_data_grid %>%
        rename(region_sub = grid_1deg_id) %>%
        mutate(region_sub = as.character(region_sub))
    }
    
   cmip_proj <- cmip_data %>%
     filter(region == region_tmp)  %>%
     filter(model %in% gcm_list) %>%
     filter_at(vlist_climate, all_vars(!is.na(.)))
   
   cmip_proj$pred_DM <- predict(model_tmp, newdata = cmip_proj)
   cmip_proj <- cmip_proj %>% mutate(pred_DM = (pred_DM + outcome_mean_level) * area)  %>%
     select(model, scenario, year, region, region_sub, pred_DM, 
            narr_temp, narr_precip, narr_rhum, narr_vpd, narr_runoff, narr_wspd, nldas_soilm) %>%
     mutate(algorithm = "Linear", spatial_res=experiment, outcome=outcome, spec=spec) 
   
   history_data <- history_data %>% filter(region == region_tmp) 
   history_data$pred_DM <- predict(model_tmp, newdata = history_data)
   history_data <- history_data %>% mutate(pred_DM = (pred_DM + outcome_mean) * area) %>%
     select(year, region, region_sub, pred_DM, 
            narr_temp, narr_precip, narr_rhum, narr_vpd, narr_runoff, narr_wspd, nldas_soilm) %>%
     mutate(algorithm = "Linear", spatial_res=experiment, outcome=outcome, spec=spec, scenario="observation") 
  }

    if (outcome == "log"){
      history_data <- readRDS(paste0(result_path,"/final_models/training_data/final_train_region_FE_log.rds"))
      history_data_grid <- readRDS(paste0(result_path,"/final_models/training_data/final_train_region_FE_log_1deg_grid.rds"))
      
      if (spec == "weighted FE"){
        model_tmp <- readRDS(paste0(result_path, "/final_models/linear/",region_label,"_linear_log_weighted_final_models.rds"))[[paste0("linear_",experiment)]]
      } else{
        model_tmp <- readRDS(paste0(result_path, "/final_models/linear/",region_label,"_linear_log_final_models.rds"))[[paste0("linear_",experiment)]]
      }
    
      if (experiment == "region"){
        cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_region_full_demean_log.rds")) %>%
          mutate(region_sub = region)
        history_data <- history_data$data_region
      }else if (experiment == "eco2"){
        cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_eco2_full_demean_log.rds"))
        history_data <- history_data$data_eco2
      }else if (experiment == "eco3"){
        cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_eco3_full_demean_log.rds"))
        history_data <- history_data$data_eco3
      }else if (experiment == "grid"){
        cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_region_1deg_full_demean_log_2099.rds")) %>%
          rename(region_sub = grid_1deg_id) %>%
          mutate(region_sub = as.character(region_sub)) 
        # %>%
        #   filter(scenario != "ssp585" & year<2071) ## save memory
        # 
        history_data <- history_data_grid %>%
          rename(region_sub = grid_1deg_id) %>%
          mutate(region_sub = as.character(region_sub))
      }
      
      cmip_proj <- cmip_data %>%
        filter(region == region_tmp)  %>%
        filter(model %in% gcm_list) %>%
        filter_at(vlist_climate, all_vars(!is.na(.)))
      
      cmip_proj$pred_DM <- predict(model_tmp, newdata = cmip_proj)
      cmip_proj <- cmip_proj %>% mutate(pred_DM = exp(pred_DM + outcome_mean_log) * area) %>%
        select(model, scenario, year, region, region_sub, pred_DM, 
               narr_temp, narr_precip, narr_rhum, narr_vpd, narr_runoff, narr_wspd, nldas_soilm) %>%
        mutate(algorithm = "Linear", spatial_res=experiment, outcome=outcome, spec=spec) 
      
      ### use the real observations to predict 
      history_data <- history_data %>% filter(region == region_tmp) 
      history_data$pred_DM <- predict(model_tmp, newdata = history_data)
      history_data <- history_data %>% mutate(pred_DM = exp(pred_DM + outcome_mean) * area) %>%
        select(year, region, region_sub, pred_DM, 
               narr_temp, narr_precip, narr_rhum, narr_vpd, narr_runoff, narr_wspd, nldas_soilm) %>%
        mutate(algorithm = "Linear", spatial_res=experiment, outcome=outcome, spec=spec, scenario="observation") 
    }
    
    linear_cmip_proj <- bind_rows(linear_cmip_proj,  cmip_proj)
    linear_history_proj <- bind_rows(linear_history_proj,  history_data)
}

saveRDS(linear_cmip_proj, paste0(result_path,"/CMIP6_fire_proj/linear_final_pred_2099.rds"))
saveRDS(linear_history_proj, paste0(result_path,"/CMIP6_fire_proj/linear_final_pred_history_2099.rds"))

#---------------------------------------------------
# Project future fire emissions using LASSO models
#---------------------------------------------------
best_lasso_models <- filter(best_models, model=="LASSO") %>%
  distinct(region, model, experiment, outcome, spec)

lasso_cmip_proj <- NULL
lasso_history <- NULL
for (ii in 1:nrow(best_lasso_models)){
  
  region_tmp <-   best_lasso_models[ii, "region"] 
  region_label <- region_label_list[which(region_list == region_tmp)]
  
  experiment <-  best_lasso_models[ii, "experiment"]
  outcome <- best_lasso_models[ii, "outcome"] 
  spec <- best_lasso_models[ii, "spec"] 
  
  print(paste(region_label, experiment, outcome, spec))
  
  if (experiment == "regional"){
    experiment <- "region"
  }
  
  if (region_tmp %in%c("Canada-Alaska","Mexico")){        
    vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
                       "narr_vpd", "narr_runoff")
    } else{
    vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
                        "narr_vpd", "narr_runoff", "nldas_soilm")
    }
    
    if (outcome == "level"){
      if (spec == "weighted FE"){
        model_tmp <- readRDS(paste0(result_path, "/final_models/LASSO/",region_label,"_LASSO_level_weighted_final_models.rds"))[[paste0("lasso_",experiment)]]
      } else{
        model_tmp <- readRDS(paste0(result_path, "/final_models/LASSO/",region_label,"_LASSO_level_final_models.rds"))[[paste0("lasso_",experiment)]]
      }
      
      coef <- model_tmp$model$beta  %>% as.matrix() %>% as.data.frame() %>% filter(s0!=0)
      print(coef)
      
      if (nrow(coef)==0) {next}
      
      history_data <- readRDS(paste0(result_path,"/final_models/training_data/final_train_region_FE_level.rds"))
      history_data_grid <- readRDS(paste0(result_path,"/final_models/training_data/final_train_region_FE_level_1deg_grid.rds"))
      
      if (experiment == "region"){
        cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_region_full_demean_level.rds")) %>%
          mutate(region_sub = region)
        history_data <- history_data$data_region
      }else if (experiment == "eco2"){
        cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_eco2_full_demean_level.rds"))
        history_data <- history_data$data_eco2
      }else if (experiment == "eco3"){
        cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_eco3_full_demean_level.rds"))
        history_data <- history_data$data_eco3
      }else if (experiment == "grid"){
        cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_region_1deg_full_demean_level.rds")) %>%
          rename(region_sub = grid_1deg_id) %>%
          mutate(region_sub = as.character(region_sub)) %>%
          filter(scenario != "ssp585" & year<2071) ## save memory
        
        history_data <- history_data_grid %>%
          rename(region_sub = grid_1deg_id) %>%
          mutate(region_sub = as.character(region_sub))
      }
      
      cmip_proj <- cmip_data %>%
        filter(region == region_tmp)  %>%
      filter(model %in% gcm_list) %>%
      filter_at(vlist_climate, all_vars(!is.na(.)))

      cmip_proj$pred_DM <- model_tmp$predict_newdata(newdata = cmip_proj) %>%
        as.data.table() %>% select(response) %>% unlist() %>% as.vector()

      cmip_proj <- cmip_proj %>% mutate(pred_DM = (pred_DM + outcome_mean_level) * area)  %>%
        select(model, scenario, year, region, region_sub, pred_DM,
               narr_temp, narr_precip, narr_rhum, narr_vpd, narr_runoff, narr_wspd, nldas_soilm) %>%
        mutate(algorithm = "LASSO", spatial_res=experiment, outcome=outcome, spec=spec)
      
      ### use the real observations to predict 
      history_data <- history_data %>% filter(region == region_tmp) 
      history_data$pred_DM <- model_tmp$predict_newdata(newdata = history_data) %>%
        as.data.table() %>% select(response) %>% unlist() %>% as.vector()
      
      history_data <- history_data %>% mutate(pred_DM = (pred_DM + outcome_mean) * area)  %>%
        select(year, region, region_sub, pred_DM,
               narr_temp, narr_precip, narr_rhum, narr_vpd, narr_runoff, narr_wspd, nldas_soilm) %>%
        mutate(algorithm = "LASSO", spatial_res=experiment, outcome=outcome, spec=spec, scenario="observation")
    }
    
    if (outcome == "log"){
      if (spec == "weighted FE"){
        model_tmp <- readRDS(paste0(result_path, "/final_models/LASSO/",region_label,"_LASSO_log_weighted_final_models.rds"))[[paste0("lasso_",experiment)]]
      } else{
        model_tmp <- readRDS(paste0(result_path, "/final_models/LASSO/",region_label,"_LASSO_log_final_models.rds"))[[paste0("lasso_",experiment)]]
      }
      coef <- model_tmp$model$beta  %>% as.matrix() %>% as.data.frame() %>% filter(s0!=0)
      print(coef)
      
      if (nrow(coef)==0) {next}
    
      history_data <- readRDS(paste0(result_path,"/final_models/training_data/final_train_region_FE_log.rds"))
      history_data_grid <- readRDS(paste0(result_path,"/final_models/training_data/final_train_region_FE_log_1deg_grid.rds"))
      
      if (experiment == "region"){
        cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_region_full_demean_log.rds")) %>%
          mutate(region_sub = region)
        history_data <- history_data$data_region
      }else if (experiment == "eco2"){
        cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_eco2_full_demean_log.rds"))
        history_data <- history_data$data_eco2
      }else if (experiment == "eco3"){
        cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_eco3_full_demean_log.rds"))
        history_data <- history_data$data_eco3
      }else if (experiment == "grid"){
        cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_region_1deg_full_demean_log.rds")) %>%
          rename(region_sub = grid_1deg_id) %>%
          mutate(region_sub = as.character(region_sub)) %>%
          filter(scenario != "ssp585" & year<2071) ## save memory
        
        history_data <- history_data_grid %>%
          rename(region_sub = grid_1deg_id) %>%
          mutate(region_sub = as.character(region_sub))
      }

      cmip_proj <- cmip_data %>%
        filter(region == region_tmp)  %>%
        filter(model %in% gcm_list) %>%
        filter_at(vlist_climate, all_vars(!is.na(.)))

      cmip_proj$pred_DM <- model_tmp$predict_newdata(newdata = cmip_proj) %>%
        as.data.table() %>% select(response) %>% unlist() %>% as.vector()

      cmip_proj <- cmip_proj %>% mutate(pred_DM = exp(pred_DM + outcome_mean_log) * area) %>%
        select(model, scenario, year, region, region_sub, pred_DM,
               narr_temp, narr_precip, narr_rhum, narr_vpd, narr_runoff, narr_wspd, nldas_soilm) %>%
        mutate(algorithm = "LASSO", spatial_res=experiment, outcome=outcome, spec=spec)
      
      ### use the real observations to predict 
      history_data <- history_data %>% filter(region == region_tmp) 
      history_data$pred_DM <- model_tmp$predict_newdata(newdata = history_data) %>%
        as.data.table() %>% select(response) %>% unlist() %>% as.vector()
      
      history_data <- history_data %>% mutate(pred_DM = exp(pred_DM + outcome_mean) * area) %>%
        select(year, region, region_sub, pred_DM,
               narr_temp, narr_precip, narr_rhum, narr_vpd, narr_runoff, narr_wspd, nldas_soilm) %>%
        mutate(algorithm = "LASSO", spatial_res=experiment, outcome=outcome, spec=spec, scenario="observation")
      
    }
    lasso_cmip_proj <- bind_rows(lasso_cmip_proj,  cmip_proj)
    lasso_history <- bind_rows(lasso_history,  history_data)
  }

saveRDS(lasso_cmip_proj, paste0(result_path,"/CMIP6_fire_proj/LASSO_final_pred.rds"))
saveRDS(lasso_history, paste0(result_path,"/CMIP6_fire_proj/LASSO_final_pred_history.rds"))

#---------------------------------------------------
# Get all the coefs from the selected models
#---------------------------------------------------
region_list <- c("Western US", "Canada-Alaska", "Mexico", "Southeastern US", "Northeastern US") ##unique(gfed_data$region)
region_label_list <- c("west", "canada", "mexico", "southeast", "northeast")

best_models <- read.xlsx(paste0(result_path,"/model_eval_final/best_metrics_rmse_bias_6models_LOOCV.xlsx"))
best_models <- bind_rows(best_models) %>% 
  filter(diff_min<0.05, metric %in% c("10 years")) %>%
  distinct(region, model, experiment, outcome, spec, .keep_all = T)

best_linear_models <- filter(best_models, model=="Linear") %>%
  distinct(region, model, experiment, outcome, spec)

linear_coef <- NULL
for (ii in 1:nrow(best_linear_models)){
  region_tmp <-   best_linear_models[ii, "region"] 
  region_label <- region_label_list[which(region_list == region_tmp)]
  
  experiment <-  best_linear_models[ii, "experiment"] 
  outcome <- best_linear_models[ii, "outcome"] 
  spec <- best_linear_models[ii, "spec"] 
  
  print(paste(region_label, experiment, outcome, spec))
  
  if (experiment == "regional"){
    experiment <- "region"
  }
  
  if (outcome == "level"){
      if (spec == "weighted FE"){
      model_tmp <- readRDS(paste0(result_path, "/final_models/linear/",region_label,"_linear_level_weighted_final_models.rds"))[[paste0("linear_",experiment)]]
    } else{
      model_tmp <- readRDS(paste0(result_path, "/final_models/linear/",region_label,"_linear_level_final_models.rds"))[[paste0("linear_",experiment)]]
    }
   }
  
  if (outcome == "log"){
    if (spec == "weighted FE"){
      model_tmp <- readRDS(paste0(result_path, "/final_models/linear/",region_label,"_linear_log_weighted_final_models.rds"))[[paste0("linear_",experiment)]]
    } else{
      model_tmp <- readRDS(paste0(result_path, "/final_models/linear/",region_label,"_linear_log_final_models.rds"))[[paste0("linear_",experiment)]]
    }
  }
 coef_tmp <- data.frame(coef=coef(model_tmp), se=se(model_tmp), pval=pvalue(model_tmp), 
                        region=region_tmp, experiment=experiment, outcome=outcome, spec=spec)
 coef_tmp$var <- row.names(coef_tmp)
 linear_coef <- bind_rows(coef_tmp, linear_coef)
}
write.xlsx(linear_coef,paste0(result_path,"/linear_coef_fire_climate_models.xlsx"))

best_lasso_models <- filter(best_models, model=="LASSO") %>%
  distinct(region, model, experiment, outcome, spec)

lasso_coef <- NULL
for (ii in 1:nrow(best_lasso_models)){
  region_tmp <-   best_lasso_models[ii, "region"] 
  region_label <- region_label_list[which(region_list == region_tmp)]
  
  experiment <-  best_lasso_models[ii, "experiment"] 
  outcome <- best_lasso_models[ii, "outcome"] 
  spec <- best_lasso_models[ii, "spec"] 
  
  print(paste(region_label, experiment, outcome, spec))
  
  if (experiment == "regional"){
    experiment <- "region"
  }
  if (outcome == "level"){
    if (spec == "weighted FE"){
      model_tmp <- readRDS(paste0(result_path, "/final_models/lasso/",region_label,"_lasso_level_weighted_final_models.rds"))[[paste0("lasso_",experiment)]]
    } else{
      model_tmp <- readRDS(paste0(result_path, "/final_models/lasso/",region_label,"_lasso_level_final_models.rds"))[[paste0("lasso_",experiment)]]
    }
  }
  
  if (outcome == "log"){
    if (spec == "weighted FE"){
      model_tmp <- readRDS(paste0(result_path, "/final_models/lasso/",region_label,"_lasso_log_weighted_final_models.rds"))[[paste0("lasso_",experiment)]]
    } else{
      model_tmp <- readRDS(paste0(result_path, "/final_models/lasso/",region_label,"_lasso_log_final_models.rds"))[[paste0("lasso_",experiment)]]
    }
  }
  
  coef <- model_tmp$model$beta  %>% as.matrix() %>% as.data.frame() %>% filter(s0!=0)
  coef_tmp <- data.frame(coef=coef,
                         region=region_tmp, experiment=experiment, outcome=outcome, spec=spec)
  coef_tmp$var <- row.names(coef_tmp)
  lasso_coef <- bind_rows(coef_tmp, lasso_coef)
}
write.xlsx(lasso_coef,paste0(result_path,"/lasso_coef_fire_climate_models.xlsx"))

#---------------------------------------------------
# Project future fire emissions using Neural Net models: this is completed in the prediction stage
# predict with existing model (that is trained somewhere else) seems very challenging...
#---------------------------------------------------
# best_nnet_models <- filter(best_models, model=="Neural Net",
#                            selected ==1) %>%
#   distinct(region, model, experiment, outcome, spec) 
# 
# h2o.init(max_mem_size = "10g",
#          nthreads = 1)
# 
# nnet_cmip_proj <- NULL
# nnet_history <- NULL
# for (ii in 1:nrow(best_nnet_models)){
#   region_tmp <-   best_nnet_models[ii, "region"] 
#   region_label <- region_label_list[which(region_list == region_tmp)]
#   
#   experiment <-  best_nnet_models[ii, "experiment"] 
#   outcome <- best_nnet_models[ii, "outcome"] 
#   spec <- best_nnet_models[ii, "spec"] 
#   
#   print(paste(region_label, experiment, outcome))
#   
#   if (region_tmp %in%c("Canada-Alaska","Mexico")){        
#     vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
#                        "narr_vpd", "narr_runoff")
#   } else{
#     vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
#                        "narr_vpd", "narr_runoff", "nldas_soilm")
#   }
#   
#   model_tmp <- readRDS(paste0(result_path, "/final_models/NN/",region_label,"_NN_h2o_final_models.rds"))[[paste(experiment,outcome,spec, sep="_")]]
#   model_tmp <- h2o.import_mojo(paste0(result_path, "/final_models/NN/NN_mexico_regional_log_full.zip"))
#   data_norm <- readRDS(paste0(result_path, "/final_models/NN/",region_label,"_NN_norm_data.rds"))[[paste(experiment,outcome,spec, sep="_")]]
# 
#   if (is.null(model_tmp)) {next} ##if all weights are zero
#   
#   if (outcome == "level"){
#     history_data <- readRDS(paste0(result_path,"/final_models/training_data/final_train_region_FE_level.rds"))
#     history_data_grid <- readRDS(paste0(result_path,"/final_models/training_data/final_train_region_FE_level_1deg_grid.rds"))
#     
#     if (experiment == "regional"){
#       cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_region_full_demean_level.rds")) %>%
#         mutate(region_sub = region)
#       history_data <- history_data$data_region
#     }else if (experiment == "eco2"){
#       cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_eco2_full_demean_level.rds"))
#       history_data <- history_data$data_eco2
#     }else if (experiment == "eco3"){
#       cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_eco3_full_demean_level.rds"))
#       history_data <- history_data$data_eco3
#     }else if (experiment == "grid"){
#       cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_region_1deg_full_demean_level.rds")) %>%
#         rename(region_sub = grid_1deg_id) %>%
#         mutate(region_sub = as.character(region_sub)) %>%
#         filter(scenario != "ssp585" & year<2071) ## save memory
#       
#       history_data <- history_data_grid %>%
#         rename(region_sub = grid_1deg_id) %>%
#         mutate(region_sub = as.character(region_sub))
#     }
#     
#     #### normalize the features in CMIP data
#     for (vv in unique(data_norm$variable)){
#       mean_tmp <- filter(data_norm, variable == vv, stat=="mean")$value
#       sd_tmp <- filter(data_norm, variable == vv, stat=="sd")$value
#       cmip_data[,vv] <- (cmip_data[,vv] - mean_tmp)/sd_tmp
#       history_data[,vv] <- (history_data[,vv] - mean_tmp)/sd_tmp
#     }
#     
#     cmip_proj <- cmip_data %>%
#       filter(region == region_tmp)  %>%
#       filter(model %in% gcm_list) %>%
#       filter_at(vlist_climate, all_vars(!is.na(.)))
#     
#     
#     cmip_proj$pred_DM <-  predict(model_tmp, cmip_proj) %>% as.vector()
#     
#     cmip_proj <- cmip_proj %>% mutate(pred_DM = (pred_DM + outcome_mean_level) * area)  %>%
#       select(model, scenario, year, region, region_sub, pred_DM,
#              narr_temp, narr_precip, narr_rhum, narr_vpd, narr_runoff, narr_wspd, nldas_soilm) %>%
#       mutate(algorithm = "Neural Net", spatial_res=experiment, outcome=outcome, spec=spec)
#     
#     ### use the real observations to predict 
#     history_data <- history_data %>% filter(region == region_tmp) 
#     history_data$pred_DM <-  predict(model_tmp, history_data) %>% as.vector()
#     
#     history_data <- history_data %>% mutate(pred_DM = (pred_DM + outcome_mean) * area)  %>%
#       select(year, region, region_sub, pred_DM,
#              narr_temp, narr_precip, narr_rhum, narr_vpd, narr_runoff, narr_wspd, nldas_soilm) %>%
#       mutate(algorithm = "Neural Net", spatial_res=experiment, outcome=outcome, spec=spec, scenario="observation")
#   }
#   
#   if (outcome == "log"){
#     history_data <- readRDS(paste0(result_path,"/final_models/training_data/final_train_region_FE_log.rds"))
#     history_data_grid <- readRDS(paste0(result_path,"/final_models/training_data/final_train_region_FE_log_1deg_grid.rds"))
#     
#     if (experiment == "regional"){
#       cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_region_full_demean_log.rds")) %>%
#         mutate(region_sub = region)
#       history_data <- history_data$data_region
#     }else if (experiment == "eco2"){
#       cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_eco2_full_demean_log.rds"))
#       history_data <- history_data$data_eco2
#     }else if (experiment == "eco3"){
#       cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_eco3_full_demean_log.rds"))
#       history_data <- history_data$data_eco3
#     }else if (experiment == "grid"){
#       cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_region_1deg_full_demean_log.rds")) %>%
#         rename(region_sub = grid_1deg_id) %>%
#         mutate(region_sub = as.character(region_sub)) %>%
#         filter(scenario != "ssp585" & year<2071) ## save memory
#       
#       history_data <- history_data_grid %>%
#         rename(region_sub = grid_1deg_id) %>%
#         mutate(region_sub = as.character(region_sub))
#     }
#     
#     #### normalize the features in CMIP data
#     for (vv in unique(data_norm$variable)){
#       mean_tmp <- filter(data_norm, variable == vv, stat=="mean")$value
#       sd_tmp <- filter(data_norm, variable == vv, stat=="sd")$value
#       cmip_data[,vv] <- (cmip_data[,vv] - mean_tmp)/sd_tmp
#       history_data[,vv] <- (history_data[,vv] - mean_tmp)/sd_tmp
#     }
#     
#     cmip_proj <- cmip_data %>%
#       filter(region == region_tmp)  %>%
#       filter(model %in% gcm_list) %>%
#       filter_at(vlist_climate, all_vars(!is.na(.)))
#     
#     cmip_proj.hex <- as.h2o(cmip_proj)
#     cmip_proj$pred_DM <-  h2o.predict(object = model_tmp, newdata = cmip_proj.hex) %>% as.vector()
# 
#     
#     cmip_proj <- cmip_proj %>% mutate(pred_DM = exp(pred_DM + outcome_mean_log) * area) %>%
#       select(model, scenario, year, region, region_sub, pred_DM,
#              narr_temp, narr_precip, narr_rhum, narr_vpd, narr_runoff, narr_wspd, nldas_soilm) %>%
#       mutate(algorithm = "Neural Net", spatial_res=experiment, outcome=outcome, spec=spec)
#     
#     ### use the real observations to predict 
#     history_data <- history_data %>% filter(region == region_tmp) 
#     history_data$pred_DM <-  predict(model_tmp, history_data) %>% as.vector()
#     
#     history_data <- history_data %>% mutate(pred_DM = exp(pred_DM + outcome_mean) * area)  %>%
#       select(year, region, region_sub, pred_DM,
#              narr_temp, narr_precip, narr_rhum, narr_vpd, narr_runoff, narr_wspd, nldas_soilm) %>%
#       mutate(algorithm = "Neural Net", spatial_res=experiment, outcome=outcome, spec=spec, scenario="observation")
#   }
#   
#   nnet_cmip_proj <- bind_rows(nnet_cmip_proj,  cmip_proj)
#   nnet_history <- bind_rows(nnet_history,  history_data)
# }
# 
# saveRDS(nnet_cmip_proj, paste0(result_path,"/CMIP6_fire_proj/nnet_final_pred.rds"))
# saveRDS(nnet_history, paste0(result_path,"/CMIP6_fire_proj/nnet_final_pred_history.rds"))


# 
# 
# 
# gfed_agg_data <- readRDS(paste0(data_path, "/gfed_agg_data_veg.rds"))
# 
# data_region <- gfed_agg_data$data_region 
# data_eco3 <- gfed_agg_data$data_eco3 
# data_eco2 <- gfed_agg_data$data_eco2 
# data_grid <- gfed_agg_data$data_grid 
# 
# area_region <- data_region %>% distinct(region, area)
# area_eco3 <- data_eco3 %>% distinct(region, region_sub, area)
# 
# #---------------------------------------------------
# # Project future fire emissions for Western US
# #---------------------------------------------------
# rrr <- "Western US"
# 
# ### Linear models
# linear_models <- readRDS(paste0(result_path, "/final_models/west_linear_models_final.rds"))
# 
# linear_region_level <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_region_full_demean_level.rds")) %>%
#   filter(region == rrr)  %>% 
#   mutate(algorithm = "Linear", spatial_res="region", outcome="level") %>%
#   filter_at(vars(narr_temp:narr_wspd), all_vars(!is.na(.)))
# linear_region_level$pred_DM <- predict(linear_models$linear_region_level, newdata = linear_region_level)
# linear_region_level <- linear_region_level %>% mutate(pred_DM = (pred_DM + outcome_mean_level) * area)
# 
# linear_region_log <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_region_full_demean_log.rds")) %>%
#   filter(region == rrr)  %>% 
#   mutate(algorithm = "Linear", spatial_res="region", outcome="log") %>%
#   filter_at(vars(narr_temp:narr_wspd), all_vars(!is.na(.)))
# linear_region_log$pred_DM <- predict(linear_models$linear_region_log, newdata = linear_region_log)
# linear_region_log <- linear_region_log %>% mutate(pred_DM = exp(pred_DM + outcome_mean_log) * area)
# 
# linear_eco2_log <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_eco2_full_demean_log.rds")) %>%
#   filter(region == rrr)  %>% 
#   mutate(algorithm = "Linear", spatial_res="eco2", outcome="log") %>%
#   filter_at(vars(narr_temp:narr_wspd), all_vars(!is.na(.)))
# linear_eco2_log$pred_DM <- predict(linear_models$linear_eco2_log, newdata = linear_eco2_log)
# linear_eco2_log <- linear_eco2_log %>% mutate(pred_DM = exp(pred_DM + outcome_mean_log) * area)
# 
# linear_eco3_log <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_eco3_full_demean_log.rds")) %>%
#   filter(region == rrr)  %>% 
#   mutate(algorithm = "Linear", spatial_res="eco3", outcome="log") %>%
#   filter_at(vars(narr_temp:narr_wspd), all_vars(!is.na(.)))
# linear_eco3_log$pred_DM <- predict(linear_models$linear_eco3_log, newdata = linear_eco3_log)
# linear_eco3_log <- linear_eco3_log %>% mutate(pred_DM = exp(pred_DM + outcome_mean_log) * area)
# 
# cmip_pred <- bind_rows(linear_region_level, linear_region_log, linear_eco2_log, linear_eco3_log) %>%
#   select(algorithm, model, scenario, year, region, region_sub, pred_DM, 
#          narr_temp, narr_precip, narr_rhum, narr_vpd, narr_soilm, nldas_soilm, narr_runoff, narr_wspd) %>%
#   rename(gcm = model)
# saveRDS(cmip_pred, paste0(result_path,"/CMIP6_fire_proj/west_pred_linear_nosoil.rds"))
# 
# 
# ### LASSO models
# #lasso_level_models  <- readRDS(paste0(result_path, "/final_models/west_LASSO_level_weighted_final_models.rds"))
# lasso_log_models  <- readRDS(paste0(result_path, "/final_models/west_LASSO_log_weighted_final_models.rds"))
# 
# lasso_log_models$lasso_region$model$beta %>% as.matrix() %>% as.data.frame() %>% filter(s0!=0)
# lasso_log_models$lasso_eco2$model$beta %>% as.matrix() %>% as.data.frame() %>% filter(s0!=0)
# 
# lasso_region_log <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_region_full_demean_log.rds")) %>%
#   filter(region == rrr)  %>% 
#   mutate(algorithm = "LASSO", spatial_res="regional", outcome="log") %>%
#   filter_at(vars(narr_temp:narr_wspd), all_vars(!is.na(.)))
# lasso_region_log$pred_DM <- lasso_log_models$lasso_region$predict_newdata(newdata = lasso_region_log) %>% 
#   as.data.table() %>% select(response) %>% unlist() %>% as.vector() 
# lasso_region_log <- lasso_region_log %>% mutate(pred_DM = exp(pred_DM + outcome_mean_log) * area)
# 
# lasso_eco2_log <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_eco2_full_demean_log.rds")) %>%
#   filter(region == rrr)  %>% 
#   mutate(algorithm = "LASSO", spatial_res="eco2", outcome="log") %>%
#   filter_at(vars(narr_temp:narr_wspd), all_vars(!is.na(.)))
# lasso_eco2_log$pred_DM <- lasso_log_models$lasso_eco2$predict_newdata(newdata = lasso_eco2_log) %>% 
#   as.data.table() %>% select(response) %>% unlist() %>% as.vector() 
# lasso_eco2_log <- lasso_eco2_log %>% mutate(pred_DM = exp(pred_DM + outcome_mean_log) * area)
# 
# lasso_eco3_log <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_eco3_full_demean_log.rds")) %>%
#   filter(region == rrr)  %>% 
#   mutate(algorithm = "LASSO", spatial_res="eco3", outcome="log") %>%
#   filter_at(vars(narr_temp:narr_wspd), all_vars(!is.na(.)))
# lasso_eco3_log$pred_DM <- lasso_log_models$lasso_eco3$predict_newdata(newdata = lasso_eco3_log) %>% 
#   as.data.table() %>% select(response) %>% unlist() %>% as.vector() 
# lasso_eco3_log <- lasso_eco3_log %>% mutate(pred_DM = exp(pred_DM + outcome_mean_log) * area)
# 
# cmip_pred <- bind_rows(lasso_region_log, lasso_eco2_log, lasso_eco3_log) %>%
#   select(algorithm, model, scenario, year, region, region_sub, pred_DM, 
#          narr_temp, narr_precip, narr_rhum, narr_vpd, narr_soilm, nldas_soilm, narr_runoff, narr_wspd) %>%
#   rename(gcm = model)
# saveRDS(cmip_pred, paste0(result_path,"/CMIP6_fire_proj/west_pred_lasso.rds"))
# 
# # ## XGBOOST: eco3 level
# # xgb_cmip <- readRDS(paste0(result_path,"/CMIP6_fire_proj/final_cmip_eco3_data_normalize.rds")) %>%
# #   filter(region == rrr) %>%  as.data.frame() %>% rename(region_sub = NA_L3KEY)
# # xgb_cmip <- left_join(xgb_cmip, area_eco3) %>% mutate(algorithm = "XGBOOST")
# # xgb_cmip <- xgb_cmip[complete.cases(xgb_cmip),]
# # xgb_model <- readRDS(paste0(result_path,"/final_models/west_eco3_XGBOOST_norm_final_model.rds"))
# # 
# # xgb_cmip$pred_DM <- predict(xgb_model,
# #                             newdata = xgb.DMatrix(data.matrix(xgb_cmip[, xgb_model$feature_names])))  
# # xgb_cmip$pred_DM <- xgb_cmip$pred_DM * xgb_cmip$area
# 
# # ## single NN regional level model
# # nn_model <- readRDS(paste0(result_path,"/final_models/west_regional_NN_norm_final_model.rds"))
# # 
# # nn_cmip <- readRDS(paste0(result_path,"/CMIP6_fire_proj/final_cmip_region_data_normalize_veg.rds")) %>%
# #   filter(region == rrr) %>%  as.data.frame()
# # nn_cmip <- left_join(nn_cmip, area_region) %>% mutate(algorithm = "Neural Network")
# # nn_cmip <- nn_cmip[complete.cases(nn_cmip),]
# # 
# # nn_cmip$pred_DM <- predict(nn_model, newdata = nn_cmip) %>% exp()
# # nn_cmip$pred_DM <- nn_cmip$pred_DM * nn_cmip$area
# 
# 
# ### LASSO: regional
# lasso_level <- readRDS(paste0(result_path, "/final_models/west_LASSO_level_final_model.rds"))
# lasso_log <- readRDS(paste0(result_path, "/final_models/west_LASSO_log_final_model.rds"))
# 
# lasso_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/final_cmip_region_data.rds")) %>%
#   filter(region == rrr) %>% as.data.frame() 
# lasso_data <- left_join(lasso_data, area_region)
# lasso_data <- lasso_data[complete.cases(lasso_data),]
# 
# lasso_level_cmip <- lasso_data %>% mutate(algorithm = "LASSO (level)")
# lasso_level_cmip$pred_DM <- lasso_level$predict_newdata(newdata = lasso_level_cmip) %>% 
#   as.data.table() %>% select(response) %>% unlist() %>% as.vector() 
# lasso_level_cmip$pred_DM <- lasso_level_cmip$pred_DM * lasso_level_cmip$area
# 
# lasso_log_cmip <- lasso_data %>% mutate(algorithm = "LASSO (log)")
# lasso_log_cmip$pred_DM <- lasso_log$predict_newdata(newdata = lasso_log_cmip) %>%
#   as.data.table() %>% select(response) %>% unlist() %>% as.vector() %>% exp() 
# lasso_log_cmip$pred_DM <- lasso_log_cmip$pred_DM * lasso_log_cmip$area
# 
# cmip_pred <- bind_rows(xgb_cmip, nn_cmip, linear_level_cmip, linear_log_cmip) %>%
#   select(algorithm, model, scenario, year, region, region_sub, pred_DM, 
#          narr_temp, narr_precip, narr_rhum, narr_vpd, narr_soilm, nldas_soilm, narr_runoff, narr_wspd) %>%
#   rename(gcm = model)
# saveRDS(cmip_pred, paste0(result_path,"/CMIP6_fire_proj/west_pred_final_veg.rds"))
# 
# #------------------------------------------------------------------
# #----------      Use the regional data that weighted by veg     ----------#
# #------------------------------------------------------------------
# linear_level <- lm(DM_kg ~ narr_temp + narr_precip + narr_rhum + narr_vpd + 
#                      narr_runoff + narr_wspd + nldas_soilm,
#                    data= filter(data_region, region == rrr))
# 
# linear_log <- lm(log(DM_kg) ~ narr_temp + narr_precip + narr_rhum + narr_vpd + 
#                    narr_runoff + narr_wspd + nldas_soilm,
#                  data= filter(data_region, region == rrr))
# 
# linear_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/final_cmip_region_data_veg.rds")) %>%
#   filter(region == rrr) %>%  as.data.frame()
# linear_data <- linear_data[complete.cases(linear_data),]
# 
# linear_level_cmip <- linear_data %>% mutate(algorithm = "Linear (level), veg")
# linear_level_cmip$pred_DM <- predict(linear_level, newdata = linear_level_cmip)
# 
# linear_log_cmip <- linear_data %>% mutate(algorithm = "Linear (log), veg")
# linear_log_cmip$pred_DM <- predict(linear_log, newdata = linear_log_cmip) %>% exp()
# cmip_pred <- bind_rows(linear_level_cmip, linear_log_cmip) %>%
#   select(algorithm, model, scenario, year, region, pred_DM, 
#          narr_temp, narr_precip, narr_rhum, narr_vpd, narr_soilm, narr_runoff, narr_wspd) %>%
#   rename(gcm = model)
# saveRDS(cmip_pred, paste0(result_path,"/CMIP6_fire_proj/west_pred_final_linear_veg.rds"))
# 
# #------------------------------------------------------------------
# #----------           Drop soil moisture         ----------#
# #------------------------------------------------------------------
# linear_level <- lm(DM_kg ~ narr_temp + narr_precip + narr_rhum + narr_vpd + 
#                      narr_runoff + narr_wspd,
#                    data= filter(data_region, region == rrr))
# 
# linear_log <- lm(log(DM_kg) ~ narr_temp + narr_precip + narr_rhum + narr_vpd + 
#                    narr_runoff + narr_wspd,
#                  data= filter(data_region, region == rrr))
# 
# linear_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/final_cmip_region_data_veg.rds")) %>%
#   filter(region == rrr) %>%  as.data.frame()
# linear_data <- linear_data[complete.cases(linear_data),]
# 
# linear_level_cmip <- linear_data %>% mutate(algorithm = "Linear (level), no soil")
# linear_level_cmip$pred_DM <- predict(linear_level, newdata = linear_level_cmip)
# 
# linear_log_cmip <- linear_data %>% mutate(algorithm = "Linear (log), no soil")
# linear_log_cmip$pred_DM <- predict(linear_log, newdata = linear_log_cmip) %>% exp()
# cmip_pred <- bind_rows(linear_level_cmip, linear_log_cmip) %>%
#   select(algorithm, model, scenario, year, region, pred_DM, 
#          narr_temp, narr_precip, narr_rhum, narr_vpd, narr_soilm, narr_runoff, narr_wspd) %>%
#   rename(gcm = model)
# saveRDS(cmip_pred, paste0(result_path,"/CMIP6_fire_proj/west_pred_final_linear_nosoil_veg.rds"))
# 
# 
# #------------------------------------------------------------------
# #---------- Sensitivity for western US (with upsampling)----------#
# #------------------------------------------------------------------
# rrr <- "Western US"
# xgb_cmip <- readRDS(paste0(result_path,"/CMIP6_fire_proj/final_cmip_eco3_data_normalize.rds")) %>%
#   filter(region == rrr) %>%  as.data.frame() %>% rename(region_sub = NA_L3KEY)
# xgb_cmip <- left_join(xgb_cmip, area_eco3) %>% mutate(algorithm = "XGBOOST, 10x fire")
# xgb_cmip <- xgb_cmip[complete.cases(xgb_cmip),]
# xgb_model <- readRDS(paste0(result_path,"/final_models/west_regional_XGBOOST_norm_final_model_weight10_outcome.rds"))
# 
# xgb_cmip$pred_DM <- predict(xgb_model,
#                             newdata = xgb.DMatrix(data.matrix(xgb_cmip[, xgb_model$feature_names])))  
# xgb_cmip$pred_DM <- xgb_cmip$pred_DM * xgb_cmip$area
# 
# nn_model <- readRDS(paste0(result_path,"/final_models/west_regional_NN_norm_final_model_weight10_outcome.rds"))
# 
# nn_cmip <- readRDS(paste0(result_path,"/CMIP6_fire_proj/final_cmip_region_data_normalize.rds")) %>%
#   filter(region == rrr) %>%  as.data.frame()
# nn_cmip <- left_join(nn_cmip, area_region) %>% mutate(algorithm = "Neural Network, 10x fire")
# nn_cmip <- nn_cmip[complete.cases(nn_cmip),]
# 
# nn_cmip$pred_DM <- predict(nn_model, newdata = nn_cmip) %>% exp()
# nn_cmip$pred_DM <- nn_cmip$pred_DM * nn_cmip$area
# 
# ### Linear models: regional (with upsampling)
# weight <- 10
# data_upsample <- bind_rows(filter(data_region, region == rrr),
#                            filter(data_region, region == rrr, year %in% c(2020,2021)) %>% 
#                              dplyr::slice(rep(1:n(), each = weight)))
# 
# linear_level_weight <- lm(DM_kg ~ narr_temp + narr_precip + narr_rhum + narr_vpd +
#                      narr_runoff + narr_wspd + nldas_soilm,  data= data_upsample)
# 
# linear_log_weight <- lm(log(DM_kg) ~ narr_temp + narr_precip + narr_rhum + narr_vpd +
#                    narr_runoff + narr_wspd + nldas_soilm,  data= data_upsample)
# 
# linear_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/final_cmip_region_data.rds")) %>%
#   filter(region == rrr) %>%  as.data.frame()
# linear_data <- linear_data[complete.cases(linear_data),]
# 
# linear_level_cmip <- linear_data %>% mutate(algorithm = "Linear (level), 10x fire")
# linear_level_cmip$pred_DM <- predict(linear_level_weight, newdata = linear_level_cmip)
# 
# linear_log_cmip <- linear_data %>% mutate(algorithm = "Linear (log), 10x fire")
# linear_log_cmip$pred_DM <- predict(linear_log_weight, newdata = linear_log_cmip) %>% exp()
# 
# lasso_level_weight <- readRDS(paste0(result_path, "/final_models/west_LASSO_level_final_model_weight10_outcome.rds"))
# lasso_log_weight <- readRDS(paste0(result_path, "/final_models/west_LASSO_log_final_model_weight10_outcome.rds"))
# 
# lasso_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/final_cmip_region_data.rds")) %>%
#   filter(region == rrr) %>% as.data.frame() 
# lasso_data <- left_join(lasso_data, area_region)
# lasso_data <- lasso_data[complete.cases(lasso_data),]
# 
# lasso_level_cmip <- lasso_data %>% mutate(algorithm = "LASSO (level), 10x fire")
# lasso_level_cmip$pred_DM <- lasso_level_weight$predict_newdata(newdata = lasso_level_cmip) %>% 
#   as.data.table() %>% select(response) %>% unlist() %>% as.vector() 
# lasso_level_cmip$pred_DM <- lasso_level_cmip$pred_DM * lasso_level_cmip$area
# 
# lasso_log_cmip <- lasso_data %>% mutate(algorithm = "LASSO (log), 10x fire")
# lasso_log_cmip$pred_DM <- lasso_log_weight$predict_newdata(newdata = lasso_log_cmip) %>%
#   as.data.table() %>% select(response) %>% unlist() %>% as.vector() %>% exp() 
# lasso_log_cmip$pred_DM <- lasso_log_cmip$pred_DM * lasso_log_cmip$area
# 
# cmip_pred <- bind_rows(xgb_cmip, nn_cmip, linear_level_cmip, linear_log_cmip, lasso_level_cmip, lasso_log_cmip) %>%
#   select(algorithm, model, scenario, year, region, region_sub, pred_DM, 
#          narr_temp, narr_precip, narr_rhum, narr_vpd, narr_soilm, nldas_soilm, narr_runoff, narr_wspd) %>%
#   rename(gcm = model)
# saveRDS(cmip_pred, paste0(result_path,"/CMIP6_fire_proj/west_pred_final_10xoutcome.rds"))
