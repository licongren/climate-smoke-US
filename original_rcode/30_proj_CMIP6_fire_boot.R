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
# Project fire emissions using different models of 40 runs
# Written by Minghao
# Last edited Jan 2024
#-------------------------------------------------------------------------------
###---- Load the area of different regions to be used for predictions and training purposes
region_list <- c("Western US", "Canada-Alaska", "Mexico", "Southeastern US", "Northeastern US") ##unique(gfed_data$region)
region_label_list <- c("west", "canada", "mexico", "southeast", "northeast")

best_models <- read.xlsx(paste0(result_path,"/model_eval_final/best_metrics_rmse_bias_6models_LOOCV.xlsx"))
#best_models1 <- read.xlsx(paste0(result_path,"/model_eval_final/best_metrics_rmse_bias_6models_subregion_LOOCV.xlsx"))

best_models <- best_models %>% 
  filter(selected==1, metric %in% c("5 years", "10 years")) %>%
  distinct(region, model, experiment, outcome, spec, .keep_all = T)

#### to get the complete GCM list 
cmip_data <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_eco3_full_demean_level.rds")) %>%
  filter(year == 1997) %>%
  filter_at(vars(narr_temp:narr_wspd), all_vars(!is.na(.)))

gcm_list <- unique(cmip_data$model)

nboot <- 40
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
    
    if (experiment == "regional"){
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
        mutate(region_sub = as.character(region_sub))
      
      history_data <- history_data_grid %>%
        rename(region_sub = grid_1deg_id) %>%
        mutate(region_sub = as.character(region_sub))
    }
  }
  
  if (outcome == "log"){
    history_data <- readRDS(paste0(result_path,"/final_models/training_data/final_train_region_FE_log.rds"))
    history_data_grid <- readRDS(paste0(result_path,"/final_models/training_data/final_train_region_FE_log_1deg_grid.rds"))
    
    if (experiment == "regional"){
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
  }
    
  cmip_proj <- cmip_data %>%
    filter(region == region_tmp)  %>%
    filter(model %in% gcm_list) %>%
    filter_at(vlist_climate, all_vars(!is.na(.))) %>%
    filter(scenario != "ssp585" & year<2056) ## save memory
  
  for (bbb in 1:nboot){    
    print(bbb)
    
    if (spec == "weighted FE"){
      model_tmp <- readRDS(paste0(result_path, "/final_models/bootstrap/linear_",outcome,"_weighted_final_models_boot.rds"))[[paste(region_label,experiment,bbb,sep="_")]]
    } else{
      model_tmp <- readRDS(paste0(result_path, "/final_models/bootstrap/linear_",outcome,"_final_models_boot.rds"))[[paste(region_label,experiment,bbb,sep="_")]]
    }
    
    cmip_proj_tmp <- cmip_proj
    cmip_proj_tmp$pred_DM <- predict(model_tmp, newdata = cmip_proj_tmp)
    
    history_data_tmp <- history_data %>% filter(region == region_tmp) 
    history_data_tmp$pred_DM <- predict(model_tmp, newdata = history_data_tmp)
    
    if (outcome == "level"){
    cmip_proj_tmp <- cmip_proj_tmp %>% mutate(pred_DM = (pred_DM + outcome_mean_level) * area)  %>%
      select(model, scenario, year, region, region_sub, pred_DM) %>%
      mutate(algorithm = "Linear", spatial_res=experiment, outcome=outcome, spec=spec) 
  
    history_data_tmp <- history_data_tmp %>% mutate(pred_DM = (pred_DM + outcome_mean) * area) %>%
      select(year, region, region_sub, pred_DM) %>%
      mutate(algorithm = "Linear", spatial_res=experiment, outcome=outcome, spec=spec, scenario="observation") 
  }
  
    if (outcome == "log"){
    cmip_proj_tmp <- cmip_proj_tmp %>% mutate(pred_DM = exp(pred_DM + outcome_mean_log) * area) %>%
      select(model, scenario, year, region, region_sub, pred_DM) %>%
      mutate(algorithm = "Linear", spatial_res=experiment, outcome=outcome, spec=spec) 
    
    history_data_tmp <- history_data_tmp %>% mutate(pred_DM = exp(pred_DM + outcome_mean) * area) %>%
      select(year, region, region_sub, pred_DM) %>%
      mutate(algorithm = "Linear", spatial_res=experiment, outcome=outcome, spec=spec, scenario="observation") 
  }
    
    cmip_proj_tmp$bootid <-  bbb
    history_data_tmp$bootid <-  bbb
    
  linear_cmip_proj <- bind_rows(linear_cmip_proj,  cmip_proj_tmp) 
  linear_history_proj <- bind_rows(linear_history_proj,  history_data_tmp) 
  }
}

saveRDS(linear_cmip_proj, paste0(result_path,"/CMIP6_fire_proj/linear_final_pred_boot.rds"))
saveRDS(linear_history_proj, paste0(result_path,"/CMIP6_fire_proj/linear_final_pred_history_boot.rds"))

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
    
    if (experiment == "regional"){
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
        mutate(region_sub = as.character(region_sub)) 
      
      history_data <- history_data_grid %>%
        rename(region_sub = grid_1deg_id) %>%
        mutate(region_sub = as.character(region_sub))
    }
         }
  
  if (outcome == "log"){
    history_data <- readRDS(paste0(result_path,"/final_models/training_data/final_train_region_FE_log.rds"))
    history_data_grid <- readRDS(paste0(result_path,"/final_models/training_data/final_train_region_FE_log_1deg_grid.rds"))
    
    if (experiment == "regional"){
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
        mutate(region_sub = as.character(region_sub)) 
      
      history_data <- history_data_grid %>%
        rename(region_sub = grid_1deg_id) %>%
        mutate(region_sub = as.character(region_sub))
    }
      }
  
  cmip_proj <- cmip_data %>%
    filter(region == region_tmp)  %>%
    filter(model %in% gcm_list) %>%
    filter_at(vlist_climate, all_vars(!is.na(.)))%>%
    filter(scenario != "ssp585" & year<2056) ## save memory
  
  for (bbb in 1:nboot){    
    print(bbb)
    
    if (spec == "weighted FE"){
      model_tmp <- readRDS(paste0(result_path, "/final_models/bootstrap/LASSO_",outcome,"_weighted_final_models_boot.rds"))[[paste(region_label,experiment,bbb,sep="_")]]
    } else{
      model_tmp <- readRDS(paste0(result_path, "/final_models/bootstrap/LASSO_",outcome,"_final_models_boot.rds"))[[paste(region_label,experiment,bbb,sep="_")]]
    }
    
    coef <- model_tmp$model$beta  %>% as.matrix() %>% as.data.frame() %>% filter(s0!=0)
    #print(coef)
    if (nrow(coef)==0) {next}
    
    cmip_proj_tmp <- cmip_proj
    cmip_proj_tmp$pred_DM <- model_tmp$predict_newdata(newdata = cmip_proj_tmp) %>%
    as.data.table() %>% select(response) %>% unlist() %>% as.vector()
  
  history_data_tmp <- history_data %>% filter(region == region_tmp) 
  history_data_tmp$pred_DM <- model_tmp$predict_newdata(newdata = history_data_tmp) %>%
    as.data.table() %>% select(response) %>% unlist() %>% as.vector()
  
  if (outcome == "level"){
  cmip_proj_tmp <- cmip_proj_tmp %>% mutate(pred_DM = (pred_DM + outcome_mean_level) * area)  %>%
    select(model, scenario, year, region, region_sub, pred_DM) %>%
    mutate(algorithm = "LASSO", spatial_res=experiment, outcome=outcome, spec=spec)
  
    history_data_tmp <- history_data_tmp %>% mutate(pred_DM = (pred_DM + outcome_mean) * area)  %>%
    select(year, region, region_sub, pred_DM) %>%
    mutate(algorithm = "LASSO", spatial_res=experiment, outcome=outcome, spec=spec, scenario="observation")
  }
  
  if (outcome == "log"){
    cmip_proj_tmp <- cmip_proj_tmp %>% mutate(pred_DM = exp(pred_DM + outcome_mean_log) * area) %>%
      select(model, scenario, year, region, region_sub, pred_DM) %>%
      mutate(algorithm = "LASSO", spatial_res=experiment, outcome=outcome, spec=spec)
    
    history_data_tmp <- history_data_tmp %>% mutate(pred_DM = exp(pred_DM + outcome_mean) * area) %>%
      select(year, region, region_sub, pred_DM) %>%
      mutate(algorithm = "LASSO", spatial_res=experiment, outcome=outcome, spec=spec, scenario="observation")
  }
  
  cmip_proj_tmp$bootid <-  bbb
  history_data_tmp$bootid <-  bbb
  
  lasso_cmip_proj <- bind_rows(lasso_cmip_proj,  cmip_proj_tmp)
  lasso_history <- bind_rows(lasso_history,  history_data_tmp)
  }
}

saveRDS(lasso_cmip_proj, paste0(result_path,"/CMIP6_fire_proj/LASSO_final_pred_boot.rds"))
saveRDS(lasso_history, paste0(result_path,"/CMIP6_fire_proj/LASSO_final_pred_history_boot.rds"))


#-------------------------------------------------------
## Form the final predictions across the different bootids
#-------------------------------------------------------
pred0 <- readRDS(paste0(result_path,"/CMIP6_fire_proj/final_pred_best_models_threshold_05.rds"))
gcm_list <- unique(pred0$gcm)

### for linear/LASSO we include models with the different bootids
pred_linear <- readRDS(paste0(result_path,"/CMIP6_fire_proj/linear_final_pred_boot.rds"))
pred_lasso <- readRDS(paste0(result_path,"/CMIP6_fire_proj/LASSO_final_pred_boot.rds"))
pred_nn <- readRDS(paste0(result_path,"/CMIP6_fire_proj/nnet_final_pred.rds")) %>% mutate(year=as.numeric(as.character(year)))

#### Load the best models
best_models <- read.xlsx(paste0(result_path,"/model_eval_final/best_metrics_rmse_bias_6models_LOOCV.xlsx")) %>% 
  filter(selected==1, metric %in% c("10 years"), diff_min < 0.05) %>%
  distinct(region, model, experiment, outcome, model_outcome, spec, selected) %>%
  rename(algorithm = model, spatial_res=experiment) 

### combine with the model performance to select the best-performed models
pred_linear <- pred_linear %>%
  left_join(best_models, by=c("region", "algorithm", "outcome", "spatial_res", "spec")) %>%
  mutate(pred_DM=pred_DM*as.numeric(pred_DM >= 0)) %>%
  filter(!is.na(selected))

pred_lasso <- pred_lasso %>%
  left_join(best_models, by=c("region", "algorithm", "outcome", "spatial_res", "spec")) %>%
  mutate(pred_DM=pred_DM*as.numeric(pred_DM >= 0)) %>%
  filter(!is.na(selected))

pred_nn <- pred_nn %>%
  left_join(best_models, by=c("region", "algorithm", "outcome", "spatial_res", "spec")) %>%
  mutate(pred_DM=pred_DM*as.numeric(pred_DM >= 0)) %>%
  filter(!is.na(selected)) %>%
  filter(scenario != "ssp585" & year<2056)
  

### add the deltas to the mean of observed emissions between 2001 to 2021
threshold <- 20
gfed_region <- readRDS(paste0(result_path, "/final_models/training_data/final_train_region_FE_level.rds"))$data_region %>%
  group_by(region, region_sub) %>%
  summarise(DM_2001_2021 = mean(DM_kg, na.rm=T),
            max_threshold = max(DM_kg) * threshold) %>% ungroup()

gfed_eco2 <- readRDS(paste0(result_path, "/final_models/training_data/final_train_region_FE_level.rds"))$data_eco2 %>%
  group_by(region, region_sub) %>%
  summarise(DM_2001_2021 = mean(DM_kg, na.rm=T),
            max_threshold = max(DM_kg) * threshold) %>% ungroup()

gfed_eco3 <- readRDS(paste0(result_path, "/final_models/training_data/final_train_region_FE_level.rds"))$data_eco3 %>%
  group_by(region, region_sub) %>%
  summarise(DM_2001_2021 = mean(DM_kg, na.rm=T),
            max_threshold = max(DM_kg) * threshold) %>% ungroup()

gfed_grid <- readRDS(paste0(result_path, "/final_models/training_data/final_train_region_FE_level_1deg_grid.rds")) %>%
  mutate(region_sub = as.character(grid_1deg_id)) %>%
  group_by(region, region_sub) %>%
  summarise(DM_2001_2021 = mean(DM_kg, na.rm=T),
            max_threshold = max(DM_kg) * threshold) %>% ungroup()

gfed_combined <- bind_rows(gfed_region, gfed_eco2, gfed_eco3, gfed_grid)  

pred_combined <- NULL 
nboot <- max(pred_lasso$bootid)
for (bbb in 1:nboot){
print(bbb)  
linear_tmp <- filter(pred_linear, bootid==bbb)  
lasso_tmp <- filter(pred_lasso, bootid==bbb)  

cmip_pred <- bind_rows(linear_tmp, lasso_tmp, pred_nn) %>% 
  rename(gcm = model) %>% filter(gcm %in% gcm_list) %>%
  select(-bootid)
cmip_pred$algorithm_outcome <- paste(cmip_pred$algorithm, cmip_pred$outcome, sep=", ")

### use 2001-2021 as the normalization baselines
cmip_2001_2021 <- filter(cmip_pred, scenario%in%c("historical", "ssp370"), year %in% seq(2001,2021)) 
cmip_2001_2021_mean <- cmip_2001_2021 %>%
  group_by(gcm, region, region_sub, algorithm_outcome) %>%
  summarise(DM_2001_2021 = mean(pred_DM, na.rm=T)) %>% ungroup()

#### calculate the delta differences within gcm projections and then add onto real obs
cmip_pred <-  bind_rows(filter(cmip_pred, year>2021),
                             cmip_2001_2021 %>% mutate(scenario = "historical"),
                             cmip_2001_2021 %>% mutate(scenario = "ssp126"),
                             cmip_2001_2021 %>% mutate(scenario = "ssp245"),
                             cmip_2001_2021 %>% mutate(scenario = "ssp370")) %>%
  left_join(cmip_2001_2021_mean) %>%
  mutate(DM_anomaly = pred_DM - DM_2001_2021) %>% select(-DM_2001_2021)

cmip_pred_tmp <- left_join(cmip_pred, gfed_combined) %>%
  mutate(DM_level = DM_anomaly + DM_2001_2021) %>%
  mutate(DM_level = pmin(DM_level, max_threshold)) %>%
  rename(pred_DM_final=DM_level) %>%
  group_by(scenario, region, region_sub, gcm, algorithm_outcome) %>%
  arrange(year) %>%
  mutate(pred_final_10yr=rollapplyr(pred_DM_final, width = 10, FUN=mean, partial=T)) %>%  ## calculating the 10-yr mean
  filter(scenario=="historical" | year>2021) %>% ungroup() %>%
  mutate(pred_final_10yr_positive=pmax(0, pred_final_10yr)) 

cmip_pred_tmp <- filter(cmip_pred_tmp, year==2055) %>% mutate(bootid=bbb)
# cmip_region_tmp <- cmip_pred_tmp %>%
#   group_by(year, scenario, gcm, region, algorithm_outcome) %>%
#   summarise_at(c("pred_final_10yr_positive"), .funs=sum, na.rm=T) %>% ungroup()
# 
# cmip_region_tmp <- cmip_region_tmp %>%
#   # group_by(scenario, year, region, gcm) %>%
#   # summarise_at(c("pred_final_10yr_positive"), .funs=mean, na.rm=T) %>% ungroup() %>% ## take mean across algorithms
#   # group_by(scenario, year, region) %>% 
#   # summarise(pred_mean=mean(pred_final_10yr_positive, na.rm=T)) %>% ## take mean across GCMs
#   # ungroup() %>%
#   mutate(bootid=bbb)

pred_combined <- bind_rows(pred_combined, cmip_pred_tmp)
}
saveRDS(pred_combined, paste0(result_path,"/CMIP6_fire_proj/pred_fire_model_uncertainty_boot_2050.rds"))


#---------------------------------------------
#  Plot the predictions across different boot runs
#---------------------------------------------
pred_boot <- readRDS(paste0(result_path,"/CMIP6_fire_proj/pred_fire_model_uncertainty_boot1.rds")) %>%
  group_by(year, scenario, gcm,region, bootid) %>%
  summarise_at(c("pred_final_10yr_positive"),.funs=mean, na.rm=T) %>% ungroup() %>% ##mean across algortihms
  group_by(year, scenario, region, bootid) %>%
  summarise_at(c("pred_final_10yr_positive"),.funs=mean, na.rm=T) %>% ungroup() ##mean across GCMs
saveRDS(pred_boot, paste0(result_path,"/CMIP6_fire_proj/pred_fire_model_uncertainty_boot_summ.rds"))

  # group_by(scenario, year,region, algorithm_outcome) %>%
  # summarise(pred_median=median(pred_final_10yr_positive),
  #           pred_mean=mean(pred_final_10yr_positive),
  #           pred_10=quantile(pred_final_10yr_positive,0.1),
  #           pred_90=quantile(pred_final_10yr_positive,0.9)) %>%
  # group_by(scenario, year,region) %>%
  # summarise(pred_mean=mean(pred_mean)) %>%
  ungroup()

cmip_pred_final <- readRDS(paste0(result_path,"/CMIP6_fire_proj/final_pred_best_models_threshold_05.rds")) %>%
  group_by(year, scenario, gcm, region, algorithm_outcome) %>%
  summarise_at(c("pred_final_10yr_positive"),.funs=sum, na.rm=T) %>% ungroup() %>%
  # group_by(scenario, year, region, gcm) %>%
  # summarise_at(c("pred_final_10yr_positive"), .funs=mean, na.rm=T) %>% ungroup() %>%
  group_by(scenario, year, region, algorithm_outcome) %>% 
  summarise(pred_mean=mean(pred_final_10yr_positive, na.rm=T)) %>% ungroup() 


  ggplot() +
    geom_line(data=cmip_pred_final %>% filter(scenario=="ssp370", year<2056),
              aes(x=year, y= pred_mean*1e-9), size=1.4) +
    geom_line(data=pred_boot %>% filter(scenario=="ssp370"),
                aes(x=year, y= pred_mean*1e-9), colour="red") +
    geom_ribbon(data=pred_boot %>% filter(scenario=="ssp370"),
              aes(x=year, ymin= pred_10*1e-9, ymax=pred_90*1e-9), alpha=0.4, fill="salmon") +
    theme_classic() + theme(text = element_text(size=14),
                            axis.text.x = element_text(size=14)) +
    labs(x="", fill="", colour="", y= "Predicted emissions (MT)") +
    scale_y_continuous(limits= range) +
    scale_x_continuous(breaks=seq(2000, 2050, by=10)) +
    facet_wrap(~region+algorithm_outcome, scales = "free")
  ggsave("fire_model_uncertainty_ssp370.png", width=7, height=5)



