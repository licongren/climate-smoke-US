library(ggplot2)
library(dplyr)   
library(tidyverse)
library(MetBrewer)
library(purrr)
library(data.table)
library(mlr3verse)
library(xgboost)
library(future)
library(broom)
library(cowplot)
library(zoo)


rm(list=ls())
gc()

code_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/rcode"
data_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data"
figure_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/figures"
result_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/results/fire_model"

setwd(figure_path)
states_map <- map_data("state")

db_path <- "/Users/mhqiu/BurkeLab Dropbox/Projects/smoke-climate"

#### load the model performance in Western US
linear_results <- readRDS(paste0(result_path,"/model_eval_final/linear_model_results_LOOCV.rds"))
lasso_results <- readRDS(paste0(result_path,"/model_eval_final/lasso_model_results_LOOCV.rds"))
mlp_results <- readRDS(paste0(result_path,"/model_eval_final/mlp_model_results_LOOCV_h2o.rds")) 
results <- bind_rows(linear_results, lasso_results, mlp_results)
best_metric <- read.xlsx(paste0(result_path,"/model_eval_final/best_metrics_rmse_bias_6models_LOOCV.xlsx"))

threshold <- 0.05

results_select <- left_join(filter(results, metric == "10 years"), 
                            filter(best_metric, metric == "10 years"),
                            by=c("region", "model", "outcome", "experiment", "spec", "metric"))  %>%
  filter(!is.na(model_outcome))

##### plot the performance after aggregating over selected models
west_my_model <- results_select %>%
  filter(diff_min<threshold) %>%
  group_by(region, experiment, model, outcome, spec) %>%
  mutate(year = 1:n()) %>%
  group_by(region, year) %>%
  summarise(truth = mean(truth, na.rm=T),
            response_positive = mean(response_positive, na.rm=T)) %>%
  ungroup() %>%
  filter(region == "Western US")

west_metric <- west_my_model %>%
  mutate(r=cor(response_positive, truth),
         rmse=sqrt(mean((response_positive - truth)^2)),
         rmse_mean=rmse/mean(truth)) %>% ungroup() %>%
  arrange(desc(truth)) %>% slice_head(n = 1) %>%
  mutate(bias = (response_positive - truth)/truth) %>%
  mutate(model="Our model", evaluation="LOOCV")

#-----------------------------------------------------------------------------
## Williams and Abatzgolou 2016 (extend to 2021, with gridmet VPD) log-linear
#-----------------------------------------------------------------------------
gfed_wang_region <- readRDS(paste0(data_path, "/crosswalk_GFED_wang_regions.rds"))

gfed_data <- readRDS(file.path(db_path, "climate-fire_training_data_2001-2021_clean_varname.rds"))
gfed_data <- left_join(gfed_data, gfed_wang_region, by="gfed_cell_id")

vlist <- c("barren", "cropland", "forest", "grassland", "lichen", "shrubland", "snowice", "urban", "water","wetland",
           "elevation", "slope",  "narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
           "narr_soilm", "narr_vpd", "narr_runoff",
           "gridmet_erc", "gridmet_fm100","gridmet_fm1000", "gridmet_vpd") 

data_west <- filter(gfed_data, region=="Western US", month>4, month<11) %>% 
  mutate(area=as.numeric(area)) %>%
  group_by(year) %>%
  mutate(DM_kg=sum(DM_kg,na.rm=T)) %>%
  mutate(across(vlist, weighted.mean, w=area*forest)) %>%
  distinct(year, .keep_all = T) %>% 
  dplyr::select(c(vlist, "DM_kg", "temporal_fold")) %>%
  ungroup()

########## John and Abatzgolou 2016 (extend to 2021, with gridmet VPD) 
### in-sample
vpd_log <- feols(log(DM_kg) ~ gridmet_vpd, data=data_west)
pred1 <- data_west %>% 
  distinct(year, DM_kg) %>%
  mutate(DM_kg_log=log(DM_kg),
         DM_pred_log = predict(vpd_log),
         DM_pred = exp(DM_pred_log),
         evaluation="in-sample") 

### LOOCV
num_temporal_folds = length(unique(data_west$year))
outer_resampling = rsmp("cv", folds = num_temporal_folds)
##### Load training data #####
vlist_climate <- c("gridmet_vpd")

  data_folds <- data_west %>%
        mutate(row_ids = row_number()) %>%
        select(row_ids, year, DM_kg)
      
   data_train <- data_west[,c("year", "DM_kg", vlist_climate)] %>%
     mutate(DM_kg_log=log(DM_kg)) %>% select(-DM_kg)
      
      task = as_task_regr(data_train, target = "DM_kg_log")
      
      task$set_col_roles(
        "year", roles = "group")
    
    rr = mlr3::resample(task = task, 
                        learner = lrn("regr.lm"), 
                        resampling = outer_resampling, 
                        store_models = T)
    
    # Extract inner tuning results
    pred2 <- rr$prediction() %>% as.data.table() %>%
      left_join(data_folds, by="row_ids")

    pred2 <- pred2 %>% mutate(DM_pred=exp(response)) %>%
                       mutate(evaluation="LOOCV")
    
    pred_log_linear <- bind_rows(pred1, pred2)
  
    #### using 10-yr moving average for evaluation 
    metric_log_linear <- pred_log_linear %>% group_by(evaluation) %>%
      arrange(year) %>%
      reframe(DM_kg=rollmean(DM_kg, 10, na.rm=T),
              DM_pred=rollmean(DM_pred, 10, na.rm=T)) %>% ungroup() %>%
      group_by(evaluation) %>%
      mutate(r=cor(DM_kg, DM_pred),
             rmse=sqrt(mean((DM_kg - DM_pred)^2)),
             rmse_mean=rmse/mean(DM_kg)) %>% 
      arrange(desc(DM_kg)) %>% slice_head(n = 1) %>%
      mutate(bias = (DM_pred - DM_kg)/DM_kg) %>%
      mutate(model="log(fire) ~ VPD") %>% ungroup()
    
    ggplot(pred_log_linear %>% filter(evaluation=="LOOCV"), 
           aes(x=DM_kg/1e9, y=DM_pred/1e9)) +
      geom_point(size=2) +
      geom_abline(,0,1,linetype="dashed") +
      #labs(title="in linear scale") +
      #facet_wrap(~evaluation) +
      scale_y_continuous(breaks=c(0,50,100,150,200), limits = c(0,200)) + ##, breaks=c(25,100,500,3000)
      scale_x_continuous(breaks=c(0,50,100,150,200), limits = c(0,200)) +
      theme_bw() +
      theme(text = element_text(size=20),axis.text = element_text(size=20)) 
ggsave("pred_log_vpd_level_scale.pdf", width=6, height=4)
    
ggplot(pred_log_linear %>% filter(evaluation=="LOOCV"), 
       aes(x=DM_kg/1e9, y=DM_pred/1e9)) +
  geom_point(size=2) +
  geom_abline(,0,1,linetype="dashed") +
  #labs(title="in log scale") +
  #facet_wrap(~evaluation) +
  scale_y_continuous(trans="log", breaks=c(10,25,50,100,200), limits = c(5,250)) + ##, breaks=c(25,100,500,3000)
  scale_x_continuous(trans="log", breaks=c(10,25,50,100,200), limits = c(5,250)) +
  theme_bw() +
  theme(text = element_text(size=20),axis.text = element_text(size=20)) 
ggsave("pred_log_vpd_log_scale.pdf", width=6, height=4)


########## analyze the model from Wang et al., 2022
replica_wang <- readRDS(paste0(result_path, "/Wang_replica/wang_replica_cv_pred_optimal_hyperparameters.rds")) %>%
  mutate(evaluation="random CV")

pred_wang_loocv <- readRDS(paste0(result_path, "/Wang_replica/wang_temporal_LOOCV_pred_optimal_hyperparameters.rds")) %>%
  mutate(evaluation="LOOCV")

pred_wang <- bind_rows(replica_wang, pred_wang_loocv)
pred_wang <- left_join(pred_wang, gfed_wang_region, by="gfed_cell_id") %>%
  filter(region=="Western US") %>% group_by(evaluation, year) %>%
  summarise(truth=sum(truth),
            response=sum(response*as.numeric(response>0)))

metric_wang2022 <-  pred_wang %>%
  arrange(year) %>%
  reframe(truth=rollmean(truth, 10, na.rm=T),
          response=rollmean(response, 10, na.rm=T)) %>% ungroup() %>%
  group_by(evaluation) %>%
  mutate(r=cor(truth, response),
         rmse=sqrt(mean((truth - response)^2)),
         rmse_mean=rmse/mean(truth)) %>% 
  arrange(desc(truth)) %>% slice_head(n = 1) %>%
  mutate(bias = (response - truth)/truth) %>%
  mutate(model="XGBOOST grid") %>% ungroup()

ggplot(filter(pred_wang, evaluation=="random CV"), 
       aes(y=truth, x=response)) +
  geom_point(size=2) +
  geom_abline(,0,1,linetype="dashed") +
  theme_bw() +
  theme(text = element_text(size=20),axis.text = element_text(size=20)) 
ggsave("westUS_wang_random_CV.pdf", width=6, height=4)

ggplot(filter(metric_wang2022, evaluation=="LOOCV"), 
       aes(y=truth, x=response)) +
  geom_point(size=2) +
  geom_abline(,0,1,linetype="dashed") +
  theme_bw() +
  theme(text = element_text(size=20),axis.text = element_text(size=20)) 
ggsave("westUS_wang_LOOCV.pdf", width=6, height=4)


metric_combined <- bind_rows(west_metric, metric_log_linear, metric_wang2022) %>%
  select( r, rmse, rmse_mean,  bias, model, evaluation) 
write.xlsx(metric_combined, paste0(result_path,"/comapre_fire_models.xlsx"), replace=T)

ggplot(metric_combined)

#### grid-level regression
#pred_grid <- readRDS(paste0(db_path, "/output/GFED_temporal_folds_pred_optimal_hyperparameters.rds")) 

pred_grid <- readRDS(paste0(result_path,"/final/west_nestedCV_annual_grid_results_final.rds"))
tune <- readRDS(paste0(result_path,"/final/west_nestedCV_annual_grid_tuning_final.rds")) 
tune %>% group_by(iteration) %>% arrange(regr.mse) %>% slice_head(n=3) %>% as.data.frame()


#pred_grid <- left_join(pred_grid, gfed_wang_region, by="gfed_cell_id")

pred_region <- pred_grid %>% 
  mutate(response_positive=response*as.numeric(response >= 0)) %>%
  as.data.frame() %>% 
  group_by(region, gfed_cell_id, year, upsample_weight, type) %>%
  summarise(response=mean(response, na.rm=T), 
            truth=mean(truth, na.rm=T),
            response_positive=mean(response_positive, na.rm=T)) %>% ungroup() %>%
  group_by(region, year, upsample_weight, type) %>%
  summarise(truth=sum(truth, na.rm=T),
            response=sum(response, na.rm=T),
            response_positive=sum(response_positive, na.rm=T),) %>% ungroup() %>%
  group_by(region, upsample_weight, type) %>%
mutate(r=cor(response, truth), r_pos=cor(response_positive, truth),
       rmse=sqrt(mean((response - truth)^2)), rmse_pos=sqrt(mean((response_positive - truth)^2))) %>% ungroup()

test_metrics <- pred_region %>% distinct(region, type, upsample_weight, .keep_all = T) %>% 
  dplyr::select(-truth, -response, -response_positive) %>%
  as.data.frame() %>% filter(type=="test")

ggplot(filter(pred_region, region=="Western US", type=="test"), 
       aes(y=truth, x=response)) +
  geom_point() +
  geom_abline(,0,1,linetype="dashed") +
  facet_wrap(~upsample_weight)



# pred1 <- readRDS(paste0(db_path, "/output/temporal_folds/predictions_fold_1.rds")) %>% mutate(hold_fold=1)
# pred2 <- readRDS(paste0(db_path, "/output/temporal_folds/predictions_fold_2.rds")) %>% mutate(hold_fold=2)
# pred3 <- readRDS(paste0(db_path, "/output/temporal_folds/predictions_fold_3.rds")) %>% mutate(hold_fold=3)
# pred4 <- readRDS(paste0(db_path, "/output/temporal_folds/predictions_fold_4.rds")) %>% mutate(hold_fold=4)
# pred5 <- readRDS(paste0(db_path, "/output/temporal_folds/predictions_fold_5.rds")) %>% mutate(hold_fold=5)
# pred6 <- readRDS(paste0(db_path, "/output/temporal_folds/predictions_fold_6.rds")) %>% mutate(hold_fold=6)
# pred7 <- readRDS(paste0(db_path, "/output/temporal_folds/predictions_fold_7.rds")) %>% mutate(hold_fold=7)
# 
# pred_wang_temporal <- list(pred1, pred2, pred3, pred4, pred5, pred6, pred7) %>% bind_rows()
# saveRDS(pred_wang_temporal, paste0(db_path, "/output/wang_temporal_fold.rds"))
