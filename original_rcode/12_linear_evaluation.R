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

#-------------------------------------------------------------------------------
# Analyze the evaluation performance of simple linear models
# Written by Minghao
# Last edited July 2022
#-------------------------------------------------------------------------------
###### calculate the correlation and biases for the top three years
test1 <- readRDS(paste0(result_path,"/model_eval_final/Linear/linear_model_log_loocv.rds"))   %>% 
  mutate(outcome = "log", grid_1deg_id = as.numeric(grid_1deg_id),
         spec ="non-weighted FE")
test2 <- readRDS(paste0(result_path,"/model_eval_final/Linear/linear_model_log_weighted_loocv.rds")) %>% 
  mutate(outcome = "log", grid_1deg_id = as.numeric(grid_1deg_id),
         spec ="weighted FE")
test3 <- readRDS(paste0(result_path,"/model_eval_final/Linear/linear_model_level_loocv.rds")) %>%
  mutate(outcome = "level", grid_1deg_id = as.numeric(grid_1deg_id),
         spec ="non-weighted FE")
test4 <- readRDS(paste0(result_path,"/model_eval_final/Linear/linear_model_level_weighted_loocv.rds")) %>%
  mutate(outcome = "level", grid_1deg_id = as.numeric(grid_1deg_id),
         spec ="weighted FE")

test <- bind_rows(test1, test2, test3, test4) %>%
  mutate(model = "Linear")

annual_summ <- test %>% 
  mutate(response_positive=response*as.numeric(response >= 0)) %>%
  as.data.frame() %>% 
  group_by(region, year, experiment, model, outcome, spec) %>%
  summarise(truth=sum(truth, na.rm=T),
            response=sum(response, na.rm=T),
            response_positive=sum(response_positive, na.rm=T)) %>% ungroup() %>%
  mutate(metric="1 year")

metric_1yr <- annual_summ %>% group_by(region, experiment, model, outcome, spec) %>%
  mutate(r=cor(response, truth), r_pos=cor(response_positive, truth),
         rmse=sqrt(mean((response - truth)^2)), 
         rmse_pos=sqrt(mean((response_positive - truth)^2)),
         truth_mean=mean(truth),
         slope=lm(response_positive~truth)$coefficients[2]) %>% 
  arrange(desc(truth)) %>% slice_head(n = 1) %>%
  mutate(bias = (response_positive - truth)/truth) %>%
  mutate(metric="1 year")

## 3-yr
annual_sum_3yr <- annual_summ %>% group_by(region, experiment, model, outcome, spec) %>%
  arrange(year) %>%
  reframe(truth=rollmean(truth, 3, na.rm=T),
          response=rollmean(response, 3, na.rm=T),
          response_positive=rollmean(response_positive, 3, na.rm=T)) %>% ungroup()%>%
  mutate(metric="3 years")

metric_3yr <- annual_sum_3yr %>% group_by(region, experiment, model, outcome, spec) %>%
  mutate(r=cor(response, truth), r_pos=cor(response_positive, truth),
         rmse=sqrt(mean((response - truth)^2)), 
         rmse_pos=sqrt(mean((response_positive - truth)^2)),
         truth_mean=mean(truth),
         slope=lm(response_positive~truth)$coefficients[2]) %>% 
  arrange(desc(truth)) %>% slice_head(n = 1) %>%
  mutate(bias = (response_positive - truth)/truth) %>%
  mutate(metric="3 years")

## 5-yr
annual_sum_5yr <- annual_summ %>% group_by(region, experiment, model, outcome, spec) %>%
  arrange(year) %>%
  reframe(truth=rollmean(truth, 5, na.rm=T),
          response=rollmean(response, 5, na.rm=T),
          response_positive=rollmean(response_positive, 5, na.rm=T)) %>%
  mutate(metric="5 years")

metric_5yr <- annual_sum_5yr %>% group_by(region, experiment, model, outcome, spec) %>%
  mutate(r=cor(response, truth), r_pos=cor(response_positive, truth),
         rmse=sqrt(mean((response - truth)^2)), 
         rmse_pos=sqrt(mean((response_positive - truth)^2)),
         truth_mean=mean(truth),
         slope=lm(response_positive~truth)$coefficients[2]) %>% 
  arrange(desc(truth)) %>% slice_head(n = 1) %>%
  mutate(bias = (response_positive - truth)/truth) %>%
  mutate(metric="5 years")

## 10-yr
annual_sum_10yr <- annual_summ %>% group_by(region, experiment, model, outcome, spec) %>%
  arrange(year) %>%
  reframe(truth=rollmean(truth, 10, na.rm=T),
          response=rollmean(response, 10, na.rm=T),
          response_positive=rollmean(response_positive, 10, na.rm=T)) %>%
  mutate(metric="10 years")

metric_10yr <- annual_sum_10yr %>% group_by(region, experiment, model, outcome, spec) %>%
  mutate(r=cor(response, truth), r_pos=cor(response_positive, truth),
         rmse=sqrt(mean((response - truth)^2)),
         rmse_pos=sqrt(mean((response_positive - truth)^2)),
         truth_mean=mean(truth),
         slope=lm(response_positive~truth)$coefficients[2]) %>% 
  arrange(desc(truth)) %>% slice_head(n = 1) %>%
  mutate(bias = (response_positive - truth)/truth) %>%
  mutate(metric="10 years") %>% ungroup()

metric_combined <- bind_rows(metric_1yr, metric_3yr, metric_5yr, metric_10yr) %>%
  mutate(rmse_pos_mean = rmse_pos / truth_mean)
# metric_combined$experiment <- factor(metric_combined$experiment, 
#                                      levels=c("grid","eco3","eco2","regional")) 
metric_combined$metric <- factor(metric_combined$metric, 
                                 levels=c("1 year","3 years","5 years","10 years")) 
metric_combined <- metric_combined  %>% arrange(metric, experiment)
write.xlsx(metric_combined, 
           paste0(result_path,"/model_eval_final/linear_model_metrics_combined_LOOCV.xlsx"),
           replace = T)

annual_summ_combined <- bind_rows(annual_summ, annual_sum_3yr, annual_sum_5yr, annual_sum_10yr)
annual_summ_combined$metric <- factor(annual_summ_combined$metric, 
                                      levels=c("1 year","3 years","5 years","10 years"))
# annual_summ_combined$experiment <- factor(annual_summ_combined$experiment, 
#                                           levels=c("grid","eco3","eco2","regional")) 
saveRDS(annual_summ_combined,
        paste0(result_path,"/model_eval_final/linear_model_results_LOOCV.rds"))

rrr <- "Western US"
### understand difference across metrics


ggplot(filter(annual_summ_combined, metric == "10 years", region == rrr, response_positive<1e11,
              !model %in% c("log", "level")),
       aes(x=truth/1e9, y=response_positive/1e9)) +
  geom_point(size=4, alpha=0.6) +
  geom_abline(,0,1,linetype="dashed") +
  theme_classic() + theme(text = element_text(size=16)) + 
  labs(y="Predictions (million tons DM)", x="Observations (million tons DM)") +
  facet_grid(experiment~model) 
ggsave("westUS_linear_models_soil.png", width=11, height=7.5)

ggplot(filter(annual_summ_combined, model=="log, weighted FE", experiment == "eco3", region == rrr),
       aes(x=truth/1e9, y=response_positive/1e9)) +
  geom_point(size=4, alpha=0.6) +
  geom_abline(,0,1,linetype="dashed") +
  theme_classic() + theme(text = element_text(size=16)) +
  labs(y="Predictions (million tons DM)", x="Observations (million tons DM)") +
  facet_wrap(~metric, scales = "free") 
ggsave("westUS_linear_log_eco3_metrics.png", width=8, height=5)#+

ggplot(filter(annual_summ_combined, model=="log, weighted FE", metric == "10 years", region == rrr),
       aes(x=truth/1e9, y=response_positive/1e9)) +
  geom_point(size=4, alpha=0.6) +
  geom_abline(,0,1,linetype="dashed") +
  theme_classic() + theme(text = element_text(size=16)) + 
  labs(y="Predictions (million tons DM)", x="Observations (million tons DM)") +
  facet_wrap(~experiment) 
ggsave("westUS_linear_level_spatial_res.png", width=8, height=3)#+

ggplot(filter(annual_summ_combined, experiment=="eco2", metric == "10 years", region == rrr),
       aes(x=truth/1e9, y=response_positive/1e9)) +
  geom_point(size=4, alpha=0.6) +
  geom_abline(,0,1,linetype="dashed") +
  theme_classic() + theme(text = element_text(size=16)) + 
  labs(y="Predictions (million tons DM)", x="Observations (million tons DM)") +
  facet_wrap(~model) 
ggsave("westUS_linear_level_eco2_model_choice.png", width=8, height=5)#+

#### examine performances across subregions

region_list <-  c("Western US", "Canada-Alaska","Mexico", "Southeastern US", "Northeastern US")
rrr <- "Western US"

ggplot(filter(metric_combined), 
       aes(x = experiment, y = r_pos, fill = metric)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_classic() + theme(text = element_text(size=16))  +
  labs(x="", y="Correlation coef (R)") +
  facet_grid(model~region )
ggsave("linear_models_performance.png", width=12, height=5)

ggplot(filter(metric_combined, bias<1), 
       aes(x = experiment, y = bias, fill = metric)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_classic() + theme(text = element_text(size=16))  +
  labs(x="", y="Correlation coef (R)") +
  facet_grid(model~region )

# gfed_grid <- readRDS(paste0(data_path, "/GFED4S_NorthAmerica_025deg_grid_EPSG4326.rds"))
# 
# gfed_ecoregion_link <- readRDS(paste0(db_path, "/crosswalk_GFED_NA_CEC_Eco_Level2_3.rds")) %>%
#   mutate(intersection_area=as.numeric(intersection_area)) %>%   group_by(gfed_cell_id)  %>%
#   arrange(desc(intersection_area)) %>%
#   mutate(area_weight=intersection_area/sum(intersection_area)) %>% ungroup() # %>% distinct(gfed_cell_id, .keep_all=T)
# 
# gfed_wang_region <- readRDS(paste0(data_path, "/crosswalk_GFED_wang_regions.rds")) 
# 
# gfed_data <- readRDS(file.path(db_path, "climate-fire_training_data_2001-2021_clean_varname.rds"))
# gfed_data <- left_join(gfed_data, gfed_wang_region, by="gfed_cell_id")
# gfed_data <- left_join(gfed_data, gfed_ecoregion_link, by="gfed_cell_id")
# 
# gfed_data$area <- as.numeric(gfed_data$area)
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
# gfed_agg_data <- readRDS(paste0(data_path, "/gfed_agg_data.rds"))
# 
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
# ### CV linear model
# num_temporal_folds = length(unique(gfed_data$temporal_fold))
# outer_resampling = rsmp("cv", folds = num_temporal_folds)
# 
# region_list <- unique(gfed_data$region)
# test <- NULL 
# LOG <- F
# area_normalize <- T
# for (experiment in c("eco3", "eco2", "grid")){ 
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
#   for (rrr in region_list){
#     print(rrr)
#     data_train <- data %>% filter(region == rrr) %>% select(-DM_kg)
#     
#     if (rrr %in% c("Canada-Alaska","Mexico")){
#       data_train <- data_train %>% select(-nldas_soilm) 
#     }
#     if (LOG){
#       data_train <- data_train %>% filter(outcome>0) %>% mutate(outcome=log(outcome))
#     }else{
#       data_train <- data_train 
#     }
#     
#     task = as_task_regr(data_train, target = "outcome")
#     
#     if (experiment =="grid"){
#       data_folds <- data_train %>% 
#         mutate(row_ids = row_number()) %>% 
#         select(row_ids, year, temporal_fold, region, gfed_cell_id, area)
#       
#       if (region_FE){
#         task$set_col_roles(
#           "temporal_fold", 
#           roles = "group")$set_col_roles(
#             c("year","region","area"), 
#             remove_from = "feature")
#       } else{
#         task$set_col_roles(
#           "temporal_fold", 
#           roles = "group")$set_col_roles(
#             c("year","region", "gfed_cell_id","area"), 
#             remove_from = "feature")  
#       }
#     } else{
#       data_folds <- data_train %>% 
#         mutate(row_ids = row_number()) %>% 
#         select(row_ids, year, temporal_fold, region, region_sub, area)
#       
#       if (region_FE){
#         task$set_col_roles(
#           "temporal_fold", 
#           roles = "group")$set_col_roles(
#             c("year","region","area"), 
#             remove_from = "feature")
#       } else{
#         task$set_col_roles(
#           "temporal_fold", 
#           roles = "group")$set_col_roles(
#             c("year","region","region_sub","area"), 
#             remove_from = "feature")       
#       }
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
#     if (LOG){
#       predictions <- predictions %>% mutate(truth=exp(truth),
#                                             response=exp(response),
#                                             model="log") 
#     }else{
#       predictions$model <- "level" 
#     }
#     
#     predictions <- predictions %>% 
#       mutate(truth = truth * area,
#              response = response * area,
#              region=rrr, experiment = experiment)
#     test <- bind_rows(test, predictions)
#     
#   }
# }
# 
# saveRDS(test, paste0(result_path,"/model_eval_final/linear_model_level_FE.rds"))

