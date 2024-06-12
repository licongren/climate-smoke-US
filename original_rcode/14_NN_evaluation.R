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
# Analyze the evaluation performance of single-NN models
# Written by Minghao
# Last edited Sep 2023
#-------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
##    Analyze model performance with different averaging windows (from model that are trained on 7-fold)
#--------------------------------------------------------------------------------

#------------------
## MLP models
#------------------
mlp1 <- readRDS(paste0(result_path,"/model_eval_final/NN/mlp_log_weighted_2lyr_simple_loocv_h2o.rds")) %>%
  mutate(model = "Neural Net", outcome="log", spec="simple")
mlp2 <- readRDS(paste0(result_path,"/model_eval_final/NN/mlp_level_weighted_2lyr_simple_loocv_h2o.rds")) %>%
  mutate(model = "Neural Net", outcome="level", spec="simple")
mlp3 <- readRDS(paste0(result_path,"/model_eval_final/NN/mlp_log_weighted_2lyr_full_loocv_h2o.rds")) %>%
  mutate(model = "Neural Net", outcome="log", spec="full")
mlp4 <- readRDS(paste0(result_path,"/model_eval_final/NN/mlp_level_weighted_2lyr_full_loocv_h2o.rds")) %>%
  mutate(model = "Neural Net", outcome="level", spec="full")
mlp5 <- readRDS(paste0(result_path,"/model_eval_final/NN/mlp_level_weighted_2lyr_simple_loocv_h2o_grid.rds")) %>%
  mutate(model = "Neural Net", outcome="level", spec="simple")
mlp6 <- readRDS(paste0(result_path,"/model_eval_final/NN/mlp_log_weighted_2lyr_simple_loocv_h2o_grid.rds")) %>%
  mutate(model = "Neural Net", outcome="log", spec="simple")
mlp7 <- readRDS(paste0(result_path,"/model_eval_final/NN/mlp_level_weighted_2lyr_full_loocv_h2o_grid_west.rds")) %>%
  mutate(model = "Neural Net", outcome="level", spec="full")
mlp8 <- readRDS(paste0(result_path,"/model_eval_final/NN/mlp_level_weighted_2lyr_full_loocv_h2o_grid_other.rds")) %>%
  mutate(model = "Neural Net", outcome="level", spec="full")
mlp9 <- readRDS(paste0(result_path,"/model_eval_final/NN/mlp_level_weighted_2lyr_full_loocv_h2o_grid_can.rds")) %>%
  mutate(model = "Neural Net", outcome="level", spec="full")
mlp10 <- readRDS(paste0(result_path,"/model_eval_final/NN/mlp_log_weighted_2lyr_full_loocv_h2o_grid_west_other.rds")) %>%
  mutate(model = "Neural Net", outcome="log", spec="full")
mlp11 <- readRDS(paste0(result_path,"/model_eval_final/NN/mlp_log_weighted_2lyr_full_loocv_h2o_grid_can.rds")) %>%
  mutate(model = "Neural Net", outcome="log", spec="full")

mlp <- bind_rows(mlp1, mlp2, mlp3, mlp4, mlp5, mlp6, 
                 mlp7, mlp8, mlp9, mlp10, mlp11)

####### combine all regions
annual_summ <- mlp  %>%
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
metric_combined$experiment <- factor(metric_combined$experiment, 
                                     levels=c("grid","eco3","eco2","regional")) 
metric_combined$metric <- factor(metric_combined$metric, 
                                 levels=c("1 year","3 years","5 years","10 years")) 
metric_combined <- metric_combined  %>% arrange(metric, experiment)
write.xlsx(metric_combined, 
           paste0(result_path,"/model_eval_final/mlp_model_metrics_combined_LOOCV_h2o.xlsx"),
           replace = T)

### understand difference across metrics
annual_summ_combined <- bind_rows(annual_summ, annual_sum_3yr, annual_sum_5yr, annual_sum_10yr)
annual_summ_combined$metric <- factor(annual_summ_combined$metric, 
                                      levels=c("1 year","3 years","5 years","10 years"))
annual_summ_combined$experiment <- factor(annual_summ_combined$experiment, 
                                          levels=c("grid","eco3","eco2","regional"))

saveRDS(annual_summ_combined,
        paste0(result_path,"/model_eval_final/mlp_model_results_LOOCV_h2o.rds"))


#--------------------------------------------------------------------------------
##    Check the tuning results
#--------------------------------------------------------------------------------
tune1 <- readRDS(paste0(result_path,"/model_eval_final/NN/mlp_log_weighted_2lyr_simple_loocv_h2o_tuning.rds")) %>%
  mutate(model = "Neural Net", outcome="log", spec="simple")
tune2 <- readRDS(paste0(result_path,"/model_eval_final/NN/mlp_level_weighted_2lyr_simple_loocv_h2o_tuning.rds")) %>%
  mutate(model = "Neural Net", outcome="level", spec="simple")
tune3 <- readRDS(paste0(result_path,"/model_eval_final/NN/mlp_log_weighted_2lyr_full_loocv_h2o_tuning.rds")) %>%
  mutate(model = "Neural Net", outcome="log", spec="full")
tune4 <- readRDS(paste0(result_path,"/model_eval_final/NN/mlp_level_weighted_2lyr_full_loocv_h2o_tuning.rds")) %>%
  mutate(model = "Neural Net", outcome="level", spec="full")
tune5 <- readRDS(paste0(result_path,"/model_eval_final/NN/mlp_level_weighted_2lyr_simple_loocv_h2o_grid_tuning.rds")) %>%
  mutate(model = "Neural Net", outcome="level", spec="simple")
tune6 <- readRDS(paste0(result_path,"/model_eval_final/NN/mlp_log_weighted_2lyr_simple_loocv_h2o_grid_tuning.rds")) %>%
  mutate(model = "Neural Net", outcome="log", spec="simple")
tune7 <- readRDS(paste0(result_path,"/model_eval_final/NN/mlp_level_weighted_2lyr_full_loocv_h2o_grid_west_tuning.rds")) %>%
  mutate(model = "Neural Net", outcome="level", spec="full")
tune8 <- readRDS(paste0(result_path,"/model_eval_final/NN/mlp_level_weighted_2lyr_full_loocv_h2o_grid_other_tuning.rds")) %>%
  mutate(model = "Neural Net", outcome="level", spec="full")
tune9 <- readRDS(paste0(result_path,"/model_eval_final/NN/mlp_level_weighted_2lyr_full_loocv_h2o_grid_can_tuning.rds")) %>%
  mutate(model = "Neural Net", outcome="level", spec="full")
tune10 <- readRDS(paste0(result_path,"/model_eval_final/NN/mlp_log_weighted_2lyr_full_loocv_h2o_grid_west_other_tuning.rds")) %>%
  mutate(model = "Neural Net", outcome="log", spec="full")
tune11 <- readRDS(paste0(result_path,"/model_eval_final/NN/mlp_log_weighted_2lyr_full_loocv_h2o_grid_can_tuning.rds")) %>%
  mutate(model = "Neural Net", outcome="log", spec="full")

tune <- bind_rows(tune1, tune2, tune3, tune4, tune5, tune6, 
                  tune7, tune8, tune9, tune10, tune11)

tune %>% group_by(experiment, region, outcome, spec)  %>%
  mutate_at(1:4, list(var=var)) %>%
  select(-model, -year) %>%
  filter( hidden1_var==0)

filter(tune, experiment=="grid")$hidden2 %>% unique()
#-----------------------------------------------------
## Check results
#-----------------------------------------------------
ggplot(filter(annual_summ, 
              region=="Mexico", experiment=="eco3",
              spec=="2 layers full", outcome=="level", threshold==0.1)) +
  geom_point(aes(x=year, y=response_positive))



#-----------------------------------------------------
## Compare the training and test set performance
#-----------------------------------------------------
test1 <- readRDS(paste0(result_path,"/model_eval_final/mlp_res_level_weighted_2lyr_simple.rds")) %>%
  mutate(model = "Neural Net", outcome="level", spec="2 layers simple", set="test")
test2 <- readRDS(paste0(result_path,"/model_eval_final/mlp_res_log_weighted_2lyr_simple.rds")) %>%
  mutate(model = "Neural Net", outcome="log", spec="2 layers simple", set="test")
train1 <- readRDS(paste0(result_path,"/model_eval_final/mlp_res_level_weighted_2lyr_simple_train.rds")) %>%
  mutate(model = "Neural Net", outcome="level", spec="2 layers simple", set="train")
train2 <- readRDS(paste0(result_path,"/model_eval_final/mlp_res_log_weighted_2lyr_simple_train.rds")) %>%
  mutate(model = "Neural Net", outcome="log", spec="2 layers simple", set="train")

mlp <- bind_rows(test1, test2, train1, train2)

annual_summ <- mlp  %>%
  mutate(response_positive=response*as.numeric(response >= 0)) %>%
  as.data.frame() %>% 
  group_by(region, region_sub, year, experiment, model, outcome, spec, set) %>%
  summarise(truth=mean(truth, na.rm=T),
            response=mean(response, na.rm=T),
            response_positive=mean(response_positive, na.rm=T)) %>% ungroup() %>%
  group_by(region, year, experiment, model, outcome, spec, set) %>%
  summarise(truth=sum(truth, na.rm=T),
            response=sum(response, na.rm=T),
            response_positive=sum(response_positive, na.rm=T)) %>% ungroup() %>%
  mutate(metric="1 year")

## 10-yr
annual_sum_10yr <- annual_summ %>% group_by(region, experiment, model, outcome, spec, set) %>%
  arrange(year) %>%
  reframe(truth=rollmean(truth, 10, na.rm=T),
          response=rollmean(response, 10, na.rm=T),
          response_positive=rollmean(response_positive, 10, na.rm=T)) %>%
  mutate(metric="10 years")

metric_10yr <- annual_sum_10yr %>% group_by(region, experiment, model, outcome, spec, set) %>%
  mutate(r=cor(response, truth), r_pos=cor(response_positive, truth),
         rmse=sqrt(mean((response - truth)^2)),
         rmse_pos=sqrt(mean((response_positive - truth)^2)),
         truth_mean=mean(truth),
         slope=lm(response_positive~truth)$coefficients[2]) %>% 
  arrange(desc(truth)) %>% slice_head(n = 1) %>%
  mutate(bias = (response_positive - truth)/truth) %>%
  mutate(metric="10 years") %>% ungroup()


ggplot(filter(annual_sum_10yr, set=="test", region=="Northeastern US"), 
       aes(x=truth/1e9, y=response_positive/1e9)) +
  geom_point(size=4, alpha=0.6) +
  geom_abline(,0,1,linetype="dashed") +
  theme_classic() + theme(text = element_text(size=16)) + 
  labs(y="Predictions (million tons DM)", x="Observations (million tons DM)") +
  facet_grid(outcome~experiment) 
ggsave("westUS_NN_models.png", width=11, height=6.5)

ggplot(filter(annual_summ_combined, experiment == "regional", region == rrr),
       aes(x=truth/1e9, y=response_positive/1e9)) +
  geom_point(size=4, alpha=0.6) +
  geom_abline(,0,1,linetype="dashed") +
  theme_classic() + theme(text = element_text(size=16)) +
  labs(y="Predictions (million tons DM)", x="Observations (million tons DM)") +
  facet_wrap(~metric, scales = "free") 
ggsave("westUS_singleNN_metrics_regional.png", width=8, height=5)#+

ggplot(filter(annual_summ_combined, metric =="10 years",region == rrr),
       aes(x=truth/1e9, y=response_positive/1e9)) +
  geom_point(size=4, alpha=0.6) +
  geom_abline(,0,1,linetype="dashed") +
  theme_classic() + theme(text = element_text(size=16)) +
  labs(y="Predictions (million tons DM)", x="Observations (million tons DM)") +
  facet_wrap(~experiment) 
ggsave("westUS_NN_10yr_spatial_exp.png", width=8, height=5)#+

#### Plot the model performances by regions
ggplot(filter(metric_combined), 
       aes(x = experiment, y = r_pos, fill = metric)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_classic() + theme(text = element_text(size=16))  +
  labs(x="", y="Correlation coef (R)") +
  facet_grid(~region)
ggsave("singeNN_models_performance.png", width=12, height=5)

ggplot(filter(metric_combined, metric == "10 years"), 
       aes(x = experiment, y = bias)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_classic() + theme(text = element_text(size=16))  +
  labs(x="", y="Prediction bias") +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(~region)
ggsave("singleNN_models_10yrs_bias.png", width=5, height=3.5)

#--------------------------------------------------------------------------------
##    Analyze model performance by subregions
#--------------------------------------------------------------------------------
west_eco3 <- readRDS(paste0(result_path,"/model_eval_final/westernUS/west_nestedCV_annual_spatialres_test_NN.rds")) %>%
  filter(experiment == "eco3")

####### combine all regions
annual_summ <- west_eco3 %>% #
  filter(upsample_weight ==0) %>%
  mutate(response_positive=response*as.numeric(response >= 0)) %>%
  as.data.frame() %>% 
  group_by(region, region_sub, year, experiment) %>%
  summarise(truth=sum(truth, na.rm=T),
            response=sum(response, na.rm=T),
            response_positive=sum(response_positive, na.rm=T)) %>% 
  ungroup() %>%
  mutate(metric="1 year")

annual_summ %>% mutate(ratio=response_positive/truth) %>%
  arrange(desc(ratio))

regionsub_summ <- annual_summ %>% group_by(region_sub) %>% 
  summarise(truth=sum(truth), n=length(year)) %>% arrange(desc(truth)) %>% as.data.frame()

### only evaluate sub regions with non-zero emissions every year
annual_summ <- filter(annual_summ, region_sub %in% filter(regionsub_summ, n==21)$region_sub)

metric_1yr <- annual_summ %>% group_by(region, region_sub, experiment) %>%
  mutate(r=cor(response, truth), r_pos=cor(response_positive, truth),
         rmse=sqrt(mean((response - truth)^2)), rmse_pos=sqrt(mean((response_positive - truth)^2))) %>% 
  arrange(desc(truth)) %>% slice_head(n = 1) %>%
  mutate(bias = (response_positive - truth)/truth) %>%
  mutate(metric="1 year")

## 3-yr
annual_sum_3yr <- annual_summ %>% group_by(region, region_sub, experiment) %>%
  arrange(year) %>%
  reframe(truth=rollmean(truth, 3, na.rm=T),
          response=rollmean(response, 3, na.rm=T),
          response_positive=rollmean(response_positive, 3, na.rm=T)) %>% ungroup()%>%
  mutate(metric="3 years")

metric_3yr <- annual_sum_3yr %>% group_by(region, region_sub, experiment) %>%
  mutate(r=cor(response, truth), r_pos=cor(response_positive, truth),
         rmse=sqrt(mean((response - truth)^2)), rmse_pos=sqrt(mean((response_positive - truth)^2))) %>% 
  arrange(desc(truth)) %>% slice_head(n = 1) %>%
  mutate(bias = (response_positive - truth)/truth) %>%
  mutate(metric="3 years")

## 5-yr
annual_sum_5yr <- annual_summ %>% group_by(region, region_sub, experiment) %>%
  arrange(year) %>%
  reframe(truth=rollmean(truth, 5, na.rm=T),
          response=rollmean(response, 5, na.rm=T),
          response_positive=rollmean(response_positive, 5, na.rm=T)) %>%
  mutate(metric="5 years")

metric_5yr <- annual_sum_5yr %>% group_by(region, region_sub, experiment) %>%
  mutate(r=cor(response, truth), r_pos=cor(response_positive, truth),
         rmse=sqrt(mean((response - truth)^2)), rmse_pos=sqrt(mean((response_positive - truth)^2))) %>% 
  arrange(desc(truth)) %>% slice_head(n = 1) %>%
  mutate(bias = (response_positive - truth)/truth) %>%
  mutate(metric="5 years")

## 10-yr
annual_sum_10yr <- annual_summ %>% group_by(region, region_sub, experiment) %>%
  arrange(year) %>%
  reframe(truth=rollmean(truth, 10, na.rm=T),
          response=rollmean(response, 10, na.rm=T),
          response_positive=rollmean(response_positive, 10, na.rm=T)) %>%
  mutate(metric="10 years")

metric_10yr <- annual_sum_10yr %>% group_by(region, region_sub, experiment) %>%
  mutate(r=cor(response, truth), r_pos=cor(response_positive, truth),
         rmse=sqrt(mean((response - truth)^2)), rmse_pos=sqrt(mean((response_positive - truth)^2))) %>% 
  arrange(desc(truth)) %>% slice_head(n = 1) %>%
  mutate(bias = (response_positive - truth)/truth) %>%
  mutate(metric="10 years")

metric_combined <- bind_rows(metric_1yr, metric_3yr, metric_5yr, metric_10yr)
metric_combined$metric <- factor(metric_combined$metric, 
                                 levels=c("1 year","3 years","5 years","10 years")) 
write.xlsx(metric_combined, 
           paste0(result_path,"/model_eval_final/singleNN_model_metrics_combined_subregions.xlsx"),
           replace = T)

##
hist(metric_10yr$bias, breaks=50)

### plot the performance for the top regions
top_12regions <- filter(annual_sum_10yr, region_sub %in% regionsub_summ[1:12,"region_sub"])

hist(metric_10yr$bias, breaks=100)
filter(metric_10yr, bias>1)

top_12regions$region_sub <- factor(top_12regions$region_sub,
                                   levels=regionsub_summ[1:12,"region_sub"])


ggplot(top_12regions,
       aes(x=truth/1e9, y=response_positive/1e9)) +
  geom_point(size=4, alpha=0.6) +
  geom_abline(,0,1,linetype="dashed") +
  theme_classic() + theme(text = element_text(size=13)) +
  labs(y="Predictions (million tons DM)", x="Observations (million tons DM)") +
  facet_wrap(~region_sub, scales = "free") 
ggsave("westUS_singleNN_subregions_10yr.png", width=10, height=5)


#### Plot the model performances by regions
ggplot(filter(metric_10yr, region_sub %in% regionsub_summ[1:10,"region_sub"]), 
       aes(x = experiment, y = r_pos, fill = metric)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_classic() + theme(text = element_text(size=16))  +
  labs(x="", y="Correlation coef (R)") +
  facet_grid(~region)
ggsave("singeNN_models_performance.png", width=12, height=4)


# #------------------
# #    Western US
# #------------------
# west_res <- readRDS(paste0(result_path,"/model_eval_final/westernUS/west_mlp_res_log_weighted.rds")) %>%
#   mutate(model = "log, FE weighted")
# west_res1 <- readRDS(paste0(result_path,"/model_eval_final/westernUS/west_mlp_res_level_weighted.rds")) %>%
#   mutate(model = "level, FE weighted")
# 
# tune1 <- readRDS(paste0(result_path,"/model_eval_final/westernUS/west_mlp_res_log_weighted_tuning.rds"))
# tune2 <- readRDS(paste0(result_path,"/model_eval_final/westernUS/west_mlp_res_level_weighted_tuning.rds"))# #tune2 <- readRDS(paste0(result_path,"/model_eval_final/westernUS/west_nestedCV_annual_grid_tuning_NN.rds"))
# bind_rows(tune1, tune2)  %>% group_by(experiment, region) %>%
#   arrange(RMSE) %>% slice_head(n=3) %>% as.data.frame()
# 
# # west1 <- readRDS(paste0(result_path,"/model_eval_final/westernUS/west_NN_full_log_FE.rds")) %>%
# #   mutate(model = "log, FE")
# # west2 <- readRDS(paste0(result_path,"/model_eval_final/westernUS/west_NN_full_log_FE_weighted.rds")) %>%
# #   mutate(model = "log, FE weighted")
# # west3 <- readRDS(paste0(result_path,"/model_eval_final/westernUS/west_NN_full_level_FE.rds")) %>%
# #   mutate(model = "level, FE")
# # west4 <- readRDS(paste0(result_path,"/model_eval_final/westernUS/west_NN_full_level_FE_weighted.rds")) %>%
# #   mutate(model = "level, FE weighted")
# # 
# # tune1 <- readRDS(paste0(result_path,"/model_eval_final/westernUS/west_NN_full_log_FE_tuning.rds")) %>%
# #   mutate(model = "log, FE")
# # tune2 <- readRDS(paste0(result_path,"/model_eval_final/westernUS/west_NN_full_log_FE_weighted_tuning.rds")) %>%
# #   mutate(model = "log, FE weighted")
# # tune3 <- readRDS(paste0(result_path,"/model_eval_final/westernUS/west_NN_full_level_FE_tuning.rds")) %>%
# #   mutate(model = "level, FE")
# # tune4 <- readRDS(paste0(result_path,"/model_eval_final/westernUS/west_NN_full_level_FE_weighted_tuning.rds")) %>%
# #   mutate(model = "level, FE weighted")
# # 
# # bind_rows(tune1, tune2, tune3, tune4)  %>% group_by(experiment, region, model, iteration) %>%
# #   arrange(regr.mse) %>% slice_head(n=3) %>% as.data.frame() %>%
# #   filter(experiment=="eco3")
# 
# ##### compare across NNs of different specifications
# metric1 <- read.xlsx(paste0(result_path,"/model_eval_final/west_NN_model_metrics.xlsx")) %>%
#   filter(model == "log, FE weighted") %>% mutate(spec="1 layer, full features")
# # metric2 <- read.xlsx(paste0(result_path,"/model_eval_final/singleNN_model_metrics.xlsx")) %>%
# #   filter(model == "log, FE weighted") %>% mutate(spec="1 layer, simple features")
# metric3 <- read.xlsx(paste0(result_path,"/model_eval_final/west_mlp_model_metrics.xlsx")) %>%
#   filter(model == "log, FE weighted") %>% mutate(spec="3 layer, full features")
# metric4 <- read.xlsx(paste0(result_path,"/model_eval_final/linear_model_metrics_combined.xlsx")) %>%
#   filter(model == "log, weighted", region=="Western US") %>% mutate(spec="Linear")
# metric5 <- read.xlsx(paste0(result_path,"/model_eval_final/lasso_model_metrics_combined.xlsx")) %>%
#   filter(model == "log, FE weighted", region=="Western US") %>% mutate(spec="LASSO")
# 
# metric_tmp <- filter(bind_rows(metric1, metric3, metric4, metric5), metric == "10 years")
# metric_tmp$experiment <- factor(metric_tmp$experiment, 
#                                 levels=c("grid","eco3","eco2","regional")) 
# 
# ggplot(metric_tmp, 
#        aes(x = experiment, y = bias, fill=spec)) +
#   geom_bar(stat = "identity", position = position_dodge()) +
#   theme_classic() + theme(text = element_text(size=16))  +
#   labs(x="", y="Prediction bias") +
#   scale_y_continuous(labels = scales::percent) 
# ggsave("west_across_models_10yrs_bias.png", width=8, height=4.5)

