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
library(glmnet)
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
# Analyze the evaluation performance of LASSO models
# Written by Minghao
# Last edited Sep 2022
#-------------------------------------------------------------------------------
test1 <- readRDS(paste0(result_path,"/model_eval_final/LASSO/LASSO_full_level_loocv.rds")) %>%
  mutate(model = "LASSO", outcome="level", spec= "weighted FE", grid_1deg_id = as.numeric(grid_1deg_id))
test2 <- readRDS(paste0(result_path,"/model_eval_final/LASSO/LASSO_full_log_loocv.rds")) %>%
  mutate(model = "LASSO", outcome="log", spec= "weighted FE", grid_1deg_id = as.numeric(grid_1deg_id))
test3 <- readRDS(paste0(result_path,"/model_eval_final/LASSO/LASSO_full_level_nonweighted_loocv.rds")) %>%
  mutate(model = "LASSO", outcome="level", spec= "non-weighted FE", grid_1deg_id = as.numeric(grid_1deg_id))
test4 <- readRDS(paste0(result_path,"/model_eval_final/LASSO/LASSO_full_log_nonweighted_loocv.rds")) %>%
  mutate(model = "LASSO", outcome="log", spec= "non-weighted FE", grid_1deg_id = as.numeric(grid_1deg_id))

test <- bind_rows(test1, test2, test3, test4)

#------------ Check if the hyperparams are appropriately tuned
tune1 <- readRDS(paste0(result_path,"/model_eval_final/LASSO/LASSO_full_level_loocv_tuning.rds")) %>%
  mutate(model = "LASSO", outcome="level", spec= "weighted FE")
tune2 <- readRDS(paste0(result_path,"/model_eval_final/LASSO/LASSO_full_log_loocv_tuning.rds")) %>%
  mutate(model = "LASSO", outcome="log", spec= "weighted FE")
tune3 <- readRDS(paste0(result_path,"/model_eval_final/LASSO/LASSO_full_level_nonweighted_loocv_tuning.rds")) %>%
  mutate(model = "LASSO", outcome="level", spec= "non-weighted FE")
tune4 <- readRDS(paste0(result_path,"/model_eval_final/LASSO/LASSO_full_log_nonweighted_loocv_tuning.rds")) %>%
  mutate(model = "LASSO", outcome="log", spec= "non-weighted FE")

tune <- bind_rows(tune1, tune2, tune3, tune4) %>%
  mutate(lambda=as.numeric(lambda))

tune %>% group_by(region, iteration, experiment, model, outcome, spec) %>%
  mutate(lambda_min=min(lambda), lambda_max=max(lambda),
         flag = (lambda < lambda_max & lambda > lambda_min)) %>%
  arrange(regr.mse, lambda) %>% slice_head(n=3) %>% as.data.frame() %>%
  filter(!flag)

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
metric_combined$experiment <- factor(metric_combined$experiment, 
                                     levels=c("grid","eco3","eco2","regional")) 
metric_combined$metric <- factor(metric_combined$metric, 
                                 levels=c("1 year","3 years","5 years","10 years")) 
metric_combined <- metric_combined  %>% arrange(metric, experiment)
write.xlsx(metric_combined, 
           paste0(result_path,"/model_eval_final/lasso_model_metrics_combined_loocv.xlsx"),
           replace = T)

annual_summ_combined <- bind_rows(annual_summ, annual_sum_3yr, annual_sum_5yr, annual_sum_10yr)
annual_summ_combined$metric <- factor(annual_summ_combined$metric, 
                                      levels=c("1 year","3 years","5 years","10 years"))
annual_summ_combined$experiment <- factor(annual_summ_combined$experiment, 
                                          levels=c("grid","eco3","eco2","regional")) 
saveRDS(annual_summ_combined,
        paste0(result_path,"/model_eval_final/lasso_model_results_loocv.rds"))




### plot the model performances
region_list <-  c("Western US", "Canada-Alaska","Mexico", "Southeastern US", "Northeastern US")
rrr <- "Canada-Alaska"

ggplot(filter(annual_summ_combined, metric == "10 years", region == rrr,
              experiment=="regional"),
       aes(x=truth/1e9, y=response_positive/1e9)) +
  geom_point(size=4, alpha=0.6) +
  geom_abline(,0,1,linetype="dashed") +
  theme_classic() + theme(text = element_text(size=16)) + 
  labs(y="Predictions (million tons DM)", x="Observations (million tons DM)") +
  facet_grid(experiment~outcome) 
ggsave("Canada_lasso_models_10yr.png", width=9, height=7.5)#+


ggplot(filter(annual_summ_combined, model=="log", experiment == "regional", region == rrr),
       aes(x=truth/1e9, y=response_positive/1e9)) +
  geom_point(size=4, alpha=0.6) +
  geom_abline(,0,1,linetype="dashed") +
  theme_classic() + theme(text = element_text(size=16)) +
  labs(y="Predictions (million tons DM)", x="Observations (million tons DM)") +
  facet_wrap(~metric, scales = "free") 
ggsave("westUS_lasso_region_log_metrics.png", width=8, height=5)#+

ggplot(filter(metric_combined), 
       aes(x = experiment, y = r_pos, fill = metric)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_classic() + theme(text = element_text(size=16))  +
  labs(x="", y="Correlation coef (R)") +
  facet_grid(~model)
ggsave("lasso_models_performance.png", width=12, height=6)

ggplot(filter(metric_combined), 
       aes(x = experiment, y = bias, fill = metric)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_classic() + theme(text = element_text(size=16))  +
  labs(x="", y="Correlation coef (R)") +
  facet_grid(~model)
ggsave("lasso_models_performance.png", width=12, height=6)

