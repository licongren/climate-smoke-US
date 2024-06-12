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
# Analyze the performance of three types of models and pick the best-perfromed ones
# Written by Minghao
# Last edited July 2022
#-------------------------------------------------------------------------------
linear_metric <- read.xlsx(paste0(result_path,"/model_eval_final/linear_model_metrics_combined_LOOCV.xlsx"))
lasso_metric <- read.xlsx(paste0(result_path,"/model_eval_final/lasso_model_metrics_combined_LOOCV.xlsx"))
mlp_metric <- read.xlsx(paste0(result_path,"/model_eval_final/mlp_model_metrics_combined_LOOCV_h2o.xlsx")) 
metrics <- bind_rows(linear_metric, lasso_metric, mlp_metric) %>%
  mutate(bias_abs=abs(bias),
         rmse_bias= bias_abs + rmse_pos_mean)

linear_results <- readRDS(paste0(result_path,"/model_eval_final/linear_model_results_LOOCV.rds"))
lasso_results <- readRDS(paste0(result_path,"/model_eval_final/lasso_model_results_LOOCV.rds"))
mlp_results <- readRDS(paste0(result_path,"/model_eval_final/mlp_model_results_LOOCV_h2o.rds")) 
results <- bind_rows(linear_results, lasso_results, mlp_results)

#-------------------------------------------------------------------
# 1. Positive R, RMSE + bias , two model per each family (6 models in total)
# consider both outcomes: level and log
#-------------------------------------------------------------------
best_metric <- metrics %>% 
  filter(r_pos > 0) %>%
  group_by(region, model, metric, outcome) %>%
  arrange(rmse_bias) %>% slice_head(n = 1) %>%
  select(region, model, experiment, outcome, spec, r_pos, slope, bias, rmse_pos_mean, rmse_bias, metric) %>%
  group_by(region, metric) %>%
  mutate(diff_min=rmse_bias - min(rmse_bias)) %>%
  ungroup() %>%
  mutate(selected = as.numeric(diff_min < 0.1))

best_metric$region <-  factor(best_metric$region,
                                   levels=c("Western US", "Canada-Alaska","Mexico", "Southeastern US", "Northeastern US"))
best_metric$model <-  factor(best_metric$model,
                                  levels=c("Linear", "LASSO","Neural Net"))
best_metric$model_outcome <- paste(best_metric$model, best_metric$outcome, sep=", ")
best_metric$model_outcome <-  factor(best_metric$model_outcome,
                                  levels=c("Linear, level", "Linear, log", 
                                           "LASSO, level", "LASSO, log",
                                           "Neural Net, level", "Neural Net, log"))
write.xlsx(best_metric, paste0(result_path,"/model_eval_final/best_metrics_rmse_bias_6models_LOOCV.xlsx"))

best_metric <- read.xlsx(paste0(result_path,"/model_eval_final/best_metrics_rmse_bias_6models_LOOCV.xlsx"))
color_pal <- met.brewer(name="Monet", n=6, type="discrete")[c(1,2,4,3,5,6)]
threshold <- 0.05

p1 <- ggplot(filter(best_metric, metric=="10 years"), 
             aes(x = region, y = rmse_pos_mean, fill = model_outcome, colour = model_outcome, alpha=as.numeric(diff_min < threshold))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_classic() + theme(text = element_text(size=16), axis.text = element_text(size=16))  +
  labs(x="", y="RMSE / Mean") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values=color_pal) + 
  scale_colour_manual(values=color_pal) + 
  scale_alpha(range=c(0.1,1)) + guides(alpha="none", colour="none")

p2 <- ggplot(filter(best_metric, metric=="10 years"), 
             aes(x = region, y = bias, fill = model_outcome, colour = model_outcome, alpha=as.numeric(diff_min < threshold))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_classic() + theme(text = element_text(size=16), axis.text = element_text(size=16))  +
  labs(x="", y="bias for the top 10-years") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values=color_pal) + 
  scale_colour_manual(values=color_pal) + 
  scale_alpha(range=c(0.1,1)) + guides(alpha="none", colour="none")
plot_grid(p1, p2, ncol = 1)
ggsave("model_performance_rmsebias_10yr_LOOCV_threshold05.png", width=12, height=7)

results_select <- left_join(filter(results, metric == "10 years"), 
                            filter(best_metric, metric == "10 years"),
                            by=c("region", "model", "outcome", "experiment", "spec", "metric"))  %>%
  filter(!is.na(model_outcome))
  

for (rrr in unique(best_metric$region)){
ggplot(filter(results_select, region==rrr),
       aes(x=truth/1e9, y=response_positive/1e9, colour=factor(diff_min<threshold))) +
  geom_point(size=4, alpha=0.6) +
  geom_abline(,0,1,linetype="dashed") +
  theme_classic() + theme(text = element_text(size=16), axis.text = element_text(size=16)) + 
  labs(y="Predictions (million tons DM)", x="Observations (million tons DM)") +
  facet_wrap(~model_outcome, nrow=1) +
  scale_colour_manual(values=c("grey30", "green4")) +
  guides(colour= "none")
ggsave(paste0(rrr,"_best_models_10yr_LOOCV_threshold05.png"), width=11, height=3)
}

##### plot the performance after aggregating over selected models
results_select_mean <- results_select %>%
  filter(diff_min<threshold) %>%
  group_by(region, experiment, model, outcome, spec) %>%
  mutate(year = 1:n()) %>%
  group_by(region, year) %>%
  summarise(truth = mean(truth, na.rm=T),
            response_positive = mean(response_positive, na.rm=T)) %>%
  ungroup()

for (rrr in unique(best_metric$region)){
  axis_range <- range(c(filter(results_select_mean, region==rrr)$response_positive,
                        filter(results_select_mean, region==rrr)$truth))/1e9
  ggplot(filter(results_select_mean, region==rrr),
         aes(x=truth/1e9, y=response_positive/1e9)) +
    geom_point(size=4, alpha=0.7) +
    geom_abline(,0,1,linetype="dashed") +
    theme_classic() + theme(text = element_text(size=14), axis.text = element_text(size=14)) + 
    scale_x_continuous(limits= axis_range) +
    scale_y_continuous(limits= axis_range) +
    labs(y="Predictions (Mton DM)", x="Observations (Mton DM)") +
    #scale_colour_manual(values=c("grey30", "green4")) +
    guides(colour= "none")
  ggsave(paste0(rrr,"_best_model_threshold05.pdf"), width=3, height=3)
}

#-------------------------------------------------------------------
# 2. Positive R, RMSE + bias , two models per each family (6 models in total), 
#. limiting to non-regional
#-------------------------------------------------------------------
best_metric <- metrics %>% 
  filter(r_pos > 0, experiment!="regional") %>%
  group_by(region, model, metric, outcome) %>%
  arrange(rmse_bias) %>% slice_head(n = 1) %>%
  select(region, model, experiment, outcome, spec, r_pos, slope, bias, rmse_pos_mean, rmse_bias, metric) %>%
  group_by(region, metric) %>%
  mutate(diff_min=rmse_bias - min(rmse_bias)) %>%
  ungroup() %>%
  mutate(selected = as.numeric(diff_min < 0.1))

best_metric$region <-  factor(best_metric$region,
                              levels=c("Western US", "Canada-Alaska","Mexico", "Southeastern US", "Northeastern US"))
best_metric$model <-  factor(best_metric$model,
                             levels=c("Linear", "LASSO","Neural Net"))
best_metric$model_outcome <- paste(best_metric$model, best_metric$outcome, sep=", ")
best_metric$model_outcome <-  factor(best_metric$model_outcome,
                                     levels=c("Linear, level", "Linear, log", 
                                              "LASSO, level", "LASSO, log",
                                              "Neural Net, level", "Neural Net, log"))
write.xlsx(best_metric, paste0(result_path,"/model_eval_final/best_metrics_rmse_bias_6models_subregion.xlsx"))

color_pal <- met.brewer(name="Monet", n=6, type="discrete")[c(1,2,4,3,5,6)]

p1 <- ggplot(filter(best_metric, metric=="5 years"), 
             aes(x = region, y = rmse_pos_mean, fill = model_outcome, colour = model_outcome, alpha=selected)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_classic() + theme(text = element_text(size=16), axis.text = element_text(size=16))  +
  labs(x="", y="RMSE / Mean") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values=color_pal) + 
  scale_colour_manual(values=color_pal) + 
  scale_alpha(range=c(0.1,1)) + guides(alpha="none", colour="none")

p2 <- ggplot(filter(best_metric, metric=="5 years"), 
             aes(x = region, y = bias, fill = model_outcome, colour = model_outcome, alpha=selected)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_classic() + theme(text = element_text(size=16), axis.text = element_text(size=16))  +
  labs(x="", y="bias for the top 5-years") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values=color_pal) + 
  scale_colour_manual(values=color_pal) + 
  scale_alpha(range=c(0.1,1)) + guides(alpha="none", colour="none")
plot_grid(p1, p2, ncol = 1)
ggsave("model_performance_rmsebias_5yr_subregion.png", width=12, height=7)

results_select <- left_join(filter(results, metric == "5 years"), 
                            filter(best_metric, metric == "5 years")) 

for (rrr in unique(best_metric$region)){
  ggplot(filter(results_select, region==rrr),
         aes(x=truth/1e9, y=response_positive/1e9, colour=selected)) +
    geom_point(size=4, alpha=0.6) +
    geom_abline(,0,1,linetype="dashed") +
    theme_classic() + theme(text = element_text(size=16), axis.text = element_text(size=16)) + 
    labs(y="Predictions (million tons DM)", x="Observations (million tons DM)") +
    facet_grid(region~model_outcome, scales = "free") 
  ggsave(paste0(rrr,"_best_models_rmsebias_5yr_subregion.png"), width=8, height=6.5)
}

#-------------------------------------------------------------------
# 3. Positive R, RMSE + bias , one model per each family (3 models in total)
#-------------------------------------------------------------------
best_metric <- metrics %>% 
  filter(r_pos > 0) %>%
  group_by(region, model, metric) %>%
  arrange(rmse_bias) %>% slice_head(n = 1) %>%
  select(region, model, experiment, outcome, spec, r_pos, slope, bias, rmse_pos_mean, rmse_bias, threshold, metric) %>%
  group_by(region, metric) %>%
  mutate(diff_min=rmse_bias - min(rmse_bias)) %>%
  ungroup() %>%
  mutate(selected = as.numeric(diff_min < 0.1))

best_metric$region <-  factor(best_metric$region,
                              levels=c("Western US", "Canada-Alaska","Mexico", "Southeastern US", "Northeastern US"))
best_metric$model <-  factor(best_metric$model,
                             levels=c("Linear", "LASSO","Neural Net"))
best_metric$model_outcome <- best_metric$model
# best_metric$model_outcome <-  factor(best_metric$model_outcome,
#                                      levels=c("Linear, level", "Linear, log", 
#                                               "LASSO, level", "LASSO, log",
#                                               "Neural Net, level", "Neural Net, log"))
write.xlsx(best_metric, paste0(result_path,"/model_eval_final/best_metrics_rmse_bias_3models.xlsx"))

color_pal <- met.brewer(name="Monet", n=6, type="discrete")[c(1,4,5)]

p1 <- ggplot(filter(best_metric, metric=="10 years"), 
             aes(x = region, y = rmse_pos_mean, fill = model_outcome, colour = model_outcome, alpha=selected)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_classic() + theme(text = element_text(size=16), axis.text = element_text(size=16))  +
  labs(x="", y="RMSE / Mean") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values=color_pal) + 
  scale_colour_manual(values=color_pal) + 
  scale_alpha(range=c(0.1,1)) + guides(alpha="none", colour="none")

p2 <- ggplot(filter(best_metric, metric=="10 years"), 
             aes(x = region, y = bias, fill = model_outcome, colour = model_outcome, alpha=selected)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_classic() + theme(text = element_text(size=16), axis.text = element_text(size=16))  +
  labs(x="", y="bias for the top 10-years") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values=color_pal) + 
  scale_colour_manual(values=color_pal) + 
  scale_alpha(range=c(0.1,1)) + guides(alpha="none", colour="none")
plot_grid(p1, p2, ncol = 1)
ggsave("model_performance_rmsebias_10yr_3models.png", width=12, height=7)

results_select <- left_join(filter(results, metric == "10 years"), filter(best_metric, metric == "10 years")) 

for (rrr in unique(best_metric$region)){
  ggplot(filter(results_select, selected==1, region==rrr),
         aes(x=truth/1e9, y=response_positive/1e9)) +
    geom_point(size=4, alpha=0.6) +
    geom_abline(,0,1,linetype="dashed") +
    theme_classic() + theme(text = element_text(size=16), axis.text = element_text(size=16)) + 
    labs(y="Predictions (million tons DM)", x="Observations (million tons DM)") +
    facet_wrap(~model_outcome, ncol=2) 
  ggsave(paste0(rrr,"_best_models_rmsebias_10yr.png"), width=8, height=6.5)
}

# #-------------------------------------------------------------------
# # 3. Positive R, only RMSE, non-regional, one model per each family (6 models in total)
# #-------------------------------------------------------------------
# best_metric <- metrics %>% 
#   filter(r_pos > 0) %>%
#   filter(experiment )
#   group_by(region, model, metric) %>%
#   arrange(rmse_pos_mean) %>% slice_head(n = 1) %>%
#   select(region, experiment, model, outcome, r_pos, slope, bias, rmse_pos_mean)
# 
# #filter(best_metric, metric=="5 years")
# best_metric_10yr <- filter(best_metric, metric=="10 years") %>%
#   group_by(region) %>%
#   mutate(diff_min_rmse=rmse_pos_mean - min(rmse_pos_mean)) %>%
#   ungroup() %>%
#   mutate(selected = as.numeric(diff_min_rmse < 0.1)) 
# 
# best_metric_10yr$region <-  factor(best_metric_10yr$region,
#                                    levels=c("Western US", "Canada-Alaska","Mexico", "Southeastern US", "Northeastern US"))
# best_metric_10yr$model <-  factor(best_metric_10yr$model,
#                                   levels=c("Linear", "LASSO","MLP"),
#                                   labels=c("Linear", "LASSO","Neural Net"))
# 
# color_pal <- met.brewer(name="Juarez", n=3, type="discrete")
# p1 <- ggplot(best_metric_10yr, 
#              aes(x = region, y = rmse_pos_mean, fill = model, alpha=selected)) +
#   geom_bar(stat = "identity", position = position_dodge()) +
#   theme_classic() + theme(text = element_text(size=16))  +
#   labs(x="", y="RMSE / Mean") +
#   scale_y_continuous(labels = scales::percent) +
#   scale_fill_manual(values=color_pal) + 
#   scale_alpha(range=c(0.2,1)) + guides(alpha="none")
# 
# p2 <- ggplot(best_metric_10yr, 
#              aes(x = region, y = bias, fill = model, alpha=selected)) +
#   geom_bar(stat = "identity", position = position_dodge()) +
#   theme_classic() + theme(text = element_text(size=16))  +
#   labs(x="", y="bias for the top 10-years") +
#   scale_y_continuous(labels = scales::percent) +
#   scale_fill_manual(values=color_pal) + 
#   scale_alpha(range=c(0.2,1)) + guides(alpha="none")
# plot_grid(p1, p2, ncol = 1)
# ggsave("model_performance_rmse_only.png", width=12, height=8)
# 


