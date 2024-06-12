library(terra);library(ncdf4)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(plotly)
library(sp);library(sf)
library(dplyr)   
library(fixest)
library(tidyverse)
library(MetBrewer)
library(stringr)
library(splines2)
library(exactextractr)
library(zoo)

rm(list=ls())
gc()

code_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/rcode"
data_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data"
figure_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/figures"
result_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/results/fire_model"

db_path <- "/Users/mhqiu/BurkeLab Dropbox/Projects/smoke-climate"

setwd(figure_path)
states_map <- map_data("state")

color_map <- c("grey65",
               rev(met.brewer(name="Homer1", n=8, type="discrete"))[c(2,6,8)])

#-------------------------------------------------------------------------------
# Analyze the projections of fire DM from CMIP6
# Written by Minghao
# Last edited July 2022
#-------------------------------------------------------------------------------
## load gfed annual emissions

#-----------------------------------------------
#### Load the fire DM predictions
#-----------------------------------------------
cmip_pred1 <- readRDS(paste0(result_path,"/CMIP6_fire_proj/linear_final_pred_2099.rds"))
cmip_pred2 <- readRDS(paste0(result_path,"/CMIP6_fire_proj/LASSO_final_pred.rds"))
cmip_pred3 <- readRDS(paste0(result_path,"/CMIP6_fire_proj/nnet_final_pred.rds")) %>%
  mutate(year=as.numeric(as.character(year)))

cmip_pred <- bind_rows(cmip_pred1, cmip_pred2, cmip_pred3) %>%
  rename(gcm = model) 
cmip_pred$algorithm_outcome <- paste(cmip_pred$algorithm, cmip_pred$outcome, sep=", ")

cmip_pred %>% group_by(algorithm_outcome, region) %>%
  summarise(max_year=max(year)) %>% ungroup()

### find the gcms that can predict for the 5808 grid cell (a proxy for it can proj for all the selected subregions)
gcm_select <- filter(cmip_pred, year==2001) %>%
  group_by(gcm) %>%
  summarise(n=length(gcm)) %>% ungroup() %>%
  filter(n > 1500) %>% ##this selects 28 gcms with near-complete coverage (no eco-2/eco-3 is missing)
  select(gcm) %>% unlist() %>% as.vector()

#cmip_pred <-  filter(cmip_pred, gcm %in% gcm_select, year<2071, scenario!="ssp585")
cmip_pred <-  filter(cmip_pred, gcm %in% gcm_select)

##### select a common scope for projection (i.e. the same regions_subs)
region_sub_count <- cmip_pred %>% 
  filter(year==2001) %>%
  group_by(region, region_sub, spatial_res, algorithm_outcome, spec) %>%
  summarise(n=length(region)) %>% ungroup() %>%
  filter(n==length(gcm_select)) %>% ### choose the common scopes that are available from at least 28 gcms
  distinct(region, region_sub, n) 

cmip_pred <- cmip_pred %>%
  left_join(region_sub_count, by=c("region", "region_sub")) %>%
  filter(!is.na(n)) %>% select(-n)
# model <- readRDS(paste0(result_path,"/CMIP6_fire_proj/NN_h2o_final_models.rds"))
# data_norm <- readRDS(paste0(result_path,"/CMIP6_fire_proj/NN_norm_data.rds"))

#### Load the best models
best_models <- read.xlsx(paste0(result_path,"/model_eval_final/best_metrics_rmse_bias_6models_LOOCV.xlsx"))
#best_models <- read.xlsx(paste0(result_path,"/model_eval_final/best_metrics_rmse_bias_6models_subregion_LOOCV.xlsx"))
best_models <- best_models %>% 
  filter(selected==1, metric %in% c("10 years"),
         diff_min < 0.05) %>%
  distinct(region, model, experiment, outcome, model_outcome, spec, selected, diff_min) %>%
  rename(algorithm = model, spatial_res=experiment) 

### combine with the model performance to select the best-performed models
cmip_pred_best <- cmip_pred %>%
  mutate(spatial_res=replace(spatial_res, spatial_res=="region", "regional")) %>%
  left_join(best_models, by=c("region", "algorithm", "outcome", "spatial_res", "spec")) %>%
  mutate(pred_DM=pred_DM*as.numeric(pred_DM >= 0)) %>%
  filter(!is.na(selected))

### use 2001-2021 as the normalization baselines
cmip_2001_2021 <- filter(cmip_pred_best, scenario%in%c("historical", "ssp370"), year %in% seq(2001,2021)) 
cmip_2001_2021_mean <- cmip_2001_2021 %>%
  group_by(gcm, region, region_sub, algorithm_outcome) %>%
  summarise(DM_2001_2021 = mean(pred_DM, na.rm=T)) %>% ungroup()

#### calculate the delta differences within gcm projections and then add onto real obs
cmip_pred_best <-  bind_rows(filter(cmip_pred_best, year>2021),
                          cmip_2001_2021 %>% mutate(scenario = "historical"),
                          cmip_2001_2021 %>% mutate(scenario = "ssp126"),
                          cmip_2001_2021 %>% mutate(scenario = "ssp245"),
                          cmip_2001_2021 %>% mutate(scenario = "ssp370"),
                          cmip_2001_2021 %>% mutate(scenario = "ssp585")) %>%
  left_join(cmip_2001_2021_mean) %>%
  mutate(DM_anomaly = pred_DM - DM_2001_2021) %>% select(-DM_2001_2021)

### add the deltas to the mean of observed emissions between 2001 to 2021
threshold <- 50
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

cmip_pred_final <- left_join(cmip_pred_best, gfed_combined) %>%
  mutate(DM_level = DM_anomaly + DM_2001_2021) %>%
  mutate(DM_level = pmin(DM_level, max_threshold)) %>%
  rename(pred_DM_raw=pred_DM,
         pred_DM_final=DM_level) %>%
  group_by(scenario, region, region_sub, gcm, algorithm_outcome) %>%
  arrange(year) %>%
  mutate(pred_raw_10yr=rollapplyr(pred_DM_raw, width = 10, FUN=mean, partial=T),
         pred_final_10yr=rollapplyr(pred_DM_final, width = 10, FUN=mean, partial=T)) %>%  ## calculating the 10-yr mean
  filter(scenario=="historical" | year>2021) %>% ungroup() ### to drop the fake 2001-2021 for the SSP scenarios

### top-coded 
cmip_pred_final <- cmip_pred_final %>%
         mutate(pred_final_10yr_positive=pmax(0, pred_final_10yr))  
  
saveRDS(cmip_pred_final, paste0(result_path,"/CMIP6_fire_proj/final_pred_smoke_2001_2099_annualmax50.rds"))

### for the subregion
# cmip_pred_final <- filter(cmip_pred_final, !(region=="Western US" & algorithm_outcome=="Neural Net, log")) 
# saveRDS(cmip_pred_final, paste0(result_path,"/CMIP6_fire_proj/final_pred_best_models_subregion.rds"))

### look at the in-sample performance to make sure everything is correct 
# obs1 <- readRDS(paste0(result_path,"/CMIP6_fire_proj/linear_final_pred_history.rds")) 
# obs2 <- readRDS(paste0(result_path,"/CMIP6_fire_proj/LASSO_final_pred_history.rds"))
# obs3 <- readRDS(paste0(result_path,"/CMIP6_fire_proj/nnet_final_pred_history.rds")) %>%
#   mutate(year=as.numeric(as.character(year)))
# 
# obs <- bind_rows(obs1, obs2, obs3) %>%
#   mutate(spatial_res=replace(spatial_res, spatial_res=="region", "regional")) %>%
#   left_join(best_models, by=c("region", "algorithm", "outcome", "spatial_res", "spec"))
# 
# obs_region <- obs %>%
#   filter(selected == 1)  %>%
#   mutate(pred_DM=pred_DM*as.numeric(pred_DM >= 0)) %>%
#   group_by(algorithm, outcome, scenario, year, region) %>%
#   summarise(pred_DM = sum(pred_DM, na.rm=T)) %>% ungroup()
# 
# obs_region$algorithm_outcome <- paste(obs_region$algorithm, obs_region$outcome, sep=", ")
# 
# obs_region <- obs_region %>%
#   group_by(scenario, region, algorithm_outcome) %>%
#   arrange(year) %>%
#   mutate(pred_10yr=rollapplyr(pred_DM, width = 10, FUN=mean, partial=T))  %>%
#   ungroup() 


#---------------------------------------------------------------
### ---- First calculate 10-yr mean, then take median across GCMs
#---------------------------------------------------------------
## just to connect 2021 and 2022 (for visualization purpose)
cmip_region <- cmip_pred_final %>%
  group_by(year, scenario, gcm, region, algorithm_outcome) %>%
  summarise_at(c("pred_DM_raw", "pred_DM_final", "pred_raw_10yr", "pred_final_10yr", "pred_final_10yr_positive"),
               .funs=sum, na.rm=T) %>% ungroup()

cmip_region <- bind_rows(cmip_region,
                          filter(cmip_region, year==2021) %>% mutate(scenario = "ssp126"),
                          filter(cmip_region, year==2021) %>% mutate(scenario = "ssp245"),
                          filter(cmip_region, year==2021) %>% mutate(scenario = "ssp370"))

### calculate the projection results from each algorithm 
variable <- "pred_final_10yr_positive"
cmip_algorithm <-  cmip_region %>%
  group_by(scenario, year, region, algorithm_outcome) %>%
  summarise(pred_median=median(!!as.name(variable), na.rm=T), 
            pred_mean=mean(!!as.name(variable), na.rm=T), 
            pred_10=quantile(!!as.name(variable), 0.1, na.rm=T), pred_90=quantile(!!as.name(variable), 0.9, na.rm=T),
            pred_25=quantile(!!as.name(variable), 0.25, na.rm=T),pred_75=quantile(!!as.name(variable), 0.75, na.rm=T)) %>% 
  ungroup() 


### for plotting purposes to calculate 10-yr mean of obs data
gfed_region_summ <- readRDS(paste0(result_path, "/final_models/training_data/final_train_region_FE_level.rds"))$data_region  %>%
  group_by(region) %>%
  arrange(year) %>%
  mutate(DM_10yr=rollapplyr(DM_kg, width = 10, FUN=mean, partial=T)) %>% ungroup()

region_list <- c("Western US", "Canada-Alaska", "Mexico", "Southeastern US", "Northeastern US") ##unique(gfed_data$region)
region_label_list <- c("west", "canada", "mexico", "southeast", "northeast")
for (rrr in 1:length(region_list)){
region_tmp <- region_list[rrr]  
region_label <- region_label_list[rrr]  

range <-  range(c(filter(gfed_region_summ, region==region_tmp, year>2009)$DM_10yr,
                  #filter(obs_region, region==region_tmp, year>2009)$pred_10yr,
                  filter(cmip_algorithm, year>2000, year<2056,
                         region==region_tmp)$pred_median))*1e-9

nalgorithms <- unique(filter(cmip_algorithm, region==region_tmp)$algorithm_outcome) %>% length()
width <- 2.5*nalgorithms + 2

ggplot(filter(cmip_algorithm, year>2000,
              region==region_tmp),
       aes(x=year, y= pred_median*1e-9,
           colour=scenario, fill=scenario)) +
  geom_line(size=1.4) +
  # geom_point(data = filter(gfed_region_summ, region==region_tmp, year>2009),
  #             aes(x=year, y=DM_10yr*1e-9), colour="black", fill=NA, size=1.5, shape=1) +
  # geom_hline(aes(yintercept = 1e-9*filter(gfed_region_summ, region==region_tmp, year==2021)$DM_10yr), 
  #            linetype="dashed") +
  # geom_hline(aes(yintercept = 1e-9*mean(filter(gfed_region_summ, region==region_tmp)$DM_kg, na.rm=T)), 
  #            linetype="dashed") +
  # geom_point(data = filter(obs_region, region==region_tmp, year>2009),
  #            aes(x=year, y=pred_10yr*1e-9), colour="green4", fill=NA, size=1.5, shape=2) +
  scale_fill_manual(values = color_map) +
  scale_colour_manual(values = color_map) +
  theme_classic() + theme(text = element_text(size=14),
                          axis.text.x = element_text(size=12)) +
  labs(x="", fill="", colour="", y= "Predicted emissions (million tons)") +
  facet_wrap(~algorithm_outcome, nrow = 1) +
  #scale_y_continuous(limits= range) +
  scale_x_continuous(breaks=seq(2000, 2050, by=10))
ggsave(paste0(region_label,"_CMIP6_10yr_median_algorithm_subregion.png"), width=width, height=2.5)
 }

#-------------------------------------
##     average overall algorithms
#-------------------------------------
### calculate the ensemble mean across different algorithms
cmip_algorithm_mean <-  cmip_region %>%
  group_by(scenario, year, region, gcm) %>%
  summarise_at(c("pred_DM_raw", "pred_DM_final", "pred_raw_10yr", "pred_final_10yr", "pred_final_10yr_positive"),
               .funs=mean, na.rm=T) %>% ungroup() %>%
  group_by(scenario, year, region) %>% 
  summarise(pred_median=median(!!as.name(variable), na.rm=T), 
            pred_mean=mean(!!as.name(variable), na.rm=T), 
            pred_10=quantile(!!as.name(variable), 0.1, na.rm=T), pred_90=quantile(!!as.name(variable), 0.9, na.rm=T),
            pred_25=quantile(!!as.name(variable), 0.25, na.rm=T),pred_75=quantile(!!as.name(variable), 0.75, na.rm=T)) %>% 
  ungroup() 

region_list <- c("Western US", "Canada-Alaska", "Mexico", "Southeastern US", "Northeastern US") ##unique(gfed_data$region)
region_label_list <- c("west", "canada", "mexico", "southeast", "northeast")
for (rrr in 1:length(region_list)){
  region_tmp <- region_list[rrr]  
  region_label <- region_label_list[rrr]  
  
  range <-  range(c(mean(filter(gfed_region_summ, region==region_tmp)$DM_kg),
                    filter(cmip_algorithm_mean, year>2000, year<2056,
                           region==region_tmp)$pred_median))*1e-9
  
  ggplot(filter(cmip_algorithm_mean, year>2000, year<2056, 
                region==region_tmp),
         aes(x=year, y= pred_median*1e-9,
             colour=scenario, fill=scenario)) +
    geom_line(size=1.4) +
    # geom_point(data = filter(gfed_region_summ, region==region_tmp, year>2009),
    #            aes(x=year, y=DM_10yr*1e-9), colour="black", fill=NA, size=1.5, shape=1) +
    geom_hline(aes(yintercept = 1e-9*mean(filter(gfed_region_summ, region==region_tmp)$DM_kg, na.rm=T)), 
               linetype="dashed") +
    scale_fill_manual(values = color_map) +
    scale_colour_manual(values = color_map) +
    theme_classic() + theme(text = element_text(size=14),
                            axis.text.x = element_text(size=14)) +
    labs(x="", fill="", colour="", y= "Predicted emissions (million tons)") +
    scale_y_continuous(limits= range) +
    scale_x_continuous(breaks=seq(2000, 2050, by=10)) #+
    #facet_wrap(~region)
  ggsave(paste0(region_label,"_CMIP6_10yr_median_subregion.pdf"), width=5.5, height=3)
  }



#### combine northeast and southeast
cmip_east <-  cmip_region %>%
  filter(region %in% c("Southeastern US", "Northeastern US")) %>%
  group_by(scenario, year, region, gcm) %>%
  summarise_at(c("pred_DM_raw", "pred_DM_final", "pred_raw_10yr", "pred_final_10yr", "pred_final_10yr_positive"),
               .funs=mean, na.rm=T) %>% ungroup() %>%
  group_by(scenario, year, gcm) %>%
  summarise_at(c("pred_DM_raw", "pred_DM_final", "pred_raw_10yr", "pred_final_10yr", "pred_final_10yr_positive"),
               .funs=sum, na.rm=T) %>% ungroup() %>%
  group_by(scenario, year) %>% 
  summarise(pred_median=median(!!as.name(variable), na.rm=T), 
            pred_mean=mean(!!as.name(variable), na.rm=T), 
            pred_10=quantile(!!as.name(variable), 0.1, na.rm=T), pred_90=quantile(!!as.name(variable), 0.9, na.rm=T),
            pred_25=quantile(!!as.name(variable), 0.25, na.rm=T),pred_75=quantile(!!as.name(variable), 0.75, na.rm=T)) %>% 
  ungroup() 

gfed_east <- filter(gfed_region_summ, region %in% c("Southeastern US", "Northeastern US")) %>%
  group_by(year) %>%
  summarise(DM_kg=sum(DM_kg, na.rm=T)) %>%
  ungroup()

range <-  range(c(mean(gfed_east$DM_kg),
                  filter(cmip_east, year>2000, year<2056)$pred_median))*1e-9

ggplot(filter(cmip_east, year>2000, year<2056),
       aes(x=year, y= pred_median*1e-9,
           colour=scenario, fill=scenario)) +
  geom_line(size=1.4) +
  geom_hline(aes(yintercept = 1e-9*mean(gfed_east$DM_kg, na.rm=T)), 
             linetype="dashed") +
  scale_fill_manual(values = color_map) +
  scale_colour_manual(values = color_map) +
  theme_classic() + theme(text = element_text(size=14),
                          axis.text.x = element_text(size=14)) +
  labs(x="", fill="", colour="", y= "Predicted emissions (million tons)") +
  scale_y_continuous(limits= range) +
  scale_x_continuous(breaks=seq(2000, 2050, by=10)) 
ggsave("east_CMIP6_10yr_median.pdf", width=5.5, height=3)


ggplot(filter(cmip_combine, year>2045, year<2050, scenario =="ssp370", algorithm=="Linear (region, log)")) +
  geom_histogram(aes(x=pred_10yr, group=factor(year), fill=factor(year)), 
            size=1) +
  geom_vline(data=filter(cmip_10yr_algorithm, year>2045, year<2050, 
                         scenario =="ssp370", algorithm=="Linear (region, log)"), 
             aes(xintercept = pred_median)) +
  facet_wrap(~year, ncol=1)

### show projections from all gcms
ggplot(filter(cmip_combine, year>2000, year<2056, algorithm=="Linear (eco2, log)",
              scenario!="ssp585"),
       aes(x=year, y= pred_10yr*1e-9, colour=scenario)) +
  geom_line(size=1) +
  geom_vline(xintercept = c(2046, 2049), linetype="dashed") +
  #geom_ribbon(aes(ymin=anomaly_25, ymax=anomaly_75), alpha=0.3, colour=NA) +
  # scale_fill_manual(values = color_map) +
  scale_colour_manual(values = color_map) +
  theme_classic() + theme(text = element_text(size=16)) +
  labs(x="", fill="", colour="") +
  facet_wrap(~gcm)
ggsave("west_fire_pred_linear_eco2_gcm_20x.png", width=12, height=8)


  geom_hline(aes(yintercept = 1e-9*filter(gfed_region, region=="Western US", year==2021)$DM_10yr), linetype="dashed") +
  #geom_hline(aes(yintercept = 1e-9*filter(gfed_region, region=="Western US", year==2020)$DM_kg)) +
  geom_vline(xintercept = c(2046, 2049), linetype="dashed") +
  scale_fill_manual(values = color_map) +
  scale_colour_manual(values = color_map) +
  theme_classic() + theme(text = element_text(size=16),
                          axis.text.x = element_text(size=13)) +
  labs(x="", fill="", colour="", y= "Predicted emissions (million tons)") +
  facet_wrap(~algorithm, ncol=2) +
  #scale_y_continuous(limits=c(-10, 520)) +
  scale_x_continuous(breaks=seq(2000, 2050, by=10))
ggsave("west_fire_CMIP6_10yr_linear_sens_nosoil.png", width=8, height=4.5)

#---------------------------------------------------------------
### ---- First take median across GCMs then calculate 10-yr mean
#---------------------------------------------------------------
# variable <-  "DM_level" ##"pred_DM" "DM_anomaly"
# cmip_combine <- bind_rows(filter(cmip_region, year>2021),
#                           cmip_2001_2021 %>% mutate(scenario = "historical", DM_level=pred_DM),
#                           cmip_2001_2021 %>% mutate(scenario = "ssp126", DM_level=pred_DM),
#                           cmip_2001_2021 %>% mutate(scenario = "ssp245", DM_level=pred_DM),
#                           cmip_2001_2021 %>% mutate(scenario = "ssp370", DM_level=pred_DM),
#                           cmip_2001_2021 %>% mutate(scenario = "ssp585", DM_level=pred_DM)) %>%
#   group_by(scenario, region, algorithm, year) %>%
#   summarise(pred_median=median(!!as.name(variable), na.rm=T)) %>%
#   filter(scenario=="historical" | year>2021) %>% ungroup() ### to drop the fake 2001-2021 for the SSP scenarios
# 
# ## just to connect 2021 and 2022 (for visualization purpose)
# cmip_combine <- bind_rows(cmip_combine,
#                           filter(cmip_combine, year==2021) %>% mutate(scenario = "ssp126"),
#                           filter(cmip_combine, year==2021) %>% mutate(scenario = "ssp245"),
#                           filter(cmip_combine, year==2021) %>% mutate(scenario = "ssp370"),
#                           filter(cmip_combine, year==2021) %>% mutate(scenario = "ssp585"))
# 
# ### clauclate the projection results from each algorithm 
# cmip_10yr_algorithm <-  cmip_combine %>%
#   group_by(scenario, region, algorithm) %>%
#   arrange(year) %>%
#   mutate(pred_10yr=rollapplyr(pred_median, width = 10, FUN=mean, partial=T)) %>%  ## just for calculating the 10-yr mean
#   ungroup() 
# 
# cmip_10yr_algorithm$algorithm <- factor(cmip_10yr_algorithm$algorithm, levels=c("Neural Network","XGBOOST","LASSO (level)",
#                                                                                 "LASSO (log)", "Linear (level)", "Linear (log)"))
# 
# ggplot(filter(cmip_10yr_algorithm, year>2000, year<2056, scenario !="ssp585"), 
#        aes(x=year, y= pred_10yr*1e-9,
#            colour=scenario, fill=scenario)) +
#   geom_line(size=1.4) +
#   geom_point(data = filter(gfed_region, region=="Western US"), 
#               aes(x=year, y=DM_10yr*1e-9), colour="black", fill=NA, size=2.5) +
#   geom_hline(aes(yintercept = 1e-9*filter(gfed_region, region=="Western US", year==2021)$DM_10yr), linetype="dashed") +
#   #geom_hline(aes(yintercept = 1e-9*filter(gfed_region, region=="Western US", year==2020)$DM_kg)) +
#   scale_fill_manual(values = color_map) +
#   scale_colour_manual(values = color_map) +
#   theme_classic() + theme(text = element_text(size=16),
#                           axis.text.x = element_text(size=13)) +
#   labs(x="", fill="", colour="", y= "Predicted emissions (million tons)") +
#   facet_wrap(~algorithm, ncol=2) +
#   scale_y_continuous(limits=c(-10, 405)) +
#   scale_x_continuous(breaks=seq(2000, 2050, by=10))
# ggsave("west_fire_CMIP6_10yr_algorithm_median1st.png", width=8, height=6.5)





filter(cmip_10yr_smooth, scenario !="ssp585", year==2055) %>%
  mutate(ratio=pred_median/filter(gfed_region, region=="Western US", year==2021)$DM_10yr) %>%
  select(scenario, algorithm, ratio) 

filter(cmip_10yr_smooth, scenario !="ssp585", year==2055) %>%
  mutate(ratio=pred_75/filter(gfed_region, region=="Western US", year==2021)$DM_10yr) %>%
  select(scenario, algorithm, ratio) 

filter(cmip_10yr_smooth, scenario !="ssp585", year==2055) %>%
  mutate(ratio=pred_25/filter(gfed_region, region=="Western US", year==2021)$DM_10yr) %>%
  select(scenario, algorithm, ratio) 

ggplot(filter(cmip_10yr_smooth, year<2061, scenario !="ssp585"), 
       aes(x=year, y= pred_median, colour=scenario, fill=scenario)) +
  geom_line(size=1) +
  #geom_ribbon(aes(ymin=pred_25, ymax=pred_75), alpha=0.5, colour=NA) +
  geom_point(data = filter(gfed_region, region=="Western US"), 
             aes(x=year, y=DM_10yr), colour="black", fill=NA, size=2) +
  # geom_point(data = filter(gfed_region, region=="Western US"), 
  #            aes(x=year, y=DM_kg), colour="black", fill=NA, size=2, alpha=0.3) +
  scale_fill_manual(values = color_map) +
  scale_colour_manual(values = color_map) +
  theme_classic() + theme(text = element_text(size=16)) +
  labs(x="", fill="", colour="", y= "Predcited fire emissions (kg)") +
  facet_wrap(~region, scales = "free") +
  scale_y_continuous(limits=c(0,1.25e11))
ggsave("west_fire_CMIP6_10yr_level_linear_vpd.png", width=10, height=5)

cmip_gcm_summ <- fire_region  %>% group_by(scenario, year, region) %>%
  summarise(pred_median=median(pred_DM), 
            pred_10=quantile(pred_DM, 0.1, na.rm=T),
            pred_90=quantile(pred_DM, 0.9),
            pred_25=quantile(pred_DM, 0.25),
            pred_75=quantile(pred_DM, 0.75)) %>% ungroup()

## just to connect 2014 and 2015 (for visualization purpose)
cmip_gcm_summ <- bind_rows(cmip_gcm_summ,
                             filter(cmip_gcm_summ, year==2014) %>% mutate(scenario = "ssp126"),
                             filter(cmip_gcm_summ, year==2014) %>% mutate(scenario = "ssp245"),
                             filter(cmip_gcm_summ, year==2014) %>% mutate(scenario = "ssp370"),
                             filter(cmip_gcm_summ, year==2014) %>% mutate(scenario = "ssp585"))
cmip_gcm_summ$region <- factor(cmip_gcm_summ$region,
                                 levels=c("Canada-Alaska", "Northeastern US", "Southeastern US",
                                          "Western US", "Mexico"))

ggplot(filter(cmip_gcm_summ, scenario!="ssp585"), 
       aes(x=year, y= pred_median, colour=scenario, fill=scenario)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=pred_25, ymax=pred_75), alpha=0.5, colour=NA) +
  geom_point(data = filter(gfed_region, region=="Western US"), aes(x=year, y=DM_kg), colour="black", fill=NA, size=2) +
  scale_fill_manual(values = color_map) +
  scale_colour_manual(values = color_map) +
  theme_classic() + theme(text = element_text(size=16)) +
  labs(x="", fill="", colour="", y= "Predcited fire emissions") +
  facet_wrap(~region, scales = "free")
ggsave(paste0("west_fire_CMIP6_xgboost.png"), width=10, height=5)

#------------- Load GFED data 
gfed_grid <- readRDS(paste0(data_path, "/GFED4S_NorthAmerica_025deg_grid_EPSG4326.rds"))

gfed_coregion_link <- readRDS(paste0(db_path, "/crosswalk_GFED_NA_CEC_Eco_Level2_3.rds")) %>%
  mutate(intersection_area=as.numeric(intersection_area)) %>%   group_by(gfed_cell_id)  %>%
  arrange(desc(intersection_area)) %>% distinct(gfed_cell_id, .keep_all=T)

gfed_wang_region <- readRDS(paste0(data_path, "/crosswalk_GFED_wang_regions.rds")) %>%
  mutate(region=as.character(region)) %>%
  mutate(region=replace(region, region%in%c("mediterranean california", "southwestern US"), "western forest area"),
         region=replace(region, (region=="northeastern US" & gfed_cell_id<30000), "southeastern US")) %>%
  mutate(region=replace(region, region=="western forest area", "Western US"),
         region=replace(region, region=="southeastern US", "Southeastern US"),
         region=replace(region, region=="northeastern US", "Northeastern US"))

gfed_data <- readRDS(file.path(db_path, "climate-fire_training_data_1997-2021_clean_varname.rds"))
gfed_data <- left_join(gfed_data, gfed_wang_region, by="gfed_cell_id")
gfed_data <- left_join(gfed_data, gfed_coregion_link, by="gfed_cell_id")

crosswalk_gfed_1deg <- readRDS(paste0(data_path, "/crosswalk_gfed_1deg.rds"))
crosswalk_1deg_region <- left_join(crosswalk_gfed_1deg, gfed_wang_region) %>%
  group_by(grid_1deg_id, region) %>% summarise(coverage_area =sum(coverage_area)) %>%
  ungroup()

#------------- Load debiased CMIP6 data 
var_long_list <- c("temperature_surface", "precipitation", "RH_surface", "vpd","soil_moisture", "runoff", "wind_speed")
var_labels <- c("temperature (degree C)", 
                "precipitation (kg/m2)", 
                "RH (%)", 
                "vpd (kPa)",
                "soil moisture (kg/m2)", 
                "runoff (kg/m2)", 
                "wind speed (m/s)")
var_short_list <- c("tas", "pr", "hurs", "vpd", "mrsos", "mrro", "sfcWind")

color_map <- c("grey65",
               rev(met.brewer(name="Homer1", n=8, type="discrete"))[c(2,5,6,8)])

for (vv in 1: length(var_long_list)){
  var_long <- var_long_list[vv]  
  var_short <- var_short_list[vv]  
  print(var_short)
  
  # narr <- readRDS(paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/",
  #                          var_short,"_narr_1deg_1997_2021.rds"))
  
  cmip_debias <- readRDS(paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/",
                                var_long, "/CMIP6_annual_debias_1deg.rds")) %>%
    mutate(ratio = cmip_value / narr_value)
  print(quantile(cmip_debias$ratio, seq(0,1,0.1)))
  #grid_list <- intersect(unique(narr$grid_1deg_id), unique(cmip_debias$grid_1deg_id))
  print(length(unique(cmip_debias$model)))
  
  cmip_region <- filter(cmip_debias) %>%
    left_join(crosswalk_1deg_region) %>% filter(!is.na(region)) %>%
    group_by(region, year, model, scenario) %>%
    summarise(cmip_debias = weighted.mean(cmip_debias, w= coverage_area, na.rm=T),
              narr_value = weighted.mean(narr_value, w= coverage_area, na.rm=T)) %>%
    ungroup() %>%
    mutate(cmip_anomaly = cmip_debias - narr_value)
  
  # narr_mean <- filter(narr, grid_1deg_id %in% grid_list) %>% filter(year>1996, year<2015) %>%
  #   rename(narr_value = !!as.name(var_short)) %>%
  #   group_by(grid_1deg_id) %>%
  #   summarise(narr_value = mean(narr_value, na.rm=T)) %>% ungroup()
  # 
  # narr_region <- left_join(narr_mean, crosswalk_1deg_region) %>% filter(!is.na(region)) %>%
  #   group_by(region) %>%
  #   summarise(narr_value = weighted.mean(narr_value, w= coverage_area, na.rm=T)) %>%
  #   ungroup()
  
  # cmip_region <- left_join(cmip_region, narr_region) %>%
  #   mutate(cmip_anomaly = cmip_debias - narr_value)
  
  print(filter(cmip_region, scenario=="historical") %>% group_by(model) %>% 
          summarise(bias_history=mean(cmip_anomaly, na.rm=T)))
  
  print(filter(cmip_debias, scenario=="historical") %>% group_by(model) %>% 
          summarise(bias=mean(bias, na.rm=T)))
  
    
}
