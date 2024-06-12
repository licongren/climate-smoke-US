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

# color_scenario <- c("grey65",
#                rev(met.brewer(name="Homer1", n=8, type="discrete"))[c(2,6,8)])

# color_scenario <- c("grey65", "navy", "darkorange","#CB172F")
color_scenario <- c("grey65", "#1F2E52", "#E78928","#CB172F") # "#38A2C6", "#870A21"

#-------------------------------------------------------------------------------
# Plot the figure 1
# Written by Minghao
# Last edited Jan 2024
#-------------------------------------------------------------------------------

#------------------------------------------------
### figure 1A: the performance of the fire models
#------------------------------------------------
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
results_select_mean <- results_select %>%
  filter(diff_min<threshold) %>%
  group_by(region, experiment, model, outcome, spec) %>%
  mutate(year = 1:n()) %>%
  group_by(region, year) %>%
  summarise(truth = mean(truth, na.rm=T),
            response_positive = mean(response_positive, na.rm=T)) %>%
  ungroup()

results_select_mean %>% group_by(region) %>%
summarise(r=cor(response_positive, truth),
  rmse_pos=sqrt(mean((response_positive - truth)^2))/mean(truth)) %>% ungroup()

for (rrr in unique(best_metric$region)){
  #rrr <- "Canada-Alaska"
  rrr <- "Northeastern US"
  axis_range <- range(c(filter(results_select_mean, region==rrr)$response_positive,
                        filter(results_select_mean, region==rrr)$truth))/1e9
  ggplot(filter(results_select_mean, region==rrr),
         aes(x=truth/1e9, y=response_positive/1e9)) +
    geom_point(size=4, alpha=0.7) +
    geom_abline(,0,1,linetype="dashed") +
    theme_classic() + theme(text = element_text(size=17), axis.text = element_text(size=17)) + 
    scale_x_continuous(limits= axis_range, breaks=c(1.8,2.0,2.2)) +
    scale_y_continuous(limits= axis_range, breaks=c(1.8,2.0,2.2)) +
    labs(y="Predictions\n(DM emissions, MT)", x="Observations\n(DM emissions, MT)") +
    #scale_colour_manual(values=c("grey30", "green4")) +
    guides(colour= "none")
  ggsave(paste0(rrr,"_best_model_threshold05.pdf"), width=3.1, height=3.1)
}

#------------------------------------------------
### figure 1B: projected wildfire emissions (time series)
#------------------------------------------------
fire_region <-  readRDS(paste0(result_path,"/CMIP6_fire_proj/final_pred_best_models_threshold_05.rds")) %>%
  group_by(year, scenario, gcm, region, algorithm_outcome) %>%
  summarise_at(c("pred_DM_raw", "pred_DM_final", "pred_raw_10yr", "pred_final_10yr", "pred_final_10yr_positive"),
               .funs=sum, na.rm=T) %>% ungroup()

fire_region <- bind_rows(fire_region,
                         filter(fire_region, year==2021) %>% mutate(scenario = "ssp126"),
                         filter(fire_region, year==2021) %>% mutate(scenario = "ssp245"),
                         filter(fire_region, year==2021) %>% mutate(scenario = "ssp370"))

fire_emis_summ <- fire_region %>%
  group_by(scenario, year, region, gcm) %>%
  summarise(fire_emis=mean(pred_final_10yr_positive, na.rm=T)) %>%## take mean across algortihms
  ungroup() %>%
  mutate(region=replace(region, region%in%c("Northeastern US","Southeastern US"), "Eastern US"))  %>% 
  group_by(gcm, scenario, year, region) %>%
  summarise(fire_emis=sum(fire_emis, na.rm=T)) %>% ungroup() %>%
  group_by(scenario, year, region) %>% ## take mean or median across GCMs
  summarise(fire_emis=median(fire_emis)) %>% ungroup() %>%
  filter(year>2009, year<2056)

gfed_region_summ <- readRDS(paste0(result_path, "/final_models/training_data/final_train_region_FE_level.rds"))$data_region  %>%
  mutate(region=replace(region, region%in%c("Northeastern US","Southeastern US"), "Eastern US")) %>%
   group_by(region, year) %>%
  summarise(DM_kg=sum(DM_kg)) %>%
  arrange(year) %>%
  mutate(DM_10yr=rollapplyr(DM_kg, width = 10, FUN=mean, partial=T)) %>% ungroup()

region_list <- c("Western US", "Canada-Alaska", "Mexico", "Eastern US") ##unique(gfed_data$region)
region_label_list <- c("west", "canada", "mexico", "east")
range <- data.frame(min=c(0, 100, 35, 13),
                    max=c(352, 290, 60, 20))

for (rrr in 1:length(region_list)){
  region_tmp <- region_list[rrr]  
  region_label <- region_label_list[rrr]  
  
  range_tmp <-  unlist(range[rrr,])
  
  ggplot(filter(fire_emis_summ, region==region_tmp),
         aes(x=year, y= fire_emis*1e-9, colour=scenario, fill=scenario)) +
    geom_line(size=1.5) +
    # geom_point(data = filter(gfed_region_summ, region==region_tmp, year>2009),
    #            aes(x=year, y=DM_10yr*1e-9), colour="black", fill=NA, size=1.5, shape=1) +
    geom_hline(aes(yintercept = 1e-9*mean(filter(gfed_region_summ, region==region_tmp)$DM_kg, na.rm=T)), 
               linetype="dashed") +
    scale_fill_manual(values = color_scenario) +
    scale_colour_manual(values = color_scenario) +
    theme_classic() + theme(text = element_text(size=16),
                            axis.text.x = element_text(size=16),
                            axis.text.y = element_text(size=16)) +
    labs(x="", fill="", colour="", y= "Emissions (Mtons)") +
    scale_y_continuous(limits= range_tmp) +
    scale_x_continuous(breaks=seq(2000, 2050, by=10)) 
  ggsave(paste0(region_label,"_CMIP6_10yr_median.pdf"), width=5.5, height=3)
  
}

fire_emis_summ %>% filter(year==2055) %>%
  left_join(gfed_region_summ %>% group_by(region) %>% 
              summarise(DM_mean=sum(DM_kg*as.numeric(year%in%2011:2020)/10),
                        DM_max=max(DM_kg))) %>%
  ungroup() %>%
  mutate(mean_ratio=fire_emis/DM_mean,
         max_ratio=fire_emis/DM_max) %>%
  filter(region=="Mexico")
#------------------------------------------------
### figure 1C: plot downscaled fire emissions
#------------------------------------------------
proj_emis_grid <- readRDS(paste0(result_path,"/CMIP6_fire_proj/downscaled_emis_final_threshold_05.rds")) %>%
  group_by(scenario, year, gfed_cell_id) %>%
  summarise(pred_DM_grid=mean(pred_DM_grid)) %>% ### take mean across GCMs
  ungroup()

gfed_2001_2021 <- readRDS(file.path(db_path, "climate-fire_training_data_2001-2021_clean_varname.rds"))  %>%
  group_by(year, gfed_cell_id) %>%
  summarise(DM_kg=sum(DM_kg, na.rm=T)) %>%
  group_by(gfed_cell_id) %>%
  summarise(DM_2001_2021 = mean(DM_kg, na.rm=T)) %>% ungroup() 

gfed_grid <- readRDS(paste0(data_path, "/GFED4S_NorthAmerica_025deg_grid_EPSG4326.rds"))
gfed_wang_region <- readRDS(paste0(data_path, "/crosswalk_GFED_wang_regions.rds")) 
gfed_grid <- left_join(gfed_grid, gfed_wang_region) %>% filter(!is.na(region)) %>%
  mutate(region=replace(region, region%in%c("Northeastern US","Southeastern US"), "Eastern US")) 
  
proj_emis_2050 <-  right_join(gfed_grid, proj_emis_grid) %>% filter(year==2055) %>%
  mutate(DM_kton = pred_DM_grid / 1e6)
gfed_2001_2021 <-  right_join(gfed_grid, gfed_2001_2021) %>%
  mutate(DM_kton = DM_2001_2021 / 1e6)

#### plot the emissions by subregions 
region_list <- c("Western US", "Canada-Alaska", "Mexico", "Eastern US") ##unique(gfed_data$region)
region_label_list <- c("west", "canada", "mexico", "east")

col_disc <- c("grey70","#c3f4f6", "steelblue3","green4",
              "#f5b642","orangered3","#551f00")

#### West
region_tmp <- "Western US"
tmp_2050 <- filter(proj_emis_2050, region == region_tmp) %>% filter(scenario%in%c("ssp370", "ssp126"))
tmp_2001_2021 <- filter(gfed_2001_2021, region == region_tmp) %>% mutate(scenario="2001-2021")
tmp_emis <- bind_rows(tmp_2050, tmp_2001_2021)


max_emis <- max(tmp_2050$DM_kton, na.rm=T)
label_emis <- c(0,1,5,10,100,250,1000,max_emis)
tmp_emis$DM_cut <- cut(tmp_emis$DM_kton, label_emis)
tmp_emis$DM_cut <- factor(tmp_emis$DM_cut, levels= rev(levels(tmp_emis$DM_cut)))

ggplot(tmp_emis, aes(colour=DM_cut, fill=DM_cut)) +
  geom_sf() + 
  scale_color_manual(values=rev(col_disc)) +
  scale_fill_manual(values=rev(col_disc)) +
  theme_classic() + theme(text = element_text(size=14),
                          axis.ticks = element_blank(),
                          axis.line = element_blank(),
                          axis.text = element_blank()) +
  labs(colour="DM (1000 tons)", fill="DM (1000 tons)") +
  facet_wrap(~scenario, nrow=1)
ggsave("proj_median_DM_map_west_allscenarios.pdf", width=8.5, height=6)


#### Canada
region_tmp <- "Canada-Alaska"
tmp_2050 <- filter(proj_emis_2050, region == region_tmp) %>% filter(scenario=="ssp370")
tmp_2001_2021 <- filter(gfed_2001_2021, region == region_tmp) %>% mutate(scenario="2001-2021")
tmp_emis <- bind_rows(tmp_2050, tmp_2001_2021) 

max_emis <- max(tmp_2050$DM_kton, na.rm=T)
label_emis <-  c(0,1,5,10,50,100,500,max_emis) 
tmp_emis$DM_cut <- cut(tmp_emis$DM_kton, label_emis, include.lowest = T)
tmp_emis$DM_cut <- factor(tmp_emis$DM_cut, levels= rev(levels(tmp_emis$DM_cut)))

ggplot(tmp_emis, aes(colour=DM_cut, fill=DM_cut)) +
  geom_sf() + 
  scale_color_manual(values=rev(col_disc)) +
  scale_fill_manual(values=rev(col_disc)) +
  theme_classic() + theme(text = element_text(size=14),
                          axis.ticks = element_blank(),
                          axis.line = element_blank(),
                          axis.text = element_blank()) +
  labs(colour="DM (1000 tons)", fill="DM (1000 tons)") +
  facet_wrap(~scenario)
ggsave("proj_median_DM_map_canada.pdf", width=8, height=6)


#### Eastern US
region_tmp <- "Eastern US"
tmp_2050 <- filter(proj_emis_2050, region == region_tmp) %>% filter(scenario=="ssp370")
tmp_2001_2021 <- filter(gfed_2001_2021, region == region_tmp) %>% mutate(scenario="2001-2021")
tmp_emis <- bind_rows(tmp_2050, tmp_2001_2021)

max_emis <- max(tmp_emis$DM_kton, na.rm=T)
label_emis <-  c(0,1,5,10,20,40,60, max_emis)
tmp_emis$DM_cut <- cut(tmp_emis$DM_kton, label_emis, include.lowest = T)
tmp_emis$DM_cut <- factor(tmp_emis$DM_cut, levels= rev(levels(tmp_emis$DM_cut)))

ggplot(tmp_emis, aes(colour=DM_cut, fill=DM_cut)) +
  geom_sf() + 
  scale_color_manual(values=rev(col_disc)) +
  scale_fill_manual(values=rev(col_disc)) +
  theme_classic() + theme(text = element_text(size=14),
                          axis.ticks = element_blank(),
                          axis.line = element_blank(),
                          axis.text = element_blank()) +
  labs(colour="DM (1000 tons)", fill="DM (1000 tons)") +
  facet_wrap(~scenario)
ggsave("proj_median_DM_map_east.pdf", width=8, height=6)


#### Mexico
region_tmp <- "Mexico"
tmp_2050 <- filter(proj_emis_2050, region == region_tmp) %>% filter(scenario=="ssp370")
tmp_2001_2021 <- filter(gfed_2001_2021, region == region_tmp) %>% mutate(scenario="2001-2021")
tmp_emis <- bind_rows(tmp_2050, tmp_2001_2021)

max_emis <- max(tmp_emis$DM_kton, na.rm=T)
label_emis <- c(0,1,5,10,25,50,100,max_emis)
tmp_emis$DM_cut <- cut(tmp_emis$DM_kton, label_emis, include.lowest = T)
tmp_emis$DM_cut <- factor(tmp_emis$DM_cut, levels= rev(levels(tmp_emis$DM_cut)))

ggplot(tmp_emis, aes(colour=DM_cut, fill=DM_cut)) +
  geom_sf() + 
  scale_color_manual(values=rev(col_disc)) +
  scale_fill_manual(values=rev(col_disc)) +
  theme_classic() + theme(text = element_text(size=14),
                          axis.ticks = element_blank(),
                          axis.line = element_blank(),
                          axis.text = element_blank()) +
  labs(colour="DM (1000 tons)", fill="DM (1000 tons)") +
  facet_wrap(~scenario)
ggsave("proj_median_DM_map_mexico.pdf", width=8, height=5)


#---------------------------------------------------------------------------------------------------
# Figure Sx: performance of the stat models across spatial scales (an example of Western US)
#---------------------------------------------------------------------------------------------------
linear_results <- readRDS(paste0(result_path,"/model_eval_final/linear_model_results_LOOCV.rds"))
lasso_results <- readRDS(paste0(result_path,"/model_eval_final/lasso_model_results_LOOCV.rds"))
mlp_results <- readRDS(paste0(result_path,"/model_eval_final/mlp_model_results_LOOCV_h2o.rds")) 
results <- bind_rows(linear_results, lasso_results, mlp_results)

west_results <- filter(results, region=="Western US") %>%
  mutate(model_outcome_spec=paste(model, outcome, spec, sep="-")) %>%
  filter(model_outcome_spec %in% c("Linear-level-weighted FE", "LASSO-log-non-weighted FE", "Neural Net-log-simple"))

west_results$experiment <- factor(west_results$experiment, levels=c("grid","eco3","eco2","regional"))
west_results$model <- factor(west_results$model, levels=c("Linear","LASSO","Neural Net"))

ggplot(west_results %>% filter(metric == "10 years"),
       aes(x=truth/1e9, y=response_positive/1e9)) +
  geom_point(size=4, alpha=0.6) +
  geom_abline(,0,1,linetype="dashed") +
  theme_classic() + theme(text = element_text(size=22)) + 
  labs(y="Predictions (emissions, MT)", x="Observations (emissions, MT)") +
  facet_wrap(model~experiment, scales = "free") 
ggsave("westUS_linear_level_spatial_res.pdf", width=10.5, height=8)#+

west_results <- filter(results, region=="Western US") %>%
  mutate(model_outcome_spec_spat=paste(model, outcome, spec, experiment, sep="-")) %>%
  filter(model_outcome_spec_spat %in% c("Linear-level-weighted FE-regional", "LASSO-log-non-weighted FE-regional", "Neural Net-log-simple-eco3"))

west_results$experiment <- factor(west_results$experiment, levels=c("grid","eco3","eco2","regional"))
west_results$model <- factor(west_results$model, levels=c("Linear","LASSO","Neural Net"))

ggplot(west_results,
       aes(x=truth/1e9, y=response_positive/1e9)) +
  geom_point(size=4, alpha=0.6) +
  geom_abline(,0,1,linetype="dashed") +
  theme_classic() + theme(text = element_text(size=22)) + 
  labs(y="Predictions (emissions, MT)", x="Observations (emissions, MT)") +
  facet_wrap(~model+metric, scales = "free") 
ggsave("westUS_linear_level_temporal.pdf", width=10.5, height=8)#+


#---------------------------------------------------------------------------------------------------
# Figure Sx: Projected changes in key climate variables
#---------------------------------------------------------------------------------------------------
cmip_region <- readRDS(paste0(result_path, "/CMIP_region_climate_vars.rds"))

cmip_10yr_smooth <- bind_rows(cmip_region,
                              filter(cmip_region, year<2015) %>% mutate(scenario = "ssp126"),
                              filter(cmip_region, year<2015) %>% mutate(scenario = "ssp245"),
                              filter(cmip_region, year<2015) %>% mutate(scenario = "ssp370"),
                              filter(cmip_region, year<2015) %>% mutate(scenario = "ssp585")) %>%
  group_by(scenario, region, model, var) %>%
  arrange(year) %>%
  mutate(anomaly_10yr=rollapplyr(cmip_anomaly, width = 10, FUN=mean, partial=T)) %>%
  group_by(scenario, year, region, var) %>%
  summarise(anomaly_median=median(anomaly_10yr, na.rm=T), 
            anomaly_10=quantile(anomaly_10yr, 0.1, na.rm=T),
            anomaly_90=quantile(anomaly_10yr, 0.9, na.rm=T),
            anomaly_25=quantile(anomaly_10yr, 0.25, na.rm=T),
            anomaly_75=quantile(anomaly_10yr, 0.75, na.rm=T)) %>% ungroup() %>%
  filter(scenario=="historical" | year>2014)

## just to connect 2014 and 2015 (for visualization purpose)
cmip_10yr_smooth <- bind_rows(cmip_10yr_smooth,
                              filter(cmip_10yr_smooth, year==2014) %>% mutate(scenario = "ssp126"),
                              filter(cmip_10yr_smooth, year==2014) %>% mutate(scenario = "ssp245"),
                              filter(cmip_10yr_smooth, year==2014) %>% mutate(scenario = "ssp370"),
                              filter(cmip_10yr_smooth, year==2014) %>% mutate(scenario = "ssp585"))

cmip_10yr_smooth$region <- factor(cmip_10yr_smooth$region,
                                  levels=c("Western US", "Northeastern US", "Southeastern US", "Canada-Alaska", "Mexico"))

cmip_10yr_smooth <- bind_rows(cmip_10yr_smooth)
var_long_list <- c("temperature_surface", "precipitation", "RH_surface", "vpd","soil_moisture", "runoff", "wind_speed")
var_labels <- c("temperature (degree C)", 
                "precipitation (kg/m2)", 
                "RH (%)", 
                "vpd (kPa)",
                "soil moisture (kg/m2)", 
                "runoff (kg/m2)", 
                "wind speed (m/s)")
var_short_list <- c("tas", "pr", "hurs", "vpd", "mrsos", "mrro", "sfcWind")
cmip_10yr_smooth$var <- factor(cmip_10yr_smooth$var, levels = c("temperature_surface", "precipitation", "RH_surface", "wind_speed", "vpd","runoff", "soil_moisture"))

ggplot(filter(cmip_10yr_smooth, year>2000, year <2061, scenario!="ssp585"), 
       aes(x=year, y= anomaly_median, colour=scenario, fill=scenario)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=anomaly_25, ymax=anomaly_75), alpha=0.3, colour=NA) +
  scale_fill_manual(values = color_scenario) +
  scale_colour_manual(values = color_scenario) +
  theme_classic() + theme(text = element_text(size=25),
                          panel.spacing = unit(2, "lines")) +
  labs(x="", fill="", colour="", y= "") +
  facet_grid(var~region, scales = "free")
ggsave("allvar_CMIP6_debias_10yr_smooth.pdf", width=19, height=13)

ggplot(filter(cmip_10yr_smooth, year <2061, scenario!="ssp585", !var%in%c("temperature_surface", "precipitation", "RH_surface","wind_speed")), 
       aes(x=year, y= anomaly_median, colour=scenario, fill=scenario)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=anomaly_25, ymax=anomaly_75), alpha=0.3, colour=NA) +
  scale_fill_manual(values = color_scenario) +
  scale_colour_manual(values = color_scenario) +
  theme_classic() + theme(text = element_text(size=22),
                          panel.spacing.x = unit(2, "lines")) +
  labs(x="", fill="", colour="", y= var_labels[vv]) +
  facet_grid(var~region, scales = "free")
ggsave("droughtvar_CMIP6_debias_10yr_smooth.pdf", width=19, height=6.2)

for (vv in 1:length(var_long_list)){
ggplot(filter(cmip_10yr_smooth, year <2061, scenario!="ssp585", var==var_long_list[vv]), 
       aes(x=year, y= anomaly_median, colour=scenario, fill=scenario)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=anomaly_25, ymax=anomaly_75), alpha=0.3, colour=NA) +
  scale_fill_manual(values = color_scenario) +
  scale_colour_manual(values = color_scenario) +
  theme_classic() + theme(text = element_text(size=22)) +
  labs(x="", fill="", colour="", y= var_labels[vv]) +
  facet_wrap(~region, scales = "free", nrow=1)
ggsave(paste0(var_short_list[vv],"_CMIP6_debias_10yr_smooth.pdf"), width=11, height=3)
}


#------------------------------------------------------------
# Figure Sx: Analyze the emissions by veg types
#-----------------------------------------------------------
gfed_veg <- readRDS(paste0(data_path, "/GFED4S_monthly_2003_2021_NorthAmerica_025deg_VEG.rds")) %>%
  filter(year %in% 2001:2021) %>%
  group_by(gfed_cell_id) %>%
  summarise_at(c("DM_AGRI","DM_BORF", "DM_DEFO","DM_PEAT","DM_SAVA","DM_TEMF"),
                 .funs = sum, na.rm=T) %>% ungroup() %>%
  pivot_longer(cols = c("DM_AGRI","DM_BORF", "DM_DEFO","DM_PEAT","DM_SAVA","DM_TEMF")) %>%
  mutate(name=replace(name, name=="DM_AGRI","agriculture"),
         name=replace(name, name=="DM_DEFO","landuse change"),
         name=replace(name, name%in%c("DM_BORF","DM_TEMF"),"forest"),
         name=replace(name, name%in%c("DM_SAVA"),"savanna"),
         name=replace(name, name%in%c("DM_PEAT"),"peatland")) 

gfed_veg_ratio <- gfed_veg %>% group_by(gfed_cell_id) %>%
  mutate(ratio=value/sum(value)) %>% ungroup()

gfed_wang_region <- readRDS(paste0(data_path, "/crosswalk_GFED_wang_regions.rds")) 
gfed_veg_summ <- left_join(gfed_veg, gfed_wang_region) %>%
  group_by(region, name) %>%
  mutate(value=value/1e9/21) %>%
  summarise(value=sum(value)) %>% ungroup() %>%
  group_by(region) %>%
  mutate(percent=value/sum(value)) %>% ungroup() 

### analyze the projected emissions
proj_emis_grid <- readRDS("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/results/fire_model/CMIP6_fire_proj/downscaled_emis_final_threshold_05.rds") %>%
  filter(year=="2055", scenario=="ssp370") %>%
  group_by(gfed_cell_id) %>%
  summarise(pred_DM_grid=median(pred_DM_grid)) %>% ### take mean across GCMs
  ungroup()

proj_emis_veg <- right_join(proj_emis_grid, gfed_veg_ratio) %>%
  filter(!is.na(pred_DM_grid)) %>%
  mutate(pred_DM_ssp370=pred_DM_grid*ratio) %>%
  left_join(gfed_wang_region) %>%
  group_by(region, name) %>%
  summarise(pred_DM_ssp370=sum(pred_DM_ssp370/1e9)) %>% ungroup() 

gfed_veg_summ <- left_join(gfed_veg_summ, proj_emis_veg)

gfed_veg_summ$name <- factor(gfed_veg_summ$name, levels=c("forest", "savanna", "agriculture",
                                                          "landuse change", "peatland"))

gfed_veg_summ$region <- factor(gfed_veg_summ$region, levels=c("Western US", "Southeastern US","Northeastern US",
                                                              "Canada-Alaska", "Mexico"))

gfed_veg_summ <- gfed_veg_summ %>% arrange(region, name)
write.xlsx(gfed_veg_summ, "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/results/summ_by_veg.xlsx", replace=T)

gfed_veg_summ %>% select(region, name, percent) %>% as.data.frame() %>%
  filter(!is.na(region))





gfed_grid <- readRDS(paste0(data_path, "/GFED4S_NorthAmerica_025deg_grid_EPSG4326.rds"))
gfed_wang_region <- readRDS(paste0(data_path, "/crosswalk_GFED_wang_regions.rds")) 
gfed_grid <- left_join(gfed_grid, gfed_wang_region) %>% filter(!is.na(region)) %>%
  mutate(region=replace(region, region%in%c("Northeastern US","Southeastern US"), "Eastern US")) 

proj_emis_2050 <-  right_join(gfed_grid, proj_emis_grid) %>% filter(year==2055) %>%
  mutate(DM_kton = pred_DM_grid / 1e6) 


###### plot the ratio of emission changes
tmp_2050 <- filter(proj_emis_2050) %>% filter(scenario=="ssp370") %>% rename(emis_ssp370=DM_kton) %>% 
  as.data.frame() %>%  select(gfed_cell_id, region, emis_ssp370) 
tmp_2001_2021 <- filter(gfed_2001_2021) %>% mutate(scenario="2001-2021") %>% rename(emis_2001_2021=DM_kton)  %>% 
  as.data.frame() %>%  select(gfed_cell_id, region, emis_2001_2021) 
tmp_emis <- left_join(tmp_2050, tmp_2001_2021) %>%
  mutate(ratio = emis_ssp370 / emis_2001_2021)
tmp_emis$ratio_cut <- cut(tmp_emis$ratio, c(0,0.8,1,1.2,1.5,3,5,100), include.lowest = T)
tmp_emis <- right_join(gfed_grid, tmp_emis) 

ggplot(tmp_emis, aes(fill=ratio_cut, colour=ratio_cut)) +
  geom_sf() +
  scale_color_manual(values = col_disc) +
  scale_fill_manual(values = col_disc)  +
  theme_classic() + theme(text = element_text(size=14),
                          axis.ticks = element_blank(),
                          axis.line = element_blank(),
                          axis.text = element_blank())
ggsave("ratios_ssp370_2001_2021.png", width=6, height=4)
  

