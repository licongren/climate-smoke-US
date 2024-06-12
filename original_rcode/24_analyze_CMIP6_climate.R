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

rm(list=ls())
gc()

code_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/rcode"
data_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data"
figure_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/figures"
result_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/results/smoke_model"

db_path <- "/Users/mhqiu/BurkeLab Dropbox/Projects/smoke-climate"

setwd(figure_path)
states_map <- map_data("state")

#------------- Load GFED data 
gfed_wang_region <- readRDS(paste0(data_path, "/crosswalk_GFED_wang_regions.rds"))

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

final_proj <- readRDS(paste0(result_path, "/smoke_proj/smoke_2046_2055_debias_gcm.rds"))
model_list <- unique(final_proj$gcm)

color_map <- c("grey65",
               rev(met.brewer(name="Homer1", n=8, type="discrete"))[c(2,5,6,8)])

cmip_region_combined <- NULL 
for (vv in 1: length(var_long_list)){
  var_long <- var_long_list[vv]  
  var_short <- var_short_list[vv]  
  print(var_short)
  
# narr <- readRDS(paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/",
#                          var_short,"_narr_1deg_1997_2021.rds"))
  if (var_long == "soil_moisture"){
  cmip_debias <- readRDS(paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/",
                                var_long, "/CMIP6_annual_debias_1deg_nldas.rds")) %>%
    mutate(ratio = cmip_value / nldas_value)

  print(quantile(cmip_debias$ratio, seq(0,1,0.1)))
  #grid_list <- intersect(unique(narr$grid_1deg_id), unique(cmip_debias$grid_1deg_id))
  print(length(unique(cmip_debias$model)))
  
  cmip_region <- filter(cmip_debias) %>%
    left_join(crosswalk_1deg_region) %>% filter(!is.na(region)) %>%
    group_by(region, year, model, scenario) %>%
    summarise(cmip_debias = weighted.mean(cmip_debias, w= coverage_area, na.rm=T),
              nldas_value = weighted.mean(nldas_value, w= coverage_area, na.rm=T)) %>%
    ungroup() %>%
    mutate(cmip_anomaly = cmip_debias - nldas_value) %>%
    filter(!region %in% c("Canada-Alaska","Mexico"))
 } else{
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
 }
  
  cmip_region <- cmip_region %>% filter(model %in% model_list)
  cmip_region$var <- var_long
  cmip_region_combined <- bind_rows(cmip_region_combined, cmip_region) 
}
saveRDS(cmip_region_combined, paste0(result_path, "/CMIP_region_climate_vars.rds"))



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

cmip_model_summ <- cmip_region %>% group_by(scenario, year, region) %>%
  summarise(anomaly_median=median(cmip_anomaly), 
            anomaly_10=quantile(cmip_anomaly, 0.1),
            anomaly_90=quantile(cmip_anomaly, 0.9),
            anomaly_25=quantile(cmip_anomaly, 0.25),
            anomaly_75=quantile(cmip_anomaly, 0.75)) %>% ungroup()

## just to connect 2014 and 2015 (for visualization purpose)
cmip_model_summ <- bind_rows(cmip_model_summ,
                           filter(cmip_model_summ, year==2014) %>% mutate(scenario = "ssp126"),
                           filter(cmip_model_summ, year==2014) %>% mutate(scenario = "ssp245"),
                           filter(cmip_model_summ, year==2014) %>% mutate(scenario = "ssp370"),
                           filter(cmip_model_summ, year==2014) %>% mutate(scenario = "ssp585"))
cmip_model_summ$region <- factor(cmip_model_summ$region,
                                 levels=c("Canada-Alaska", "Northeastern US", "Southeastern US",
                                          "Western US", "Mexico"))
ggplot(filter(cmip_model_summ, year<2056, scenario!="ssp585"), 
       aes(x=year, y= anomaly_median, colour=scenario, fill=scenario)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=anomaly_25, ymax=anomaly_75), alpha=0.3, colour=NA) +
  scale_fill_manual(values = color_map) +
  scale_colour_manual(values = color_map) +
  theme_classic() + theme(text = element_text(size=16)) +
  labs(x="", fill="", colour="", y= var_labels[vv]) +
  facet_wrap(~region, scales = "free")
ggsave(paste0(var_short,"_CMIP6_debias.png"), width=9, height=5)

cmip_10yr_smooth <- bind_rows(cmip_region,
                              filter(cmip_region, year<2015) %>% mutate(scenario = "ssp126"),
                              filter(cmip_region, year<2015) %>% mutate(scenario = "ssp245"),
                              filter(cmip_region, year<2015) %>% mutate(scenario = "ssp370"),
                              filter(cmip_region, year<2015) %>% mutate(scenario = "ssp585")) %>%
  group_by(scenario, region, model) %>%
  arrange(year) %>%
  mutate(anomaly_10yr=rollapplyr(cmip_anomaly, width = 10, FUN=mean, partial=T)) %>%
    group_by(scenario, year, region) %>%
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
                                 levels=c("Canada-Alaska", "Northeastern US", "Southeastern US",
                                          "Western US", "Mexico"))

ggplot(filter(cmip_10yr_smooth, year <2056, scenario!="ssp585"), 
       aes(x=year, y= anomaly_median, colour=scenario, fill=scenario)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=anomaly_25, ymax=anomaly_75), alpha=0.3, colour=NA) +
  scale_fill_manual(values = color_map) +
  scale_colour_manual(values = color_map) +
  theme_classic() + theme(text = element_text(size=16)) +
  labs(x="", fill="", colour="", y= var_labels[vv]) +
  facet_wrap(~region, scales = "free")
ggsave(paste0(var_short,"_CMIP6_debias_10yr_smooth.png"), width=9, height=5)

}



