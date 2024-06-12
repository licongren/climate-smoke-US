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

color_map <- c("grey65",
               rev(met.brewer(name="Homer1", n=8, type="discrete"))[c(2,6,7,8)])
#-------------------------------------------------------------------------------
# Analyze the fire CO2 emissions from CMIP6 models
# Written by Minghao
# Last edited July 2022
#-------------------------------------------------------------------------------
# gfed_ecoregion_link <- readRDS(paste0(db_path, "/crosswalk_GFED_NA_CEC_Eco_Level2_3.rds")) %>%
#   mutate(intersection_area=as.numeric(intersection_area)) %>%   group_by(gfed_cell_id)  %>%
#   arrange(desc(intersection_area)) %>%
#   mutate(area_weight=intersection_area/sum(intersection_area)) %>% ungroup() # %>% distinct(gfed_cell_id, .keep_all=T)
# 
# gfed_wang_region <- readRDS(paste0(data_path, "/crosswalk_GFED_wang_regions.rds")) 
# 
# ### create crosswalk between ecoregions and 1deg cells of CMIP6 
# crosswalk_gfed_1deg <- readRDS(paste0(data_path, "/crosswalk_gfed_1deg.rds"))
# crosswalk_ecoregion <- left_join(crosswalk_gfed_1deg, gfed_ecoregion_link) %>%
#   left_join(gfed_wang_region)
# 
# cmip6_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/fireCO2/North_America_1deg"
# file_list <- list.files(path=cmip6_path, pattern = "*.nc", full.names = T)
# nc_open(file_list[1])

#---------------------
## Aggregate CMIP6 fire CO2 to regional level
#---------------------
# var_short <- "fFire"
# var_long <- "fireCO2"
# file_list <- list.files(path=paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/",
#                                     var_long),
#                         pattern = "*.rds", full.names = T)
# 
# cmip_files <- data.frame(filename=file_list) %>%
#   mutate(model=str_match(filename, paste0(".*", var_short, "\\_(.*?)\\_.*"))[,2]) %>%
#   mutate(scenario=str_match(filename, paste0(".*",model,"\\_(.*?)\\_.*"))[,2])
# 
# model_list <- unique(cmip_files$model)
# 
# month_seq <-  data.frame(date=seq(as.Date("1997-01-01"), as.Date("2100-12-01"), by="months")) %>%
#   mutate(year=as.numeric(substr(date,1,4)), month=as.numeric(substr(date,6,7)),
#          ndays=lubridate::days_in_month(date)) %>% select(-date)
# 
# fire_region <- NULL
# for (ii in 1:nrow(cmip_files)){
#     print(ii)
#     cmip_data <-  readRDS(cmip_files[ii,"filename"]) %>% 
#       filter(year>1996, year< 2101)
#     
#     cmip_data <- left_join(cmip_data, month_seq, by=c("year", "month")) %>%
#       mutate(fFire = fFire * ndays *24 * 3600) %>%
#       left_join(crosswalk_ecoregion) %>%
#       filter(!is.na(region)) %>%
#       group_by(year, region) %>%
#       summarise(fire_co2_kg = sum(fFire * intersection_area, na.rm=T))%>%
#       mutate(scenario = cmip_files[ii, "scenario"], model = cmip_files[ii, "model"])
#     
#     fire_region <- bind_rows(fire_region, cmip_data)
#   }
#   
# saveRDS(fire_region, paste0(result_path,"/CMIP6_fire_proj/cmip_fireCO2_region.rds"))
# 

fire_cmip6 <- readRDS(paste0(result_path,"/CMIP6_fire_proj/cmip_fireCO2_region.rds")) %>%
  rename(gcm = model)

#%>%
#   filter(model %in% c("CESM2-WACCM","GFDL-ESM4","CNRM-ESM2-1"))
##################
gfed_region <- readRDS(paste0(result_path, "/regional_annual_gfed_DM.rds")) 

# gfed_2014 <- filter(gfed_region, year<2015) %>% group_by(region) %>%
#   summarise(DM_2001_2014 = mean(DM_kg, na.rm=T)) %>% ungroup()

gfed_2021 <- gfed_region %>% group_by(region) %>%
  summarise(DM_2001_2021 = mean(DM_kg, na.rm=T)) %>% ungroup()

gfed_region <- left_join(gfed_region, gfed_2021) %>%
  mutate(DM_anomaly = DM_kg - DM_2001_2021) %>%
  group_by(region) %>%
  arrange(year) %>%
  mutate(DM_10yr=rollapplyr(DM_kg, width = 10, FUN=mean, partial=T),
         DM_anomaly_10yr=rollapplyr(DM_anomaly, width = 10, FUN=mean, partial=T)) %>% ungroup()

fire_region_mean <- filter(fire_cmip6, scenario%in%c("historical", "ssp370"), year %in% seq(2001,2021)) %>% 
  group_by(gcm, region) %>%
  summarise(fireco2_2001_2021 = mean(fire_co2_kg, na.rm=T)) %>% ungroup() %>%
  # group_by(region) %>%
  # summarise(fireco2_2001_2021 = median(fireco2_2001_2021, na.rm=T)) %>% ungroup() %>%
  left_join(gfed_2021) %>%
  mutate(dm_co2=DM_2001_2021/fireco2_2001_2021) %>%
  select(region, gcm, dm_co2)
#---------------------
## calculate a relationship between DM and fire CO2 using historical data 
#---------------------
fire_region <- left_join(fire_cmip6, fire_region_mean) %>%
  mutate(pred_DM = fire_co2_kg * dm_co2)

fire_2001_2021 <- filter(fire_region, scenario%in%c("historical", "ssp370"), year %in% seq(2001,2021)) 
fire_2001_2021_mean <- fire_2001_2021 %>%
  group_by(gcm, region) %>%
  summarise(DM_2001_2021 = mean(pred_DM, na.rm=T)) %>% ungroup() %>%
  select(gcm, region, DM_2001_2021)

#### calculate the delta differences within gcm projections and then add onto real obs
fire_region <-  bind_rows(filter(fire_region, year>2021),
                          fire_2001_2021 %>% mutate(scenario = "historical"),
                          fire_2001_2021 %>% mutate(scenario = "ssp126"),
                          fire_2001_2021 %>% mutate(scenario = "ssp245"),
                          fire_2001_2021 %>% mutate(scenario = "ssp370"),
                          fire_2001_2021 %>% mutate(scenario = "ssp585")) %>%
  left_join(fire_2001_2021_mean) %>%
  mutate(DM_anomaly = pred_DM - DM_2001_2021) %>% select(-DM_2001_2021)

fire_region <- left_join(fire_region, gfed_2021) %>%
  mutate(DM_level = DM_anomaly + DM_2001_2021)

variable <-  "DM_level" ##"pred_fireco2" "fireco2_anomaly"
cmip_combine <- fire_region %>%
  group_by(scenario, region, gcm) %>%
  arrange(year) %>%
  mutate(pred_10yr=rollapplyr(!!as.name(variable), width = 10, FUN=mean, partial=T)) %>%  ## just for calculating the 10-yr mean
  filter(scenario=="historical" | year>2021) %>% ungroup() ### to drop the fake 2001-2021 for the SSP scenarios

## just to connect 2021 and 2022 (for visualization purpose)
cmip_combine <- bind_rows(cmip_combine,
                          filter(cmip_combine, year==2021) %>% mutate(scenario = "ssp126"),
                          filter(cmip_combine, year==2021) %>% mutate(scenario = "ssp245"),
                          filter(cmip_combine, year==2021) %>% mutate(scenario = "ssp370"),
                          filter(cmip_combine, year==2021) %>% mutate(scenario = "ssp585"))

ggplot(filter(cmip_combine,
              region=="Western US", gcm=="CESM2-WACCM"),
       aes(x=year, y= pred_10yr*1e-9,
           colour=scenario, fill=scenario)) +
  geom_line(size=1.4) +
  geom_point(data = filter(gfed_region, region=="Western US"),
             aes(x=year, y=DM_10yr*1e-9), colour="black", fill=NA, size=1.5, shape=1) +
  geom_hline(aes(yintercept = 1e-9*filter(gfed_region, region=="Western US", year==2021)$DM_10yr), linetype="dashed") +
  #geom_hline(aes(yintercept = 1e-9*filter(gfed_region, region=="Western US", year==2020)$DM_kg)) +
  scale_fill_manual(values = color_map) +
  scale_colour_manual(values = color_map) +
  theme_classic() + theme(text = element_text(size=16),
                          axis.text.x = element_text(size=12),
                          panel.spacing = unit(1.5, "lines")) +
  labs(x="", fill="", colour="", y= "Predicted DM (MT)") +
  facet_wrap(~gcm, ncol=3) +
  scale_y_continuous(limits=c(-10, 160)) +
  scale_x_continuous(breaks=seq(2000, 2100, by=25))
ggsave("west_fireCO2_CMIP6_10yr_cesm2.png", width=5.5, height=3.5)



