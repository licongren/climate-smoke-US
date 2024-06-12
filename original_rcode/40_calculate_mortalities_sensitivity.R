library(foreach)
library(doParallel)
library(doSNOW)
library(dplyr)
library(tidyr)
library(ncdf4)
library(terra)
library(ncdf4.helpers)
library(sf)
library(exactextractr)
library(ggplot2)
library(lubridate)
library(stringr)
library(MetBrewer)
library(openxlsx)
library(tigris)

rm(list=ls())
gc()

code_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/rcode"
data_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data"
figure_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/figures"
result_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/results/smoke_model"
health_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/results/health"

setwd(figure_path)
states_map <- map_data("state")

db_path <- "/Users/mhqiu/BurkeLab Dropbox/Projects/smoke-climate"
db_proj_path <- "/Users/mhqiu/BurkeLab Dropbox/Projects"

color_map <- rev(met.brewer(name="Homer1", n=8, type="discrete"))[c(2,6,8)]
#-------------------------------------------------------------------------------
# Calculate moralities associated smoke exposure for the sensitivity scenarios
# Written by Minghao
# Last edited Dec 2023
#-------------------------------------------------------------------------------
#### convert smoke data to county 
grid_10km <- st_read(paste0(db_proj_path, "/smoke_PM_prediction/data/1_grids/grid_10km_wgs84/grid_10km_wgs84.shp")) %>%
  rename(grid_id_10km=ID)

crosswalk_county_10km <- readRDS(paste0(data_path, "/health/county_grid_area_crosswalk.rds")) %>%
  mutate(area=as.numeric(area))
pop_10km_grid <- readRDS(paste0(data_path, "/health/us_population_grid_2020_2022_2050.rds")) %>%
  as.data.frame() %>% select(grid_id_10km, pop_2020)
crosswalk_county_10km <- left_join(crosswalk_county_10km, pop_10km_grid) %>%
  group_by(grid_id_10km) %>%
  mutate(area_weight=area/sum(area)) %>% ungroup()

county_pop <- readRDS(paste0(data_path, "/health/US_county_pop_2019_2050_age.rds"))  %>% filter(age_group=="all_ages")

####---------------------------------------------------
### convert smoke to county-level
####---------------------------------------------------
file_list <- list.files(path=paste0(result_path, "/smoke_proj/sensitivity"),
                        pattern = "*_smoke_bootsame.rds", full.names = T)

smoke_county_combined <- NULL
for (ii in file_list){
smoke_future <-  readRDS(ii) %>%
  rename(smokePM = smoke_grid)  %>%
  select(year, scenario, bootid, grid_id_10km, smokePM)

smoke_county <- left_join(smoke_future, crosswalk_county_10km) %>%
  filter(!is.na(GEOID)) %>%
  group_by(year, scenario, GEOID, bootid) %>%
  summarise(smokePM = weighted.mean(smokePM, w=area_weight *pop_2020, na.rm=T)) %>%
  ungroup() %>%
  mutate(GEOID = as.numeric(GEOID)) %>%  rename(fips=GEOID) %>%
  left_join(county_pop) %>% filter(!is.na(pop_2050))

smoke_county_combined <- bind_rows(smoke_county_combined, smoke_county)
}

saveRDS(smoke_county_combined, paste0(result_path, "/smoke_proj/sensitivity/US_smoke_county_smoke_bootsame.rds"))

#### for smoke boot whcih the dataset os too big
# file_list <- list.files(path=paste0(result_path, "/smoke_proj/sensitivity"),
#                         pattern = "*_smoke_boot_gcm_mean.rds", full.names = T)
# 
# smoke_future <-  readRDS(file_list[1]) %>%
#   rename(smokePM = smoke_grid)  %>%
#   select(year, scenario, bootid, grid_id_10km, smokePM)
# 
# gcm_list <- unique(smoke_future$gcm)
# 
# smoke_county_combined <- NULL
# for (ii in gcm_list){
#   print(ii)
#   smoke_future_tmp <- filter(smoke_future, gcm==ii)
#   
#   smoke_county <- left_join(smoke_future_tmp, crosswalk_county_10km) %>%
#     filter(!is.na(GEOID)) %>%
#     group_by(year, scenario, GEOID, bootid) %>%
#     summarise(smokePM = weighted.mean(smokePM, w=area_weight *pop_2020, na.rm=T),
#               smokePM = weighted.mean(smokePM, w=area_weight *pop_2020, na.rm=T)) %>%
#     ungroup() %>%
#     mutate(GEOID = as.numeric(GEOID)) %>%  rename(fips=GEOID) %>%
#     left_join(county_pop) %>% filter(!is.na(pop_2050))
#   
#   smoke_county_combined <- bind_rows(smoke_county_combined, smoke_county)
# }
# 
# saveRDS(smoke_county_combined, paste0(result_path, "/smoke_proj/sensitivity/US_smoke_county_smoke_boot.rds"))


#-----------------------------------------------------------
### Using poisson bin regression to calculate the mortality
#----------------------------------------------------------
coef <- readRDS(paste0(health_path, "/coef_poisson_bins.rds")) %>% filter(age_group=="all_ages",  model=="poisson bins") %>%
  group_by(bins, age_group, model) %>%
  summarise(coef=mean(exp(coef)-1)) %>% ungroup()

cutoff <- c(-0.01, 0.1, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6, 10, 20)

file_list <- list.files(path=paste0(result_path, "/smoke_proj/sensitivity"),
                        pattern = "US_smoke_county_smoke_bootsame*", full.names = T)

smoke_death_tmp  <- readRDS(file_list[1]) %>% 
    mutate(smokePM=replace(smokePM, smokePM>10, 10)) %>%
    mutate(bins=cut(smokePM, cutoff)) %>%
    left_join(coef) %>%
    mutate(coef=replace(coef, is.na(coef), 0)) %>%
    mutate(death_smoke =  (exp(coef)-1)*death_rate_avg*pop_2050)
saveRDS(smoke_death_tmp, paste0(health_path, "/US_smoke_death_smoke_bootsame.rds"))
  
