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
library(ggridges)
library(purrr)

rm(list=ls())
gc()

code_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/rcode"
data_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data"
figure_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/figures"
result_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/results/smoke_model"
health_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/results/health"

db_path <- "/Users/mhqiu/BurkeLab Dropbox/Projects/smoke-climate"
db_proj_path <- "/Users/mhqiu/BurkeLab Dropbox/Projects"

setwd(figure_path)
states_map <- map_data("state")

col_scenario <- c("grey65",rev(met.brewer(name="Homer1", n=8, type="discrete"))[c(2,6,8)])
#-------------------------------------------------------------------------------
# Plot the uncertainty of the causal chain
# Written by Minghao
# Last edited May 2024
#-------------------------------------------------------------------------------
#------ load the projected smoke from MC simulations (combined uncertainty of GCM and fire models)
smoke_mc <- list.files(paste0(result_path, "/smoke_proj"), pattern = "proj_smoke_ssp370_MC*", full.names = T) %>%
    purrr::map_df(readRDS) %>% bind_rows()

#------ first de-bias the smoke using historical smoke simulations ------#
pred_smoke_gfed <- readRDS(paste0(result_path, "/smoke_proj/gfedDM_final_smoke_pred_grid_wind_9region_90cone_gridmet.rds")) %>%
  filter(year %in% 2011:2020) %>%
  filter(model %in% c("wind_all", "obs")) %>%
  group_by(grid_id_10km, model) %>%
  summarise(pred_smoke=mean(pred_smoke, na.rm=T)) %>%
  ungroup() %>%
  pivot_wider(values_from=pred_smoke, names_from = model, names_prefix = "smoke_") 

total_smoke_pm <- readRDS(paste0(data_path, "/totalPM/totalPM_smoke_combined_2006_2020.rds")) 

grid_list <- intersect(unique(pred_smoke_gfed$grid_id_10km), 
                       unique(total_smoke_pm$grid_id_10km))

smoke_mc_annual <- smoke_mc %>%
  filter(grid_id_10km %in% grid_list) %>%
  filter(!is.na(pred_smoke)) %>%
  group_by(grid_id_10km, year, scenario, mc_id) %>%
  summarise(pred_smoke=mean(pred_smoke)) %>% ungroup()   ## get the annual mean

##de-bias based on annual mean in 2011-2020
smoke_2011_2020 <- total_smoke_pm %>% 
  filter(year %in% 2011:2020) %>%
  group_by(grid_id_10km) %>%
  summarise_at(c("smokePM_mean"),.funs = mean, na.rm=T) %>%
  rename(smoke_2011_2020=smokePM_mean) %>% ungroup()

smoke_debias <- left_join(smoke_mc_annual, pred_smoke_gfed) %>%
  left_join(smoke_2011_2020) %>%
  mutate(pred_smoke_debias=pred_smoke - smoke_wind_all + smoke_2011_2020) %>%
  rename(smoke_grid = pred_smoke_debias) %>% select(grid_id_10km, year, scenario, mc_id, smoke_grid)

smoke_debias[smoke_debias$smoke_grid<0, "smoke_grid"] <- 0
saveRDS(smoke_debias, paste0(result_path, "/smoke_proj/smoke_debias_MC_1_100.rds"))

#---- convert to county-level smoke -------#
crosswalk_county_10km <- readRDS(paste0(data_path, "/health/county_grid_area_crosswalk.rds")) %>%
  mutate(area=as.numeric(area))
pop_10km_grid <- readRDS(paste0(data_path, "/health/us_population_grid_2020_2022_2050.rds")) %>%
  as.data.frame() %>% select(grid_id_10km, pop_2020)
crosswalk_county_10km <- left_join(crosswalk_county_10km, pop_10km_grid) %>%
  group_by(grid_id_10km) %>%
  mutate(area_weight=area/sum(area)) %>% ungroup()

county_pop <- readRDS(paste0(data_path, "/health/US_county_pop_2019_2050_age.rds"))  %>%
  filter(age_group!="under_5")

smoke_county <- left_join(smoke_debias, crosswalk_county_10km) %>%
  filter(!is.na(GEOID)) %>%
  group_by(year, scenario, mc_id, GEOID) %>%
  summarise(smokePM = weighted.mean(smoke_grid, w=area_weight *pop_2020, na.rm=T)) %>%
  ungroup() %>%
  mutate(GEOID = as.numeric(GEOID)) %>%  rename(fips=GEOID) %>%
  left_join(county_pop) %>% filter(!is.na(pop_2050))
saveRDS(smoke_county, paste0(data_path, "/health/US_smoke_county_2050_MC_1_100.rds"))

smoke_pop_weight <- left_join(smoke_debias, crosswalk_county_10km) %>%
  group_by(mc_id) %>%
  summarise(smokePM = weighted.mean(smoke_grid, w=area_weight *pop_2020, na.rm=T))%>%
  ungroup()

###### randomly match with the health dose-response function 
# set.seed(1000)
# bootid_list <-  sample(1:500, size=100)
#coef <- readRDS(paste0(health_path, "/coef_poisson_bins.rds")) %>%
  # filter(age_group!="under_5",  model=="poisson bins", 
  #        bootid %in% bootid_list) %>%
  # rename(mc_id=bootid)

coef <- readRDS(paste0(health_path, "/coef_poisson_bins.rds")) %>%
filter(age_group!="under_5",  model=="poisson bins") %>%
rename(mc_id=bootid)

cutoff <- c(-0.01, 0.1, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6, 10, 20)

smoke_death_combined <- NULL
for (iii in 1:5){
coef_tmp <- filter(coef, mc_id %in%seq(1+100*(iii-1), 100*iii)) %>%
  mutate(mc_id=mc_id-100*(iii-1))
  
smoke_death <- smoke_county %>%
  mutate(smokePM=replace(smokePM, smokePM>10, 10)) %>%
  mutate(bins=cut(smokePM, cutoff)) %>%
  left_join(coef_tmp) %>%
  mutate(coef=replace(coef, is.na(coef), 0)) %>%
  mutate(death_smoke = (exp(coef)-1)*death_rate_avg*pop_2019) %>%
  group_by(age_group, mc_id) %>%
  summarise(death_smoke=sum(death_smoke)) %>% ungroup()

smoke_death_combined <- bind_rows(smoke_death_combined, smoke_death)
}

saveRDS(smoke_death_combined, paste0(health_path, "/ssp370_smoke_death_MC_1_500.rds"))

summ_death <- smoke_death_combined  %>%
  group_by(age_group) %>% 
  summarise(death_median=median(death_smoke),  ### get uncertainty across boot ids
            death_mean=mean(death_smoke),  ### get uncertainty across boot ids
            death_p025=quantile(death_smoke, 0.025),
            death_p975=quantile(death_smoke, 0.975),
            death_p25=quantile(death_smoke, 0.25),
            death_p75=quantile(death_smoke, 0.75)) %>% ungroup() 


