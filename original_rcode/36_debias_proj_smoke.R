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

rm(list=ls())
gc()

code_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/rcode"
data_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data"
figure_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/figures"
result_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/results/smoke_model"

setwd(figure_path)
states_map <- map_data("state")

db_path <- "/Users/mhqiu/BurkeLab Dropbox/Projects/smoke-climate"
db_proj_path <- "/Users/mhqiu/BurkeLab Dropbox/Projects"

color_map <- rev(met.brewer(name="Homer1", n=8, type="discrete"))[c(2,6,8)]
#-------------------------------------------------------------------------------
# Plot the future projected smoke and calculate its contributiosn based on 2016-2020 PM
# Written by Minghao
# Last edited Oct 2022
#-------------------------------------------------------------------------------
# read the smoke grid shape file
# grid_10km <- st_read(paste0(db_proj_path, "/smoke_PM_prediction/data/1_grids/grid_10km_wgs84/grid_10km_wgs84.shp")) %>%
#   rename(grid_id_10km=ID)
# 
# pop_raster <- rast("~/BurkeLab Dropbox/Data/population/gpw_v4_population_count_rev11_2020_2pt5_min.tif")
# # Crop population raster to CONUS
# #population = raster::crop(population, extent(-124.848974, -66.885444, 24.396308, 49.384358))
#   
# crosswalk_pop_10yr <- exact_extract(pop_raster, grid_10km,
#                                     coverage_area = T,
#                                     include_cell = T,
#                                     include_cols="grid_id_10km",
#                                     include_xy = T) %>%
#   bind_rows() %>% 
#   rename(pop = value) %>% select(-cell)
# 
# pop_grid_10km <- crosswalk_pop_10yr %>% group_by(x,y) %>%
#   mutate(weight=coverage_area/sum(coverage_area, na.rm=T)) %>%
#   mutate(pop=pop*weight) %>% ungroup() %>%
#   group_by(grid_id_10km) %>%
#   summarise(pop_2020=sum(pop, na.rm=T)) %>%
#   ungroup()
# 
# pop_grid <- left_join(grid_10km, pop_grid_10km) 
# pop_2020_model <- sum(pop_grid$pop_2020)
# ### take the 2020 and 2022 population from US census to scale
# pop_2020_obs <- read.xlsx(paste0(data_path, "/health/us_population_2020_2022.xlsx"))[5,"X3"] %>% as.numeric()
# pop_2022_obs <- read.xlsx(paste0(data_path, "/health/us_population_2020_2022.xlsx"))[5,"X5"] %>% as.numeric()
# 
# pop_grid <- pop_grid %>% mutate(pop_2020 = pop_2020 / pop_2020_model*pop_2020_obs,
#                                 pop_2022 = pop_2020 / pop_2020_obs * pop_2022_obs)
# ### take the 2050 population and scale the 2020 population
# future_pop <- read.csv(paste0(data_path, "/health/US_census_population_proj.csv")) %>%
#   filter(SEX==0, ORIGIN==0, RACE==0, YEAR%in%c(2022, 2050)) %>%
#   select(YEAR, TOTAL_POP)
# 
# future_pop_2022 <- filter(future_pop, YEAR==2022)$TOTAL_POP %>% unlist()
# future_pop_2050 <- filter(future_pop, YEAR==2050)$TOTAL_POP %>% unlist()
# pop_grid <- pop_grid %>% mutate(pop_2050 = pop_2022 / future_pop_2022 * future_pop_2050)
# 
# saveRDS(pop_grid, paste0(paste0(data_path, "/health/us_population_grid_2020_2022_2050.rds")))
#   
# ggplot(pop_grid, aes(fill=pop_2050, colour=pop_2050)) +
#   geom_sf() +
#   scale_fill_continuous(trans="sqrt") +
#   scale_colour_continuous(trans="sqrt")

### read observed smoke between 2006 and 2020
total_smoke_pm <- readRDS(paste0(data_path, "/totalPM/totalPM_smoke_combined_2006_2020.rds")) 

smoke_2011_2020 <- total_smoke_pm %>% 
  filter(year %in% 2011:2020) %>%
  group_by(grid_id_10km) %>%
  summarise_at(c("smokePM_mean"),.funs = mean, na.rm=T) %>%
  rename(smoke_2011_2020=smokePM_mean) %>% ungroup()

### read the predicted smoke using the same reg model  
pred_smoke_gfed <- readRDS(paste0(result_path, "/smoke_proj/gfedDM_final_smoke_pred_grid_wind_9region_90cone_gridmet.rds")) %>%
  filter(year %in% 2011:2020) %>%
  filter(model %in% c("wind_all", "obs")) %>%
  group_by(grid_id_10km, model) %>%
  summarise(pred_smoke=mean(pred_smoke, na.rm=T)) %>%
  ungroup() %>%
  pivot_wider(values_from=pred_smoke, names_from = model, names_prefix = "smoke_") 

grid_list <- intersect(unique(pred_smoke_gfed$grid_id_10km), 
                       unique(smoke_2011_2020$grid_id_10km))

#### get the non-smoke conc from 2016-2020
smoke_2016_2020 <- total_smoke_pm %>% 
  filter(year %in% 2016:2020) %>%
  group_by(grid_id_10km) %>%
  summarise_at(c("smokePM_mean", "totalPM"),.funs = mean, na.rm=T) %>%
  mutate(nonsmoke_2016_2020=totalPM-smokePM_mean) %>% ungroup()

smoke_2011_2020 <- smoke_2011_2020 %>% filter(grid_id_10km %in% grid_list)
pred_smoke_gfed <- pred_smoke_gfed %>% filter(grid_id_10km %in% grid_list)
smoke_2016_2020 <- smoke_2016_2020 %>% filter(grid_id_10km %in% grid_list) %>%
  select(grid_id_10km, nonsmoke_2016_2020)

#----------------------------------------------------------------
## Preprocess and debias the smoke for 2025, 2035, 2045, and 2055
#----------------------------------------------------------------
file_list <- list.files(path=paste0(result_path, "/smoke_proj/mean_10yrs/"),
                        pattern="proj_smoke_gridmetobs*",
                        full.names = T)

smoke_debias_combined <- NULL
for (ii in file_list){
  print(ii)

smoke_tmp <- readRDS(ii) %>%
  filter(grid_id_10km %in% grid_list) %>%
  filter(!is.na(pred_smoke)) %>%
  group_by(grid_id_10km, year, scenario, gcm) %>%
  summarise(pred_smoke=mean(pred_smoke)) %>% ungroup()   ## get the annual mean

##de-bias based on annual mean in 2011-2020
smoke_tmp_debias <- left_join(smoke_tmp, pred_smoke_gfed) %>%
  left_join(smoke_2011_2020) %>%
  mutate(pred_smoke_debias=pred_smoke - smoke_wind_all + smoke_2011_2020) %>%
  rename(smoke_grid = pred_smoke_debias) %>% select(grid_id_10km, gcm, year, scenario, smoke_grid)
smoke_debias_combined <- bind_rows(smoke_debias_combined, smoke_tmp_debias)
}

smoke_debias_combined[smoke_debias_combined$smoke_grid<0, "smoke_grid"] <- 0

smoke_debias_combined <- left_join(smoke_debias_combined, smoke_2016_2020) %>%
  mutate(totalPM=smoke_grid+nonsmoke_2016_2020)
saveRDS(smoke_debias_combined, paste0(result_path, "/smoke_proj/smoke_10yrs_debias_gcm.rds"))


#----------------------------------------------------------------
## Preprocess and debias the smoke for 2046-2055, individual years
#----------------------------------------------------------------
file_list <- list.files(path=paste0(result_path, "/smoke_proj/2046-2055"),
                        pattern="proj_smoke_gridmetobs*",
                        full.names = T)

smoke_debias_combined <- NULL
for (ii in file_list){
  print(ii)

  smoke_tmp <- readRDS(ii) %>%
    filter(grid_id_10km %in% grid_list) %>%
    filter(!is.na(pred_smoke)) %>%
    group_by(grid_id_10km, year, scenario, gcm) %>%
    summarise(pred_smoke=mean(pred_smoke)) %>% ungroup()   ## get the annual mean
  
  ##de-bias based on annual mean in 2011-2020
  smoke_tmp_debias <- left_join(smoke_tmp, pred_smoke_gfed) %>%
    left_join(smoke_2011_2020) %>%
    mutate(pred_smoke_debias=pred_smoke - smoke_wind_all + smoke_2011_2020) %>%
    rename(smoke_grid = pred_smoke_debias) %>% select(grid_id_10km, gcm, year, scenario, smoke_grid)
  smoke_debias_combined <- bind_rows(smoke_debias_combined, smoke_tmp_debias)
}

smoke_debias_combined[smoke_debias_combined$smoke_grid<0, "smoke_grid"] <- 0

moke_debias_combined <- left_join(smoke_debias_combined, smoke_2016_2020) %>%
  mutate(totalPM=smoke_grid+nonsmoke_2016_2020)

saveRDS(smoke_debias_combined, paste0(result_path, "/smoke_proj/smoke_2046_2055_debias_gcm.rds"))

#----------------------------------------------------------------
## Preprocess and debias the smoke for subregion sensitivity
#----------------------------------------------------------------
smoke_proj <- readRDS(paste0(result_path, 
                             "/smoke_proj/sensitivity/subregion/proj_smoke_gridmetobs_2050_subregion_ssp370.rds")) %>%
  ungroup() %>% filter(grid_id_10km %in% grid_list) %>% filter(grid_id_10km %in% grid_list) %>%
  filter(!is.na(pred_smoke)) %>%
  group_by(grid_id_10km, year, scenario, gcm) %>%
  summarise(pred_smoke=mean(pred_smoke)) %>% ungroup()   ## get the annual mean

smoke_debias_combined <- left_join(smoke_proj, pred_smoke_gfed) %>%
    left_join(smoke_2011_2020) %>%
    mutate(pred_smoke_debias=pred_smoke - smoke_wind_all + smoke_2011_2020) %>%
    rename(smoke_grid = pred_smoke_debias) %>% select(grid_id_10km, year, scenario, smoke_grid, gcm)

smoke_debias_combined[smoke_debias_combined$smoke_grid<0, "smoke_grid"] <- 0

saveRDS(smoke_debias_combined, paste0(result_path, "/smoke_proj/sensitivity/subregion/smoke_2050_10yr_ssp370_subregion.rds"))

#----------------------------------------------------------------
## Preprocess and debias the smoke for the fire uncertainty runs (across many boot)
#----------------------------------------------------------------
smoke_proj <- readRDS(paste0(result_path, 
                             "/smoke_proj/sensitivity/fire_model_boot/proj_smoke_gridmetobs_2050_ssp370_fire_boot.rds")) %>%
  ungroup() %>% filter(grid_id_10km %in% grid_list) %>%
  filter(!is.na(pred_smoke)) %>%
  group_by(grid_id_10km, year, scenario, bootid) %>%
  summarise(pred_smoke=mean(pred_smoke)) %>% ungroup()   ## get the annual mean

nboot <- max(smoke_proj$bootid)
smoke_debias_combined <- NULL
for (bbb in 1:nboot){
  print(bbb)
  
  smoke_tmp <- smoke_proj %>% filter(bootid == bbb)  

  smoke_tmp_debias <- left_join(smoke_tmp, pred_smoke_gfed) %>%
    left_join(smoke_2011_2020) %>%
    mutate(pred_smoke_debias=pred_smoke - smoke_wind_all + smoke_2011_2020) %>%
    rename(smoke_grid = pred_smoke_debias) %>% select(grid_id_10km, year, scenario, smoke_grid, bootid)
  smoke_debias_combined <- bind_rows(smoke_debias_combined, smoke_tmp_debias)
}

smoke_debias_combined[smoke_debias_combined$smoke_grid<0, "smoke_grid"] <- 0

saveRDS(smoke_debias_combined, paste0(result_path, "/smoke_proj/sensitivity/fire_model_boot/smoke_2050_10yr_ssp370_fire_boot.rds"))

#----------------------------------------------------------------
## Preprocess and debias the smoke for the smoke uncertainty runs (across many boot)
#----------------------------------------------------------------
total_smoke_pm <- readRDS(paste0(data_path, "/totalPM/totalPM_smoke_combined_2006_2020.rds")) 

smoke_2011_2020 <- total_smoke_pm %>% 
  filter(year %in% 2011:2020) %>%
  group_by(grid_id_10km) %>%
  summarise_at(c("smokePM_mean"),.funs = mean, na.rm=T) %>%
  rename(smoke_2011_2020=smokePM_mean) %>% ungroup()

pred_smoke_gfed <- readRDS(paste0(result_path, 
                                  "/smoke_proj/sensitivity/smoke_model_boot/gfedDM_final_smoke_pred_bootsame_ssp370.rds")) %>%
  ungroup()   %>%
  filter(year %in% 2011:2020) %>%
  filter(model %in% c("wind_all", "obs")) %>%
  group_by(grid_id_10km, model, bootid) %>%
  summarise(pred_smoke=mean(pred_smoke, na.rm=T)) %>%
  ungroup() %>%
  pivot_wider(values_from=pred_smoke, names_from = model, names_prefix = "smoke_") 

grid_list <- intersect(unique(pred_smoke_gfed$grid_id_10km), 
                       unique(smoke_2011_2020$grid_id_10km))

smoke_2011_2020 <- smoke_2011_2020 %>% filter(grid_id_10km %in% grid_list)
pred_smoke_gfed <- pred_smoke_gfed %>% filter(grid_id_10km %in% grid_list)

smoke_proj <- readRDS(paste0(result_path, 
                             "/smoke_proj/sensitivity/smoke_model_boot/proj_smoke_2050_10yrs_ssp370_smoke_bootsame_gcm_mean.rds")) %>%
  ungroup() %>% filter(grid_id_10km %in% grid_list) %>%
  filter(!is.na(pred_smoke))

nboot <- max(smoke_proj$bootid)
smoke_debias_combined <- NULL
for (bbb in 1:nboot){
  print(bbb)
  
  smoke_tmp <- smoke_proj %>% filter(bootid == bbb)
  
  pred_gfed_tmp <- pred_smoke_gfed %>% filter(bootid == bbb) 
  ##de-bias based on annual mean in 2011-2020
  smoke_tmp_debias <- left_join(smoke_tmp, pred_gfed_tmp) %>%
    left_join(smoke_2011_2020) %>%
    mutate(pred_smoke_debias=pred_smoke - smoke_wind_all + smoke_2011_2020) %>%
    rename(smoke_grid = pred_smoke_debias) %>% select(grid_id_10km, year, scenario, smoke_grid, bootid)
  smoke_debias_combined <- bind_rows(smoke_debias_combined, smoke_tmp_debias)
}

#smoke_debias_combined[smoke_debias_combined$smoke_grid<0, "smoke_grid"] <- 0
smoke_debias_combined <- smoke_debias_combined %>%
  mutate(smoke_grid_pos = smoke_grid,
         smoke_grid_pos=replace(smoke_grid_pos, smoke_grid_pos<0, 0))

saveRDS(smoke_debias_combined, paste0(result_path, "/smoke_proj/sensitivity/smoke_model_boot/smoke_2050_10yr_ssp370_smoke_bootsame.rds"))

# ##----------------------------------------------------
# ## Calculate the pop-weighted PM data
# ##----------------------------------------------------
# ### read the projected smoke in both future and history
# proj_smoke <- readRDS(paste0(result_path, "/smoke_proj/smoke_10yrs_debias_percent.rds")) 
# proj_smoke_pop <- left_join(proj_smoke, pop_df)
# 
# proj_smoke_pop_weight <- proj_smoke_pop %>% group_by(year, scenario) %>%
#   summarise(pop_smoke_2020 = weighted.mean(smoke_grid, w=pop_2020),
#             pop_smoke_2050 = weighted.mean(smoke_grid, w=pop_2050)) %>%
#   ungroup() %>%
#   mutate(year = as.numeric(year) - 5)
# 
# ### just for visualization 
# proj_smoke_pop_weight <- bind_rows(proj_smoke_pop_weight,
#                                    filter(proj_smoke_pop_weight, scenario == "Observation") %>% mutate(year=2018, scenario="ssp126"),
#                                    filter(proj_smoke_pop_weight, scenario == "Observation") %>% mutate(year=2018, scenario="ssp245"),
#                                    filter(proj_smoke_pop_weight, scenario == "Observation") %>% mutate(year=2018, scenario="ssp370"))
# 
# ggplot() +
#   geom_point(data = proj_smoke_pop_weight %>% filter(scenario %in% c("ssp126", "ssp245", "ssp370"), year>2020), 
#              aes(x=year, y=pop_smoke_2050, colour=scenario), size=3, position = position_dodge(width=0.1)) +
#   geom_line(data = proj_smoke_pop_weight %>% filter(scenario %in% c("ssp126", "ssp245", "ssp370"), year!=2020), 
#              aes(x=year, y=pop_smoke_2050, colour=scenario), size=1) +
#   geom_point(data =total_smoke_year, aes(x=year, y=pop_smoke_2020), colour="grey65", size=2, shape=1) +
#   geom_point(data = proj_smoke_pop_weight %>% filter(scenario == "ssp126", year==2018), 
#             aes(x=year, y=pop_smoke_2050), colour="grey65", size=3) +
#   # geom_hline(data = proj_smoke_pop_weight %>% filter(scenario == "Observation"), 
#   #                                                   aes(yintercept=pop_smoke_2020), linetype="dashed", size=1) +
#   scale_colour_manual(values=color_map) +
#   scale_x_continuous(breaks=c(2018,2030,2040,2050),
#                      labels=c("2016-2020",2030,2040,2050)) +
#   theme_classic() +
#   theme(text = element_text(size=16)) +
#   labs(x= "", y="Pop-weighted smoke PM2.5")
# ggsave("pop_smoke_timeseries.png", width=6.5, height=3.5)
# 
#   
# 
# ##----------------------------------------------------
# ## Calculate the pop-weighted PM data across different GCMs
# ##----------------------------------------------------
# proj_smoke_year <- readRDS(paste0(result_path, "/threshold_05/proj_smoke_gridmetobs_2050_10yrs_annual.rds"))
# 
# 
# 
