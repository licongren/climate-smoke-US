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
# Calculate moralities associated smoke exposure
# Written by Minghao
# Last edited Dec 2023
#-------------------------------------------------------------------------------
#### process the county mortality data 
# county_mortality <- read.xlsx(paste0(data_path, "/health/US_county_mortality_1999_2016.xlsx")) %>%
#   filter(Year==2016) %>%
#   mutate(death_rate=Deaths/Population) %>%
#   select(-Notes)
# 
# mortality_2016 <- sum(county_mortality$Deaths) 
# 
# # future_pop <- read.csv(paste0(data_path, "/health/US_census_population_proj.csv")) %>%
# #   filter(SEX==0, ORIGIN==0, RACE==0) %>%
# #    select(YEAR, TOTAL_POP)
# 
# future_death <- read.csv(paste0(data_path, "/health/US_census_death_proj.csv")) %>%
#   filter(SEX==0, RACE_HISP==0) %>%
#   select(YEAR, TOTAL_DEATHS)
# 
# # future_pop_death <- left_join(future_pop, future_death) %>%
# #   mutate(death_rate = TOTAL_DEATHS / TOTAL_POP)
# 
# #future_death_2023 <- filter(future_death, YEAR==2023)$TOTAL_DEATHS %>% unlist()
# future_death_2050 <- filter(future_death, YEAR==2050)$TOTAL_DEATHS %>% unlist()
# 
# county_mortality <- county_mortality %>% mutate(Deaths_2050 = Deaths / mortality_2016 * future_death_2050) %>%
#   rename(Deaths_2016 = Deaths,
#          Population_2016 =  Population)
# 
# colnames(county_mortality)[2] <- "GEOID"
# write.xlsx(county_mortality, paste0(data_path, "/health/county_motality_2016_2050.xlsx"))

###----------------------------------------------------
#### merge with the death rate data at the county level
# mort_data <- readRDS("/Users/mhqiu/BurkeLab Dropbox/Projects/adaptation/data/mortality/processed/usa/age_standardized_rates/age_standardized_rates-county_month-all_sex-all_marital-all_cause.RDS") %>%
#   filter(year>2005, year<2020, race_eth=="all") %>%
#   rename(fips = fipsihme) %>%
#   mutate(fips = as.numeric(fips)) %>%
#   group_by(death_type, age_group, fips, year, pop) %>%
#   summarise(n_deaths=sum(n_deaths, na.rm=T)) %>%
#   mutate(death_rate=n_deaths/pop) %>%
#   ungroup() %>%
#   group_by(death_type, age_group, fips) %>%
#   summarise(pop_2019=sum(pop*as.numeric(year==2019)),
#             death_rate_avg=mean(death_rate)) %>% ungroup()
# 
# future_pop <- read.csv(paste0(data_path, "/health/US_census_population_proj.csv")) %>%
#   filter(SEX==0, ORIGIN==0, RACE==0) %>%
#   rowwise() %>% 
#   mutate(under_5 = sum(c_across(POP_0: POP_4), na.rm = T),
#          under_65 = sum(c_across(POP_0: POP_64), na.rm = T),
#          up_65 = sum(c_across(POP_65: POP_100), na.rm = T)) %>%
#   ungroup() %>%
#   select(YEAR, TOTAL_POP, under_5, under_65, up_65)
# 
# pop_change <- future_pop %>% filter(YEAR %in% c(2050,2022)) %>%
#   pivot_wider(values_from = TOTAL_POP:up_65, 
#               names_from=YEAR) %>%
#   mutate(change_all_ages=TOTAL_POP_2050/TOTAL_POP_2022,
#          change_under_5=under_5_2050/under_5_2022,
#          change_under_65=under_65_2050/under_65_2022,
#          change_65_and_up= up_65_2050 / up_65_2022) %>% 
#   select(change_all_ages, change_under_5, change_under_65,
#          change_65_and_up) 
# 
# data_under5 <- mort_data %>% filter(age_group=="under_5") %>%
#   mutate(pop_2050=pop_2019*pop_change$change_under_5)
# 
# data_under65 <- mort_data %>% filter(age_group=="under_65") %>%
#   mutate(pop_2050=pop_2019*pop_change$change_under_65)
# 
# data_65_and_up <- mort_data %>% filter(age_group=="65_and_up") %>%
#   mutate(pop_2050=pop_2019*pop_change$change_65_and_up)
# 
# data_all_age <- bind_rows(data_under65, data_65_and_up) %>%
#   group_by(death_type, fips) %>% summarise(pop_2050=sum(pop_2050)) %>%
#   ungroup() %>%
#   left_join( mort_data %>% filter(age_group=="all_ages") )
# 
# county_data_2050 <- bind_rows(data_under5, data_under65, data_65_and_up, data_all_age)
# saveRDS(county_data_2050, paste0(data_path, "/health/US_county_pop_2019_2050_age.rds"))


#### convert smoke data to county 
grid_10km <- st_read(paste0(db_proj_path, "/smoke_PM_prediction/data/1_grids/grid_10km_wgs84/grid_10km_wgs84.shp")) %>%
  rename(grid_id_10km=ID)

# county_shp <- counties() %>%
#   st_transform(st_crs(grid_10km))
# 
# crosswalk_county_10km <- st_intersection(grid_10km, county_shp) %>%
#   mutate(intersect_area=st_area(.))

crosswalk_county_10km <- readRDS(paste0(data_path, "/health/county_grid_area_crosswalk.rds")) %>%
  mutate(area=as.numeric(area))
pop_10km_grid <- readRDS(paste0(data_path, "/health/us_population_grid_2020_2022_2050.rds")) %>%
  as.data.frame() %>% select(grid_id_10km, pop_2020)
crosswalk_county_10km <- left_join(crosswalk_county_10km, pop_10km_grid) %>%
  group_by(grid_id_10km) %>%
  mutate(area_weight=area/sum(area)) %>% ungroup()

county_pop <- readRDS(paste0(data_path, "/health/US_county_pop_2019_2050_age.rds"))  %>%
      filter(age_group!="under_5")
####---------------------------------------------------
### read the projected smoke in both future and present
####---------------------------------------------------
smoke_future <-  readRDS(paste0(result_path, "/smoke_proj/smoke_2046_2055_debias_gcm.rds")) %>%
  rename(smokePM = smoke_grid)  %>%
  select(year, scenario, gcm, grid_id_10km, smokePM)

smoke_county <- left_join(smoke_future, crosswalk_county_10km) %>%
  filter(!is.na(GEOID)) %>%
  group_by(year, gcm, scenario, GEOID) %>%
  summarise(smokePM = weighted.mean(smokePM, w=area_weight *pop_2020, na.rm=T)) %>%
  ungroup() %>%
  mutate(GEOID = as.numeric(GEOID)) %>%  rename(fips=GEOID) %>%
  left_join(county_pop) %>% filter(!is.na(pop_2050))
saveRDS(smoke_county, paste0(data_path, "/health/US_smoke_county_2026_2045_gcm.rds"))

smoke_obs <-  readRDS(paste0(data_path, "/totalPM/totalPM_smoke_combined_2006_2020.rds"))  %>%
  filter(year %in% 2011:2020) %>%
  rename(smokePM = smokePM_mean) %>%
  filter(grid_id_10km %in% unique(smoke_future$grid_id_10km))

smoke_county <- left_join(smoke_obs, crosswalk_county_10km) %>%
  filter(!is.na(GEOID)) %>%
  group_by(year, GEOID) %>%
  summarise(smokePM = weighted.mean(smokePM, w=area_weight *pop_2020, na.rm=T)) %>%
  ungroup() %>%
  mutate(GEOID = as.numeric(GEOID)) %>%  rename(fips=GEOID) %>%
  left_join(county_pop) %>% filter(!is.na(pop_2050)) %>%
  mutate(scenario="Observation")
saveRDS(smoke_county, paste0(data_path, "/health/US_smoke_county_2011_2020.rds"))

totalPM_county <- left_join(smoke_obs, crosswalk_county_10km) %>%
  filter(!is.na(GEOID)) %>%
  group_by(year, GEOID) %>%
  summarise(smokePM = weighted.mean(smokePM, w=area_weight *pop_2020, na.rm=T),
            totalPM = weighted.mean(totalPM, w=area_weight *pop_2020, na.rm=T)) %>%
  ungroup() %>%
  mutate(GEOID = as.numeric(GEOID)) %>%  rename(fips=GEOID) %>%
  left_join(county_pop) %>% filter(!is.na(pop_2050)) %>%
  mutate(scenario="Observation")
saveRDS(totalPM_county, paste0(data_path, "/health/US_totalPM_county_2011_2020.rds"))

######## also process other years (10-year mean or individual years)
smoke_future <-  readRDS(paste0(result_path, "/smoke_proj/smoke_10yrs_debias_gcm.rds")) %>%
  rename(smokePM = smoke_grid)  %>%
  select(year, scenario, gcm, grid_id_10km, smokePM)

smoke_county <- left_join(smoke_future, crosswalk_county_10km) %>%
  filter(!is.na(GEOID)) %>%
  group_by(year, gcm, scenario, GEOID) %>%
  summarise(smokePM = weighted.mean(smokePM, w=area_weight *pop_2020, na.rm=T)) %>%
  ungroup() %>%
  mutate(GEOID = as.numeric(GEOID)) %>%  rename(fips=GEOID) %>%
  left_join(county_pop) %>% filter(!is.na(pop_2050))
saveRDS(smoke_county, paste0(data_path, "/health/US_smoke_county_2025_2055_10yr_gcm.rds"))



smoke_future <- readRDS(paste0(data_path, "/health/US_smoke_county_2046_2055_gcm.rds"))
smoke_obs <- readRDS(paste0(data_path, "/health/US_smoke_county_2011_2020.rds")) %>%
  mutate(gcm="none") ## just to be able to use the same function
gcm_list <- smoke_future$gcm %>% unique()
nboot <- 500
###---------------------------------------------------------------------------------
## function: get county and total summ of deaths given bootstrap results 
###---------------------------------------------------------------------------------
county_death <- function(df=df){
  county_death <- df %>%
    group_by(year, scenario, gcm, fips, age_group, smokePM) %>%
    summarise(death_smoke=median(death_smoke)) %>% ungroup() ## take the median across bootrap runs
 
  return(county_death) 
}

summ_death <- function(df=df){
  summ_death <- df  %>%
    ## commneted out because this step is perfromed in the for loop already
    # group_by(age_group, year, scenario, bootid, gcm) %>%
    # summarise(death_smoke=sum(death_smoke)) %>% ungroup() %>% ### calculate total mortalities in US
    group_by(age_group, scenario, bootid) %>%
    summarise(death_smoke=mean(death_smoke)) %>% ungroup() %>% ### take mean across years and GCMs
    group_by(age_group, scenario) %>% 
    summarise(death_median=median(death_smoke),  ### get uncertainty across boot ids
              death_p025=quantile(death_smoke, 0.025),
              death_p975=quantile(death_smoke, 0.975),
              death_p25=quantile(death_smoke, 0.25),
              death_p75=quantile(death_smoke, 0.75)) %>% ungroup() 
  
  return(summ_death)
}

#----------------------------------------
##. calculate smoke deaths: annual linear
#---------------------------------------
coef <- readRDS(paste0(health_path, "/coef_linear.rds")) %>% filter(age_group!="under_5")

smoke_death_2020  <- smoke_obs %>% ungroup() %>%
  slice(rep(1:n(), each = nboot)) %>%
  mutate(bootid=rep(1:nboot, nrow(smoke_obs))) %>%
  left_join(coef) %>%
  mutate(death_smoke = (smokePM*coef) *pop_2019/1e5)

county_death_2020  <- county_death(smoke_death_2020)
summ_death_2020  <- smoke_death_2020 %>% 
           group_by(age_group, year, scenario, bootid, gcm) %>%
                   summarise(death_smoke=sum(death_smoke)) %>% ungroup() %>%
                       summ_death()

rm(smoke_death_2020)

### use for-loop instead of join to save memory
county_death_2050  <- NULL
summ_death_2050 <- NULL
for (gg in gcm_list){
  for (ss in c("ssp126","ssp245","ssp370")){
  print(gg)
  smoke_pop_tmp <-   smoke_future %>% filter(gcm == gg, scenario==ss) %>% ungroup()
  
  smoke_death_tmp  <- smoke_pop_tmp %>% 
    slice(rep(1:n(), each = nboot)) %>%
    mutate(bootid=rep(1:nboot, nrow(smoke_pop_tmp))) %>%
    left_join(coef) %>%
    mutate(death_smoke = (smokePM*coef) *pop_2050/1e5)
  
  county_death_tmp  <- county_death(smoke_death_tmp)
  county_death_2050  <- bind_rows(county_death_2050, county_death_tmp)
  
  summ_death_tmp  <- smoke_death_tmp %>% 
    group_by(age_group, year, scenario, bootid, gcm) %>%
    summarise(death_smoke=sum(death_smoke)) %>% ungroup()
  summ_death_2050  <- bind_rows(summ_death_2050, summ_death_tmp)
    
  rm(smoke_death_tmp)
  }
}

summ_death_2050   <- summ_death(summ_death_2050)

county_linear <- bind_rows(county_death_2020, county_death_2050) %>% mutate(CRF="linear")
summ_linear <- bind_rows(summ_death_2020, summ_death_2050) %>% mutate(CRF="linear")

#----------------------------------------
##. calculate smoke deaths: annual bins
#---------------------------------------
### coarse bins
coef <- readRDS(paste0(health_path, "/coef_linear_bins.rds")) %>%
  filter(age_group!="under_5",  model=="linear bins (coarse)")

cutoff <- c(-0.01,0.1, 0.25, 0.5, 0.75,1, 2.5, 5, 10, 20)
smoke_death_2020 <- smoke_obs %>%
  mutate(smokePM=replace(smokePM, smokePM>10, 10)) %>%
  mutate(bins=cut(smokePM, cutoff)) %>%
  slice(rep(1:n(), each = nboot)) %>%
  mutate(bootid=rep(1:nboot, nrow(smoke_obs))) %>%
  left_join(coef) %>%
  mutate(coef=replace(coef, is.na(coef), 0)) %>%
  mutate(death_smoke = coef*pop_2019/1e5)

county_death_2020  <- county_death(smoke_death_2020)
summ_death_2020  <- summ_death(smoke_death_2020)
rm(smoke_death_2020)

county_death_2050  <- NULL
summ_death_2050 <- NULL
for (ss in c("ssp126","ssp245","ssp370")){
  smoke_pop_tmp <-   smoke_future %>% filter(scenario == ss) 
  
  smoke_death_tmp  <- smoke_pop_tmp %>% 
    mutate(smokePM=replace(smokePM, smokePM>10, 10)) %>%
    mutate(bins=cut(smokePM, cutoff)) %>%
    slice(rep(1:n(), each = nboot)) %>%
    mutate(bootid=rep(1:nboot, nrow(smoke_pop_tmp))) %>%
    left_join(coef) %>%
    mutate(coef=replace(coef, is.na(coef), 0)) %>%
    mutate(death_smoke = coef*pop_2050/1e5)
  
  county_death_tmp  <- county_death(smoke_death_tmp)
  summ_death_tmp  <- summ_death(smoke_death_tmp)
  
  county_death_2050  <- bind_rows(county_death_2050, county_death_tmp)
  summ_death_2050    <- bind_rows(summ_death_2050, summ_death_tmp)
  rm(smoke_death_tmp)
}

county_linear_bin_coarse <- bind_rows(county_death_2020, county_death_2050) %>% mutate(CRF="linear bins (coarse)")
summ_linear_bin_coarse <- bind_rows(summ_death_2020, summ_death_2050) %>% mutate(CRF="linear bins (coarse)")

### finer bins
coef <- readRDS(paste0(health_path, "/coef_linear_bins.rds")) %>%
  filter(age_group!="under_5",  model=="linear bins")

cutoff <- c(-0.01, 0.1, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6, 10, 20)
smoke_death_2020 <- smoke_obs %>%
  mutate(smokePM=replace(smokePM, smokePM>10, 10)) %>%
  mutate(bins=cut(smokePM, cutoff)) %>%
  slice(rep(1:n(), each = nboot)) %>%
  mutate(bootid=rep(1:nboot, nrow(smoke_obs))) %>%
  left_join(coef) %>%
  mutate(coef=replace(coef, is.na(coef), 0)) %>%
  mutate(death_smoke = coef*pop_2019/1e5)

county_death_2020  <- county_death(smoke_death_2020)
summ_death_2020  <- summ_death(smoke_death_2020)
rm(smoke_death_2020)

county_death_2050  <- NULL
summ_death_2050 <- NULL
for (ss in c("ssp126","ssp245","ssp370")){
  smoke_pop_tmp <-   smoke_future %>% filter(scenario == ss) 
  
  smoke_death_tmp  <- smoke_pop_tmp %>% 
    mutate(smokePM=replace(smokePM, smokePM>10, 10)) %>%
    mutate(bins=cut(smokePM, cutoff)) %>%
    slice(rep(1:n(), each = nboot)) %>%
    mutate(bootid=rep(1:nboot, nrow(smoke_pop_tmp))) %>%
    left_join(coef) %>%
    mutate(coef=replace(coef, is.na(coef), 0)) %>%
    mutate(death_smoke = coef*pop_2050/1e5)
  
  county_death_tmp  <- county_death(smoke_death_tmp)
  summ_death_tmp  <- summ_death(smoke_death_tmp)
  
  county_death_2050  <- bind_rows(county_death_2050, county_death_tmp)
  summ_death_2050    <- bind_rows(summ_death_2050, summ_death_tmp)
  rm(smoke_death_tmp)
}

county_linear_bin <- bind_rows(county_death_2020, county_death_2050) %>% mutate(CRF="linear bins")
summ_linear_bin <- bind_rows(summ_death_2020, summ_death_2050) %>% mutate(CRF="linear bins")

#----------------------------------------
##. calculate smoke deaths: annual quadratic
#---------------------------------------
coef <- readRDS(paste0(health_path, "/coef_quadratic.rds")) %>%
  rename(beta1=smokePM_mean,
         beta2=smokePM_mean_sqr) %>%
  filter(age_group!="under_5")

smoke_death_2020  <- smoke_obs %>% 
  slice(rep(1:n(), each = nboot)) %>%
  mutate(bootid=rep(1:nboot, nrow(smoke_obs))) %>%
  left_join(coef) %>%
  mutate(death_smoke = (smokePM*beta1 + smokePM^2*beta2) *pop_2019/1e5) 

county_death_2020  <- county_death(smoke_death_2020)
summ_death_2020  <- summ_death(smoke_death_2020)
rm(smoke_death_2020)

county_death_2050  <- NULL
summ_death_2050 <- NULL
for (ss in c("ssp126","ssp245","ssp370")){
  smoke_pop_tmp <-   smoke_future %>% filter(scenario == ss) 
  
  smoke_death_tmp  <- smoke_pop_tmp %>% 
    slice(rep(1:n(), each = nboot)) %>%
    mutate(bootid=rep(1:nboot, nrow(smoke_pop_tmp))) %>%
    left_join(coef) %>%
    mutate(death_smoke = (smokePM*beta1 + smokePM^2*beta2) *pop_2050/1e5) 
  
  county_death_tmp  <- county_death(smoke_death_tmp)
  summ_death_tmp  <- summ_death(smoke_death_tmp)
  
  county_death_2050  <- bind_rows(county_death_2050, county_death_tmp)
  summ_death_2050    <- bind_rows(summ_death_2050, summ_death_tmp)
  rm(smoke_death_tmp)
}

county_quadratic <- bind_rows(county_death_2020, county_death_2050) %>% mutate(CRF="quadratic")
summ_quadratic <- bind_rows(summ_death_2020, summ_death_2050) %>% mutate(CRF="quadratic")

#----------------------------------------
##. calculate smoke deaths: annual cubic
#---------------------------------------
coef <- readRDS(paste0(health_path, "/coef_cubic.rds")) %>%
  rename(beta1=smokePM_mean,
         beta2=smokePM_mean_sqr,
         beta3=smokePM_mean_cub) %>%
  filter(age_group!="under_5")

smoke_death_2020  <- smoke_obs %>% 
  slice(rep(1:n(), each = nboot)) %>%
  mutate(bootid=rep(1:nboot, nrow(smoke_obs))) %>%
  left_join(coef) %>%
  mutate(death_smoke = (smokePM*beta1 + smokePM^2*beta2 + smokePM^3*beta3) *pop_2019/1e5) 

county_death_2020  <- county_death(smoke_death_2020)
summ_death_2020  <- summ_death(smoke_death_2020)
rm(smoke_death_2020)

county_death_2050  <- NULL
summ_death_2050 <- NULL
for (ss in c("ssp126","ssp245","ssp370")){
  smoke_pop_tmp <-   smoke_future %>% filter(scenario == ss) 
  
  smoke_death_tmp  <- smoke_pop_tmp %>% 
    slice(rep(1:n(), each = nboot)) %>%
    mutate(bootid=rep(1:nboot, nrow(smoke_pop_tmp))) %>%
    left_join(coef) %>%
    mutate(death_smoke = (smokePM*beta1 + smokePM^2*beta2 + smokePM^3*beta3) *pop_2050/1e5) 
  
  county_death_tmp  <- county_death(smoke_death_tmp)
  summ_death_tmp  <- summ_death(smoke_death_tmp)
  
  county_death_2050  <- bind_rows(county_death_2050, county_death_tmp)
  summ_death_2050    <- bind_rows(summ_death_2050, summ_death_tmp)
  rm(smoke_death_tmp)
}

county_cubic <- bind_rows(county_death_2020, county_death_2050) %>% mutate(CRF="cubic")
summ_cubic <- bind_rows(summ_death_2020, summ_death_2050) %>% mutate(CRF="cubic")

#----------------------------------------
##. calculate smoke deaths: annual poisson
#---------------------------------------
coef <- readRDS(paste0(health_path, "/coef_poisson.rds")) %>% filter(age_group!="under_5")

smoke_death_2020  <- smoke_obs %>% 
  slice(rep(1:n(), each = nboot)) %>%
  mutate(bootid=rep(1:nboot, nrow(smoke_obs))) %>%
  left_join(coef) %>%
  mutate(death_smoke = (exp(smokePM*coef)-1)*death_rate_avg*pop_2019) 

county_death_2020  <- county_death(smoke_death_2020)
summ_death_2020  <- summ_death(smoke_death_2020)
rm(smoke_death_2020)

### use for-loop instead of join to save memory
county_death_2050  <- NULL
summ_death_2050 <- NULL
for (ss in c("ssp126","ssp245","ssp370")){
  smoke_pop_tmp <-   smoke_future %>% filter(scenario == ss) 
  
  smoke_death_tmp  <- smoke_pop_tmp %>% 
    slice(rep(1:n(), each = nboot)) %>%
    mutate(bootid=rep(1:nboot, nrow(smoke_pop_tmp))) %>%
    left_join(coef) %>%
    mutate(death_smoke = (exp(smokePM*coef)-1)*death_rate_avg*pop_2050) 
  
  county_death_tmp  <- county_death(smoke_death_tmp)
  summ_death_tmp  <- summ_death(smoke_death_tmp)
  
  county_death_2050  <- bind_rows(county_death_2050, county_death_tmp)
  summ_death_2050    <- bind_rows(summ_death_2050, summ_death_tmp)
  rm(smoke_death_tmp)
}

county_poisson <- bind_rows(county_death_2020, county_death_2050) %>% mutate(CRF="poisson")
summ_poisson <- bind_rows(summ_death_2020, summ_death_2050) %>% mutate(CRF="poisson")

#----------------------------------------
##. calculate smoke deaths: annual poisson bins
#---------------------------------------
### coarse bins
coef <- readRDS(paste0(health_path, "/coef_poisson_bins.rds")) %>%
  filter(age_group!="under_5",  model=="poisson bins (coarse)")

cutoff <- c(-0.01,0.1, 0.25, 0.5, 0.75,1, 2.5, 5, 10, 20)
smoke_death_2020 <- smoke_obs %>%
  mutate(smokePM=replace(smokePM, smokePM>10, 10)) %>%
  mutate(bins=cut(smokePM, cutoff)) %>%
  slice(rep(1:n(), each = nboot)) %>%
  mutate(bootid=rep(1:nboot, nrow(smoke_obs))) %>%
  left_join(coef) %>%
  mutate(coef=replace(coef, is.na(coef), 0)) %>%
  mutate(death_smoke = (exp(coef)-1)*death_rate_avg*pop_2019)

county_death_2020  <- county_death(smoke_death_2020))
summ_death_2020  <- summ_death(smoke_death_2020))
rm(smoke_death_2020)

county_death_2050  <- NULL
summ_death_2050 <- NULL
for (ss in c("ssp126","ssp245","ssp370")){
  smoke_pop_tmp <-   smoke_future %>% filter(scenario == ss) 
  
  smoke_death_tmp  <- smoke_pop_tmp %>% 
    mutate(smokePM=replace(smokePM, smokePM>10, 10)) %>%
    mutate(bins=cut(smokePM, cutoff)) %>%
    slice(rep(1:n(), each = nboot)) %>%
    mutate(bootid=rep(1:nboot, nrow(smoke_pop_tmp))) %>%
    left_join(coef) %>%
    mutate(coef=replace(coef, is.na(coef), 0)) %>%
    mutate(death_smoke =  (exp(coef)-1)*death_rate_avg*pop_2050)
  
  county_death_tmp  <- county_death(smoke_death_tmp)
  summ_death_tmp  <- summ_death(smoke_death_tmp)
  
  county_death_2050  <- bind_rows(county_death_2050, county_death_tmp)
  summ_death_2050    <- bind_rows(summ_death_2050, summ_death_tmp)
  rm(smoke_death_tmp)
}

county_poisson_bin_coarse <- bind_rows(county_death_2020, county_death_2050) %>% mutate(CRF="poisson bins (coarse)")
summ_poisson_bin_coarse <- bind_rows(summ_death_2020, summ_death_2050) %>% mutate(CRF="poisson bins (coarse)")

### finer bins
coef <- readRDS(paste0(health_path, "/coef_poisson_bins.rds")) %>%
  filter(age_group!="under_5",  model=="poisson bins")

cutoff <- c(-0.01, 0.1, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6, 10, 20)
smoke_death_2020 <- smoke_obs %>%
  mutate(smokePM=replace(smokePM, smokePM>10, 10)) %>%
  mutate(bins=cut(smokePM, cutoff)) %>%
  slice(rep(1:n(), each = nboot)) %>%
  mutate(bootid=rep(1:nboot, nrow(smoke_obs))) %>%
  left_join(coef) %>%
  mutate(coef=replace(coef, is.na(coef), 0)) %>%
  mutate(death_smoke = (exp(coef)-1)*death_rate_avg*pop_2019)

county_death_2020  <- county_death(smoke_death_2020))
summ_death_2020  <- summ_death(smoke_death_2020))
rm(smoke_death_2020)

county_death_2050  <- NULL
summ_death_2050 <- NULL
for (ss in c("ssp126","ssp245","ssp370")){
  smoke_pop_tmp <-   smoke_future %>% filter(scenario == ss) 
  
  smoke_death_tmp  <- smoke_pop_tmp %>% 
    mutate(smokePM=replace(smokePM, smokePM>10, 10)) %>%
    mutate(bins=cut(smokePM, cutoff)) %>%
    slice(rep(1:n(), each = nboot)) %>%
    mutate(bootid=rep(1:nboot, nrow(smoke_pop_tmp))) %>%
    left_join(coef) %>%
    mutate(coef=replace(coef, is.na(coef), 0)) %>%
    mutate(death_smoke =  (exp(coef)-1)*death_rate_avg*pop_2050)
  
  county_death_tmp  <- county_death(smoke_death_tmp)
  summ_death_tmp  <- summ_death(smoke_death_tmp)
  
  county_death_2050  <- bind_rows(county_death_2050, county_death_tmp)
  summ_death_2050    <- bind_rows(summ_death_2050, summ_death_tmp)
  rm(smoke_death_tmp)
}

county_poisson_bin <- bind_rows(county_death_2020, county_death_2050) %>% mutate(CRF="poisson bins")
summ_poisson_bin <- bind_rows(summ_death_2020, summ_death_2050) %>% mutate(CRF="poisson bins")


county_combined <- bind_rows(county_linear,county_linear_bin_coarse, county_linear_bin,
                             county_quadratic, county_cubic,
                             county_poisson, county_poisson_bin_coarse, county_poisson_bin)
saveRDS(county_combined, paste0(health_path, "/county_smoke_death_smokePM_CRF.rds"))

summ_combined <- bind_rows(summ_linear,summ_linear_bin_coarse, summ_linear_bin,
                             summ_quadratic, summ_cubic,
                             summ_poisson, summ_poisson_bin_coarse, summ_poisson_bin)
saveRDS(summ_combined, paste0(health_path, "/summ_smoke_death_smokePM_CRF.rds"))


##### Plot
summ_combined$CRF <- factor(summ_combined$CRF, levels=
                              c("linear","linear bins", "linear bins (coarse)",
                                "poisson","poisson bins", "poisson bins (coarse)",
                                "cubic"))
ggplot(summ_combined %>% filter(!age_group%in%c("under_5", "under_65"), CRF!="quadratic"),
       aes(x=scenario, y=death_median, fill=CRF)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  geom_errorbar(aes(ymin=death_p25, ymax=death_p75), stat = "identity", 
                position = position_dodge2(.9, padding = .9), linewidth=0.7) +
  scale_fill_manual(values=c(met.brewer(name="Veronese", n=10, type="continuous")[c(1:3,5,7,9)],
                    "blue4"))+
  scale_colour_manual(values=c(met.brewer(name="Veronese", n=10, type="continuous")[c(1:3,5,7,9)],
                             "blue4"))+
  theme_classic() +
  theme(text = element_text(size=20)) +
  labs(x="", y = "Annual mortalities\ndue to smoke PM2.5") +
  theme(text = element_text(size=20)) +
  facet_wrap(~age_group, ncol=1, scale="free_y")
ggsave("smoke_deaths_crf_25_75.pdf", width=10, height=6)


##### examine the difference across the bins 
smoke_crf <- smoke_deaths_combined %>% 
  select(year, scenario, fips, smokePM, age_group, death_smoke, CRF) %>%
  pivot_wider(id_cols = year:age_group, names_from = CRF, values_from = death_smoke) %>%
  filter(age_group=="all_ages") 

smoke_crf %>% arrange(desc(linear))

colnames(smoke_crf)[9] <- c("linear_bins_coarse")

aa <- smoke_crf %>% #filter(scenario=="observation") %>%
  mutate(diff= linear - quadratic) %>%
  arrange(desc(diff))

cutoff <-  c(-0.01, 0.1, 0.5, 1, 2, 3, 4, 5, 6, 10,20)
smoke_crf <- smoke_deaths_combined %>% 
  mutate(bins=cut(smokePM, cutoff)) %>%
  group_by(year, scenario, bins, age_group, CRF) %>%
  summarise(death_smoke=sum(death_smoke)) %>%
  pivot_wider(id_cols = year:age_group, names_from = CRF, values_from = death_smoke) %>%
  filter(age_group=="all_ages") 


  
  
  
rr_10 <- 1.08
smoke_county_mortality <- smoke_county_mortality %>%
  mutate(mort_loglinear_2016=Deaths_2016*(rr_10^(smokePM/10)-1),
         mort_loglinear_2050=Deaths_2050*(rr_10^(smokePM/10)-1),
         mort_smoke_2016 = Population_2016*smokePM*2.434*1e-5)
         #mort_smoke_2050 = Population_2016*smokePM*2.434*1e-5,)

mortality_summ <- smoke_county_mortality %>% group_by(year, scenario) %>%
  summarise_at(c("mort_loglinear_2016","mort_loglinear_2050","mort_smoke_2016"), .funs = sum, na.rm=T) %>%
  ungroup() %>%
  mutate(year= as.numeric(year) - 5)


### Plotting the results 
mortality_tmp <- filter(mortality_summ, year==2050) %>% bind_rows(
  filter(mortality_summ, is.na(year)))

ggplot(mortality_tmp) +
  geom_bar(aes(x=NA, y=mort_loglinear_2050, fill=scenario), stat = "identity",
           position = position_dodge()) +
  scale_fill_manual(values=c("grey60",color_map)) + 
  theme_classic() +
  labs(x="", y = "Annual mortalities\ndue to smoke PM2.5") +
  theme(text = element_text(size=14))
ggsave("total_mortality.png", width=5, height=3.5)

mortality_plot_data <- bind_rows(mortality_summ ,
                                   filter(mortality_summ, is.na(year)) %>% mutate(year=2018, scenario="ssp126"),
                                   filter(mortality_summ, is.na(year)) %>% mutate(year=2018, scenario="ssp245"),
                                   filter(mortality_summ, is.na(year)) %>% mutate(year=2018, scenario="ssp370"))

mortality_plot_data
ggplot() +
  geom_point(data = mortality_plot_data %>% filter(scenario %in% c("ssp126", "ssp245", "ssp370"), year>2020), 
             aes(x=year, y=mort_loglinear_2050, colour=scenario), size=3, position = position_dodge(width=0.1)) +
  # geom_line(data = mortality_plot_data %>% filter(scenario %in% c("ssp126", "ssp245", "ssp370"), year!=2020), 
  #           aes(x=year, y=mort_loglinear_2050, colour=scenario), size=1) +
  geom_point(data = mortality_plot_data %>% filter(scenario == "Observation"),
             aes(x=year, y=mort_loglinear_2016), colour="grey65", size=2, shape=1) +
  geom_point(data = mortality_plot_data %>% filter(scenario == "ssp126", year==2018), 
             aes(x=year, y=mort_loglinear_2016), colour="grey65", size=3) +
  # geom_hline(data = proj_smoke_pop_weight %>% filter(scenario == "Observation"), 
  #                                                   aes(yintercept=pop_smoke_2020), linetype="dashed", size=1) +
  scale_colour_manual(values=color_map) +
  scale_x_continuous(breaks=c(2018,2030,2040,2050),
                     labels=c("2016-2020",2030,2040,2050)) +
  theme_classic() +
  theme(text = element_text(size=16)) +
  labs(x= "", y="Deaths due to smoke PM2.5")
ggsave("death_smoke_timeseries.png", width=6.5, height=3.5)
  

