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
library(RColorBrewer)

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

#color_scenario <- c("grey65",rev(met.brewer(name="Homer1", n=8, type="discrete"))[c(2,6,8)])
color_scenario <- c("grey65", "#1F2E52", "#E78928","#CB172F") # "#38A2C6", "#870A21"

#-------------------------------------------------------------------------------
# Process the carleton et al  
# Written by Minghao
# Last edited April 2024
#-------------------------------------------------------------------------------
mort_data <- readRDS(paste0(db_proj_path, "/adaptation/data/joined_panels/US_mortality_all-cause_pwtemperature_panel_complete-fips.rds")) %>%
  filter(year>2000) %>%
  rename(fips = fipsihme) %>%
  mutate(fips = as.numeric(fips)) %>%
  filter(fips!=48301, death_type=="all_cause")  %>% # remove a county with only one year of obs (fips==48301, very small pop<100)
  group_by(state_fips, fips, year, pop) %>%
  summarise(n_deaths=sum(n_deaths, na.rm=T)) %>%
  mutate(death_rate=n_deaths/pop) %>%
  ungroup() %>%
  group_by(fips) %>%
  mutate(pop_2019=sum(pop*as.numeric(year==2019))) %>% ungroup() %>%
  filter(year %in% 2001:2010) %>%
  group_by(fips, state_fips, pop_2019) %>%
  summarise(pop_2001_2010=mean(pop),
            death_rate_2001_2010=mean(death_rate)) %>% ungroup() 


future_pop <- read.csv(paste0(data_path, "/health/US_census_population_proj.csv")) %>%
  filter(SEX==0, ORIGIN==0, RACE==0) %>%
  select(YEAR, TOTAL_POP)

pop_change <- future_pop %>% filter(YEAR %in% c(2050,2090,2022)) %>%
  pivot_wider(values_from = TOTAL_POP,
              names_from=YEAR, names_prefix = "pop_") 

mort_pop_data <- mort_data %>% 
  mutate(pop_2050=pop_2019*pop_change$pop_2050/pop_change$pop_2022,
         pop_2090=pop_2019*pop_change$pop_2090/pop_change$pop_2022)


##### temperature deaths data from Carleton 
state_shp <- tigris::states() %>% as.data.frame() %>%
  rename(state_fips=STATEFP,
         state_abbrev = NAME) %>%
  select(state_fips, state_abbrev)

rcp45 <- read.csv(paste0(data_path, "/Carleton_QJE_data/unit_change_in_deathrate_geography_US_states_years_averaged_rcp45_SSP3_quantiles_mean.csv"))
rcp45 <- left_join(rcp45, state_shp)

mort_rcp45 <- left_join(mort_pop_data, rcp45) %>%
  mutate(temp_death_2040_2059=pop_2050*(years_2040_2059/1e5),
         temp_death_2080_2099=pop_2090*(years_2080_2099/1e5)) %>%
 select(fips, state_fips, state_abbrev, 
         pop_2050, pop_2090,
         years_2040_2059, years_2080_2099, 
         temp_death_2040_2059, temp_death_2080_2099)
saveRDS(mort_rcp45, paste0(data_path, "/Carleton_QJE_data/county_temp_deaths_rcp45.rds"))

rcp85 <- read.csv(paste0(data_path, "/Carleton_QJE_data/unit_change_in_deathrate_geography_US_states_years_averaged_rcp85_SSP3_quantiles_mean.csv"))
rcp85 <- left_join(rcp85, state_shp)

mort_rcp85 <- left_join(mort_pop_data, rcp85) %>%
  mutate(temp_death_2040_2059=pop_2050*(years_2040_2059/1e5),
         temp_death_2080_2099=pop_2090*(years_2080_2099/1e5)) %>%
  select(fips, state_fips, state_abbrev, 
         pop_2050, pop_2090,
         years_2040_2059, years_2080_2099, 
         temp_death_2040_2059, temp_death_2080_2099)
saveRDS(mort_rcp85, paste0(data_path, "/Carleton_QJE_data/county_temp_deaths_rcp85.rds"))

