library(dplyr)   
library(fixest)
library(tidyverse)
library(foreach)
library(doParallel)
library(doSNOW)

rm(list=ls())
gc()

input_path <- "~/oak_space/mhqiu/smoke_climate/input_data"
output_path <- "~/oak_space/mhqiu/smoke_climate/output/health"

#-------------------------------------------------------------------------------
# Calculate moralities associated smoke exposure with all years between 2025 and 2055
# Written by Minghao
# Last edited Dec 2023
#-------------------------------------------------------------------------------
year_list <- c(2025,2035,2045,2055)
smoke_future <- readRDS(paste0(output_path, "/US_smoke_county_2025_2055_10yr_gcm.rds")) %>%
  filter(year %in% year_list) %>% ungroup()

gcm_list <- smoke_future$gcm %>% unique()
nboot <- 500

###---------------------------------------------------------------------------------
## function: get county and total summ of deaths given bootstrap results 
###---------------------------------------------------------------------------------
county_death <- function(df=df){
  county_death <- df %>%
    group_by(year, scenario, gcm, fips, age_group, smokePM) %>%
    summarise(death_smoke=mean(death_smoke)) %>% ungroup() ## take the mean/median across bootrap runs
  return(county_death) 
}

summ_death <- function(df=df){
  summ_death <- df  %>%
    ## commneted out because this step is performed in the for loop already
    # group_by(age_group, year, scenario, bootid, gcm) %>%
    # summarise(death_smoke=sum(death_smoke)) %>% ungroup() %>% ### calculate total mortalities in US
    group_by(year, age_group, scenario, bootid, gcm) %>%
    summarise(death_smoke=mean(death_smoke)) %>% ungroup() %>% ### take mean across years 
    group_by(year, age_group, scenario, bootid) %>%
    summarise(death_smoke=median(death_smoke)) %>% ungroup() %>% ### take median across GCMs 
    group_by(year, age_group, scenario) %>% 
    summarise(death_median=median(death_smoke),  ### get uncertainty across boot ids
              death_mean=mean(death_smoke),
              death_p025=quantile(death_smoke, 0.025),
              death_p975=quantile(death_smoke, 0.975),
              death_p10=quantile(death_smoke, 0.1),
              death_p90=quantile(death_smoke, 0.9),
              death_p25=quantile(death_smoke, 0.25),
              death_p75=quantile(death_smoke, 0.75)) %>% ungroup() 
  
  return(summ_death)
}

#----------------------------------------
##. calculate smoke deaths: annual poisson bins
#---------------------------------------
coef <- readRDS(paste0(output_path, "/coef_poisson_bins.rds")) %>%
  filter(age_group!="under_5",  model=="poisson bins")

cutoff <- c(-0.01, 0.1, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6, 10, 20)

county_death_2050  <- NULL
summ_death_2050 <- NULL
for (gg in gcm_list){
  for (ss in c("ssp126","ssp245","ssp370")){
    print(gg)
    smoke_pop_tmp <-   smoke_future %>% filter(gcm == gg, scenario==ss) %>% ungroup()
    
    smoke_pop_tmp <- smoke_pop_tmp %>%
      mutate(pop_scale_factor=1) %>%
      mutate(pop_scale_factor=replace(pop_scale_factor, year%in%2025:2035, 0.956),
             pop_scale_factor=replace(pop_scale_factor, year%in%2036:2045, 0.986)) ## consider the ratios between pop_2050 to other years
    
    smoke_death_tmp  <- smoke_pop_tmp %>% 
      mutate(smokePM=replace(smokePM, smokePM>10, 10)) %>%
      mutate(bins=cut(smokePM, cutoff)) %>%
      slice(rep(1:n(), each = nboot)) %>%
      mutate(bootid=rep(1:nboot, nrow(smoke_pop_tmp))) %>%
      left_join(coef) %>%
      mutate(coef=replace(coef, is.na(coef), 0)) %>%
      mutate(death_smoke =  (exp(coef)-1)*death_rate_avg*pop_2050*pop_scale_factor)
    
    county_death_tmp  <- county_death(smoke_death_tmp)
    county_death_2050  <- bind_rows(county_death_2050, county_death_tmp)
    
    summ_death_tmp  <- smoke_death_tmp %>% 
      group_by(age_group, year, scenario, bootid, gcm) %>%
      summarise(death_smoke=sum(death_smoke)) %>% ungroup()
    summ_death_2050  <- bind_rows(summ_death_2050, summ_death_tmp)
    rm(smoke_death_tmp)
  }
}

summ_death_2050   <- summ_death(summ_death_2050) %>% mutate(CRF="poisson bins")
county_death_2050 <- county_death_2050 %>% mutate(CRF="poisson bins")

saveRDS(summ_death_2050, paste0(output_path, "/summ_smoke_death_2025_2055_10yrs.rds"))
saveRDS(county_death_2050, paste0(output_path, "/county_smoke_death_2025_2055_10yrs.rds"))



