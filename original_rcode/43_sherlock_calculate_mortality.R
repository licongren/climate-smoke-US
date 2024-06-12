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
# Calculate moralities associated smoke exposure
# Written by Minghao
# Last edited Dec 2023
#-------------------------------------------------------------------------------
smoke_future <- readRDS(paste0(output_path, "/US_smoke_county_2046_2055_gcm.rds")) %>% ungroup()
smoke_obs <- readRDS(paste0(output_path, "/US_smoke_county_2011_2020.rds")) %>%
  mutate(gcm="none") %>% ungroup() ## just to be able to use the same function
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
    group_by(age_group, scenario, bootid, gcm) %>%
    summarise(death_smoke=mean(death_smoke)) %>% ungroup() %>% ### take mean across years 
    group_by(age_group, scenario, bootid) %>%
    summarise(death_smoke=median(death_smoke)) %>% ungroup() %>% ### take median across GCMs 
    group_by(age_group, scenario) %>% 
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
##. calculate smoke deaths: annual linear
#---------------------------------------
coef <- readRDS(paste0(output_path, "/coef_linear.rds")) %>% filter(age_group!="under_5")

smoke_death_2020  <- smoke_obs %>%
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
coef <- readRDS(paste0(output_path, "/coef_linear_bins.rds")) %>%
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
summ_death_2020  <- smoke_death_2020 %>% 
  group_by(age_group, year, scenario, bootid, gcm) %>%
  summarise(death_smoke=sum(death_smoke)) %>% ungroup() %>%
  summ_death()

rm(smoke_death_2020)

county_death_2050  <- NULL
summ_death_2050 <- NULL
for (gg in gcm_list){
  for (ss in c("ssp126","ssp245","ssp370")){
    print(gg)
    smoke_pop_tmp <-   smoke_future %>% filter(gcm == gg, scenario==ss) %>% ungroup()
  
  smoke_death_tmp  <- smoke_pop_tmp %>% 
    mutate(smokePM=replace(smokePM, smokePM>10, 10)) %>%
    mutate(bins=cut(smokePM, cutoff)) %>%
    slice(rep(1:n(), each = nboot)) %>%
    mutate(bootid=rep(1:nboot, nrow(smoke_pop_tmp))) %>%
    left_join(coef) %>%
    mutate(coef=replace(coef, is.na(coef), 0)) %>%
    mutate(death_smoke = coef*pop_2050/1e5)
  
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

county_linear_bin_coarse <- bind_rows(county_death_2020, county_death_2050) %>% mutate(CRF="linear bins (coarse)")
summ_linear_bin_coarse <- bind_rows(summ_death_2020, summ_death_2050) %>% mutate(CRF="linear bins (coarse)")

### finer bins
coef <- readRDS(paste0(output_path, "/coef_linear_bins.rds")) %>%
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
summ_death_2020  <- smoke_death_2020 %>%    group_by(age_group, year, scenario, bootid, gcm) %>%   summarise(death_smoke=sum(death_smoke)) %>% ungroup() %>%   summ_death()
rm(smoke_death_2020)

county_death_2050  <- NULL
summ_death_2050 <- NULL
for (gg in gcm_list){
  for (ss in c("ssp126","ssp245","ssp370")){
    print(gg)
    smoke_pop_tmp <-   smoke_future %>% filter(gcm == gg, scenario==ss) %>% ungroup()
    
  smoke_death_tmp  <- smoke_pop_tmp %>% 
    mutate(smokePM=replace(smokePM, smokePM>10, 10)) %>%
    mutate(bins=cut(smokePM, cutoff)) %>%
    slice(rep(1:n(), each = nboot)) %>%
    mutate(bootid=rep(1:nboot, nrow(smoke_pop_tmp))) %>%
    left_join(coef) %>%
    mutate(coef=replace(coef, is.na(coef), 0)) %>%
    mutate(death_smoke = coef*pop_2050/1e5)
  
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

county_linear_bin <- bind_rows(county_death_2020, county_death_2050) %>% mutate(CRF="linear bins")
summ_linear_bin <- bind_rows(summ_death_2020, summ_death_2050) %>% mutate(CRF="linear bins")

#----------------------------------------
##. calculate smoke deaths: annual quadratic
#---------------------------------------
coef <- readRDS(paste0(output_path, "/coef_quadratic.rds")) %>%
  rename(beta1=smokePM_mean,
         beta2=smokePM_mean_sqr) %>%
  filter(age_group!="under_5")

smoke_death_2020  <- smoke_obs %>% 
  slice(rep(1:n(), each = nboot)) %>%
  mutate(bootid=rep(1:nboot, nrow(smoke_obs))) %>%
  left_join(coef) %>%
  mutate(death_smoke = (smokePM*beta1 + smokePM^2*beta2) *pop_2019/1e5) 

county_death_2020  <- county_death(smoke_death_2020)
summ_death_2020  <- smoke_death_2020 %>% group_by(age_group, year, scenario, bootid, gcm) %>% 
  summarise(death_smoke=sum(death_smoke)) %>% ungroup() %>%  summ_death()
rm(smoke_death_2020)

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
    mutate(death_smoke = (smokePM*beta1 + smokePM^2*beta2) *pop_2050/1e5) 
  
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
county_quadratic <- bind_rows(county_death_2020, county_death_2050) %>% mutate(CRF="quadratic")
summ_quadratic <- bind_rows(summ_death_2020, summ_death_2050) %>% mutate(CRF="quadratic")

#----------------------------------------
##. calculate smoke deaths: annual cubic
#---------------------------------------
coef <- readRDS(paste0(output_path, "/coef_cubic.rds")) %>%
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
summ_death_2020  <- smoke_death_2020 %>%
  group_by(age_group, year, scenario, bootid, gcm) %>%  
  summarise(death_smoke=sum(death_smoke)) %>% ungroup() %>%  summ_death()
rm(smoke_death_2020)

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
    mutate(death_smoke = (smokePM*beta1 + smokePM^2*beta2 + smokePM^3*beta3) *pop_2050/1e5) 
  
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

county_cubic <- bind_rows(county_death_2020, county_death_2050) %>% mutate(CRF="cubic")
summ_cubic <- bind_rows(summ_death_2020, summ_death_2050) %>% mutate(CRF="cubic")

#----------------------------------------
##. calculate smoke deaths: annual poisson
#---------------------------------------
coef <- readRDS(paste0(output_path, "/coef_poisson.rds")) %>% filter(age_group!="under_5")

smoke_death_2020  <- smoke_obs %>% 
  slice(rep(1:n(), each = nboot)) %>%
  mutate(bootid=rep(1:nboot, nrow(smoke_obs))) %>%
  left_join(coef) %>%
  mutate(death_smoke = (exp(smokePM*coef)-1)*death_rate_avg*pop_2019) 

county_death_2020  <- county_death(smoke_death_2020)
summ_death_2020  <- smoke_death_2020 %>%    group_by(age_group, year, scenario, bootid, gcm) %>%   summarise(death_smoke=sum(death_smoke)) %>% ungroup() %>%   summ_death()
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
    mutate(death_smoke = (exp(smokePM*coef)-1)*death_rate_avg*pop_2050) 
  
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

county_poisson <- bind_rows(county_death_2020, county_death_2050) %>% mutate(CRF="poisson")
summ_poisson <- bind_rows(summ_death_2020, summ_death_2050) %>% mutate(CRF="poisson")

#----------------------------------------
##. calculate smoke deaths: annual poisson bins
#---------------------------------------
### coarse bins
coef <- readRDS(paste0(output_path, "/coef_poisson_bins.rds")) %>%
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

county_death_2020  <- county_death(smoke_death_2020)
summ_death_2020  <- smoke_death_2020 %>%  group_by(age_group, year, scenario, bootid, gcm) %>% 
  summarise(death_smoke=sum(death_smoke)) %>% ungroup() %>%  summ_death()
rm(smoke_death_2020)

county_death_2050  <- NULL
summ_death_2050 <- NULL
for (gg in gcm_list){
  for (ss in c("ssp126","ssp245","ssp370")){
    print(gg)
    smoke_pop_tmp <-   smoke_future %>% filter(gcm == gg, scenario==ss) %>% ungroup()
    
  smoke_death_tmp  <- smoke_pop_tmp %>% 
    mutate(smokePM=replace(smokePM, smokePM>10, 10)) %>%
    mutate(bins=cut(smokePM, cutoff)) %>%
    slice(rep(1:n(), each = nboot)) %>%
    mutate(bootid=rep(1:nboot, nrow(smoke_pop_tmp))) %>%
    left_join(coef) %>%
    mutate(coef=replace(coef, is.na(coef), 0)) %>%
    mutate(death_smoke =  (exp(coef)-1)*death_rate_avg*pop_2050)
  
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

county_poisson_bin_coarse <- bind_rows(county_death_2020, county_death_2050) %>% mutate(CRF="poisson bins (coarse)")
summ_poisson_bin_coarse <- bind_rows(summ_death_2020, summ_death_2050) %>% mutate(CRF="poisson bins (coarse)")

### finer bins
coef <- readRDS(paste0(output_path, "/coef_poisson_bins.rds")) %>%
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

county_death_2020  <- county_death(smoke_death_2020)
summ_death_2020  <- smoke_death_2020 %>% group_by(age_group, year, scenario, bootid, gcm) %>%  
  summarise(death_smoke=sum(death_smoke)) %>% ungroup() %>% summ_death()
rm(smoke_death_2020)

county_death_2050  <- NULL
summ_death_2050 <- NULL
for (gg in gcm_list){
  for (ss in c("ssp126","ssp245","ssp370")){
    print(gg)
    smoke_pop_tmp <-   smoke_future %>% filter(gcm == gg, scenario==ss) %>% ungroup()
    
  smoke_death_tmp  <- smoke_pop_tmp %>% 
    mutate(smokePM=replace(smokePM, smokePM>10, 10)) %>%
    mutate(bins=cut(smokePM, cutoff)) %>%
    slice(rep(1:n(), each = nboot)) %>%
    mutate(bootid=rep(1:nboot, nrow(smoke_pop_tmp))) %>%
    left_join(coef) %>%
    mutate(coef=replace(coef, is.na(coef), 0)) %>%
    mutate(death_smoke =  (exp(coef)-1)*death_rate_avg*pop_2050)
  
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
county_poisson_bin <- bind_rows(county_death_2020, county_death_2050) %>% mutate(CRF="poisson bins")
summ_poisson_bin <- bind_rows(summ_death_2020, summ_death_2050) %>% mutate(CRF="poisson bins")


# county_combined <- bind_rows(county_linear,county_linear_bin_coarse, county_linear_bin,
#                              county_quadratic, county_cubic,
#                              county_poisson, county_poisson_bin_coarse, county_poisson_bin)
# saveRDS(county_combined, paste0(output_path, "/county_smoke_death_smokePM_CRF.rds"))

summ_combined <- bind_rows(summ_linear,summ_linear_bin_coarse, summ_linear_bin,
                           summ_quadratic, summ_cubic,
                           summ_poisson, summ_poisson_bin_coarse, summ_poisson_bin)
saveRDS(summ_combined, paste0(output_path, "/summ_smoke_death_smokePM_CRF_mediangcm.rds"))



