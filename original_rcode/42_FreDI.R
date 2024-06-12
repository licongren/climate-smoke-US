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
library(cowplot)

rm(list=ls())
gc()

code_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/rcode"
data_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data"
figure_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/figures"
result_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/results/fire_model"
health_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/results/health"


setwd(figure_path)
states_map <- map_data("state")

db_path <- "/Users/mhqiu/BurkeLab Dropbox/Projects/smoke-climate"

#-------------------------------------------------------------------------------
# Compare our estimates with FrEDI
# Written by Minghao
# Last edited Feb 2024
#-------------------------------------------------------------------------------
##convert to 2019 dollars 
inflation_19_15 <- 1.08

#### calculate smoke 
summ_death <- readRDS(paste0(health_path, "/summ_smoke_death_smokePM_CRF_mediangcm.rds")) %>%
  filter(age_group=="all_ages", CRF=="poisson bins") %>%
  select(scenario, death_mean) %>%
  pivot_wider(names_from = scenario, values_from = death_mean)

vsl <- 10.95*1e6
our_estimate <- as.numeric(summ_death$ssp370 - summ_death$Observation)*vsl*(1.02^31)

###--- FrEDI
load(paste0(data_path,"/FrEDI/defaultResults.rda"))

default_2050 <- defaultResults %>% filter(year==2050, region=="National Total",
                                          model%in%c("Average")) %>%
  filter(!sector %in% c("ATS Extreme Temperature", "Extreme Temperature")) %>%
  filter(variant %in% c("2040 Emissions", "N/A", "Without CO2 Fertilization",
                        "Median","No Additional Adaptation"))

sum(default_2050$annual_impacts)*inflation_19_15/1e9

unique(default_2050$gdp_usd)*inflation_19_15/1e9 * c(0.4, 0.8)*1e-2




