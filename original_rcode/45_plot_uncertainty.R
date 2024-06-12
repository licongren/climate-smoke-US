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

###### Fig 5A: quanitfy the uncertainty in ssp370 mortality across four components
county_death <- readRDS(paste0(health_path, "/county_smoke_death_smokePM_CRF_bootmean.rds")) %>%
  filter(age_group=="all_ages", CRF=="poisson bins")

gcm_uncertainty <- county_death  %>%
  group_by(scenario, gcm, year) %>%
  summarise(death_smoke=sum(death_smoke)) %>% ungroup() %>% ##sum over CONUS
  group_by(scenario, gcm) %>%
  summarise(death_smoke=mean(death_smoke)) %>% ungroup() %>% ##mean across years
  group_by(scenario) %>%
  summarise(death_mean=mean(death_smoke),
            death_median=median(death_smoke),
            death_p025=quantile(death_smoke, 0.025),
            death_p975=quantile(death_smoke, 0.975),
            death_p10=quantile(death_smoke, 0.1),
            death_p90=quantile(death_smoke, 0.9),
            death_p25=quantile(death_smoke, 0.25),
            death_p75=quantile(death_smoke, 0.75)) %>% ungroup() %>%
  mutate(type="GCM uncertainty") ## median over GCMs 

crf_uncertainty <- readRDS(paste0(health_path, "/summ_smoke_death_smokePM_CRF_mediangcm.rds")) %>%
  filter(age_group=="all_ages", CRF=="poisson bins") %>%
  mutate(type="CRF uncertainty")

fire_uncertainty <- readRDS(paste0(health_path, "/US_smoke_death_fire_boot.rds")) %>%
  group_by(scenario, year, bootid) %>%
  summarise(death_smoke=sum(death_smoke)) %>% ungroup() %>% ##sum over CONUS
  group_by(scenario, bootid) %>%
  summarise(death_smoke=mean(death_smoke)) %>% ungroup() %>% ##mean across years
  group_by(scenario) %>%
  summarise(death_mean=mean(death_smoke),
            death_median=median(death_smoke),
            death_p025=quantile(death_smoke, 0.025),
            death_p975=quantile(death_smoke, 0.975),
            death_p10=quantile(death_smoke, 0.1),
            death_p90=quantile(death_smoke, 0.9),
            death_p25=quantile(death_smoke, 0.25),
            death_p75=quantile(death_smoke, 0.75)) %>% ungroup() %>%
  mutate(type="fire uncertainty") ## median over GCMs 


# smoke_uncertainty <- readRDS(paste0(health_path, "/US_smoke_death_smoke_boot_gcm_mean.rds")) %>%
#   group_by(scenario, year, bootid, gcm) %>%
#   summarise(death_smoke=sum(death_smoke)) %>% ungroup() %>% ##sum over CONUS
#   group_by(scenario, bootid, gcm) %>%
#   summarise(death_smoke=mean(death_smoke)) %>% ungroup() %>% ##mean across years 
#   group_by(scenario, bootid) %>%
#   summarise(death_smoke=mean(death_smoke)) %>% ungroup() %>% ##median across GCMs
#   group_by(scenario) %>%
#   summarise(death_mean=mean(death_smoke),
#             death_median=median(death_smoke),
#             death_p025=quantile(death_smoke, 0.025),
#             death_p975=quantile(death_smoke, 0.975),
#             death_p10=quantile(death_smoke, 0.1),
#             death_p90=quantile(death_smoke, 0.9),
#             death_p25=quantile(death_smoke, 0.25),
#             death_p75=quantile(death_smoke, 0.75)) %>% ungroup() %>%
#   mutate(type="smoke uncertainty") ## median over GCMs 

smoke_uncertainty <- readRDS(paste0(health_path, "/US_smoke_death_smoke_bootsame.rds")) %>%
  group_by(scenario, year, bootid) %>%
  summarise(death_smoke=sum(death_smoke)) %>% ungroup() %>% ##sum over CONUS
  group_by(scenario, bootid) %>%
  summarise(death_smoke=mean(death_smoke)) %>% ungroup() %>% ##mean across years 
  group_by(scenario) %>%
  summarise(death_mean=mean(death_smoke),
            death_median=median(death_smoke),
            death_p025=quantile(death_smoke, 0.025),
            death_p975=quantile(death_smoke, 0.975),
            death_p10=quantile(death_smoke, 0.1),
            death_p90=quantile(death_smoke, 0.9),
            death_p25=quantile(death_smoke, 0.25),
            death_p75=quantile(death_smoke, 0.75)) %>% ungroup() %>%
  mutate(type="smoke uncertainty")

#----------  analyze the combined uncertainty with Monte-Carlo simulations 
smoke_death_mc <-  readRDS(paste0(health_path, "/ssp370_smoke_death_MC_1_500.rds"))

mc_uncertainty <- smoke_death_mc  %>%
  group_by(age_group) %>% 
  summarise(death_median=median(death_smoke),  ### get uncertainty across boot ids
            death_mean=mean(death_smoke),  ### get uncertainty across boot ids
            death_p025=quantile(death_smoke, 0.025),
            death_p975=quantile(death_smoke, 0.975),
            death_p10=quantile(death_smoke, 0.1),
            death_p90=quantile(death_smoke, 0.9),
            death_p25=quantile(death_smoke, 0.25),
            death_p75=quantile(death_smoke, 0.75)) %>% ungroup() %>%
  filter(age_group=="all_ages") %>% mutate(scenario="ssp370", type="combined uncertainty")

uncertainty_combined <- bind_rows(fire_uncertainty, gcm_uncertainty, crf_uncertainty, mc_uncertainty)

ggplot(uncertainty_combined %>% filter(scenario == "ssp370")) +
  geom_errorbar(aes(x=type, ymin=death_p10, ymax=death_p90), size=5, width=0.0001) +
  geom_errorbar(aes(x=type, ymin=death_p025, ymax=death_p975), width=0.0001) +
  coord_flip() +
   theme_classic() +
  theme(text = element_text(size=18)) + labs(x="") +
  geom_hline(aes(yintercept=filter(crf_uncertainty, scenario=="ssp370")$death_mean), linetype="dashed", colour="red")
ggsave("uncertainty_combined.pdf", width=7, height = 4)
