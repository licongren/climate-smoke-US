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
# Calculate the mortality using smoke related effects and 
# Written by Minghao
# Last edited May 2024
#-------------------------------------------------------------------------------

#--------------------------------------------------------------
## Compare binned CRF w/ Pope et al CRF
#--------------------------------------------------------------
smoke_conc_vec <- seq(0,7, by=0.01)
mort_smoke_year <- readRDS(paste0(data_path, "/county_mort_smoke_year_2006_2020_age.rds"))
death_rate <- mort_smoke_year %>% filter(year%in%2006:2019, age_group == "all_ages") %>%
  group_by(year) %>%
  summarise(pop_total=sum(pop), death_total=sum(pop*rate/1e5)) %>%
  ungroup() %>%
  mutate(death_rate=death_total/pop_total*1e5)
death_rate_avg <- mean(death_rate$death_rate)

coef_poisson_bin <- readRDS(paste0(health_path, "/coef_poisson_bins.rds"))  %>%
  filter(panel=="year", model=="poisson bins", age_group=="all_ages") %>%
  group_by(age_group, bins) %>%
  summarise(mean=mean(exp(coef)-1),
            p025=quantile(exp(coef)-1,0.025),
            p975=quantile(exp(coef)-1,0.975)) %>% ungroup() %>%
  mutate(smoke=NA)
coef_poisson_bin <- coef_poisson_bin %>%
  mutate(smoke=replace(smoke, bins=="(0.1,0.25]", 0.175),
         smoke=replace(smoke, bins=="(0.25,0.5]", 0.375),
         smoke=replace(smoke, bins=="(0.5,0.75]", 0.625),
         smoke=replace(smoke, bins=="(0.75,1]", 0.875),
         smoke=replace(smoke, bins=="(1,2]", 1.5),
         smoke=replace(smoke, bins=="(2,3]", 2.5),
         smoke=replace(smoke, bins=="(3,4]", 3.5),
         smoke=replace(smoke, bins=="(4,5]", 4.5),
         smoke=replace(smoke, bins=="(5,6]", 5.5),
         smoke=replace(smoke, bins=="(6,10]", 6.5))

pope_response <- data.frame(smoke=smoke_conc_vec) %>%
  mutate(mean=1.08^(smoke/10)-1,
         p025=1.05^(smoke/10)-1,
         p975=1.11^(smoke/10)-1)

ggplot(pope_response, aes(x=smoke, y=mean)) +
  geom_line(aes(colour="red3", fill="red3")) +
  geom_ribbon(aes(ymin=p025, ymax=p975,  fill="red3"), alpha=0.5) +
  geom_point(data=coef_poisson_bin, size=2)+
  geom_errorbar(data=coef_poisson_bin, 
                aes(ymin=p025,ymax=p975), width=0.001, linewidth=1) +
  geom_hline(aes(yintercept=0), linetype="dashed") +
  theme_classic() +
  theme(text = element_text(size=22)) +
  guides(fill=F, colour=F) +
  scale_y_continuous(labels = scales::percent, breaks=seq(0,0.08,0.02)) 
ggsave("response_curve_pope_compare.pdf", width=7, height=3.5)

#--------------------------------------------------------------
## compute the smoke mortalities using total PM CRF
#--------------------------------------------------------------
rr_df <- data.frame(
  study= c("Pope et al. 2020", "Krewski et al. 2009","Di et al. 2017","Lepeule et al. 2012", "Hoek et al. 2013"),
  rr_main=c(1.08, 1.056, 1.073, 1.14, 1.06),
  rr_min=c(1.05, 1.035, 1.071, 1.07, 1.04),
  rr_max=c(1.11, 1.078, 1.075, 1.22, 1.08)) 

smoke_future <- readRDS(paste0(data_path, "/health/US_smoke_county_2046_2055_gcm.rds")) %>%
  filter(age_group=="all_ages")
smoke_obs <- readRDS(paste0(data_path, "/health/US_smoke_county_2011_2020.rds")) %>%
  filter(age_group=="all_ages") %>%
  mutate(gcm="none") ## just to be able to use the same function
smoke_combined <- bind_rows(smoke_future, smoke_obs)

smoke_all <- smoke_combined %>% 
  slice(rep(1:n(), each = nrow(rr_df))) %>%
  mutate(study=rep(c("Pope et al. 2020", "Krewski et al. 2009","Di et al. 2017","Lepeule et al. 2012", "Hoek et al. 2013"), 
                   nrow(smoke_combined)))

totalPM <- readRDS(paste0(data_path, "/health/US_totalPM_county_2011_2020.rds")) %>%
  filter(age_group=="all_ages") %>%
  mutate(nonsmokePM=totalPM - smokePM) %>%
  group_by(fips) %>%
  summarise(nonsmokePM=mean(nonsmokePM)) %>%
  ungroup()
  
mortality_crf <- left_join(smoke_all, rr_df) %>%
  left_join(totalPM) %>%
  mutate(totalPM = smokePM + nonsmokePM,
         death_mean = pop_2050*death_rate_avg* ((rr_main^(totalPM/10)-1) - (rr_main^(nonsmokePM/10)-1)),
         death_p025 = pop_2050*death_rate_avg* ((rr_min^(totalPM/10)-1) - (rr_min^(nonsmokePM/10)-1)),
         death_p975 = pop_2050*death_rate_avg* ((rr_max^(totalPM/10)-1) - (rr_max^(nonsmokePM/10)-1))) %>%
  group_by(year, scenario, gcm, study) %>%
  summarise_at(c("death_mean", "death_p025", "death_p975"), .funs = sum, na.rm=T) %>%
  group_by(scenario, gcm, study) %>%
  summarise_at(c("death_mean", "death_p025", "death_p975"), .funs = mean, na.rm=T) %>% ### take mean across years 
  group_by(scenario, study) %>%
  summarise_at(c("death_mean", "death_p025", "death_p975"), .funs = median, na.rm=T) %>% ### take mean across years 
  ungroup()

smoke_death_summ <-  readRDS(paste0(health_path, "/summ_smoke_death_smokePM_CRF_mediangcm.rds")) %>%
  filter(age_group=="all_ages", CRF=="poisson bins") %>%
  mutate(study="our study") %>% select(-CRF, -age_group)

mortality_crf <- bind_rows(mortality_crf, smoke_death_summ)
mortality_crf$study <- factor(mortality_crf$study, levels=c("our study","Pope et al. 2020", "Krewski et al. 2009","Lepeule et al. 2012", "Hoek et al. 2013","Di et al. 2017"))

mortality_crf <- mortality_crf %>% arrange(scenario) 
write.xlsx(mortality_crf, paste0(health_path, "/smoke_mortality_totalPM_CRF.xlsx"))

ggplot(mortality_crf %>% filter(study!="Lepeule et al. 2012"),  # 
       aes(x=scenario, y=death_mean, fill=study))+
  geom_bar(stat = "identity", position = position_dodge()) +
  # geom_errorbar(aes(ymin=death_p025, ymax=death_p975,group=study), stat = "identity",
  #               position = position_dodge2(.9, padding = .9), linewidth=0.9, width=0.01) +
  scale_fill_manual(values=c("grey40",met.brewer(name="Egypt", n=11, type="continuous")[c(1,5,7,9,11)]))+
  scale_colour_manual(values=c("grey40",met.brewer(name="Egypt", n=11, type="continuous")[c(1,5,7,9,11)]))+
  theme_classic() +
  theme(text = element_text(size=20)) +
  labs(x="", y = "Annual smoke deaths", fill="") +
  theme(text = element_text(size=20)) 
ggsave("SI_CRF_totalPM_no_lep.pdf", width = 9, height=4)

##### just weatrn US
geoid_west <- c(16, 35, 6, 41, 31, 53, 49, 48, 8, 40, 56, 38, 32, 30, 20, 46 ,4)

mortality_west_crf <- left_join(smoke_all, rr_df) %>%
  left_join(totalPM) %>%
  mutate(totalPM = smokePM + nonsmokePM,
         death_mean = pop_2050*death_rate_avg* ((rr_main^(totalPM/10)-1) - (rr_main^(nonsmokePM/10)-1)),
         death_p025 = pop_2050*death_rate_avg* ((rr_min^(totalPM/10)-1) - (rr_min^(nonsmokePM/10)-1)),
         death_p975 = pop_2050*death_rate_avg* ((rr_max^(totalPM/10)-1) - (rr_max^(nonsmokePM/10)-1))) %>%
  mutate(state_id=floor(fips/1000)) %>%
  filter(state_id %in% geoid_west) %>%
  group_by(year, scenario, gcm, study) %>%
  summarise_at(c("death_mean", "death_p025", "death_p975"), .funs = sum, na.rm=T) %>%
  group_by(scenario, gcm, study) %>%
  summarise_at(c("death_mean", "death_p025", "death_p975"), .funs = mean, na.rm=T) %>% ### take mean across years 
  group_by(scenario, study) %>%
  summarise_at(c("death_mean", "death_p025", "death_p975"), .funs = median, na.rm=T) %>% ### take mean across years 
  ungroup()

smoke_death_summ_west <- readRDS(paste0(health_path, "/county_smoke_death_smokePM_CRF_bootmean.rds")) %>%
  filter(age_group=="all_ages", CRF=="poisson bins") %>%
  mutate(state_id=floor(fips/1000)) %>%
  filter(state_id %in% geoid_west) %>%
  group_by(year, scenario, gcm) %>%
  summarise_at(c("death_smoke"), .funs = sum, na.rm=T) %>%
  group_by(scenario, gcm) %>%
  summarise_at(c("death_smoke"), .funs = mean, na.rm=T) %>% ### take mean across years 
  group_by(scenario) %>%
  summarise_at(c("death_smoke"), .funs = median, na.rm=T) %>% ### take mean across years 
  ungroup() %>%
  mutate(study="our study") 

mortality_crf_west <- bind_rows(mortality_crf_west, smoke_death_summ_west)
mortality_crf_west$study <- factor(mortality_crf_west$study, levels=c("our study","Pope et al. 2020", "Krewski et al. 2009","Lepeule et al. 2012", "Hoek et al. 2013","Di et al. 2017"))

mortality_crf_west <- mortality_crf_west %>% arrange(scenario) 

#--------------------------------------------------------------
## Calculate percentage of smoke mortality in the  
#--------------------------------------------------------------
rr_df <- data.frame(
  study= c("Pope et al. 2020", "Krewski et al. 2009","Di et al. 2017","Lepeule et al. 2012", "Hoek et al. 2013"),
  rr_main=c(1.08, 1.056, 1.073, 1.14, 1.06),
  rr_min=c(1.05, 1.035, 1.071, 1.07, 1.04),
  rr_max=c(1.11, 1.078, 1.075, 1.22, 1.08)) 

smoke_future <- readRDS(paste0(data_path, "/health/US_smoke_county_2046_2055_gcm.rds")) %>%
  filter(age_group=="all_ages")
smoke_obs <- readRDS(paste0(data_path, "/health/US_smoke_county_2011_2020.rds")) %>%
  filter(age_group=="all_ages") %>%
  mutate(gcm="none") ## just to be able to use the same function
smoke_combined <- bind_rows(smoke_future, smoke_obs)

smoke_all <- smoke_combined %>% 
  slice(rep(1:n(), each = nrow(rr_df))) %>%
  mutate(study=rep(c("Pope et al. 2020", "Krewski et al. 2009","Di et al. 2017","Lepeule et al. 2012", "Hoek et al. 2013"), 
                   nrow(smoke_combined)))

totalPM <- readRDS(paste0(data_path, "/health/US_totalPM_county_2011_2020.rds")) %>%
  filter(age_group=="all_ages") %>%
  mutate(nonsmokePM=totalPM - smokePM) %>%
  group_by(fips) %>%
  summarise(nonsmokePM=mean(nonsmokePM)) %>%
  ungroup()

mortality_crf <- left_join(smoke_all, rr_df) %>%
  left_join(totalPM) %>%
  mutate(totalPM = smokePM + nonsmokePM,
         death_mean = pop_2050*death_rate_avg* ((rr_main^(totalPM/10)-1)),
         death_p025 = pop_2050*death_rate_avg* ((rr_min^(totalPM/10)-1)),
         death_p975 = pop_2050*death_rate_avg* ((rr_max^(totalPM/10)-1))) %>%
  group_by(year, scenario, gcm, study) %>%
  summarise_at(c("death_mean", "death_p025", "death_p975","totalPM"), .funs = sum, na.rm=T) %>%
  group_by(scenario, gcm, study) %>%
  summarise_at(c("death_mean", "death_p025", "death_p975","totalPM"), .funs = mean, na.rm=T) %>% ### take mean across years 
  group_by(scenario, study) %>%
  summarise_at(c("death_mean", "death_p025", "death_p975","totalPM"), .funs = median, na.rm=T) %>% ### take mean across years 
  ungroup()

mortality_summ <- mortality_crf %>% filter(study!="Lepeule et al. 2012") %>%
  group_by(scenario) %>%
  summarise(totalPM_death=mean(death_mean)) %>%
  ungroup()

smoke_death_summ <-  readRDS(paste0(health_path, "/summ_smoke_death_smokePM_CRF_mediangcm.rds")) %>%
  filter(age_group=="all_ages", CRF=="poisson bins") %>%
  mutate(study="our study") %>% select(-CRF, -age_group)

smoke_death_percent <- left_join(smoke_death_summ, mortality_summ) %>%
  mutate(smoke_percent=death_median/totalPM_death)

