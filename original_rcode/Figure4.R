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
# Plot the figure 4
# Written by Minghao
# Last edited Jan 2024
#-------------------------------------------------------------------------------
### Fig4A: response function
coef <- readRDS(paste0(health_path, "/coef_poisson_bins.rds")) %>% filter(age_group!="under_5")

coef_plot <- coef %>% filter(panel=="year", model=="poisson bins") %>%
  group_by(age_group, bins) %>%
  summarise(median=quantile(exp(coef)-1,0.5),
            p025=quantile(exp(coef)-1,0.025),
            p975=quantile(exp(coef)-1,0.975)) %>% ungroup()

coef_plot$age_group <- factor(coef_plot$age_group, levels= c("65_and_up", "under_65", "all_ages"))

ggplot(coef_plot %>% filter(age_group=="all_ages"),
       aes(x=bins, y=median)) +
  geom_point(size=3)+
  geom_errorbar(aes(ymin=p025,
                    ymax=p975), width=0.001, linewidth=1) +
  facet_wrap(~age_group, ncol=1, scales="free") +
  geom_hline(aes(yintercept=0), linetype="dashed") +
  theme_classic() +
  theme(text = element_text(size=22)) +
  scale_y_continuous(labels = scales::percent, breaks=seq(0,0.08,0.02)) +
  labs(x="", y="mortalities rate change\nper 1ug smoke PM2.5")
ggsave("poisson_year_bins_boot.pdf", width=6.5, height=3.5)

ggplot(coef_plot,
       aes(x=bins, y=median, colour=age_group)) +
  geom_point(size=3, position = position_dodge(0.45))+
  geom_errorbar(aes(ymin=p025,
                    ymax=p975), width=0.001, linewidth=1,
                position = position_dodge(0.45)) +
  #facet_wrap(~age_group, ncol=1, scales="free") +
  geom_hline(aes(yintercept=0), linetype="dashed") +
  theme_classic() +
  theme(text = element_text(size=22)) +
  scale_y_continuous(labels = scales::percent, breaks=seq(0,0.12,0.03)) +
  labs(x="", y="mortalities rate change\nper 1ug smoke PM2.5") +
  scale_colour_manual(values=c("salmon","black","dodgerblue3"))
ggsave("poisson_year_bins_boot_age.pdf", width=8.5, height=4)

### Fig4B: mortality summ by different scenario with 95 CI%
# summ_death <- readRDS(paste0(health_path, "/summ_smoke_death_smokePM_CRF.rds"))
# 
# ggplot(summ_death %>% filter(age_group=="all_ages", CRF%in%c("poisson bins")),
#        aes(x=scenario, y=death_mean, fill=scenario)) +
#   geom_bar(stat = "identity", position = position_dodge()) +
#  # geom_errorbar(aes(ymin=death_p025, ymax=death_p975), stat = "identity", width=0.001, linewidth=0.8) +
#   scale_fill_manual(values=col_scenario) +
#   theme_classic() +
#   theme(text = element_text(size=20)) +
#   labs(x="", y = "Annual mortalities\ndue to smoke PM2.5") +
#   theme(text = element_text(size=20)) +
#   facet_wrap(~age_group, ncol=1, scale="free_y")
# ggsave("smoke_deaths_poisson_bins.pdf", width=7.5, height=3.5)

### Fig4B: mortality summ by different smoke conditions
county_death <- readRDS(paste0(health_path, "/county_smoke_death_smokePM_CRF_bootmean.rds"))

death_summ <- county_death %>% 
  filter(age_group=="all_ages", CRF=="poisson bins") %>%
  group_by(scenario, year, gcm) %>%
  summarise(death_smoke=sum(death_smoke)) %>% ungroup() %>%
  group_by(scenario, gcm) %>%
  summarise(death_smoke=mean(death_smoke)) %>% ungroup() %>%
  group_by(scenario) %>%
  summarise(death_smoke_mean=mean(death_smoke),
            death_smoke_median=median(death_smoke)) %>% ungroup()


cutoff <-  c(-0.01,0.1, 0.5, 1, 2, 5, 20)
death_conc_bin <- county_death %>% 
  mutate(smokePM=replace(smokePM, smokePM>20, 20)) %>%
  mutate(bins=cut(smokePM, cutoff)) %>%
  group_by(year, scenario, gcm, bins, age_group, CRF) %>%
  summarise(death_smoke=sum(death_smoke)) %>%
  ungroup() %>%
  mutate(n_gcm_yr=if_else(scenario=="Observation",10,10*28)) %>% ## assign #of years*gcms in each row to perform the correct avg
  group_by(scenario, bins, age_group, CRF) %>%
  summarise(death_smoke=sum(death_smoke/n_gcm_yr),
            death_smoke_median=median(death_smoke/n_gcm_yr)) %>% ungroup()

death_conc_bin %>%
  filter(age_group=="all_ages", CRF=="poisson bins") %>%
  group_by(scenario) %>% summarise(death_smoke=sum(death_smoke_median)) %>% ungroup()

col_smoke <- c("#c3f4f6","#fff179","#f5b642","orangered3","#551f00")
death_conc_bin$bins <- as.character(death_conc_bin$bins)
death_conc_bin$bins <- factor(death_conc_bin$bins,level=rev(c("(0.1,0.5]", "(0.5,1]", "(1,2]","(2,5]","(5,20]")))

ggplot(death_conc_bin %>% filter(bins!="(-0.01,0.1]", age_group=="all_ages", CRF=="poisson bins"),
       aes(x=scenario)) +
  geom_bar(aes(y=death_smoke, fill=bins, colour=bins), stat = "identity") +
  # geom_errorbar(data=summ_death %>% filter(age_group=="all_ages", CRF%in%c("poisson bins")),
  #               aes(ymin=death_p025, ymax=death_p975), width=0.1, linewidth=1) +
  scale_fill_manual(values=rev(col_smoke)) + 
  scale_colour_manual(values=rev(col_smoke)) + 
  theme_classic() +
  theme(text = element_text(size=20)) +
  labs(x="", y = "Annual mortalities\ndue to smoke PM2.5") 
ggsave("deaths_poisson_bins_smoke_conc.pdf", width=7, height=3.5)


###### Fig 4C: contributions from different bins vs # counties
cutoff <-  c(-0.01,0.1, 0.25,0.5, 0.75,1, 2, 3, 4,5,6, 20)
county_num <-  county_death %>% 
  filter(age_group=="all_ages", CRF=="poisson bins") %>%
  mutate(smokePM=replace(smokePM, smokePM>20, 20)) %>%
  mutate(bins=cut(smokePM, cutoff)) %>%
  group_by(scenario, gcm, bins) %>%
  summarise(ncounty=length(bins)) %>%
  ungroup() %>%
  mutate(n_gcm=if_else(scenario=="Observation",1,28)) %>% ## assign #of years*gcms in each row to perform the correct avg
  group_by(scenario, bins) %>%
  summarise(ncounty=sum(ncounty/n_gcm)) %>%
  ungroup() %>% pivot_wider(names_from = "scenario", values_from = "ncounty")

county_num_percent <- county_num
county_num_percent[,2:5] <- county_num_percent[,2:5]/30850

county_num_summ <- left_join(county_num, county_num_percent, by="bins")
county_num_summ <- county_num_summ[,c(1,2,6,3,7,4,8,5,9)]
write.xlsx(county_num_summ, paste0(health_path, "/county_smoke_conc_bins_num.xlsx"))

death_bin <- county_death %>% 
  filter(age_group=="all_ages", CRF=="poisson bins") %>%
  mutate(smokePM=replace(smokePM, smokePM>20, 20)) %>%
  mutate(bins=cut(smokePM, cutoff)) %>%
  group_by(year, scenario, gcm, bins) %>%
  summarise(death_smoke=sum(death_smoke)) %>%
  ungroup() %>%
  mutate(n_gcm_yr=if_else(scenario=="Observation",10,10*28)) %>% ## assign #of years*gcms in each row to perform the correct avg
  group_by(scenario, bins) %>%
  summarise(death_smoke=sum(death_smoke/n_gcm_yr)) %>% ungroup() %>% 
  pivot_wider(names_from = "scenario", values_from = "death_smoke")

death_percent_bin <- county_death %>% 
  filter(age_group=="all_ages", CRF=="poisson bins") %>%
  mutate(smokePM=replace(smokePM, smokePM>20, 20)) %>%
  mutate(bins=cut(smokePM, cutoff)) %>%
  group_by(year, scenario, gcm, bins) %>%
  summarise(death_smoke=sum(death_smoke)) %>%
  ungroup() %>%
  mutate(n_gcm_yr=if_else(scenario=="Observation",10,10*28)) %>% ## assign #of years*gcms in each row to perform the correct avg
  group_by(scenario, bins) %>%
  summarise(death_smoke=sum(death_smoke/n_gcm_yr)) %>% ungroup() %>% 
  group_by(scenario) %>%
  mutate(death_percent=death_smoke/sum(death_smoke)) %>% select(-death_smoke) %>%
  pivot_wider(names_from = "scenario", values_from = "death_percent")
death_bin <- left_join(death_bin, death_percent_bin , by="bins")
death_bin <- death_bin[,c(1,2,6,3,7,4,8,5,9)]
write.xlsx(death_bin, paste0(health_path, "/death_percent_smoke_conc_bins.xlsx"))

#### combine county percent with death percent
county_percent <- county_num_percent %>% 
  pivot_longer(cols=Observation:ssp370, names_to = "scenario", values_to = "county_percent")

death_percent%>%
  left_join(death_percent_bin %>% pivot_longer(cols=Observation:ssp370, names_to = "scenario", values_to = "death_percent"),
            by=c("bins","scenario")) 
# %>%
#   group_by(scenario) %>%
#   mutate(county_percent_sum=cumsum(county_percent),
#          death_percent_sum=cumsum(death_percent)) %>% ungroup()

ggplot(county_percent %>% filter(bins!="(-0.01,0.1]")) +
  geom_bar(stat="identity", aes(x=bins, y=county_percent, fill=scenario),
           position = position_dodge()) +
  scale_fill_manual(values=col_scenario) +
  theme_classic() +
  theme(text = element_text(size=25)) +
  labs(x="", y="") +
  scale_y_continuous(labels = scales::percent) 
ggsave("county_num_percent.pdf", width=11, height=2.5)



## Fig4c: mortality maps
county_death <- readRDS(paste0(health_path, "/county_smoke_death_smokePM_CRF_bootmean.rds"))

county_death_summ <- county_death %>% 
  filter(age_group=="all_ages", CRF=="poisson bins") %>%
  group_by(fips, scenario) %>%
  summarise(death_smoke=mean(death_smoke)) %>% ungroup() %>%
  pivot_wider(names_from=scenario, values_from=death_smoke) %>%
  mutate(delta_death=ssp245-Observation)

library(tigris)
county_sf <- counties() %>% mutate(fips=as.numeric(GEOID))

county_mortality_sf <- right_join(county_sf, county_death_summ) %>%
  filter(!is.na(delta_death)) 

county_mortality_sf <- county_mortality_sf %>% 
  mutate(delta_death_cut=cut(delta_death, c(1, 10, 50, 100, 200,2000))) #mutate(pm25_delta=replace(pm25_delta, pm25_delta>0.1, 0.1))

#col_map <- c("#c3f4f6","steelblue1","#f5b642","darkorange","#551f00")
col_map <- c("#ECD9B5", rev(met.brewer("Greek", type = "continuous", n=4)))

ggplot() + 
  geom_sf(data=county_mortality_sf, 
          aes(fill=delta_death_cut, colour=delta_death_cut),linewidth=0.001) +
  geom_polygon(data=states_map,aes(long, lat, group = group),fill=NA,color ="black", linewidth=0.3) +
  scale_fill_manual(values=col_map, na.value = "white") +
  scale_colour_manual(values=col_map, na.value = "white") +
  coord_sf(xlim = c(-125,-65), ylim=c(25,51)) + 
  labs(x="",y="",) + 
  theme_void() +
  theme(axis.ticks=element_blank(), axis.text=element_blank(),
        text = element_text(size=18)) +
  labs(fill="Deaths (SSP2-4.5 relative to 2011-2020)") +
  guides(colour="none")
ggsave("delta_death_ssp245_obs.png", width = 8, height = 4)

## Fig4d: comparing with temperature deaths
temp_death_county <- readRDS(paste0(data_path, "/Carleton_QJE_data/county_temp_deaths_rcp45.rds"))
temp_death_state <- temp_death_county %>% group_by(state_abbrev, GEOID_state) %>%
  summarise(death_2040_2059=sum(death_2040_2059)) %>% ungroup()

county_death_summ <- readRDS(paste0(health_path, "/county_smoke_death_smokePM_CRF_bootmean.rds")) %>% 
  filter(age_group=="all_ages", CRF=="poisson bins") %>%
  group_by(fips, scenario) %>%
  summarise(death_smoke=mean(death_smoke)) %>% ungroup() 

county_df <- tigris::counties() %>% as.data.frame() %>%
  distinct(STATEFP, GEOID) %>%
  rename(fips=GEOID, GEOID_state=STATEFP) %>%
  mutate(fips=as.numeric(fips), GEOID_state=as.numeric(GEOID_state))

state_death_summ <- left_join(county_death_summ, county_df) %>%
  group_by(GEOID_state, scenario) %>%
  summarise(death_smoke=sum(death_smoke)) %>%
  pivot_wider(names_from=scenario, values_from=death_smoke) %>%
  mutate(delta_death_smoke=ssp245-Observation)

temp_smoke_death_state <- left_join(state_death_summ, temp_death_state) %>% ungroup() %>%
  select(state_abbrev, delta_death_smoke,  death_2040_2059) %>%
  rename(death_smoke=delta_death_smoke,
         death_temp=death_2040_2059) %>% arrange(desc(death_smoke))

state_df <- tigris::states() %>%
  as.data.frame() %>% rename(state_abbrev=NAME) %>% select(STUSPS, state_abbrev)

temp_smoke_death_state <- left_join(temp_smoke_death_state, state_df)
write.csv(temp_smoke_death_state, paste0(health_path, "/temp_smoke_death_state.csv"))

temp_smoke_death_state <- read.csv(paste0(health_path, "/temp_smoke_death_state.csv"))[1:49,]
ggplot(temp_smoke_death_state %>% filter(death_smoke>75), 
       aes(x=death_smoke, y=death_temp)) +
  annotate("rect", xmin = 55, xmax = 6500, ymin = 0, ymax = 12000, alpha = .5,fill = "salmon") +
  annotate("rect", xmin = 55, xmax = 6500, ymin = -7000, ymax = 0, alpha = .5,fill = "steelblue") +
  geom_point(size=3.5) +
  geom_text_repel(aes(label=STUSPS), size=6, box.padding=0.32) +
  scale_x_continuous(trans='log2', breaks=c(25,100,500,3000)) +
  geom_hline(yintercept = 0,linetype="dashed", size=1) +
  theme_classic() +
  theme(text = element_text(size=25)) +
  coord_cartesian(xlim = c(55,6000), ylim=c(-7000,12000), expand = F)
ggsave("state_death_temp_smoke.pdf", width=6.5, height=5)  

### create a bar chart 
total_death <- data.frame(type=c("smoke", "NA","heat","cold","net"),
                          death=c(sum(temp_smoke_death_state$death_smoke),0,
                                  sum(temp_smoke_death_state$death_temp*as.numeric(temp_smoke_death_state$death_temp>0)),
                                  sum(temp_smoke_death_state$death_temp*as.numeric(temp_smoke_death_state$death_temp<0)),
                                  sum(temp_smoke_death_state$death_temp)))
total_death$type <- factor(total_death$type, levels = c("smoke","NA","heat","cold","net"))

ggplot() +
  geom_bar(data=total_death %>% filter(type%in%c("smoke","net")), 
           stat = "identity",
           aes(x=NA,y=death, fill=type), alpha=0.9, colour="black",
           position=position_dodge()) +
  theme_classic() +
  theme(text = element_text(size=25),
        axis.line.x = element_blank(), 
        axis.text.x = element_blank()) +
  #scale_fill_manual(values=c("#551f00","steelblue")) +
  scale_fill_manual(values=c("grey90","grey90")) +
  scale_y_continuous(breaks=seq(-16000, 12000, by=4000)) +
  geom_hline(aes(yintercept=0), linetype="dashed", linewidth=1)
ggsave("bar_temp_smoke.pdf", width=5, height=5)





### Fig SX: across different CRFs
summ_death <- readRDS(paste0(health_path, "/summ_smoke_death_smokePM_CRF.rds"))
summ_death$CRF <- factor(summ_death$CRF, levels=
                           c("linear","linear bins", "linear bins (coarse)",
                             "poisson","poisson bins", "poisson bins (coarse)",
                             "quadratic", "cubic"))

ggplot(summ_death %>% filter(!age_group%in%c("under_5", "under_65", "65_and_up"), CRF!="cubic"),
       aes(x=scenario, y=death_median, fill=CRF)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  geom_errorbar(aes(ymin=death_p025, ymax=death_p975), stat = "identity",
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
ggsave("smoke_deaths_crf_95ci_all_age.pdf", width=10, height=4)

#### Comapring across different CRFs

county_death <- readRDS(paste0(health_path, "/county_smoke_death_smokePM_CRF_bootmean.rds"))
cutoff <-  c(-0.01,0.1, 0.5, 1, 2, 5, 20)
death_bin <- county_death %>% 
  filter(age_group=="all_ages", CRF%in%c("poisson","poisson bins","quadratic")) %>%
  mutate(smokePM=replace(smokePM, smokePM>20, 20)) %>%
  mutate(bins=cut(smokePM, cutoff)) %>%
  group_by(year, scenario, gcm, bins, CRF) %>%
  summarise(death_smoke=sum(death_smoke)) %>%
  ungroup() %>%
  mutate(n_gcm_yr=if_else(scenario=="Observation",10,10*28)) %>% ## assign #of years*gcms in each row to perform the correct avg
  group_by(scenario, bins, CRF) %>%
  summarise(death_smoke=sum(death_smoke/n_gcm_yr)) %>% ungroup() %>% 
  pivot_wider(names_from = "CRF", values_from = "death_smoke")
