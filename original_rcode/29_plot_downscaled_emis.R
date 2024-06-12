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

rm(list=ls())
gc()

code_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/rcode"
data_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data"
figure_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/figures"
result_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/results/fire_model"

db_path <- "/Users/mhqiu/BurkeLab Dropbox/Projects/smoke-climate"

setwd(figure_path)
states_map <- map_data("state")

color_map <- c("grey65",
               rev(met.brewer(name="Homer1", n=8, type="discrete"))[c(2,6,8)])

#-------------------------------------------------------------------------------
# Plot the projected downscale emissions
# Written by Minghao
# Last edited October 2022
#-------------------------------------------------------------------------------
## load historical gfed emissions
gfed_grid <- readRDS(paste0(data_path, "/GFED4S_NorthAmerica_025deg_grid_EPSG4326.rds"))
gfed_wang_region <- readRDS(paste0(data_path, "/crosswalk_GFED_wang_regions.rds")) 
gfed_grid <- left_join(gfed_grid, gfed_wang_region) %>% filter(!is.na(region))
  
fire_emis_grid <- readRDS(paste0(result_path,"/CMIP6_fire_proj/downscaled_emis_final_best_models.rds"))

fire_emis_summ1 <- fire_emis_grid %>%
  group_by(year, scenario, gcm, region, algorithm_outcome) %>%
  summarise(pred_DM=sum(pred_DM_grid, na.rm=T)) %>%
  group_by(year, scenario, gcm, region) %>%
  summarise(pred_DM_grid=mean(pred_DM, na.rm=T)) %>% ## mean across stat algorithms
  ungroup() %>%
  group_by(year, scenario, region) %>%
  summarise(DM_grid_mean = mean(pred_DM_grid, na.rm=T),
            DM_grid_median = median(pred_DM_grid, na.rm=T)) %>%  ## mean/median across gcms
  ungroup()

fire_emis_summ <- fire_emis_grid %>%
  group_by(year, scenario, gcm, region, gfed_cell_id, algorithm_outcome) %>%
  summarise(pred_DM_grid=sum(pred_DM_grid, na.rm=T)) %>% ## sum over grid cells that span multiple subregions
  group_by(year, scenario, gcm, region, gfed_cell_id) %>%
  summarise(pred_DM_grid=mean(pred_DM_grid, na.rm=T)) %>% ## mean across stat algorithms
  ungroup() %>%
  group_by(year, scenario, region, gfed_cell_id) %>%
  summarise(DM_grid_mean = mean(pred_DM_grid, na.rm=T),
         DM_grid_median = median(pred_DM_grid, na.rm=T)) %>%  ## mean/median across gcms
  ungroup() %>%
  group_by(year, scenario, gfed_cell_id) %>%
  summarise(DM_grid_mean = sum(DM_grid_mean, na.rm=T),
            DM_grid_median = sum(DM_grid_median, na.rm=T)) %>% ## sum over grid cells that span multiple regions
  ungroup()
saveRDS(fire_emis_summ, paste0(result_path,"/CMIP6_fire_proj/downscaled_emis_final_best_models_agg.rds"))


#### plot the projected and historical emissions over map ####
fire_emis_summ <- readRDS(paste0(result_path,"/CMIP6_fire_proj/downscaled_emis_final_best_models_agg.rds"))
gfed_data <- readRDS(file.path(db_path, "climate-fire_training_data_2001-2021_clean_varname.rds")) %>%
  group_by(year, gfed_cell_id) %>%
  summarise(DM_kg=sum(DM_kg, na.rm=T)) %>%
  group_by(gfed_cell_id) %>%
  summarise(DM_2001_2021 = mean(DM_kg, na.rm=T)) %>% ungroup() 

plot_var <- "DM_grid_median" ## select the variable for plotting in the CMIP6 runs
gfed_data <- gfed_data %>%
  mutate(scenario="Observation",
         year="2001-2021") %>%
  mutate(DM_kton=DM_2001_2021/1e6)

fire_emis_summ <- fire_emis_summ %>%
  mutate(year=as.character(year)) %>%
  mutate(DM_kton= !!as.name(plot_var)/1e6) 

fire_emis_plot <- bind_rows(fire_emis_summ, gfed_data)
emis_grid <-  right_join(gfed_grid, fire_emis_plot) %>% filter(!is.na(scenario))

emis_grid %>% filter(year %in% c("2001-2021", "2055")) %>%
  group_by(region, year, scenario) %>% 
  summarise(DM_summ=sum(DM_kton)/1e3) %>%
  ungroup()

### just for plotting to top-code the large values
col_disc <- c("grey70","#c3f4f6", "steelblue3","green4",
              "#f5b642","orangered3","#551f00")

#"#551f00" "#a62f00" "#df7700" "#f5b642" "#fff179" "#c3f4f6" "#6ad5e8" "#32b2da"
### for the entire North America
emis_grid_tmp <- emis_grid 
max_emis <- max(emis_grid_tmp$DM_kton, na.rm=T)
label_emis <- c(0,1,5,10,100,500,1000, max_emis)
emis_grid_tmp$DM_cut <- cut(emis_grid_tmp$DM_kton, label_emis)
emis_grid_tmp$DM_cut <- factor(emis_grid_tmp$DM_cut, levels= rev(levels(emis_grid_tmp$DM_cut)))
#levels(emis_grid_tmp$DM_cut) <- rev(levels(emis_grid_tmp$DM_cut)) 

ggplot(emis_grid_tmp %>% filter(year %in% c("2001-2021", "2055")), 
       aes(colour=DM_cut, fill=DM_cut)) +
  geom_sf() + 
  # scale_color_stepsn(colors = col_disc, limits=c(0, max_emis), breaks = label_emis, trans="pseudo_log") +
  # scale_fill_stepsn(colors = col_disc, limits=c(0, max_emis), breaks = label_emis,trans="pseudo_log") +
  scale_color_manual(values=rev(col_disc)) +
  scale_fill_manual(values=rev(col_disc)) +
  theme_classic() + theme(text = element_text(size=14),
                          axis.ticks = element_blank(),
                          axis.line = element_blank(),
                          axis.text = element_blank()) +
  labs(colour="DM (1000 tons)", fill="DM (1000 tons)") +
  facet_wrap(~year+scenario)
ggsave("proj_2055_median_DM_map.png", width=12, height=8)


##### plot by subregions 
emis_grid_tmp <- filter(emis_grid, region=="Western US")
max_emis <- max(emis_grid_tmp$DM_kton, na.rm=T)
label_emis <- c(0,1,5,10,100,250,1000,max_emis)
emis_grid_tmp$DM_cut <- cut(emis_grid_tmp$DM_kton, label_emis)
emis_grid_tmp$DM_cut <- factor(emis_grid_tmp$DM_cut, levels= rev(levels(emis_grid_tmp$DM_cut)))

ggplot(emis_grid_tmp %>% filter(year %in% c("2001-2021", "2055")), 
       aes(colour=DM_cut, fill=DM_cut)) +
  geom_sf() + 
  # scale_color_stepsn(colors = col_disc, limits=c(0, max_emis), breaks = label_emis, trans="pseudo_log") +
  # scale_fill_stepsn(colors = col_disc, limits=c(0, max_emis), breaks = label_emis,trans="pseudo_log") +
  scale_color_manual(values=rev(col_disc)) +
  scale_fill_manual(values=rev(col_disc)) +
  theme_classic() + theme(text = element_text(size=14),
                          axis.ticks = element_blank(),
                          axis.line = element_blank(),
                          axis.text = element_blank()) +
  labs(colour="DM (1000 tons)", fill="DM (1000 tons)") +
  facet_wrap(~year+scenario)
ggsave("proj_2055_median_DM_map_west.png", width=6, height=6)

### Canada
emis_grid_tmp <- filter(emis_grid, region=="Canada-Alaska")
emis_grid_tmp[emis_grid_tmp$DM_kton>3000, "DM_kton"] <- 3000
max_emis <- max(emis_grid_tmp$DM_kton, na.rm=T)
label_emis <- c(0,1,5,10,50,100,500,max_emis)

emis_grid_tmp$DM_cut <- cut(emis_grid_tmp$DM_kton, label_emis)
emis_grid_tmp$DM_cut <- factor(emis_grid_tmp$DM_cut, levels= rev(levels(emis_grid_tmp$DM_cut)))

ggplot(emis_grid_tmp %>% filter(year %in% c("2001-2021", "2055")), 
       aes(colour=DM_cut, fill=DM_cut)) +
  geom_sf() + 
  # scale_color_stepsn(colors = col_disc, limits=c(0, max_emis), breaks = label_emis, trans="pseudo_log") +
  # scale_fill_stepsn(colors = col_disc, limits=c(0, max_emis), breaks = label_emis,trans="pseudo_log") +
  scale_color_manual(values=rev(col_disc)) +
  scale_fill_manual(values=rev(col_disc)) +
  theme_classic() + theme(text = element_text(size=14),
                          axis.ticks = element_blank(),
                          axis.line = element_blank(),
                          axis.text = element_blank()) +
  labs(colour="DM (1000 tons)", fill="DM (1000 tons)") +
  facet_wrap(~year+scenario)
ggsave("proj_2055_median_DM_map_canada.png", width=8.5, height=6)

  
### Mexico
emis_grid_tmp <- filter(emis_grid, region=="Mexico")
#emis_grid_tmp[emis_grid_tmp$DM_kton>3000, "DM_kton"] <- 3000
max_emis <- max(emis_grid_tmp$DM_kton, na.rm=T)
label_emis <- c(0,1,5,10,25,50,100,max_emis)

emis_grid_tmp$DM_cut <- cut(emis_grid_tmp$DM_kton, label_emis)
emis_grid_tmp$DM_cut <- factor(emis_grid_tmp$DM_cut, levels= rev(levels(emis_grid_tmp$DM_cut)))

ggplot(emis_grid_tmp %>% filter(year %in% c("2001-2021", "2055")), 
       aes(colour=DM_cut, fill=DM_cut)) +
  geom_sf() + 
  # scale_color_stepsn(colors = col_disc, limits=c(0, max_emis), breaks = label_emis, trans="pseudo_log") +
  # scale_fill_stepsn(colors = col_disc, limits=c(0, max_emis), breaks = label_emis,trans="pseudo_log") +
  scale_color_manual(values=rev(col_disc)) +
  scale_fill_manual(values=rev(col_disc)) +
  theme_classic() + theme(text = element_text(size=14),
                          axis.ticks = element_blank(),
                          axis.line = element_blank(),
                          axis.text = element_blank()) +
  labs(colour="DM (1000 tons)", fill="DM (1000 tons)") +
  facet_wrap(~year+scenario)
ggsave("proj_2055_median_DM_map_mexico.png", width=7, height=6)


### Eastern US
emis_grid_tmp <- filter(emis_grid, region%in%c("Southeastern US", "Northeastern US"))
#emis_grid_tmp[emis_grid_tmp$DM_kton>3000, "DM_kton"] <- 3000
max_emis <- max(emis_grid_tmp$DM_kton, na.rm=T)
label_emis <- c(0,1,5,10,20,40,60, max_emis)

emis_grid_tmp$DM_cut <- cut(emis_grid_tmp$DM_kton, label_emis)
emis_grid_tmp$DM_cut <- factor(emis_grid_tmp$DM_cut, levels= rev(levels(emis_grid_tmp$DM_cut)))

ggplot(emis_grid_tmp %>% filter(year %in% c("2001-2021", "2055")), 
       aes(colour=DM_cut, fill=DM_cut)) +
  geom_sf() + 
  # scale_color_stepsn(colors = col_disc, limits=c(0, max_emis), breaks = label_emis, trans="pseudo_log") +
  # scale_fill_stepsn(colors = col_disc, limits=c(0, max_emis), breaks = label_emis,trans="pseudo_log") +
  scale_color_manual(values=rev(col_disc)) +
  scale_fill_manual(values=rev(col_disc)) +
  theme_classic() + theme(text = element_text(size=14),
                          axis.ticks = element_blank(),
                          axis.line = element_blank(),
                          axis.text = element_blank()) +
  labs(colour="DM (1000 tons)", fill="DM (1000 tons)") +
  facet_wrap(~year+scenario)
ggsave("proj_2055_median_DM_map_eastUS.png", width=7, height=6)



