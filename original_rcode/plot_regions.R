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
library(exactextractr)

rm(list=ls())
gc()

code_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/rcode"
data_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data"
figure_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/figures"
result_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/results/fire_model"

setwd(figure_path)
states_map <- map_data("state")

db_path <- "/Users/mhqiu/BurkeLab Dropbox/Projects/smoke-climate"

#-------------------------------------------------------------------------------
# Analyze the evaluation performance of simple linear models
# Written by Minghao
# Last edited July 2022
#-------------------------------------------------------------------------------
gfed_grid <- readRDS(paste0(data_path, "/GFED4S_NorthAmerica_025deg_grid_EPSG4326.rds"))

### create crosswalk between gfed grid and ecoregions
#ecoregion <- readRDS(paste0(data_path, "/ecoregion/NA_Eco_level3_MULTIPOLYGON_WGS84.rds"))
#gfed_ecoregion_link <- exact_extract(gfed_grid, ecoregion)

gfed_ecoregion_link <- readRDS(paste0(db_path, "/crosswalk_GFED_NA_CEC_Eco_Level2_3.rds")) %>%
  mutate(intersection_area=as.numeric(intersection_area)) %>% group_by(gfed_cell_id)  %>%
  arrange(desc(intersection_area)) %>%
  mutate(area_weight=intersection_area/sum(intersection_area)) %>% ungroup() # %>% distinct(gfed_cell_id, .keep_all=T)

# gfed_wang_region <- readRDS(paste0(data_path, "/crosswalk_GFED_wang_regions.rds")) %>%
#   mutate(region=as.character(region)) %>%
#   mutate(region=replace(region, region%in%c("mediterranean california", "southwestern US"), "western forest area"),
#          region=replace(region, (region=="northeastern US" & gfed_cell_id<30000), "southeastern US")) %>%
#   mutate(region=replace(region, region=="western forest area", "Western US"),
#          region=replace(region, region=="southeastern US", "Southeastern US"),
#          region=replace(region, region=="northeastern US", "Northeastern US")) %>%
#   mutate(region=replace(region, gfed_cell_id == 27036, "Southeastern US"))
# saveRDS(gfed_wang_region, paste0(data_path, "/crosswalk_GFED_wang_regions.rds"))

gfed_wang_region <- readRDS(paste0(data_path, "/crosswalk_GFED_wang_regions.rds")) 
gfed_grid <- left_join(gfed_grid, gfed_wang_region, by="gfed_cell_id")
gfed_grid <- left_join(gfed_grid, gfed_ecoregion_link, by="gfed_cell_id")

gfed_data <- readRDS(file.path(db_path, "climate-fire_training_data_2001-2021_clean_varname.rds"))
gfed_data <- left_join(gfed_data, gfed_wang_region, by="gfed_cell_id")
gfed_data <- left_join(gfed_data, gfed_ecoregion_link, by="gfed_cell_id")

ggplot()+
  geom_sf(data=filter(gfed_grid, gfed_cell_id %in% unique(gfed_data$gfed_cell_id)), 
          aes(fill=region, colour=region)) +
  #geom_polygon(data=states_map,aes(long, lat, group = group),fill=NA,color ="black") +
  #geom_sf(data=ecoregion, fill=NA) +
  guides(fill=F, colour=F) +
  theme_classic() +
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave("regions_fire_model.pdf", width=5, height=3)


ggplot(gfed_grid, aes(fill=NA_L3KEY, colour=NA_L3KEY))+
  geom_sf(size=0.001) +
  guides(fill=F, colour=F)

