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

rm(list=ls())
gc()

code_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/rcode"
data_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data"
figure_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/figures"
result_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/results/fire_model"

setwd(figure_path)
states_map <- map_data("state")

db_path <- "/Users/mhqiu/BurkeLab Dropbox/Projects/smoke-climate"
db_proj_path <- "/Users/mhqiu/BurkeLab Dropbox/Projects"
#-------------------------------------------------------------------------------
# Load the WUSTL total PM data
# Written by Minghao
# Last edited Oct 2022
#-------------------------------------------------------------------------------

# read the smoke grid shape file
smoke_pm <- readRDS(paste0(db_proj_path, "/smoke_PM_prediction/output/smokePM_predictions_20060101_20201231.rds"))
smokepm_annual <- smoke_pm %>%
  mutate(year=as.numeric(substr(date,1,4))) %>%
  group_by(grid_id_10km, year) %>%
  summarise(smokePM_pred=sum(smokePM_pred)) %>% ungroup() %>%
  mutate(n=if_else(leap_year(year),366,365),
         smokePM_mean=smokePM_pred/n) 
  
grid_10km <- st_read(paste0(db_proj_path, "/smoke_PM_prediction/data/1_grids/grid_10km_wgs84/grid_10km_wgs84.shp"))

file_list <- list.files(paste0(data_path, "/totalPM"), pattern = "*.nc", full.names = T)

### create a corssswalk between the smoke PM grid and the totalPM grid
file_sample <- nc_open(file_list[1])
lon <- ncvar_get(file_sample, "lon")
lat <- ncvar_get(file_sample, "lat")
data_sample <- ncvar_get(file_sample, "GWRPM25")
var_raster <- rast(t(data_sample[,length(lat):1]), crs="WGS84")
ext(var_raster) <- c(range(lon)[1]-0.05,range(lon)[2]+0.05,range(lat)[1]-0.05,range(lat)[2]+0.05)

crosswalk_pm_10kmgrid <- exact_extract(var_raster, grid_10km,
                                       coverage_area = T,
                                       include_cell = T,
                                       include_cols="ID",
                                       include_xy = T) %>%
  bind_rows() %>% select(-value, -cell)
saveRDS(crosswalk_pm_10kmgrid, paste0(data_path, "/totalPM/crosswalk_totalPM_smoke_10kmgrid.rds"))

crosswalk_pm_10kmgrid <- readRDS(paste0(data_path, "/totalPM/crosswalk_totalPM_smoke_10kmgrid.rds"))

# ggplot() +
#   geom_sf(data= grid_10km) +
#   geom_raster(data=fortify(var_raster))

pm_combined <- NULL
for (ii in file_list){
year_tmp <- str_extract(ii, "\\d{4}") %>% as.numeric()

 file_tmp <- nc_open(ii)
 lon <- ncvar_get(file_tmp, "lon")
 lat <- ncvar_get(file_tmp, "lat")
 data_tmp <- ncvar_get(file_tmp, "GWRPM25")
 var_raster <- rast(t(data_tmp[,length(lat):1]))
 ext(var_raster) <- c(range(lon)[1]-0.05,range(lon)[2]+0.05,range(lat)[1]-0.05,range(lat)[2]+0.05)
 
 pm_df <- as.data.frame(var_raster, xy=T)
 colnames(pm_df)[3] <- "totalPM" 
 
 pm_df <- left_join(pm_df, crosswalk_pm_10kmgrid, by=c("x","y")) %>%
   filter(!is.na(ID)) %>%
   group_by(ID) %>%
   summarise(totalPM = weighted.mean(totalPM, w=coverage_area)) %>%
   mutate(year = year_tmp) %>%
   rename(grid_id_10km = ID) %>%
   ungroup()
 
 pm_combined <- bind_rows(pm_combined, pm_df)
}

total_smoke_pm <- left_join(pm_combined, smokepm_annual[,c("grid_id_10km", "year", "smokePM_mean")]) %>%
  filter(year %in% 2006:2020) %>%
  mutate(smokePM_mean=replace(smokePM_mean, is.na(smokePM_mean), 0))
saveRDS(total_smoke_pm, paste0(data_path, "/totalPM/totalPM_smoke_combined_2006_2020.rds"))

# total_smoke_pm  %>% filter(smokePM_mean > totalPM) ### 631 rows out of 489800 have totalPM < smokePM
# total_smoke_pm  %>% mutate(nonsmokePM = totalPM - smokePM_mean) %>%
#   filter(nonsmokePM  < 0) %>% arrange(nonsmokePM) 


######### calculate and plot the smoke contribution to annual mean PM
total_smoke_pm <- readRDS(paste0(data_path, "/totalPM/totalPM_smoke_combined_2016_2020.rds")) 

smoke_percent <- total_smoke_pm  %>%
  #filter(year<2020) %>%
  group_by(grid_id_10km ) %>%
  summarise(smokePM=mean(smokePM_mean, na.rm=T),
            totalPM=mean(totalPM, na.rm=T)) %>%
  mutate(smoke_percent = smokePM/totalPM,
         nonsmokePM = totalPM - smokePM)

smoke_percent  <- right_join(grid_10km, smoke_percent, by=c("ID" = "grid_id_10km"))

col_map <- c(rep("white",2), rev(met.brewer(name="Demuth", n=10, type="continuous")[1:6]))
max_pm <- 15
smoke_percent[smoke_percent$nonsmokePM>max_pm, "nonsmokePM"] <- max_pm
ggplot() +
  geom_sf(data=smoke_percent, aes(colour=nonsmokePM , fill=nonsmokePM ), linewidth=0.0001) +
  #geom_point(data=station_mean, aes(x=lon, y=lat, fill=station_smokePM), shape=21, colour="black") +
  geom_polygon(data=states_map,aes(long, lat, group = group),
               fill=NA,color ="black",size=0.3) +
  scale_color_gradientn(colors = col_map, limits=c(0, max_pm)) +
  scale_fill_gradientn(colors= col_map, limits=c(0, max_pm)) +
  theme_classic() + theme(text = element_text(size=12)) 
ggsave("nonsmoke_PM_2016_2020.png", width = 8, height = 4)

max_percent <- 0.5
smoke_percent[smoke_percent$smoke_percent>max_percent, "smoke_percent"] <- max_percent
ggplot() +
  geom_sf(data=smoke_percent, aes(colour=smoke_percent, fill=smoke_percent), linewidth=0.0001) +
  #geom_point(data=station_mean, aes(x=lon, y=lat, fill=station_smokePM), shape=21, colour="black") +
  geom_polygon(data=states_map,aes(long, lat, group = group),
               fill=NA,color ="black",size=0.3) +
  scale_color_gradientn(colors = col_map, limits=c(0, max_percent), labels=scales::percent) +
  scale_fill_gradientn(colors= col_map, limits=c(0, max_percent), labels=scales::percent) +
  theme_classic() + theme(text = element_text(size=12)) 
ggsave("smoke_percent_2016_2020.png", width = 8, height = 4)



