library(terra);library(ncdf4)
library(raster)
library(gstat)
library(openxlsx)
library(ggplot2)
library(sf)
library(dplyr)   
library(tidyverse)
library(MetBrewer)
library(dplyr)
library(lubridate)

rm(list=ls())
gc()

setwd("~/BurkeLab Dropbox/Projects/smokePM-contribution/output")

states_map <- map_data("state")
us_shp <- vect("~/BurkeLab Dropbox/Projects/smoke_PM_prediction/data/1_grids/grid_10km_wgs84/grid_10km_wgs84.shp")
population = rast("~/BurkeLab Dropbox/Data/population/gpw_v4_population_count_rev11_2020_2pt5_min.tif")
# Crop population raster to CONUS
population = raster::crop(population, extent(-124.848974, -66.885444, 24.396308, 49.384358))

#smoke_pm_raw <- readRDS("~/BurkeLab Dropbox/Projects/smokePM-contribution/data/epa_station_day_totalPM_smokePM_panel_20000101-20230705.rds")
smoke_pm_raw <- readRDS("~/BurkeLab Dropbox/projects/smokePM-prediction/data/EPA/station_smokePM_panel/epa_station_day_totalPM_smokePM_panel_20000101-20230914.rds")

# date_range <- c(seq(as.Date("2007-04-01"), as.Date("2007-05-30"), by="days"),
#                 seq(as.Date("2018-11-01"), as.Date("2018-11-30"), by="days"),
#                 seq(as.Date("2020-08-01"), as.Date("2020-09-30"), by="days"),
#                 seq(as.Date("2021-08-01"), as.Date("2021-09-30"), by="days"),
#                 seq(as.Date("2023-06-01"), as.Date("2023-07-05"), by="days"))

date_range <- seq(as.Date("2023-01-01"), as.Date("2023-08-31"), by="days")

#pick_date <- as.Date("2023-06-09")
data_summ_combined <- NULL
for (pick_date in date_range){
  print(as.Date(pick_date))
  smoke_pm =   smoke_pm_raw %>% 
    filter(date == pick_date, !is.na(smokePM))
  
  smoke_data_1 <- smoke_pm %>%
    rename(x=lon, y=lat, z=smokePM) %>%
    dplyr::select(x, y, z)
  
  idw <- terra::interpIDW(population, as.matrix(smoke_data_1), 
                          radius=2.5, power=5, smooth=0.5) 
  names(idw) <- "interp_smokePM"
  
  ### mask to US range
  idw <- terra::mask(idw, us_shp)
  population <- terra::mask(population, us_shp)
  writeRaster(idw, paste0("~/BurkeLab Dropbox/Projects/smokePM-contribution/output/interpolation/interpolated_smokePM_",
              as.Date(pick_date),".tif"), overwrite=T)
  
  df_smoke_pop <- data.frame(as.matrix(idw), as.matrix(population)) %>%
    rename(population=gpw_v4_population_count_rev11_2020_2pt5_min) %>%
    filter(!is.na(population), !is.na(interp_smokePM))
  
  data_summ <- df_smoke_pop %>%
    summarise( 
      ppl_times_smoke=sum(population*interp_smokePM),
      pop_weight_smoke=weighted.mean(interp_smokePM, w=population),
      ppl_above_50=sum(population*as.numeric(interp_smokePM>50)),
      ppl_above_75=sum(population*as.numeric(interp_smokePM>75)),
      ppl_above_100=sum(population*as.numeric(interp_smokePM>100))) %>%
    mutate(date=as.Date(pick_date))
  
  data_summ_combined <- bind_rows(data_summ_combined, data_summ)
}

data_summ_combined$date <- as.Date(data_summ_combined$date)
saveRDS(data_summ_combined, 
        "~/BurkeLab Dropbox/Projects/smokePM-contribution/output/df_ppl_smoke_20230101_20230831.rds")

data_summ_combined <-
        readRDS("~/BurkeLab Dropbox/Projects/smokePM-contribution/output/df_ppl_smoke_20230101_20230831.rds")

################## plotting ##############
data_sort <- data_summ_combined %>% arrange(desc(pop_weight_smoke)) 
plot_date <- data_sort[1:5,"date"]. ## only plot the top 5 polluted days

for (pick_date in plot_date){
  smoke_pm =   smoke_pm_raw %>% 
    filter(date == pick_date, !is.na(smokePM))
  
  smoke_data_1 <- smoke_pm %>%
    rename(x=lon, y=lat, z=smokePM) %>%
    dplyr::select(x, y, z)
  
  idw <- terra::interpIDW(population, as.matrix(smoke_data_1), 
                          radius=2.5, power=5, smooth=0.5) 
  names(idw) <- "interp_smokePM"
  idw_us <- terra::mask(idw, us_shp)
  values(idw_us)[values(idw_us)>400] <- 400
  values(idw_us)[is.na(values(idw_us))] <- 0
  
  smoke_data_1[smoke_data_1$z>400,"z"] <- 400
  
  library(tidyterra)
  ggplot() +
    geom_spatraster(data=idw_us, aes(fill=interp_smokePM)) +
    geom_polygon(data=states_map,aes(long, lat, group = group),fill=NA,color ="black") +
    geom_point(data=smoke_data_1, aes(x=x, y=y, fill=z), colour="black", shape=21) +
    scale_colour_gradientn(colors=c(rep("white",2),rev(met.brewer(name="Homer1", n=14, type="continuous"))),
                           limits=c(0,400), trans="pseudo_log", breaks=c(0,50,100,200,400)) +
    scale_fill_gradientn(colors=c(rep("white",2),rev(met.brewer(name="Homer1", n=14, type="continuous"))),
                         limits=c(0,400), trans="pseudo_log", breaks=c(0,50,100,200,400)) +
    theme_bw() + labs(x="", y="", fill="Smoke PM2.5 (ug)") +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          text = element_text(size=14)) +
    ggtitle(as.Date(pick_date))
  ggsave(paste0(as.Date(pick_date),"_idw_250.pdf"), width=10, height=5)
}

### calculate cumulated exposure by number of days in a year
data_cumsum <- data_summ_combined %>%
  mutate(year=as.numeric(substr(date,1,4)),
         nday=yday(date)) %>%
  group_by(year) %>%
  arrange(nday) %>%
  mutate(cum_pop_smoke=cumsum(pop_weight_smoke)) %>%
  ungroup()

ggplot(data_cumsum, aes(x=nday, y=cum_pop_smoke, colour=factor(year), group=factor(year))) +
  geom_line(linewidth=0.8) +
  theme_bw() +    
  labs(x= "", y="PM2.5 (ug/m3)\n", colour="",
       title="Cumulative population average smoke PM2.5") +
  scale_x_continuous(breaks=cumsum(c(31,28,31,30,31,30,31,31,30,31,30,31)),
                     labels=c("Jan"," Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep"," Oct","Nov","Dec")) +
  guides(colour=guide_legend(ncol=2)) +
  scale_colour_manual(values=c(rep("grey80",14), 
                               met.brewer(name="Homer1", n=7, type="discrete")[4:1])) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(linewidth=0.5), 
        text = element_text(size=18)) 
ggsave("cum_pop_smoke_0705_bydays.pdf", width=9, height=4.5)

### calculate cumulated exposure by year
# data_annual_summ <- data_summ_combined %>%
#   mutate(year=as.numeric(substr(date,1,4)),
#          month=as.numeric(substr(date,6,7))) %>%
#   group_by(year) %>%
#   summarise(cum_pop_smoke=sum(pop_weight_smoke),
#             cum_pop_smoke_july=sum(pop_weight_smoke*as.numeric(month<7))) %>%
#   ungroup() %>%
#   pivot_longer(cols=cum_pop_smoke:cum_pop_smoke_july)
# 
# data_annual_summ$name <- factor(data_annual_summ$name,
#                                 labels=c("Annual total", "Up to July.1st"))
# 
# ggplot(data_annual_summ, aes(x=year, y=value, colour=name)) +
#   geom_point(size=2) + geom_line(linewidth=0.8) +
#   theme_bw() +    
#   labs(x= "", y="PM2.5 (ug/m3)\n", colour="",
#        title="Cumulative population average smoke PM2.5") +
#   scale_x_continuous(breaks=seq(2006,2023, by=2)) +
#   scale_colour_manual(values=met.brewer(name="Cassatt2", n=2, type="discrete")) +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         axis.line = element_line(linewidth=0.5), 
#         text = element_text(size=18)) 
# ggsave("cum_pop_smoke_0705.pdf", width=8.5, height=4.5)
  

# data_summ_combined <- readRDS("~/BurkeLab Dropbox/Projects/smokePM-contribution/output/df_ppl_totalPM_2006_20230705.rds")
# 
# data_annual_summ <- data_summ_combined %>%
#   mutate(year=as.numeric(substr(date,1,4)),
#          month=as.numeric(substr(date,6,7))) %>%
#   group_by(year) %>%
#   summarise(cum_pop_smoke=sum(pop_weight_total, na.rm=T)) %>%
#   ungroup()
#   
# 
# #####
# plot_date <- c(
# seq(as.Date("2020-08-01"), as.Date("2020-09-30"), by="days"),
# seq(as.Date("2021-08-01"), as.Date("2021-09-30"), by="days"),
# seq(as.Date("2023-06-01"), as.Date("2023-06-09"), by="days"))
# 
# for (pick_date in plot_date){
#   smoke_pm =   smoke_pm_raw %>% 
#     filter(date == pick_date, !is.na(smokePM))
#   
#   smoke_data_1 <- smoke_pm %>%
#     rename(x=lon, y=lat, z=smokePM) %>%
#     dplyr::select(x, y, z)
#   
#   idw_us <- rast(paste0("~/BurkeLab Dropbox/Projects/smokePM-contribution/output/interpolation/interpolated_smokePM_",
#                      as.Date(pick_date),".tif"))
# 
#   #idw_us <- terra::mask(idw, us_shp)
#   values(idw_us)[values(idw_us)>400] <- 400
#   values(idw_us)[is.na(values(idw_us))] <- 0
#   
#   smoke_data_1[smoke_data_1$z>400,"z"] <- 400
#   
#   library(tidyterra)
#   ggplot() +
#     geom_spatraster(data=idw_us, aes(fill=interp_smokePM)) +
#     geom_polygon(data=states_map,aes(long, lat, group = group),fill=NA,color ="black") +
#     geom_point(data=smoke_data_1, aes(x=x, y=y, fill=z), colour="black", shape=21) +
#     scale_colour_gradientn(colors=c(rep("white",2),rev(met.brewer(name="Homer1", n=14, type="continuous"))),
#                            limits=c(0,400), trans="pseudo_log", breaks=c(0,50,100,200,400)) +
#     scale_fill_gradientn(colors=c(rep("white",2),rev(met.brewer(name="Homer1", n=14, type="continuous"))),
#                          limits=c(0,400), trans="pseudo_log", breaks=c(0,50,100,200,400)) +
#     theme_bw() + labs(x="", y="", fill="Smoke PM2.5 (ug)") +
#     theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#           axis.text = element_blank(), axis.ticks = element_blank(),
#           text = element_text(size=14)) +
#     ggtitle(as.Date(pick_date))
#   ggsave(paste0(as.Date(pick_date),"_idw_250.pdf"), width=10, height=5)
#   
# }
# 


