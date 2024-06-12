library(maps);library(mapdata);library(fields)
library(terra);library(ncdf4)
library(openxlsx)
library(ggplot2)
library(plotly)
library(maptools);library(rgdal);library(rgeos)
library(dplyr)   
library(sf)
library(tidyverse)
library(stargazer)
library(splines2)
library(MetBrewer)
library(rhdf5)

code_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/rcode"
data_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data"
figure_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/figures"
result_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/results"
function_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/general_data/general_code/r_functions/"
db_path <- "/Users/mhqiu/BurkeLab Dropbox/Projects/smoke-climate"

source(paste(function_path,"aggregate_matrix_nx_ny.r",sep=""))
source(paste(function_path,"used_palette.r",sep=""))
source(paste(function_path,"us_states_name.r",sep=""))
source(paste(function_path,"calculate_grid_area.r",sep=""))

setwd(figure_path)
states_map <- map_data("state")

# global_region <- readOGR("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/mortality_climate/data/IPCC-WGI-regions-v4/IPCC-WGI-reference-regions-v4.shp")
# global_region_vec <- vect(global_region)

world <- readOGR("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/mortality_climate/data/world_shp/ne_110m_admin_0_countries.shp")
#world_shp <- gSimplify(world, tol=50, topologyPreserve=TRUE)
world_shp_data <- fortify(world)

us_shp <- readOGR("/Users/mhqiu/Dropbox_Personal/Dropbox/general_data/shapefiles/us_nation/cb_2018_us_nation_20m.shp")
us_shp_data <- fortify(us_shp)

globfire <- readRDS(paste0(data_path,"/globfire_daily_burned_area_2001_2021.rds"))
gfed <- readRDS(paste0(data_path,"/GFED4S_monthly_1997_2021_NorthAmerica_025deg_df.rds"))

####### aggregate to the annual level ##########
globfire_year <- globfire  %>% filter(lat>32, lat<49, lon> -125, lon< -102) %>%
  mutate(year=as.numeric(substr(date,1,4)), month=as.numeric(substr(date,6,7)))%>% 
  group_by(year) %>%
  summarise(ba_globfire=sum(as.numeric(burned_area_m2), na.rm=T)*1e-6) %>% ungroup()

gfed_year <- gfed %>% group_by(year) %>% filter(lat>32, lat<49, lon> -125, lon< -102) %>%
  summarise(ba_gfed=sum(burned_area_km2, na.rm=T), DM=sum(DM_kg, na.rm=T)*1e-7) %>% ungroup() %>%
  mutate(ba_gfed=replace(ba_gfed, ba_gfed==0, NA))

globfire_gfed <- full_join(globfire_year, gfed_year, by="year") %>%
  pivot_longer(cols=ba_globfire:DM)

ggplot(filter(globfire_gfed), aes(x=year, y=value, colour=name)) +
  geom_line() + theme(text=element_text(size=16)) + labs(y="Burned area (km2)") 
ggsave("gfed_globfire_ba_west_us.png", width=10, height=6)

####### aggregate to the grid annual level ##########
# globfire_year <- globfire  %>%
#   mutate(year=as.numeric(substr(date,1,4)), month=as.numeric(substr(date,6,7)))%>% 
#   group_by(year, gfed_cell_id) %>%
#   summarise(ba_globfire=sum(as.numeric(burned_area_m2), na.rm=T)*1e-6) %>% ungroup()

gfed_year <- gfed %>% group_by(gfed_cell_id, lon, lat) %>% filter(year>1997, year<2017) %>%
  summarise(ba_gfed=sum(burned_area_km2, na.rm=T), DM=sum(DM_kg, na.rm=T)*1e-6) %>% ungroup() %>%
  mutate(ba_gfed=replace(ba_gfed, ba_gfed==0, NA)) %>%
  mutate(dm_ba=DM/ba_gfed) %>% filter(!is.na(ba_gfed), !is.na(DM))

ggplot(gfed_year, aes(x=lon, y=lat, fill=dm_ba)) +
  geom_raster(alpha=1) + theme(text=element_text(size=16), axis.text = element_blank()) + labs(fill="DM / Burned area\n(1000 ton per km2)") +
  scale_fill_gradientn(colors=rev(met.brewer(name="Homer1", n=15, type="continuous")), limit=c(0,25), breaks=c(0,2,10,25), trans = "pseudo_log") #+
  #facet_wrap(~year)
ggsave("gfed_dm_ba_ratio.png", width=12, height=7)


gfed_month <- filter(gfed, lon >-125, lon< -115, lat>45, lat<50) %>% group_by(year) %>% filter(year>1997, year<2017) %>%
  summarise(ba_gfed=sum(burned_area_km2, na.rm=T), DM=sum(DM_kg, na.rm=T)*1e-6) %>% ungroup() %>%
  mutate(ba_gfed=replace(ba_gfed, ba_gfed==0, NA)) %>%
  mutate(dm_ba=DM/ba_gfed) %>% filter(!is.na(ba_gfed), !is.na(DM))

ggplot(gfed_month, aes(x=year, y=dm_ba)) + geom_line(size=1) + 
  labs(y="DM / Burned area\n(1000 ton per km2)") +
  theme(text=element_text(size=14)) + scale_x_continuous(breaks=1998:2016,labels =1998:2016)
ggsave("gfed_dm_ba_ratio_year.png", width=10, height=5)


####################  focusing on western US forest #####################
path_project <- "~/BurkeLab Dropbox/Projects/smoke-climate/"

aba = read.csv("~/Downloads/pnas.1607171113.sd01.csv", skip = 10) %>%
  mutate(area.burned = log(area.burned)) %>%
  select(year, area.burned, Z.VPD)

vpd_west <-  readRDS(paste0(path_project,"gridmet_vpd_1979_2022_1deg_forest_cover.rds")) %>%
  filter(west==1, month %in% 3:9) %>% group_by(year) %>% summarise(vpd=weighted.mean(vpd, forest)) %>% ungroup()

globfire_year <- globfire %>% filter(lon > -125, lon < -102.5, lat> 32, lat<50) %>%
  mutate(year=as.numeric(substr(date,1,4)), month=as.numeric(substr(date,6,7)))%>% 
  group_by(year) %>% filter( month %in% 3:9) %>%
  summarise(ba_globfire=sum(as.numeric(burned_area_m2), na.rm=T)*1e-6) %>% ungroup()

gfed_year <- gfed %>% group_by(year) %>% 
  filter(lon > -125, lon < -102.5, lat> 32, lat<50, month %in% 3:9) %>%
  summarise(ba_gfed=sum(burned_area_km2, na.rm=T), DM=sum(DM_kg, na.rm=T)*1e-7) %>% ungroup() %>%
  mutate(ba_gfed=replace(ba_gfed, ba_gfed==0, NA)) 

globfire_gfed <- full_join(globfire_year, gfed_year, by="year") %>% 
  full_join(vpd_west, by="year") %>%
  left_join(aba , by="year") %>%
  filter(!is.na(ba_globfire), !is.na(ba_gfed), !is.na(area.burned)) %>%
  mutate(ba_globfire=log(ba_globfire), ba_gfed=log(ba_gfed), DM=log(DM))

ggplot(globfire_gfed) + 
  geom_point(aes(Z.VPD, DM)) + 
  geom_smooth(aes(Z.VPD, DM), method = "lm") + 
  geom_text(aes(-Inf, Inf), data = data.frame(x=1,y=2), vjust = 3.5, hjust = -0.1, 
            label = paste("R2 =", round(summary(lm(DM ~ Z.VPD, globfire_gfed))$r.squared, 2))) + 
  theme_light() + 
  labs(y = "log(DM)", title = "GFED DM and VPD, march-sep")
ggsave("gfed_dm_vpd_2000-2015.png",  width = 8, height = 8)

p = ggplot(aba_summer) + 
  geom_point(aes(narr_vpd, DM_kg)) + 
  geom_smooth(aes(narr_vpd, DM_kg), method = "lm") + 
  geom_text(aes(-Inf, Inf), data = data.frame(x=1,y=2), vjust = 3.5, hjust = -0.1, 
            label = paste("R2 =", round(summary(lm(DM_kg ~ narr_vpd, aba_summer))$r.squared, 2))) + 
  theme_light() + 
  labs(x = "std_vpd", y = "log(DM emissions)", title = "gfed DM emissions, narr vpd, mar-sep")
ggsave(file.path(path_github, "figures", "gfed_narr_2000-2015.png"), 
       width = 8, height = 8)



