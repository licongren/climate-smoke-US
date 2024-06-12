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

#########  read monthly burned area and emissions from GFED4  ########
gfed_path <- paste(data_path, "/GFED4",sep="")
gfed_files <- list.files(path=gfed_path, full.names = T, pattern="*.hdf5")

veg_region <- as.numeric(h5read(gfed_files[1], "/emissions/01/partitioning"))
veg_region <- rast(t(veg_region ))
ext(veg_region ) <- c(-180,180,-90,90)

#h5ls(gfed_files[11], all = T)

grid_area <- rast(t(h5read(gfed_files[1], "/ancill/grid_cell_area")))
ext(grid_area) <- c(-180,180,-90,90)
month_list <- c("01","02","03","04","05","06","07","08","09","10","11","12")

gfed_combined <- NULL
for (ii in 1:length(gfed_files)){
year <- as.numeric(substr(gfed_files[ii], 90,  93))
print(year)

for (mm in month_list){
  print(mm)
month <- as.numeric(mm)

tmp_emis_DM  <- terra::rast(t(h5read(gfed_files[ii], paste("/emissions/",mm,"/DM",sep="")))) # unit: kg per m2 per month
ext(tmp_emis_DM) <- c(-180,180,-90,90)
values(tmp_emis_DM) <- values(tmp_emis_DM) * values(grid_area)  ## unit: kg per month

tmp_emis_C  <- terra::rast(t(h5read(gfed_files[ii], paste("/emissions/",mm,"/C",sep="")))) # unit: g C per m2 per month
ext(tmp_emis_C) <- c(-180,180,-90,90)
values(tmp_emis_C) <- values(tmp_emis_C) * values(grid_area) *1e-3  ## unit: kg per month

gfed_emis1 <- as.data.frame(tmp_emis_DM, xy=T)
colnames(gfed_emis1) <- c("lon","lat","DM_kg")
gfed_emis2 <- as.data.frame(tmp_emis_C, xy=T)
colnames(gfed_emis2) <- c("lon","lat","carbon_emis_kg")

gfed_emis <- inner_join(gfed_emis1, gfed_emis2, by=c("lon","lat")) %>% filter(lat > -60) 

if(year < 2017){   ### from 2017 onwards GFED4s does not have burned area data
tmp_frac  <- terra::rast(t(h5read(gfed_files[ii], paste("/burned_area/",mm,"/burned_fraction",sep=""))))
ext(tmp_frac) <- c(-180,180,-90,90)
values(tmp_frac) <- values(tmp_frac) * values(grid_area) *1e-6

area_raster <- tmp_frac
gfed_area <- as.data.frame(area_raster, xy=T)
colnames(gfed_area) <- c("lon","lat","burned_area_km2")
gfed_area <- gfed_area  %>% filter(lat > -60) 

gfed_df <- inner_join(gfed_emis, gfed_area, by=c("lon","lat"))
}
else{
  gfed_df <- gfed_emis
}
  gfed_df <- gfed_df %>% mutate(year=year, month=month)
  gfed_combined <- bind_rows(gfed_combined, gfed_df) 
}
}

saveRDS(gfed_combined, paste(data_path,"/GFED4S_monthly_1997_2021_global_025deg_df.rds", sep=""))



############################################################################# 
################   get the emissions from each partition veg type  ###############
#############################################################################
gfed_path <- paste(data_path, "/GFED4",sep="")
gfed_files <- list.files(path=gfed_path, full.names = T, pattern="*.hdf5")

#h5ls(gfed_files[11], all = T)

# grid_area <- rast(t(h5read(gfed_files[1], "/ancill/grid_cell_area")))
# ext(grid_area) <- c(-180,180,-90,90)
month_list <- c("01","02","03","04","05","06","07","08","09","10","11","12")
veg_list <- c("AGRI","BORF","DEFO","PEAT","SAVA","TEMF")
  
gfed_combined <- NULL
for (ii in 1:length(gfed_files)){
  year <- as.numeric(substr(gfed_files[ii], 90,  93))
  print(year)
  
  gfed_month <- NULL
  for (mm in month_list){
    print(mm)
    month <- as.numeric(mm)
    
    gfed_emis <- NULL
    for (vv in veg_list){
    tmp_emis_DM  <- terra::rast(t(h5read(gfed_files[ii], paste("/emissions/",mm,"/partitioning/DM_",vv,sep="")))) # unit: kg per m2 per month
    ext(tmp_emis_DM) <- c(-180,180,-90,90)

    tmp_emis_DM <- as.data.frame(tmp_emis_DM, xy=T)
    colnames(tmp_emis_DM) <- c("lon","lat",vv)
     
    tmp_emis_DM <-  tmp_emis_DM %>% filter(lat > -60) 
    
    if (is.null(gfed_emis)){
      gfed_emis <- tmp_emis_DM
    }
    else{
    gfed_emis <- inner_join(gfed_emis, tmp_emis_DM, by=c("lon","lat")) 
    }
    }
    gfed_emis <- gfed_emis %>% mutate(year=year, month=month, total=AGRI+BORF+DEFO+PEAT+SAVA+TEMF) %>% filter(total>0) %>% select(-total)
    gfed_month <- bind_rows(gfed_month, gfed_emis) 
  
  # gfed_df <- gfed_month %>% group_by(lon, lat) %>% summarise_all(.funs = sum, na.rm=T) %>% ungroup() %>%
  #   mutate(year=year)
  }
  gfed_combined <- bind_rows(gfed_combined, gfed_month)
}

saveRDS(gfed_combined, paste(data_path,"/GFED4S_monthly_1997_2021_global_025deg_veg_partition.rds", sep=""))


######################### get monthly and daily emissions by veg types ################
veg_part <- readRDS(paste(data_path,"/GFED4S_monthly_1997_2021_global_025deg_veg_partition.rds", sep=""))
veg_part %>% rowwise() %>% mutate(total=AGRI+BORF+DEFO+PEAT+SAVA+TEMF) %>% ungroup()

gfed_daily <- readRDS(paste(data_path,"/GFED4S_daily_2003_2021_NorthAmerica_025deg_df.rds", sep=""))
gfed_monthly <- readRDS(paste(data_path,"/GFED4S_monthly_1997_2021_NorthAmerica_025deg_df.rds", sep=""))

gfed_daily_veg <- left_join(gfed_daily, veg_part, by=c("lon","lat","year","month")) %>%
  replace(is.na(.), 0) %>%
  mutate(DM_AGRI=AGRI*daily_DM_kg, DM_BORF=BORF*daily_DM_kg, DM_DEFO=DEFO*daily_DM_kg,
         DM_PEAT=PEAT*daily_DM_kg, DM_SAVA=SAVA*daily_DM_kg, DM_TEMF=TEMF*daily_DM_kg) %>% select(-AGRI, -BORF, -DEFO, -TEMF, -SAVA, -PEAT)
saveRDS(gfed_daily_veg, paste(data_path,"/GFED4S_daily_2003_2021_NorthAmerica_025deg_VEG.rds", sep=""))


gfed_monthly_veg <- left_join(gfed_monthly, veg_part, by=c("lon","lat","year","month")) %>%
  replace(is.na(.), 0) %>%
  mutate(DM_AGRI=AGRI*DM_kg, DM_BORF=BORF*DM_kg, DM_DEFO=DEFO*DM_kg,
         DM_PEAT=PEAT*DM_kg, DM_SAVA=SAVA*DM_kg, DM_TEMF=TEMF*DM_kg) %>% select(-AGRI, -BORF, -DEFO, -TEMF, -SAVA, -PEAT)
  
saveRDS(gfed_monthly_veg, paste(data_path,"/GFED4S_monthly_2003_2021_NorthAmerica_025deg_VEG.rds", sep=""))


gfed_veg <- gfed_combined %>% select(-year) %>% group_by(lon, lat) %>% 
  summarise_all(.funs = sum, na.rm=T) %>% ungroup() 

gfed_veg1 <- gfed_veg %>% rowwise() %>% mutate(total=AGRI+BORF+DEFO+PEAT+SAVA+TEMF) %>%
  mutate_at(3:8, ~. / total) %>% ungroup()
gfed_veg_long <- gfed_veg1 %>% filter(total>0) %>% select(-total) %>% 
  pivot_longer(cols=AGRI:TEMF) %>% group_by(lon, lat) %>% arrange(desc(value)) %>% top_n(1) %>% ungroup()
  
gfed_coord <- readRDS(paste0(db_path, "/Fire_emis/GFED4S_monthly_1997_2021_NorthAmerica_025deg_df.rds")) %>% distinct(lat, lon, gfed_cell_id)
gfed_veg_na <- inner_join(gfed_coord, gfed_veg_long, by=c("lat","lon")) %>% mutate(veg=name) %>% select(-name, -value)  
saveRDS(gfed_veg_na, paste0(data_path, "/GFED_veg_type_northamerica.rds"))

### to simplify and only use three veg types
gfed_veg1 <- gfed_veg %>% rowwise() %>% mutate(total=AGRI+BORF+DEFO+PEAT+SAVA+TEMF) %>%
  mutate(TEMF=TEMF + DEFO + BORF) %>% select(-DEFO, -BORF) %>%
  mutate_at(3:6, ~. / total) %>% ungroup()

gfed_veg_long <- gfed_veg1 %>% filter(total>0) %>% select(-total) %>% 
  pivot_longer(cols=AGRI:TEMF) %>% group_by(lon, lat) %>% arrange(desc(value)) %>% top_n(1) %>% ungroup()

ggplot(filter(gfed_veg_long, name!="PEAT")) + geom_tile(aes(x=lon,y=lat, fill=name, alpha=value), colour=NA, size=0.00001) +
  labs(x="",y="", fill="Vegetation", alpha="Percent") + theme(text = element_text(size=16)) +
  coord_cartesian(xlim=c(-170, -55), ylim=c(15,65))
ggsave("gfed_veg_3type_northAmerica.png", width=9, height=4)


#################   compare emissions from veg and non-veg ################# 
gfed_path <- paste(data_path, "/GFED4",sep="")
gfed_files <- list.files(path=gfed_path, full.names = T, pattern="*.hdf5")

h5ls(gfed_files[10], all = T)[300:500,]

mm <- "07"
total <- terra::rast(t(h5read(gfed_files[10], paste("/emissions/",mm,"/DM",sep="")))) # unit: kg per m2 per month
AGRI  <- terra::rast(t(h5read(gfed_files[10], paste("/emissions/",mm,"/partitioning/DM_AGRI",sep="")))) # unit: kg per m2 per month
BORF  <- terra::rast(t(h5read(gfed_files[10], paste("/emissions/",mm,"/partitioning/DM_BORF",sep="")))) # unit: kg per m2 per month
DEFO  <- terra::rast(t(h5read(gfed_files[10], paste("/emissions/",mm,"/partitioning/DM_DEFO",sep="")))) # unit: kg per m2 per month
PEAT  <- terra::rast(t(h5read(gfed_files[10], paste("/emissions/",mm,"/partitioning/DM_PEAT",sep="")))) # unit: kg per m2 per month
SAVA  <- terra::rast(t(h5read(gfed_files[10], paste("/emissions/",mm,"/partitioning/DM_SAVA",sep="")))) # unit: kg per m2 per month
TEMF  <- terra::rast(t(h5read(gfed_files[10], paste("/emissions/",mm,"/partitioning/DM_TEMF",sep="")))) # unit: kg per m2 per month

combined <- AGRI
values(combined) <- values(AGRI) + values(BORF) + values(DEFO) + values(PEAT) + values(SAVA) + values(TEMF)
diff <- total
values(diff) <- values(total) - values(combined)

plot(diff)
  
  
########################################################### 
################   merge with US Ecoregions ###############
########################################################### 
# gfed1 <- readRDS(paste(data_path,"/GFED4S_monthly_1997_2005_global_025deg_df.rds", sep=""))
# gfed2 <- readRDS(paste(data_path,"/GFED4S_monthly_2006_2021_global_025deg_df.rds", sep=""))
# gfed_1997_2021 <- bind_rows(gfed1, gfed2)
# saveRDS(gfed_1997_2021, paste(data_path,"/GFED4S_monthly_1997_2021_global_025deg_df.rds", sep=""))

#gfed_1997_2021 <- readRDS(paste(data_path,"/GFED4S_monthly_1997_2021_global_025deg_df.rds", sep=""))
#gfed_glob_summ <- gfed_1997_2021 %>% group_by(year, lat, lon) %>% summarise(carbon_emis_kg=sum(carbon_emis_kg,na.rm=T)) %>% ungroup()
# gfed_glob_summ <- gfed_1997_2021 %>% group_by(year) %>% summarise(carbon_emis_kg=sum(carbon_emis_kg,na.rm=T)) %>% ungroup()
# 
# #gfed_glob_summ$date <- as.Date(paste(gfed_glob_summ$year, gfed_glob_summ$month, "15", sep="-")) 
# ggplot(gfed_glob_summ, aes(x=year, y=carbon_emis_kg)) + geom_line()
# ggsave("GFED_global_year.png",width=10, height=4.5)

# gfed_glob_summ$carbon_emis_log <- log(gfed_glob_summ$carbon_emis_kg)
# gfed_glob_summ <- gfed_glob_summ %>% mutate(carbon_emis_log=replace(carbon_emis_log, carbon_emis_log<0 ,0))
# ggplot(filter(gfed_glob_summ, year<2010)) +  #geom_polygon(data=world_shp_data, aes(x=long, y=lat)) + 
#   geom_tile(aes(x=lon, y=lat, fill=carbon_emis_log)) +
#   scale_fill_gradientn(colors=c("white",rev(met.brewer(name="Greek", n=16, type="continuous"))))  + 
#   facet_wrap(~year, ncol=5)
# ggsave("global_1997_2021_gfed_emis.png", width=13, height=7)
dropbox_path <- "/Users/mhqiu/BurkeLab Dropbox/Projects/smoke-climate/"

gfed_1997_2021 <- readRDS(paste0(data_path,"/GFED4S_monthly_1997_2021_global_025deg_df.rds"))
gfed_na <- filter(gfed_1997_2021 , lon > -170, lon < -55, lat > 10, lat <75)
gfed_na_coord <- readRDS(paste0(data_path,"/GFED4S_monthly_1997_2021_NorthAmerica_025deg_df.rds")) %>%
filter(year==1997, month==1)

gfed_na <- left_join(gfed_na, gfed_na_coord[,c("lon","lat","gfed_cell_id")], by=c("lon","lat"))
saveRDS(gfed_na, paste0(data_path,"/GFED4S_monthly_1997_2021_NorthAmerica_025deg_df.rds"))

gfed_na <- readRDS(paste0(data_path,"/GFED4S_monthly_1997_2021_NorthAmerica_025deg_df.rds"))
gfed_na_2020 <-  filter(gfed_na, year==2020) %>% group_by(year, month, lat, lon) %>% summarise(carbon_emis=sum(carbon_emis_kg,na.rm=T),
                                                              carbon_emis_log=log10(carbon_emis),
                                                              DM=sum(DM_kg,na.rm=T),
                                                              DM_log=log10(DM_kg)) %>% ungroup()

ggplot() +  geom_tile(data=filter(gfed_na_2020, year==2020, month==9),aes(x=lon, y=lat, fill=DM*1e-3)) +
  geom_polygon(data=filter(us_shp_data, long>-130, long< -50, lat>20),aes(long, lat, group = group),fill=NA,color ="blue",size=0.7)  +#geom_polygon(data=world_shp_data, aes(x=long, y=lat)) +
  scale_fill_gradientn(colors=c(rep("white",2),rev(met.brewer(name="Greek", n=16, type="continuous"))),
                      trans = "pseudo_log", breaks=c(100,1e4, 1e6, 1e8)) +
  labs(x="",y="",fill="DM, tons") + coord_sf(xlim=c(-130,-110), ylim=c(10,65)) + theme_bw() +
  theme(panel.border = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), text = element_text(size=15))
ggsave("northAmerica_2020_sep_gfed_DM_west.png", width=11, height=7)


gfed_na_summ <- gfed_na %>% group_by(lat, lon) %>% summarise(carbon_emis=sum(carbon_emis_kg,na.rm=T),
                                                             carbon_emis_log=log10(carbon_emis),
                                                             DM=sum(DM_kg,na.rm=T),
                                                             DM_log=log10(DM_kg)) %>% ungroup()
# gfed_na_summ <- gfed_na_summ %>% mutate(carbon_emis_log=replace(carbon_emis_log, carbon_emis_log<4, 0),
#                                         DM_log=replace(DM_log, DM<4, 0),)
ggplot() +  geom_tile(data=filter(gfed_na_summ, DM_log>6.9),aes(x=lon, y=lat, fill=DM_log)) +
    geom_polygon(data=filter(us_shp_data, long>-130, long< -50, lat>20),aes(long, lat, group = group),fill=NA,color ="blue",size=0.7)  +#geom_polygon(data=world_shp_data, aes(x=long, y=lat)) +
   scale_fill_gradientn(colors=c(rep("white",2),rev(met.brewer(name="Greek", n=16, type="continuous"))), breaks=7:11, labels=c("1e4","1e5","1e6","1e7","1e8"))  +
  labs(x="",y="",fill="DM, tons") + coord_sf(xlim=c(-162,-55))
ggsave("northAmerica_1997_2021_gfed_DM.png", width=11, height=7)

# gfed_na <- readRDS(paste(data_path,"/GFED4S_monthly_1997_2021_NorthAmerica_025deg_df.rds", sep=""))
# gfed_coord <- filter(gfed_na, year==1997, month==1)
# gfed_coord_raster <- rast(gfed_coord[,c("lon","lat","gfed_cell_id")])
# gfed_poly <- st_as_sf(as.polygons(gfed_coord_raster), crs=4326)
# saveRDS(gfed_poly, paste(data_path,"/GFED4S_NorthAmerica_025deg_grid.rds", sep=""))

# lonlat <- st_as_sf(gfed_coord[,c("lon","lat","gfed_cell_id")], coords = c("lon", "lat"), crs = 4326) # picking a basic crs, this depends on data. ideally if comes as spatial data, then crs is implicit already
# saveRDS(lonlat, paste(data_path,"/GFED4S_NorthAmerica_025deg_grid.rds", sep=""))


##########  merge with the EPA-level3 ecoregion shp files #############
ecoregion <- readRDS(paste0(data_path, "/ecoregion/NA_Eco_level3_MULTIPOLYGON_WGS84.rds"))
gfed_ecoregion <-  st_join(lonlat, st_transform(ecoregion, st_crs(lonlat)), largest = TRUE)
gfed_ecoregion1 <- as.data.frame(gfed_ecoregion) %>%  distinct(gfed_cell_id, .keep_all=T) %>% ungroup()

ecoregion_center1 <- st_centroid(st_transform(ecoregion, st_crs(lonlat)))
ecoregion_center1$region_center_lon <- as.numeric(gsub("^c\\((.*),.*$", "\\1", as.character(ecoregion_center1$geometry)))
ecoregion_center1$region_center_lat <- as.numeric(gsub(".*[, ](.*)\\)$", "\\1", as.character(ecoregion_center1$geometry)))
ecoregion_center1 <- as.data.frame(ecoregion_center1) %>% select(L3_index, region_center_lon,region_center_lat)
gfed_ecoregion2 <- left_join(gfed_ecoregion1, ecoregion_center1, by="L3_index")

# ecoregion_list <- unique(filter(gfed_ecoregion2,  region_center_lat<58,  region_center_lat>19.9)$L3_index)
# #gfed_ecoregion2 %>% filter(L3_index %in% ecoregion_list) %>% distinct(NA_L3CODE, NA_L2CODE)
# gfed_ecoregion3 <- gfed_ecoregion2  %>% mutate(insample=0) %>% mutate(insample=replace(insample, L3_index %in% ecoregion_list, 1)) %>%
#       mutate(L3_index_agg=L3_index) %>% mutate(L3_index_agg=replace(L3_index_agg, NA_L3CODE%in%c("14.3.1","14.3.2", "15.5.1","15.5.2"), 200),
#                                                L3_index_agg=replace(L3_index_agg, NA_L3CODE%in%c("12.2.1","12.1.2"), 201),
#                                                L3_index_agg=replace(L3_index_agg, NA_L3CODE%in%c("15.2.1","15.2.2","14.2.1"), 202),
#                                                L3_index_agg=replace(L3_index_agg, NA_L3CODE%in%c("14.1.1","14.1.2","13.3.1"), 203),
#                                                L3_index_agg=replace(L3_index_agg, NA_L3CODE%in%c("6.2.6","9.3.1"), 204),
#                                                L3_index_agg=replace(L3_index_agg, NA_L3CODE%in%c("5.2.1","5.2.2"), 205),
#                                                L3_index_agg=replace(L3_index_agg, NA_L3CODE%in%c("8.1.8","8.1.9"), 206),
#                                                L3_index_agg=replace(L3_index_agg, NA_L3CODE%in%c("5.4.2","5.4.1"), 207),
#                                                L3_index_agg=replace(L3_index_agg, NA_L3CODE%in%c("4.1.1","4.1.2"), 208),
#                                                L3_index_agg=replace(L3_index_agg, NA_L3CODE%in%c("10.2.3","14.6.1","14.6.2"), 209),
#                                                L3_index_agg=replace(L3_index_agg, NA_L3CODE%in%c("3.4.3","3.4.4"), 210))
saveRDS(gfed_ecoregion2, paste0(data_path,"/GFED4S_NorthAmerica_ecoregion_link.rds"))
# 
# gfed_na_summ <- gfed_na %>% group_by(lat, lon) %>% summarise(carbon_emis=sum(carbon_emis_kg,na.rm=T),
#                                                              carbon_emis_log=log10(carbon_emis)) %>% ungroup()
# gfed_na_summ <- gfed_na_summ %>% mutate(carbon_emis_log=replace(carbon_emis_log, carbon_emis_log<4, 0))
# 
# ecoregion_df_agg <- filter(gfed_ecoregion3, insample==1) %>% distinct(NA_L3CODE,L3_index_agg)
# ecoregion_agg <- inner_join(ecoregion, ecoregion_df_agg, by="NA_L3CODE") %>% group_by(L3_index_agg) %>% summarise() 
#   
# ggplot() +  geom_tile(data=filter(gfed_na_summ, carbon_emis_log>6.9),aes(x=lon, y=lat, fill=carbon_emis_log)) +
#   geom_sf(data=filter(ecoregion_agg),fill=NA,color ="black",size=0.5)  +
#   geom_polygon(data=filter(us_shp_data, long>-130, long< -50, lat>20),aes(long, lat, group = group),fill=NA,color ="blue",size=0.7)  +#geom_polygon(data=world_shp_data, aes(x=long, y=lat)) +
#   #geom_polygon(data=states_map,aes(long, lat, group = group),fill=NA,color ="black")  +#geom_polygon(data=world_shp_data, aes(x=long, y=lat)) +
#   scale_fill_gradientn(colors=c(rep("white",2),rev(met.brewer(name="Greek", n=16, type="continuous"))), breaks=7:10, labels=c("1e4","1e5","1e6","1e7"))  +
#   labs(x="",y="",fill="Log (Fire carbon), tons") + coord_sf(xlim=c(-162,-55))
# ggsave("northAmerica_1997_2021_gfed_emis_ecoregion_full.png", width=11, height=7)

#############################################################################################
##########  get the daily emissions for 2003-2021 (no daily frac before 2003)  ##############
#############################################################################################
gfed_na <- readRDS(paste0(data_path,"/GFED4S_monthly_1997_2021_NorthAmerica_025deg_df.rds"))
gfed_daily_frac <- readRDS(paste0(data_path,"/GFED4S_daily_frac_2003_2021_NorthAmerica_025deg_df.rds"))
gfed_2003_2021 <- filter(gfed_na, year>2002)
gfed_daily <- inner_join(gfed_2003_2021, gfed_daily_frac, by=c("year","month","lon","lat")) 
gfed_daily <- gfed_daily %>% mutate(daily_C_kg = carbon_emis_kg * daily_frac,
                                    daily_DM_kg = DM_kg * daily_frac)
#aa <- filter(gfed_daily, year==2003, month==1) %>% group_by(year, lon, lat) %>% summarise(dm_sum=sum(daily_DM_kg, na.rm=T)) %>% ungroup() %>% arrange(desc(dm_sum))
saveRDS(gfed_daily, paste0(data_path,"/GFED4S_daily_2003_2021_NorthAmerica_025deg_df.rds"))

# gfed_us_daily  <- readRDS(paste(data_path,"/GFED4S_daily_2006_2021_US_025deg_df.rds", sep=""))
# gfed_ecoregion <- readRDS(paste(data_path,"/GFED4S_US_025deg_ecoregion_cellid.rds", sep=""))

# gfed_us_daily  <- left_join(gfed_us_daily, gfed_ecoregion, by="gfed_cell_id")
# gfed_us_daily1 <- filter(gfed_us_daily, !is.na(L3_KEY)) %>% select(-carbon_emis_kg ,-geometry, -region_center_lon, -region_center_lat, -burned_area_km2, -daily_frac) 
# gfed_us_daily1 <- gfed_us_daily1 %>% mutate(gfed_lat=lat, gfed_lon=lon) %>% select(-lat, -lon)
# gfed_us_daily1$date <- as.Date(paste(gfed_us_daily1$year, gfed_us_daily1$month, gfed_us_daily1$day, sep="-"))
# gfedcells <- unique(gfed_us_daily1$gfed_cell_id)
# date_df <- data.frame(date=rep(seq(as.Date("2006-01-01"), as.Date("2021-12-31"), by="days"), each=length(gfedcells)), gfed_cell_id=gfedcells)
# gfed_us_daily2 <- full_join(gfed_us_daily1, date_df, by=c("gfed_cell_id","date"))
# gfed_us_daily2 <- gfed_us_daily2 %>% group_by(gfed_cell_id) %>% mutate(L3_KEY=replace(L3_KEY, is.na(L3_KEY), unique(L3_KEY[!is.na(L3_KEY)])), L3_index=replace(L3_index, is.na(L3_index), unique(L3_index[!is.na(L3_index)])),
#                                                     gfed_lat=replace(gfed_lat, is.na(gfed_lat), unique(gfed_lat[!is.na(gfed_lat)])),
#                                                     gfed_lon=replace(gfed_lon, is.na(gfed_lon), unique(gfed_lon[!is.na(gfed_lon)]))) %>% ungroup()
# gfed_us_daily2 <- gfed_us_daily2 %>% mutate(year=as.numeric(substr(date,1,4)), month=as.numeric(substr(date,6,7)), day=as.numeric(substr(date,9,10))) %>% 
#                                      mutate(daily_C_kg=replace(daily_C_kg, is.na(daily_C_kg), 0))  %>% ungroup()
# 
# saveRDS(gfed_us_daily1, paste(data_path,"/GFED4S_daily_2006_2021_US_ecoregion_gfedcell.rds", sep=""))
# 
# ##### aggregate at monthly, but keep gfed cell id ########
# gfed_eco_monthly <-  gfed_us_daily %>% group_by(year, month, gfed_cell_id, lon, lat, L3_KEY, L3_index) %>% summarise(carbon_emis_kg=sum(daily_C_kg, na.rm=T)) %>% ungroup()
# gfed_eco_monthly$date <- as.Date(paste(gfed_eco_monthly$year, gfed_eco_monthly$month, "01", sep="-"))
# gfed_eco_monthly <- filter(gfed_eco_monthly, !is.na(L3_KEY))
# 
# gfed_cells <- unique(gfed_eco_monthly$gfed_cell_id)
# date_df <- data.frame(date=seq(as.Date("2006-01-01"), as.Date("2021-12-01"), by="months"))
# 
# gfed_eco_monthly_combined <- NULL
# for (ii in gfed_cells){
#   data_tmp <- filter(gfed_eco_monthly, gfed_cell_id==ii)
#   
#   date_df$gfed_cell_id <- ii
#   date_df$lon <- unique(data_tmp$lon)
#   date_df$lat <- unique(data_tmp$lat)
#   date_df$L3_KEY <- unique(data_tmp$L3_KEY)
#   date_df$L3_index <- unique(data_tmp$L3_index)
#   
#   data_tmp <- merge(data_tmp, date_df, by=c("date","lon","lat","L3_KEY","L3_index","gfed_cell_id"), all=T) %>% 
#               mutate(year=as.numeric(substr(date,1,4)),  month=as.numeric(substr(date,6,7))) %>% arrange(year, month)
#   
#   data_tmp <- data_tmp %>% mutate(carbon_emis_kg=replace(carbon_emis_kg, is.na(carbon_emis_kg), 0))
#   gfed_eco_monthly_combined <- bind_rows(gfed_eco_monthly_combined, data_tmp)
# }
# gfed_eco_monthly_combined <- gfed_eco_monthly_combined %>% mutate(gfed_lat=lat, gfed_lon=lon) %>% select(-lat, -lon)
# saveRDS(gfed_eco_monthly_combined, paste(data_path,"/GFED4S_monthly_2006_2021_US_ecoregion_gfedcell.rds", sep=""))
# 
# gfed_eco_monthly1 <- gfed_eco_monthly %>% group_by(year, month, L3_KEY, L3_index, date) %>% summarise(carbon_emis_kg=sum(carbon_emis_kg, na.rm=T)) %>% ungroup()
# ggplot(filter(gfed_eco_monthly1, L3_index<10), aes(x=date, y=carbon_emis_kg)) + geom_line() +
#   facet_wrap(~L3_KEY)
# ggsave("ecoregion_gfed4_fire_carbon_emis_monthly1.png", width=12, height=8)
# 
# 
# # ecoregion <- readRDS(paste(data_path, "/ecoregions/us_eco_l3_MULTIPOLYGON.rds",sep=""))
# # ecoregion_simpl <- st_simplify(ecoregion, preserveTopology = FALSE, dTolerance = 1000)
# 
# gfed_ecoregion_summ <- gfed_us_daily  %>% group_by(L3_KEY, L3_index, region_center_lon, region_center_lat, year, month, day) %>% summarise(carbon_emis_kg=sum(daily_C_kg, na.rm=T)) %>% ungroup() %>% filter(!is.na(L3_index))
# gfed_ecoregion_summ$date <- as.Date(paste(gfed_ecoregion_summ$year, gfed_ecoregion_summ$month, gfed_ecoregion_summ$day, sep="-"))
# 
# date_df <- data.frame(date=rep(seq(as.Date("2006-01-01"), as.Date("2021-12-31"), by="days"), each=85), L3_index=1:85)
# l3_key_list <- distinct(gfed_ecoregion, L3_index, L3_KEY,.keep_all=F)
# date_df <- left_join(date_df, l3_key_list, by="L3_index")
# gfed_ecoregion_summ1 <- full_join(gfed_ecoregion_summ, date_df, by=c("date","L3_index","L3_KEY"))
# gfed_ecoregion_summ1 <- gfed_ecoregion_summ1 %>% mutate(year=as.numeric(substr(date,1,4)), month=as.numeric(substr(date,6,7)), day=as.numeric(substr(date,9,10)))
# gfed_ecoregion_summ1 <- gfed_ecoregion_summ1 %>% mutate(carbon_emis_kg=replace(carbon_emis_kg, is.na(carbon_emis_kg), 0))
# saveRDS(gfed_ecoregion_summ1, paste(data_path,"/GFED4S_daily_2006_2021_US_ecoregion.rds", sep=""))
# 
# #l3_list<- unqiue(gfed_ecoregion_summ$L3_KEY)
# ggplot(filter(gfed_ecoregion_summ1, L3_index<25), aes(x=date, y=carbon_emis_kg)) + geom_line() +
#   facet_wrap(~L3_KEY)
# ggsave("ecoregion_gfed4_fire_carbon_emis_monthly1.png", width=12, height=8)
