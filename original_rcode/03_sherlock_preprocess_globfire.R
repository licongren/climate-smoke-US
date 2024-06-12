library(foreach)
library(doParallel)
library(doSNOW)
library(raster)
library(terra)
library(sf)
library(exactextractr)
library(dplyr)
library(tidyr)

rm(list=ls())
gc()

globfire_path <- "~/oak_space/mhqiu/Globfire_2000_2021"
gfed_path <- "~/oak_space/mhqiu/smoke_climate/input_data"

sf::sf_use_s2(TRUE)

gfed_na <- readRDS(paste0(gfed_path, "/GFED4S_monthly_1997_2021_NorthAmerica_025deg_df.rds"))
gfed_coord <- filter(gfed_na, year==1997, month==1)
gfed_coord_raster <- rast(gfed_coord[,c("lon","lat","gfed_cell_id")])

gfed_grid <- readRDS(paste0(gfed_path, "/GFED4S_NorthAmerica_025deg_grid_EPSG4326.rds"))
gfed_grid$area <- st_area(gfed_grid)

##################################################
#            read globfire datasets              # 
##################################################
globfire_list <- list.files(globfire_path, pattern = "*.shp", full.names = T)

get_cell_coverage <- function(ddd){
  data_tmp <- filter(globfire_agg, IDate==ddd)
  gfed_coord_raster <- rast(gfed_coord[,c("lon","lat","gfed_cell_id")])
  coverage_tmp <- coverage_fraction(gfed_coord_raster, data_tmp)[[1]]
  coverage_tmp <- as.data.frame(coverage_tmp, xy=T) %>% mutate(date=ddd)
  colnames(coverage_tmp) <- c("lon","lat","coverage","date")
  coverage_tmp <- filter(coverage_tmp, coverage>0)
  return(coverage_tmp)
}

cell_coverage_combined <- NULL
for (ii in globfire_list){
  print(ii)
  globfire_tmp <- st_read(ii) %>% filter(Type == "ActiveArea") 
  
  globfire_na <- st_intersection(globfire_tmp, 
                                 st_set_crs(st_as_sf(as(raster::extent(-170, -55, 10, 75), "SpatialPolygons")), st_crs(globfire_tmp)))
  
  globfire_agg <- globfire_na %>% group_by(IDate) %>%
    summarize(geometry = st_union(geometry)) %>% ungroup()
  
  date_unique <- unique(globfire_agg$IDate) 
  coverage_combine <- NULL
  cl <- makeCluster(8)  #### change the number of clusters based on the configuration of the computing environment
  registerDoSNOW(cl)
  pb <- txtProgressBar(min = 1, max = length(date_unique), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  coverage_combine  <- foreach(ii=1:length(date_unique),.options.snow = opts,.packages=c("dplyr", "sf", "exactextractr", "terra"),.combine=bind_rows) %dopar% 
    get_cell_coverage(date_unique[ii])
  stopCluster(cl)
  
  cell_coverage_combined <- bind_rows(cell_coverage_combined, coverage_combine)
}

coverage_final <- left_join(cell_coverage_combined, gfed_coord[,c("lat","lon","gfed_cell_id")], by = c("lon", "lat"))
coverage_final <- left_join(coverage_final, gfed_grid, by="gfed_cell_id") %>%
  select(-geometry) %>% mutate(burned_area_m2=coverage*area)

saveRDS(coverage_final, paste0(gfed_path, "/globfire_daily_burned_area_2001_2021.rds"))

#globfire_tmp1 <- globfire_tmp %>% st_simplify(preserveTopology = FALSE, dTolerance = 1000) ### simplify to 1km
# globfire_tmp1$geometry <- globfire_tmp1$geometry %>%
#   s2::s2_rebuild() %>% sf::st_as_sfc()

#globfire_tmp1 <- globfire_tmp %>% filter(Type == "ActiveArea") 


