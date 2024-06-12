library(stringr)
library(lubridate)
library(ncdf4)
library(raster)
library(sf)
library(exactextractr)
library(dplyr)
library(tidyr)

rm(list=ls())
gc()

code_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/rcode"
data_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data"
figure_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/figures"
result_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/results/smoke_model"

setwd(figure_path)

db_proj_path <- "/Users/mhqiu/BurkeLab Dropbox/Projects"
smoke_pm <- readRDS(paste0(db_proj_path, "/smoke_PM_prediction/output/smokePM_predictions_20060101_20201231.rds"))

##################################################
# get monthly gridmet vars at the 10km smoke cell
##################################################
grid_smoke <- st_read(paste0(db_proj_path, "/smoke_PM_prediction/data/1_grids/grid_10km_wgs84/grid_10km_wgs84.shp")) 

path_gridmet <- file.path("~/BurkeLab Dropbox/Data/gridMET")
grid_gridmet <- terra::rast(file.path(path_gridmet, sprintf("tmmx_%s.nc", years[1])), lyrs=1) 
grid_smoke <- grid_smoke %>% st_transform(crs(grid_gridmet))

crosswalk_gridmet_smoke <- exact_extract(
    grid_gridmet, 
    grid_smoke, 
    coverage_area = T, 
    include_cell = T,
    include_cols = "ID") %>% bind_rows()

crosswalk_gridmet_smoke1 <- right_join(crosswalk_gridmet_smoke ,st_drop_geometry(grid_smoke), by = "ID") %>% 
    select(grid_id_10km = ID, grid_id_gridmet = cell, coverage_area)
  
saveRDS(crosswalk_gridmet_smoke1, file.path(data_path, "crosswalk_gridMET_smoke10km_grid.rds"))



# Aggregate to monthly smoke 10 km grid and save by year
# Takes ~10 minutes per year
crosswalk_gridmet_smoke <-  readRDS(file.path(data_path, "crosswalk_gridMET_smoke10km_grid.rds"))

years <- 2009:2011

origin_date = ymd("1900-01-01")
direct_variables = c(
  "pr", # Precipitation
  "sph", # Near-Surface Specific Humidity
  "vs", # Wind speed at 10 m
  "vpd", # Mean Vapor Pressure Deficit
  "erc", # Energy Release Component (model-G)
  "fm100", # 100-hour dead fuel moisture
  "fm1000" # 1000-hour dead fuel moisture
)

  derived_variables = data.frame(
    derived_variable = c("tmav", "tmav",  # Mean Near-Surface Air Temperature
                         "ravg", "ravg"), # Mean Near-Surface Relative Humidity
    variable = c("tmmn", # Minimum Near-Surface Air Temperature
                 "tmmx", # Maximum Near-Surface Air Temperature
                 "rmin", # Minimum Near-Surface Relative Humidity
                 "rmax") # Maximum Near-Surface Relative Humidity
  )
  
  variables = c(direct_variables, unique(derived_variables$derived_variable))
  variables <- c("vs")
  # Aggregate to monthly GFED 25 km grid and save by year
  # Takes ~10 minutes per year
  
  ## potential issues with vs_2010
  for (vbl in variables) {
    for (y in years) {
      out = vector("list", 12)
      if (vbl %in% direct_variables) {
        df_y = stack(file.path(path_gridmet, sprintf("%s_%s.nc", vbl, y)))
      } else {
        vbls = derived_variables %>% filter(derived_variable == vbl) %>% pull(variable)
        df_y = stack(file.path(path_gridmet, sprintf("%s_%s.nc", vbls[1], y)))
        df_y2 = stack(file.path(path_gridmet, sprintf("%s_%s.nc", vbls[2], y)))
      }
      dates_y = as.Date(as.numeric(gsub("^X", "", names(df_y))), origin = origin_date)
      for (m in 1:12) {
        d = grep(sprintf("^%s-%s", y, str_pad(m, 2, "left", 0)), dates_y)
        df = subset(df_y, d) %>% 
          values() %>% 
          as.data.frame() %>% 
          mutate(grid_id_gridmet = row_number()) %>% 
          pivot_longer(cols = starts_with("X"),
                       names_to = "date",
                       names_prefix = "X")
        if (vbl %in% derived_variables$derived_variable) {
          df2 = subset(df_y2, d) %>% 
            values() %>% 
            as.data.frame()%>% 
            mutate(grid_id_gridmet = row_number()) %>% 
            pivot_longer(cols = starts_with("X"),
                         names_to = "date",
                         names_prefix = "X")
          df = df %>% mutate(value = (value + df2$value)/2)
        }
        df = df %>% 
          mutate(date = origin_date + days(date),
                 year = as.character(y),
                 month = str_pad(m, 2, "left", 0)) %>% 
          group_by(grid_id_gridmet, year, month) %>% 
          summarize(value = mean(value, na.rm = T)) %>%
          ungroup() %>% 
          right_join(crosswalk_gridmet_smoke, by = "grid_id_gridmet") %>% 
          mutate(year = as.character(y), 
                 month = str_pad(m, 2, "left", 0)) %>% 
          group_by(grid_id_10km, year, month) %>% 
          summarize(!!vbl := weighted.mean(value, w = coverage_area, na.rm = T)) %>% 
          ungroup()
        out[[m]] = df
      }
      out = bind_rows(out)
      saveRDS(out, file.path(path_gridmet, "smoke_10km_grid", "monthly", sprintf("%s_%s.rds", vbl, y)))
    }
  }
  
  ################# combine gridmet at smoke grid level into a wide data.frame   #################
  vbl_list <- c("pr", "vs", "tmav","ravg")
  
  path_gridmet <- "~/BurkeLab Dropbox/Data/gridMET/smoke_10km_grid/monthly"
  gridmet_combined <- NULL
  for (vv in vbl_list){
    file_list <- list.files(path=path_gridmet, pattern=paste0("^",vv), full.names = T)
  data_combined <- NULL  
  for (ii in file_list){
    print(ii)
    data_tmp <- readRDS(ii)
    data_combined <- bind_rows(data_combined, data_tmp)
  }  
  if (is.null(gridmet_combined)){
    gridmet_combined <- data_combined
  }
  else{
    gridmet_combined <- left_join(gridmet_combined, data_combined, by=c("year","month","grid_id_10km")) 
  } 
  }
  
  saveRDS(gridmet_combined, file.path(data_path, "gridMET_smoke10km_monthly_2006_2020.rds"))
  