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
library(fst)
library(foreach)
library(doParallel)
library(doSNOW)
library(stringr)
library(splines2)
library(exactextractr)

rm(list=ls())
gc()

code_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/rcode"
data_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data"
figure_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/figures"
result_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/results/smoke_model"

db_path <- "/Users/mhqiu/BurkeLab Dropbox/Projects/smoke-climate"

setwd(figure_path)
states_map <- map_data("state")

#------------- Load GFED data 
gfed_grid <- readRDS(paste0(data_path, "/GFED4S_NorthAmerica_025deg_grid_EPSG4326.rds"))

gfed_coregion_link <- readRDS(paste0(db_path, "/crosswalk_GFED_NA_CEC_Eco_Level2_3.rds")) %>%
  mutate(intersection_area=as.numeric(intersection_area)) %>%   group_by(gfed_cell_id)  %>%
  arrange(desc(intersection_area)) %>% distinct(gfed_cell_id, .keep_all=T)

gfed_wang_region <- readRDS(paste0(data_path, "/crosswalk_GFED_wang_regions.rds")) %>%
  mutate(region=as.character(region)) %>%
  mutate(region=replace(region, region%in%c("mediterranean california", "southwestern US"), "western forest area"),
         region=replace(region, (region=="northeastern US" & gfed_cell_id<30000), "southeastern US")) %>%
  mutate(region=replace(region, region=="western forest area", "Western US"),
         region=replace(region, region=="southeastern US", "Southeastern US"),
         region=replace(region, region=="northeastern US", "Northeastern US"))

gfed_data <- readRDS(file.path(db_path, "climate-fire_training_data_1997-2021_clean_varname.rds"))
gfed_data <- left_join(gfed_data, gfed_wang_region, by="gfed_cell_id")
gfed_data <- left_join(gfed_data, gfed_coregion_link, by="gfed_cell_id")


#------------- regrid NARR variables to the 1X1 grid to match with CMIP6  -------------
crosswalk_gfed_1deg <- readRDS(paste0(data_path, "/crosswalk_gfed_1deg.rds"))
sample_data <- readRDS("/Users/mhqiu/BurkeLab Dropbox/Data/NARR/GFED_grid/monthly/NARR_vpd_1997.rds") %>%
  filter(!is.na(vpd))

ggplot(filter(crosswalk_gfed_1deg, gfed_cell_id %in% sample_data$gfed_cell_id), aes(x=lon, y=lat)) +
  geom_raster()

path_narr <- "/Users/mhqiu/BurkeLab Dropbox/Data/NARR/GFED_grid/monthly"

narr_var_list <- c("air.sfc", "apcp", "bgrun", "ssrun", "rhum.2m", "soilm", "vpd", "wspd.10m")
var_names <- c("tas", "pr", "bgrun", "ssrun", "hurs", "mrsos", "vpd", "sfcWind")
years <- 1997:2021

for (ii in 1:length(narr_var_list)){
  spec <- narr_var_list[ii]
  spec_name <- var_names[ii]
  print(spec)
  file_list <- list.files(path =path_narr, pattern=paste0("NARR_",spec,"*"), full.names = T)
  
  data_combined <- NULL
  for (ff in file_list){
    data_tmp <- readRDS(ff) %>% filter(gfed_cell_id %in% sample_data$gfed_cell_id)
    colnames(data_tmp)[4] <- "value"
    data_agg <- left_join(data_tmp, crosswalk_gfed_1deg) %>%
      group_by(year, month, grid_1deg_id) %>%
      summarise(value=weighted.mean(value, w=coverage_area, na.rm=T)) %>% ungroup() %>%
      mutate(year = as.numeric(year), month = as.numeric(month))
    
    colnames(data_agg)[4] <- spec_name
    data_combined <- bind_rows(data_combined, data_agg)
  }
  
  saveRDS(data_combined, 
          paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/",spec_name,"narr_1deg_1997_2021.rds"))
}


### generate NARR runoff by combining surface and sub-surface runoff
bgrun <- readRDS(paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/bgrun_narr_1deg_1997_2021.rds")) %>%
   filter(!is.na(bgrun))
ssrun <- readRDS(paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/ssrun_narr_1deg_1997_2021.rds")) %>%
  filter(!is.na(ssrun))
runoff <- full_join(bgrun, ssrun ) %>% mutate(mrro = bgrun + ssrun) %>%
  select(-bgrun, -ssrun) 
saveRDS(runoff, 
        paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/mrro_narr_1deg_1997_2021.rds"))


### check NARR units
narr <- list.files(path = "/Users/mhqiu/BurkeLab Dropbox/Data/NARR/monthly", pattern="*.nc", full.names = T)
nc_open(narr[8])


#------------- De-bias climate variables  -------------
var_long_list <- c("temperature_surface", "precipitation", "RH_surface", "vpd","soil_moisture", "runoff", "wind_speed")
var_short_list <- c("tas", "pr", "hurs", "vpd", "mrsos", "mrro", "sfcWind")

month_seq <-  data.frame(date=seq(as.Date("1997-01-01"), as.Date("2100-12-01"), by="months")) %>%
  mutate(year=as.numeric(substr(date,1,4)), month=as.numeric(substr(date,6,7)),
                                  ndays=lubridate::days_in_month(date)) %>% select(-date)

for (vv in 1: length(var_long_list)){ ##1: length(var_long_list)
var_long <- var_long_list[vv]  
var_short <- var_short_list[vv]  
print(var_short)
# ### check cmip6 units
# nc_list <- list.files(path=paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/",
#                                     var_long, "/sample_netcdf"),  pattern = "*.nc", full.names = T)
# 
# for (ii in 1:length(nc_list)){
# print(ncatt_get(nc_open(nc_list[ii]), var_short ,"units")$value)
# }
file_list <- list.files(path=paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/",
                        var_long, "/NorthAmerica_1deg"),
                       pattern = "*.rds", full.names = T)

cmip_files <- data.frame(filename=file_list) %>%
  mutate(model=str_match(filename, paste0(".*", var_short, "\\_(.*?)\\_.*"))[,2]) %>%
  mutate(scenario=str_match(filename, paste0(".*",model,"\\_(.*?)\\_.*"))[,2])

model_list <- unique(cmip_files$model)
## load the historical data and calculate monthly-mean for each grid cell
narr <- readRDS(paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/",
                       var_short,"_narr_1deg_1997_2021.rds"))

narr_mean <- narr %>% filter(year>1996, year<2015) %>%
  rename(narr_value = !!as.name(var_short)) %>%
  group_by(grid_1deg_id) %>%
  summarise(narr_value = mean(narr_value, na.rm=T)) %>% ungroup()

cmip_debias <- NULL
model_bias_summ <- NULL  
for (mm in model_list){
  print(mm)
  file_list_tmp <- filter(cmip_files, model==mm)
  cmip_history <- readRDS(filter(file_list_tmp, scenario=="historical")$filename) %>%
    filter(year > 1996, year < 2015) %>% 
    rename(cmip_value = !!as.name(var_short))
  
  if (var_short %in% c("mrro", "pr")){
    cmip_history <- left_join(cmip_history, month_seq) %>%
      mutate(cmip_value = cmip_value * 24 *3600)
  }
  
  history_mean <-  cmip_history  %>%
    group_by(grid_1deg_id) %>%
    summarise(cmip_value = mean(cmip_value, na.rm=T)) %>% ungroup()
  
  history_mean <- left_join(history_mean, narr_mean) %>%
    mutate(bias = cmip_value - narr_value)

  model_bias_summ <- bind_rows(model_bias_summ,  
                               data.frame(history_mean[,c("grid_1deg_id", "bias","narr_value")],
                                          model=mm))
  
  ssp_list <- filter(file_list_tmp)
  for (ii in 1:nrow(ssp_list)){
    ssp_tmp <- readRDS(ssp_list[ii,"filename"]) %>% 
      rename(cmip_value = !!as.name(var_short)) %>%
      filter(year<2100)
    if (ssp_list[ii,"scenario"]=="historical"){
      ssp_tmp <- filter(ssp_tmp, year>1996, year<2015)
    } 
    
    if (var_short %in% c("mrro", "pr")){
      ssp_tmp <- left_join(ssp_tmp, month_seq) %>%
        mutate(cmip_value = cmip_value * 24 *3600)
    }
    
    ssp_annual <-  ssp_tmp  %>%
      group_by(grid_1deg_id, year) %>%
      summarise(cmip_value = mean(cmip_value, na.rm=T)) %>% ungroup()
    
    ssp_annual <- left_join(ssp_annual, history_mean[,c("grid_1deg_id", "bias", "narr_value")]) %>%
      filter(!is.na(cmip_value), !is.na(narr_value)) %>%
      mutate(cmip_debias = cmip_value - bias)

    cmip_debias <- bind_rows(cmip_debias,
                            data.frame(model=ssp_list[ii,"model"], scenario=ssp_list[ii,"scenario"], ssp_annual))
  }
}
saveRDS(cmip_debias, paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/",
                            var_long, "/CMIP6_annual_debias_1deg.rds"))

saveRDS(model_bias_summ, paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/",
                            var_long, "/CMIP6_narr_bias_1deg.rds"))
}



### for soil moisture, use NLDAS-2 observation to de-bias CMIP6
path_nldas <- "/Users/mhqiu/BurkeLab Dropbox/Data/NLDAS-2/GFED_grid/monthly"

years <- 1997:2021

file_list <- list.files(path =path_nldas, pattern="NLDAS_VIC0125_M.A_SoilM*", full.names = T)

data_combined <- NULL
for (ff in file_list){
  data_tmp <- readRDS(ff)
  colnames(data_tmp)[4] <- "value"
  data_agg <- left_join(data_tmp, crosswalk_gfed_1deg) %>%
    group_by(year, month, grid_1deg_id) %>%
    summarise(value=weighted.mean(value, w=coverage_area, na.rm=T)) %>% ungroup() %>%
    mutate(year = as.numeric(year), month = as.numeric(month)) %>%
    filter(!is.na(value))
  
  colnames(data_agg)[4] <- "mrsos"
  data_combined <- bind_rows(data_combined, data_agg)
}

saveRDS(data_combined, 
        "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/soilm_nldas2_1deg_1997_2021.rds")

### for runoff, also use NLDAS-2 observation to de-bias CMIP6
path_nldas <- "/Users/mhqiu/BurkeLab Dropbox/Data/NLDAS-2/GFED_grid/monthly"

years <- 1997:2021

file_list <- list.files(path =path_nldas, pattern="NLDAS_VIC0125_M.A_SoilM*", full.names = T)

data_combined <- NULL
for (ff in file_list){
  data_tmp <- readRDS(ff)
  colnames(data_tmp)[4] <- "value"
  data_agg <- left_join(data_tmp, crosswalk_gfed_1deg) %>%
    group_by(year, month, grid_1deg_id) %>%
    summarise(value=weighted.mean(value, w=coverage_area, na.rm=T)) %>% ungroup() %>%
    mutate(year = as.numeric(year), month = as.numeric(month)) %>%
    filter(!is.na(value))
  
  colnames(data_agg)[4] <- "mrsos"
  data_combined <- bind_rows(data_combined, data_agg)
}

saveRDS(data_combined, 
        "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/soilm_nldas2_1deg_1997_2021.rds")


### check nLDAS-2 units
nc_open("/Users/mhqiu/BurkeLab Dropbox/Data/NLDAS-2/NLDAS_VIC0125_M.A202110.020.nc")
# nc_list <- list.files(path=paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/runoff/sample_netcdf"),  pattern = "*.nc", full.names = T)
# print(ncatt_get(nc_open(nc_list[1]), "mrro" ,"units")$value)


  file_list <- list.files(path="/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/soil_moisture/NorthAmerica_1deg",
                          pattern = "*.rds", full.names = T)
  
  cmip_files <- data.frame(filename=file_list) %>%
    mutate(model=str_match(filename, paste0(".*", var_short, "\\_(.*?)\\_.*"))[,2]) %>%
    mutate(scenario=str_match(filename, paste0(".*",model,"\\_(.*?)\\_.*"))[,2])
  
  model_list <- unique(cmip_files$model)
  ## load the historical data and calculate monthly-mean for each grid cell
  nldas <- readRDS("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/soilm_nldas2_1deg_1997_2021.rds")
  
  nldas_mean <- nldas %>% filter(year>1996, year<2015) %>%
    rename(nldas_value = !!as.name(var_short)) %>%
    group_by(grid_1deg_id) %>%
    summarise(nldas_value = mean(nldas_value, na.rm=T)) %>% ungroup()
  
  cmip_debias <- NULL
  model_bias_summ <- NULL  
  for (mm in model_list){
    print(mm)
    file_list_tmp <- filter(cmip_files, model==mm)
    cmip_history <- readRDS(filter(file_list_tmp, scenario=="historical")$filename) %>%
      filter(year > 1996, year < 2015) %>% 
      rename(cmip_value = !!as.name(var_short))
    
    history_mean <-  cmip_history  %>%
      group_by(grid_1deg_id) %>%
      summarise(cmip_value = mean(cmip_value, na.rm=T)) %>% ungroup()
    
    history_mean <- left_join(history_mean, nldas_mean) %>%
      mutate(bias = cmip_value - nldas_value)
    
    model_bias_summ <- bind_rows(model_bias_summ,  
                                 data.frame(history_mean[,c("grid_1deg_id", "bias","nldas_value")],
                                            model=mm))
    
    ssp_list <- filter(file_list_tmp)
    for (ii in 1:nrow(ssp_list)){
      ssp_tmp <- readRDS(ssp_list[ii,"filename"]) %>% 
        rename(cmip_value = !!as.name(var_short)) %>%
        filter(year<2100)
      if (ssp_list[ii,"scenario"]=="historical"){
        ssp_tmp <- filter(ssp_tmp, year>1996, year<2015)
      } 
      
      if (var_short %in% c("mrro", "pr")){
        ssp_tmp <- left_join(ssp_tmp, month_seq) %>%
          mutate(cmip_value = cmip_value * 24 *3600)
      }
      
      ssp_annual <-  ssp_tmp  %>%
        group_by(grid_1deg_id, year) %>%
        summarise(cmip_value = mean(cmip_value, na.rm=T)) %>% ungroup()
      
      ssp_annual <- left_join(ssp_annual, history_mean[,c("grid_1deg_id", "bias", "nldas_value")]) %>%
        filter(!is.na(cmip_value), !is.na(nldas_value)) %>%
        mutate(cmip_debias = cmip_value - bias)
      
      cmip_debias <- bind_rows(cmip_debias,
                               data.frame(model=ssp_list[ii,"model"], scenario=ssp_list[ii,"scenario"], ssp_annual))
    }
  }
  saveRDS(cmip_debias, paste0("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/",
                              var_long, "/CMIP6_annual_debias_1deg_nldas.rds"))
  


  ############ check biases for each variable
