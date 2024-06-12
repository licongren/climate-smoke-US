library(foreach)
library(doParallel)
library(doSNOW)
library(dplyr)
library(tidyr)
library(ncdf4)
library(terra)
library(ncdf4.helpers)
library(sf)

#library(exactextractr)


cmip6_path <- "~/oak_space/CMIP6/fireCO2/North_America_1deg"
#cmip6_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/CMIP6/fireCO2/North_America_1deg"
variable_name <- "fFire"
variable_name_alternate <- "fFire"
file_list <- list.files(path=cmip6_path, pattern = "*.nc", full.names = T)

######################################################################################################################
####################     regrid CMIP6 data to 1x1 cells (thne could be linked to GFED cells)    ######################
######################################################################################################################
file_list_df <- NULL

for (ii in file_list){
  print(ii)
  data_tmp <- nc_open(ii)
  experiment <- ncatt_get(data_tmp, varid=0, attname="experiment_id")$value
  model <- ncatt_get(data_tmp, varid=0, attname="source_id")$value
  variant <- ncatt_get(data_tmp, varid=0, attname="variant_label")$value
  nc_close(data_tmp)
  file_list_df <- bind_rows(file_list_df, data.frame(filename=ii, model, experiment, variant))
}

file_list_df$model_exp <- paste(file_list_df$model, file_list_df$experiment,file_list_df$variant,sep="_") 
model_experiemnt_unique <- unique(file_list_df$model_exp)

raster_1deg <- rast("/home/users/mhqiu/oak_space/mhqiu/smoke_climate/input_data/NorthAmerica_1deg_raster_empty.tif")
#raster_1deg <- rast("/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data/NorthAmerica_1deg_raster_empty.tif")
polygon_1deg <- st_as_sf(as.polygons(raster_1deg), crs=4326)

cmip_to_1deg_cell <- function(date=date, index=index){
  terra_1deg <- terra::vect(polygon_1deg)
  var_raster <- rast(t(data_value[,length(lat):1,index]))
  ext(var_raster) <- c(range(lon)[1]-0.5,range(lon)[2]+0.5,range(lat)[1]-0.5,range(lat)[2]+0.5)
  # ext_namerica <- c(-170, -55, 10, 75)
  # var_raster <- crop(var_raster, ext_namerica)
 extract_1deg <- terra::extract(var_raster, terra_1deg, fun=mean, na.rm=T)
 colnames(extract_1deg) <- c("grid_1deg_id",variable_name)
 extract_1deg$year <- as.numeric(substr(date[index], 1, 4))
 extract_1deg$month <- as.numeric(substr(date[index], 6, 7))
  return(extract_1deg)
}

#model_experiemnt_unique <- model_experiemnt_unique[13:length(model_experiemnt_unique)]
for (mm in model_experiemnt_unique){
  print(mm)
  file_list_tmp <- filter(file_list_df, model_exp==mm)
  
  value_1deg_combined <- NULL
  for (ii in 1:nrow(file_list_tmp)){
  data_tmp <- nc_open(file_list_tmp[ii,"filename"])
  dim_list <- nc.get.dim.names(data_tmp)
  variable_list <- names(data_tmp$var)
  
  if (variable_name %in% variable_list){
    var_name <- variable_name
  }
  else if (variable_name_alternate %in% variable_list){
    var_name <- variable_name_alternate
  }
    
  if ("lat"%in%dim_list){
    lat <-  ncvar_get(data_tmp, "lat")
    lon <-  ncvar_get(data_tmp, "lon")
  }
  else if("latitude"%in%dim_list){
    lat <-  ncvar_get(data_tmp, "latitude")
    lon <-  ncvar_get(data_tmp, "longitude")
  }

 lon <- lon -360  ## shift to -140 ~ -60 (approximately)

  experiment <- ncatt_get(data_tmp, varid=0, attname="experiment_id")$value
  model <- ncatt_get(data_tmp, varid=0, attname="source_id")$value
  variant <- ncatt_get(data_tmp, varid=0, attname="variant_label")$value
  time_tmp <- as.character(nc.get.time.series(f = data_tmp, time.dim.name = "time"))
  nmonths <- length(time_tmp) 
  if (nmonths < 2) {next}
  
  data_value <- ncvar_get(data_tmp, varid=var_name)
  #print(range(data_value))
  #### below is for RH values
  # data_value[data_value<0] <- NA
  # data_value[data_value>100] <- 100
  nc_close(data_tmp)
  
  # var_raster <- rast(t(data_value[,length(lat):1,2]))
  # ext(var_raster) <- c(range(lon)[1]-0.5,range(lon)[2]+0.5,range(lat)[1]-0.5,range(lat)[2]+0.5)
  # plot(var_raster)

  value_tmp_1deg <- NULL
  cl <- makeCluster(15)  #### change the number of clusters based on the configuration of the computing environment
  registerDoSNOW(cl)
  pb <- txtProgressBar(min = 1, max = nmonths, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  value_tmp_1deg  <- foreach(dd=1:nmonths,.options.snow = opts,.packages=c("dplyr","terra","ncdf4"),.combine=rbind) %dopar%
    cmip_to_1deg_cell(date=time_tmp, index=dd)
  stopCluster(cl)

  value_1deg_combined <- bind_rows(value_1deg_combined, value_tmp_1deg)

  }
  output_name <- paste(variable_name,model,experiment,variant,sep="_")
  saveRDS(value_1deg_combined, paste0(cmip6_path,"/", output_name, ".rds"))
}


  
