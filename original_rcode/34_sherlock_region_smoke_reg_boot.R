library(terra);library(ncdf4)
library(sf)
library(dplyr)   
library(fixest)
library(tidyverse)
library(foreach)
library(doParallel)
library(doSNOW)
library(splines2)

rm(list=ls())
gc()

input_path <- "~/oak_space/mhqiu/smoke_climate/input_data"
output_path <- "~/oak_space/mhqiu/smoke_climate/output"

smoke_pm <- readRDS(paste0(input_path, "/smokePM_predictions_20060101_20201231.rds"))
smokepm_month <- smoke_pm %>%
  mutate(year=lubridate::year(date), month=lubridate::month(date)) %>%
  group_by(grid_id_10km, year, month) %>% summarise(smoke_month_sum=sum(smokePM_pred, na.rm=T)) %>% ungroup()
month_seq <- data.frame(date=seq(as.Date("2006-01-01"), as.Date("2020-12-01"), by="months"))
month_seq <- month_seq %>% mutate(year=as.numeric(substr(date,1,4)), month=as.numeric(substr(date,6,7)),
                                  ndays=lubridate::days_in_month(date))

smokepm_month <- right_join(smokepm_month, month_seq, by=c("year","month")) %>%
  mutate(smoke_month=smoke_month_sum/ndays) %>% select(-smoke_month_sum, -date, -ndays)

#### if use gfed DM as the fire emission input (main method) The emissions used in the regression models
fire <- readRDS(paste0(output_path, "/smokepm_gfedDM_dailywind_90cone.rds"))  %>% select(-name) # %>% 

#### if use gfed pm2.5 as the fire emission input
#fire <- readRDS(paste0(output_path, "/smokepm_gfedPM25_dailywind_90cone.rds"))  %>% select(-name) # %>% 

#### if use globfire BA as the fire emission input
#fire <- readRDS(paste0(output_path, "/smokepm_globfire_dailywind_90cone.rds"))  %>% select(-name) # %>% 

smokegrid_climate <- readRDS(paste0(input_path, "/crosswalk_10km_grid_US_9region.rds"))
region_list <- unique(smokegrid_climate$region)

################################    pooled regression by regions   ################################
#met_formula <- "+ bSpline(temp_2m, df=4) + bSpline(precep_max, df=4) + bSpline(dewpoint_2m, df=4) + bSpline(PBLH, df=4) + bSpline(surface_pressure, df=4)"
met_formula <- "+ bSpline(tmav, df=4) + bSpline(pr, df=4) + bSpline(ravg, df=4) + bSpline(vs, df=4) "
#met_formula <- "+ poly(tmav, degree=3, raw=TRUE) + poly(pr,  degree=3, raw=TRUE) + poly(ravg,  degree=3, raw=TRUE) + poly(vs,  degree=3, raw=TRUE) "

fml_allwind <- paste("smoke_month ~", 
                     paste(paste0(rep(c("within50km", "dist50_100km", "dist100_200km","dist200_350km", "dist350_500km", 
                                        "dist500_750km", "dist750_1000km","dist1000_1500km","dist1500_2000km","above_2000km"), each=3), c("_up", "_other", "_down")), collapse = " + "))

###### Using boostrap to estimate the uncertainty of the smoke models
nboot <- 40

coef_combined <- NULL
pred_combined <- NULL

  fire_proj_list <- data.frame(
    filename=list.files(
      path = paste0(output_path, "/downscaled_emission/main_10yrs"),pattern="emis_proj_dist_wind_2050_10yrs_ssp370*",
      #path = paste0(output_path, "/downscaled_emission/2046-2055"),pattern="emis_proj_dist_wind*",
      full.names = T))
  
  fire_proj_list$output_name <- str_extract(fire_proj_list$filename, "(?<=dist_wind_).*?(?=\\.rds)")

for (iii in 1:nrow(fire_proj_list)){
    fire_proj <- readRDS(fire_proj_list[iii, "filename"])
    output_name <- fire_proj_list[iii, "output_name"]
    print(output_name)

    ### take GCM over emissions to save space
    fire_proj <- fire_proj %>% group_by(grid_id_10km, year, month, scenario, upwind) %>%
      summarise_at(vars(within50km:above_2000km), .funs = mean, na.rm=T) %>% ungroup() 
      
smoke_proj_combined <- NULL

for (rrr in region_list){
  print(rrr)
  grid_10km_list <- filter(smokegrid_climate, region==rrr)$grid_id_10km#[1:200]
  
  fire_tmp <- filter(fire, grid_id_10km %in% grid_10km_list)
  smokepm_tmp <- filter(smokepm_month, grid_id_10km %in% grid_10km_list)
  
  smoke_fire <- right_join(smokepm_tmp, fire_tmp, by=c("grid_id_10km","year","month")) %>%
    mutate(smoke_month=replace(smoke_month, is.na(smoke_month), 0)) %>% arrange(grid_id_10km,year,month) %>%
    mutate(within200km=within50km+dist50_100km+dist100_200km,
           dist50_200km=dist50_100km+dist100_200km,
           dist200_500km=dist200_350km+dist350_500km,
           dist500_1000km=dist500_750km+dist750_1000km)
  
  smoke_fire$wind <- factor(smoke_fire$upwind, levels = c(1,2,0), labels=c("up","other","down")) 
  
  smoke_fire_nowind <- smoke_fire %>% group_by(grid_id_10km,  year, month, smoke_month) %>%
    summarise_at(vars(within50km:dist500_1000km), .funs = sum, na.rm=T) %>% ungroup() 
  
  smoke_fire_wide <-   smoke_fire %>% 
    pivot_wider(id_cols=grid_id_10km:smoke_month, values_from=within50km: dist500_1000km, names_from=wind) %>% replace(is.na(.), 0)

  ### if using GRIDMET
  gridmet_month <-  readRDS(paste0(input_path,"/gridMET_smoke10km_monthly_2006_2020.rds")) %>%
    filter(grid_id_10km %in% grid_10km_list) %>%
    mutate(year=as.numeric(year), month=as.numeric(month))
  
  smoke_fire_nowind <-  left_join(smoke_fire_nowind, gridmet_month, by=c("grid_id_10km","year","month")) %>%
    filter(!is.na(smoke_month), !is.na(tmav))
  
  smoke_fire_wide <-  left_join(smoke_fire_wide, gridmet_month, by=c("grid_id_10km","year","month")) %>%
    filter(!is.na(smoke_month), !is.na(tmav))
  
  smoke_fire_wide <- smoke_fire_wide %>% mutate(dist1000_1500km=dist1000_1500km_up +  dist1000_1500km_other + dist1000_1500km_down,
                                                dist1500_2000km=dist1500_2000km_up +  dist1500_2000km_other + dist1500_2000km_down,
                                                above_2000km=above_2000km_up +  above_2000km_other + above_2000km_down)
  
  for (bbb in 1:nboot){
   set.seed(bbb+50)
  reg_data_tmp <- smoke_fire_wide#[sample(1:nrow(smoke_fire_wide), replace = T), ]
  
  wind_all <- feols(data = reg_data_tmp, fml = as.formula(paste0(fml_allwind, met_formula, "+ year | month + grid_id_10km")), cluster = "year^month")
  
  coef_tmp <- bind_rows(data.frame(model="wind_all", beta=coef(wind_all), se=se(wind_all), var= rownames(wind_all$coeftable))) %>%
    mutate(var_sub=substr(var,1,4)) %>% filter(!var_sub%in%c("bSpl")) %>% select(-var_sub) %>%
    mutate(region=rrr, bootid=bbb)
  
  coef_combined <- bind_rows(coef_combined, coef_tmp)
  
  pred <- smoke_fire_wide  %>% ungroup() %>%
    mutate(pred_smoke=predict(wind_all, newdata=smoke_fire_wide)) %>% 
    group_by(grid_id_10km, year) %>%
    summarise(pred_smoke = mean(pred_smoke)) %>%
    mutate(model="wind_all",region=rrr, bootid=bbb) ## get predictions with the boot model and full data
  pred_combined <- bind_rows(pred_combined, pred)
  
  ####  applying the model to new projected fire emissions
    era5_tmp <-   readRDS(paste0(input_path,"/era5_10km_grid_monthly_2006_2020.rds")) %>%
      filter(grid_id_10km %in% grid_10km_list) %>%
      group_by(month, grid_id_10km) %>%
      summarise_at(2:10, .funs=mean, na.rm=T) %>%
      ungroup()
    
    gridmet_tmp <-   readRDS(paste0(input_path,"/gridMET_smoke10km_monthly_2006_2020.rds")) %>%
      filter(grid_id_10km %in% grid_10km_list) %>%
      mutate(year=as.numeric(year), month=as.numeric(month)) %>%
      group_by(month, grid_id_10km) %>%
      summarise_at(2:5, .funs=mean, na.rm=T) %>%
      ungroup()
    
    fire_proj_tmp <- filter(fire_proj, grid_id_10km %in% grid_10km_list) 
    
    fire_proj_tmp$wind <- factor(fire_proj_tmp$upwind, levels = c(1,2,0), 
                                 labels=c("up","other","down")) 
    
    fire_proj_tmp <-   fire_proj_tmp %>% 
      pivot_wider(id_cols=grid_id_10km:scenario, values_from=within50km:above_2000km, names_from=wind) %>% 
      replace(is.na(.), 0)
    
    ### if using gridmet
    fire_proj_tmp <- left_join(fire_proj_tmp, gridmet_tmp, by=c("grid_id_10km", "month")) %>%
      filter(!is.na(pr))
    
    smoke_proj_tmp <- fire_proj_tmp %>% ungroup() %>%
      mutate(pred_smoke=predict(wind_all, newdata=fire_proj_tmp)) %>%  
      group_by(grid_id_10km, year, scenario, gcm) %>%
      summarise(pred_smoke = mean(pred_smoke)) %>%
      mutate(model="wind_all", bootid=bbb) 

    smoke_proj_combined <- bind_rows(smoke_proj_combined, smoke_proj_tmp)
  }
  rm("smoke_fire", "smoke_fire_nowind", "smoke_fire_wide")
}
  saveRDS(smoke_proj_combined, paste0(output_path, "/proj_smoke_",output_name, "_smoke_bootsame.rds"))
}


# names(model_list) <- paste(rep(region_list, 2), c("wind_all", "wind_main"), sep="_")
# saveRDS(model_list, paste0(output_path,"/gfedDM_final_model_wind_9region_90cone_era5.rds"))
# saveRDS(r2_combined, paste0(output_path,"/gfedDM_final_r2_wind_9region_90cone_gridmet.rds"))
saveRDS(coef_combined, paste0(output_path,"/gfedDM_final_coef_bootsame_ssp370.rds"))
saveRDS(pred_combined, paste0(output_path,"/gfedDM_final_smoke_pred_bootsame_ssp370.rds"))

