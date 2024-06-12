library(terra);library(ncdf4)
library(sf)
library(dplyr)   
library(fixest)
library(tidyverse)
library(purrr)
library(data.table)
library(mlr3verse)
library(xgboost)
library(nnet)
library(ranger)
library(future)

rm(list=ls())
gc()

input_path <- "~/oak_space/mhqiu/smoke_climate/input_data"
output_path <- "~/oak_space/mhqiu/smoke_climate/output"

nthread <- 23
#-------------------------------------------------------------------------------
# Use MLR3 models to predict GFED DM emissions at the regional levels
# Written by Minghao
# Last edited Mar 2022
#-------------------------------------------------------------------------------
MODEL <- "NN"
full_feature <- F
fire_season <- F
region_FE <- F
LOG <- T
upweight <- T
feature_normalize <- T
area_normalize <- T

gfed_grid <- readRDS(paste0(input_path, "/GFED4S_NorthAmerica_025deg_grid_EPSG4326.rds"))

gfed_ecoregion_link <- readRDS(paste0(input_path, "/crosswalk_GFED_NA_CEC_Eco_Level2_3.rds")) %>%
 mutate(intersection_area=as.numeric(intersection_area)) %>%   group_by(gfed_cell_id)  %>%
  arrange(desc(intersection_area)) %>%
  mutate(area_weight=intersection_area/sum(intersection_area)) %>% ungroup() # %>% distinct(gfed_cell_id, .keep_all=T)

gfed_region <- readRDS(paste0(input_path, "/crosswalk_GFED_wang_regions.rds"))

gfed_data <- readRDS(file.path(input_path, "climate-fire_training_data_2001-2021_clean_varname.rds")) 

gfed_data <- left_join(gfed_data, gfed_region, by="gfed_cell_id")
gfed_data <- left_join(gfed_data, gfed_ecoregion_link, by="gfed_cell_id")

#### when use the NLDAS feature, remove the following two regions as a sensitivity case due to extreme runoff values
#gfed_data <- filter(gfed_data, (!NA_L3KEY %in% c("6.2.9  Blue Mountains","10.1.2  Columbia Plateau") | year!=2017))

 if (fire_season){
  gfed_data <- filter(gfed_data, month > 4, month < 11)
 }
  
#### choose the options for training the model ####
if (full_feature){
  vlist <- c("barren", "cropland", "forest", "grassland", "lichen", "shrubland", "snowice", "urban", "water","wetland",
             "elevation", "slope",  "narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
             "narr_soilm", "narr_vpd", "narr_runoff",
             "gridmet_erc", "gridmet_fm100","gridmet_fm1000") 
  vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
                     "narr_soilm", "narr_vpd", "narr_runoff","gridmet_erc", "gridmet_fm100","gridmet_fm1000")
} else{
  vlist <- c("barren", "cropland", "forest", "grassland", "lichen", "shrubland", "snowice", "urban", "water","wetland",
             "elevation", "slope",  
             "narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
             "narr_vpd", "narr_runoff",
             "nldas_soilm") 
  vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
                     "narr_vpd", "narr_runoff", 
                     "nldas_soilm")
}

#########  summarise data at different levels ############
data_region <- gfed_data %>%
  mutate(veg_frac = forest * as.numeric(region %in% c("Western US", "Canada-Alaska")) +
                    (forest + grassland + cropland) * as.numeric(!region %in% c("Western US", "Canada-Alaska"))) %>%
  group_by(year, region) %>%
  mutate(DM_kg=sum(DM_kg*area_weight,na.rm=T)) %>%
  mutate(area = sum(intersection_area, na.rm=T)/12) %>%
  mutate(across(vlist_climate, weighted.mean, w=intersection_area*veg_frac)) %>%
  distinct(year, region, .keep_all = T) %>% select(c(vlist_climate, "DM_kg", "temporal_fold", "area")) %>%
  ungroup() %>% mutate(region_sub=region) 

data_eco2 <- gfed_data %>% group_by(year, region, NA_L2KEY) %>%
  mutate(DM_kg=sum(DM_kg*area_weight,na.rm=T)) %>%
  mutate(area = sum(intersection_area, na.rm=T)/12) %>%
  mutate(across(vlist, weighted.mean, w=intersection_area)) %>%
  distinct(year, region, NA_L2KEY, .keep_all = T) %>% select(c(vlist, "DM_kg", "temporal_fold", "area")) %>%
  ungroup() %>% mutate(region_sub=factor(NA_L2KEY)) %>% select(-NA_L2KEY)

data_eco3 <- gfed_data %>% group_by(year, region, NA_L3KEY) %>% 
  mutate(DM_kg=sum(DM_kg*area_weight,na.rm=T)) %>%
  mutate(area = sum(intersection_area, na.rm=T)/12) %>%
  mutate(across(vlist, weighted.mean, w=intersection_area)) %>%
  distinct(year, region, NA_L3KEY, .keep_all = T) %>% select(c(vlist, "DM_kg", "temporal_fold", "area")) %>% 
  ungroup() %>% mutate(region_sub=factor(NA_L3KEY)) %>% select(-NA_L3KEY)

data_grid <- gfed_data %>% group_by(year, gfed_cell_id, region) %>%
  mutate(DM_kg=sum(DM_kg*area_weight,na.rm=T)) %>%
  mutate(area = as.numeric(area)*area_weight) %>%
  mutate(across(vlist, weighted.mean, w=intersection_area)) %>%
  distinct(year, gfed_cell_id, region, .keep_all = T) %>% 
  select(c(vlist, "DM_kg", "temporal_fold","lat","lon","area")) %>%
  ungroup() 


gfed_agg_data <- list(data_region = data_region, data_eco2 = data_eco2, data_eco3 = data_eco3, data_grid = data_grid)
saveRDS(gfed_agg_data, paste0(output_path, "/gfed_agg_data_may_oct.rds"))

################### agg the GFED data to 1deg resolution ####################
crosswalk_1deg_gfed <- readRDS(paste0(input_path, "/crosswalk_gfed_1deg.rds"))
gfed_agg_data <- readRDS(paste0(output_path, "/gfed_agg_data.rds"))
data_grid <- gfed_agg_data$data_grid 
data_grid_1deg <- left_join(data_grid, crosswalk_1deg_gfed[,c("gfed_cell_id", "coverage_area", "grid_1deg_id")])

data_1deg <- data_grid_1deg %>% 
  group_by(year, grid_1deg_id, region) %>%
  mutate(DM_kg=sum(DM_kg,na.rm=T)) %>%
  mutate(across(vlist, weighted.mean, w=area, na.rm=T)) %>%
  mutate(area =sum(area)) %>%
  distinct(year, grid_1deg_id, region, .keep_all = T) %>% 
  select(c(vlist, "DM_kg", "temporal_fold","lat","lon","area")) %>%
  ungroup() 

saveRDS(data_1deg, paste0(output_path, "/gfed_1deg_data.rds"))



##### drop the eco-3 and eco-2 regions with very few grid cells
eco3_summ <- gfed_data %>% 
  group_by(region, NA_L3KEY) %>% summarise(n=length(region)/252, DM_mton=sum(DM_kg)*1e-9) %>%
  arrange(n) %>% ungroup() %>% mutate(region_sub=factor(NA_L3KEY)) %>% select(-NA_L3KEY)
  
eco2_summ <- gfed_data %>% 
  group_by(region, NA_L2KEY) %>% summarise(n=length(region)/252, DM_mton=sum(DM_kg)*1e-9) %>%
  arrange(n) %>% ungroup() %>% mutate(region_sub=factor(NA_L2KEY)) %>% select(-NA_L2KEY)

### only select regions with more than 10 grid cells in the region
data_eco3 <- left_join(data_eco3, eco3_summ) %>% filter(n>10) %>%
  select(-n, -DM_mton)

data_eco2 <- left_join(data_eco2, eco2_summ) %>% filter(n>10) %>%
  select(-n, -DM_mton)

# Set up model
if (MODEL =="XGBOOST"){
  full_learner = lrn("regr.xgboost",
                     eta = to_tune(c(0.01, 0.1, 0.4)),  #0.01, 0.1, 0.4
                     max_depth = to_tune(c(2, 4, 8)), #2, 4, 8
                     nrounds = to_tune(c(50, 100, 200, 400)), #50, 100, 200, 400
                     subsample = to_tune(c(0.2, 0.5, 0.8)), #0.2, 0.5, 0.8
                     lambda   = to_tune(c(3, 5, 10)), #3, 5, 10
                     predict_sets=c("train","test")) #, 200, 400             
} else if (MODEL =="NN"){
  full_learner = lrn("regr.nnet",
                     size = to_tune(c(5, 10, 20, 50, 100)),  #5, 10, 20, 50, 100
                     maxit = to_tune(c(10, 20, 50, 100, 200, 400)), #10, 20, 50, 100, 200, 400
                     MaxNWts = 50000,
                     predict_sets=c("train","test")) #, 200, 400             
} else if (MODEL =="RF"){
  full_learner = lrn("regr.ranger",
                     mtry.ratio = to_tune(seq(0.1, 1, by=0.1)),
                     num.trees = to_tune(c(500, 1000, 2000)),
                     min.node.size = to_tune(c(1,2,4)),
                     regularization.factor = to_tune(c(0.5, 1)),
                     predict_sets=c("train","test"))
}
set_threads(full_learner, nthread)

# Set up measure
measure = msr("regr.mse")
# Set up terminator
terminator = trm("none")
# Set up tuner
tuner = tnr("grid_search")

# Set up resampling
num_temporal_folds = length(unique(gfed_data$temporal_fold))
outer_resampling = rsmp("cv", folds = num_temporal_folds)
inner_resampling = rsmp("cv", folds = num_temporal_folds - 1)

extreme_year <- readRDS(paste0(input_path, "/extreme_vpd_year_upsample.rds"))

#region_list <- unique(gfed_data$region)
region_list <- c("Western US") ##, "Canada-Alaska" "Mexico" "Northeastern US" "Southeastern US" "Western US"     

test <- NULL; train <- NULL   
tuning_combined <- NULL
for (experiment in c("eco3", "eco2", "regional")){ 
  print(paste("exp:", experiment))
  
  if (experiment == "regional"){
    data <- data_region
  }else if (experiment == "eco3"){
    data <- data_eco3
  }else if (experiment == "eco2"){
    data <- data_eco2
  }
  
  if (area_normalize){
    data <- data %>% mutate(outcome = DM_kg / area)
  }
  else{
    data <- data %>% mutate(outcome = DM_kg)
  }
  
  for (rrr in region_list){
    print(rrr)
    for (ww in c(0)){
      print(paste0("weight:",ww))
      data_tmp <- filter(data, region == rrr) %>% select(-DM_kg)
      
      if (LOG){
        data_tmp <- data_tmp %>% filter(outcome>0) %>% mutate(outcome=log(outcome))
      }else{
        data_tmp <- data_tmp %>% filter(outcome>0)
      }
      ###### if normalizing ######
      if (feature_normalize){
        if (experiment == "regional"){
          data_tmp <- data_tmp %>% 
            mutate_at(vlist_climate, .fun=function(x){return((x-mean(x, na.rm=T))/sd(x, na.rm=T))})
        } else{
          data_tmp <- data_tmp %>% 
          mutate_at(vlist, .fun=function(x){return((x-mean(x, na.rm=T))/sd(x, na.rm=T))})
        }
      }
      
      data_tmp <- data_tmp %>% select_if(colSums(!is.na(.)) > 0)
      if (upweight){
        weight <- ww ###how many copies being added
        weight_yr <- filter(extreme_year, region==rrr)$year
        if (ww > 0){
          data_upweight <- filter(data_tmp, year %in% weight_yr) %>% 
            dplyr::slice(rep(1:n(), each = weight)) 
          data_tmp <- bind_rows(data_tmp, data_upweight)
        }
      }

        task = as_task_regr(data_tmp, target = "outcome", id = experiment)
  
        task$set_col_roles(
          "temporal_fold", 
          roles = "group")$set_col_roles(
            c("area","year","region","region_sub"), 
            remove_from = "feature")
        
      data_folds = data_tmp %>% 
        mutate(row_ids = row_number()) %>% 
        select(row_ids, year, temporal_fold, region, region_sub, area)
      
      at = auto_tuner(method = tuner,
                      learner = full_learner,
                      resampling = inner_resampling,
                      measure = measure,
                      terminator = terminator,
                      store_models = T)
      
      # Run nested resampling
      lgr::get_logger("mlr3")$set_threshold("warn")  ## only show warnings while training
      rr = mlr3::resample(task = task,
                          learner = at,
                          resampling = outer_resampling,
                          store_models = T)
      
      # Extract inner tuning archives
      archives = extract_inner_tuning_archives(rr)
      hyperparam <- data.frame(archives[,1:4]) %>% #nn:4 xgboost:7
        mutate(upsample_weight=ww, region=rrr, experiment = experiment)
      tuning_combined <- bind_rows(tuning_combined, hyperparam)
      
      # Save predictions given by inner model
      test_tmp <- lapply(rr$predictions(predict_sets="test"),
                         function(x){as.data.table(x)}) %>%
        rbindlist() %>%
        left_join(data_folds)  %>%
        mutate(upsample_weight=ww, experiment = experiment)
      test <- bind_rows(test_tmp, test)
      
      train_tmp <- lapply(rr$predictions(predict_sets="train"),
                          function(x){as.data.table(x)}) %>%
        rbindlist() %>%
        left_join(data_folds) %>%
        mutate(upsample_weight=ww, experiment = experiment)
      train <- bind_rows(train_tmp, train)
    }
  }
}

test_summ <- test %>% 
  mutate(response = exp(response), truth = exp(truth)) %>%
  mutate(response = response*area, truth = truth*area) %>% 
  mutate(response_positive=response*as.numeric(response >= 0)) %>%
  as.data.frame() %>% 
  group_by(region, region_sub, year, upsample_weight, experiment) %>%
  summarise(response=mean(response, na.rm=T), 
            truth=mean(truth, na.rm=T),
            response_positive=mean(response_positive, na.rm=T)) %>% ungroup() %>%
  mutate(type="test")
saveRDS(test_summ, paste0(output_path, "/west_nestedCV_annual_spatialres_test_", MODEL,"extremerunoff.rds"))

train_summ <- train %>% 
  mutate(response= exp(response), truth = exp(truth)) %>%
  mutate(response = response*area, truth = truth*area) %>% 
  mutate(response_positive=response*as.numeric(response >= 0)) %>%
  as.data.frame() %>% 
  group_by(region, region_sub, year, upsample_weight, experiment) %>%
  summarise(response=mean(response, na.rm=T), 
            truth=mean(truth, na.rm=T),
            response_positive=mean(response_positive, na.rm=T)) %>% ungroup() %>%
  mutate(type="train")
saveRDS(train_summ, paste0(output_path, "/west_nestedCV_annual_spatialres_train_", MODEL,"extremerunoff.rds"))

saveRDS(tuning_combined, paste0(output_path, "/west_nestedCV_annual_spatialres_tuning_", MODEL,"extremerunoff.rds"))


#final_result <- bind_rows(test, train)
#saveRDS(final_result, paste0(output_path, "/west_annual_spatialres_full_results.rds"))
#saveRDS(tuning_combined, paste0(output_path, "/west_nestedCV_annual_spatialres_tuning_final1.rds"))


############ archive: this is to only get the "best model" from the inner-loop
# # Set up autotuner
# at = auto_tuner(method = tuner,
#                 learner = learner,
#                 resampling = inner_resampling,
#                 measure = measure,
#                 terminator = terminator,
#                 store_models = T)
# 
# # Run nested resampling
# set.seed(111)
# lgr::get_logger("mlr3")$set_threshold("warn")  ## only show warnings while training
# rr = mlr3::resample(task = task,
#                     learner = at,
#                     resampling = outer_resampling,
#                     store_models = T)
# 
# # Extract inner tuning archives
# archives = extract_inner_tuning_archives(rr)
# hyperparam <- data.frame(archives[,1:8]) %>%
#   mutate(upsample_weight=ww, region=rrr, experiment=experiment)
# tuning_combined <- bind_rows(tuning_combined, hyperparam)
# 
# # Extract inner tuning results
# # results = extract_inner_tuning_results(rr) %>%
# #   mutate(upsample_weight=ww, region=rrr)
# # tuning_combined <- bind_rows(tuning_combined, results)
# 
# # Link iteration number to outer fold
# crosswalk_iteration_outer = 1:num_temporal_folds %>% map_dfr(function(iter) {
#   inner_row_ids = archives[iteration == iter][1, resample_result][[1]]$prediction() %>%
#     as.data.table() %>%
#     pull(row_ids)
#   outer_fold = data_folds %>%
#     filter(!(row_id %in% inner_row_ids)) %>%
#     pull(temporal_fold) %>%
#     unique()
#   out = data.frame(iteration = iter,
#                    outer_fold = outer_fold)
#   return(out)
# })
# 
# # Save predictions given by inner model
# test_tmp <- lapply(rr$predictions(predict_sets="test"), 
#                    function(x){as.data.table(x)}) %>% 
#   rbindlist() %>%
#   left_join(data_folds, by = c("row_ids" = "row_id"))  %>% 
#   mutate(upsample_weight=ww, type="test", experiment=experiment) 
# 
# test <- bind_rows(test_tmp, test)
# 
# train_tmp <- lapply(rr$predictions(predict_sets="train"), 
#                     function(x){as.data.table(x)}) %>% 
#   rbindlist() %>%
#   left_join(data_folds, by = c("row_ids" = "row_id")) %>% 
#   mutate(upsample_weight=ww, type="train", experiment=experiment) 
# 
# train <- bind_rows(train_tmp, train)

#### manually get all test-set predictions from the inner nested CV model 
### to enable evaluations under differnt metrics   
# for (eta in c(0.01, 0.1, 0.4, 0.8)){  ##, 0.2, 0.5
#   for (max_depth in c(1, 2 ,4, 8)){
#     for (nrounds in c(50, 100, 200, 400)){  #200, 400
#       for (subsample in c(0.25, 0.5, 0.75)){  ##, 0.75, 1
#         for (lambda in c(2, 5, 10)){  #200, 400
#           full_learner = lrn("regr.xgboost",
#                              eta = eta,
#                              max_depth = max_depth,
#                              nrounds = nrounds,
#                              subsample = subsample,
#                              lambda   = lambda,
#                              predict_sets=c("train","test"))   
#           set_threads(full_learner, nthread)
#           
#           # Tune final model
#           inner_model <- resample(task, full_learner, outer_resampling, store_models = T)
#           test_tmp <- lapply(inner_model$predictions(predict_sets="test"), 
#                              function(x){as.data.table(x)}) %>% 
#             rbindlist() %>%
#             left_join(data_folds, by = c("row_ids" = "row_id")) %>%
#             mutate(eta = eta,
#                    max_depth = max_depth,
#                    subsample = subsample,
#                    lambda   = lambda,
#                    nrounds = nrounds) %>% 
#             mutate(upsample_weight=ww, experiment = experiment) 
#           
#           train_tmp <- lapply(inner_model$predictions(predict_sets="train"), 
#                               function(x){as.data.table(x)}) %>% 
#             rbindlist() %>%
#             left_join(data_folds, by = c("row_ids" = "row_id")) %>%
#             mutate(eta = eta,
#                    max_depth = max_depth,
#                    subsample = subsample,
#                    lambda   = lambda,
#                    nrounds = nrounds) %>% 
#             mutate(upsample_weight=ww, experiment = experiment)
#           
#           test <- bind_rows(test_tmp, test)
#           train <- bind_rows(train_tmp, train)
#           
#         }}}}} ## end of the loop to iterate over all hyper-params     