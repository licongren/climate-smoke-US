library(terra);library(ncdf4)
library(sf)
library(dplyr)   
library(fixest)
library(tidyverse)
library(purrr)
library(data.table)
library(mlr3verse)
library(future)
library(caret)
library(neuralnet)
library(h2o)
library(RSNNS)
library(doSNOW)

rm(list=ls())
gc()

input_path <- "~/oak_space/mhqiu/smoke_climate/input_data"
output_path <- "~/oak_space/mhqiu/smoke_climate/output"

nthread <- 1
#-------------------------------------------------------------------------------
# Use the H2O framework to perfrom multi-layer NN
# Written by Minghao
# Last edited Mar 2022
#-------------------------------------------------------------------------------
gfed_grid <- readRDS(paste0(input_path, "/GFED4S_NorthAmerica_025deg_grid_EPSG4326.rds"))

gfed_ecoregion_link <- readRDS(paste0(input_path, "/crosswalk_GFED_NA_CEC_Eco_Level2_3.rds")) %>%
  mutate(intersection_area=as.numeric(intersection_area)) %>%   group_by(gfed_cell_id)  %>%
  arrange(desc(intersection_area)) %>%
  mutate(area_weight=intersection_area/sum(intersection_area)) %>% ungroup() # %>% distinct(gfed_cell_id, .keep_all=T)

gfed_region <- readRDS(paste0(input_path, "/crosswalk_GFED_wang_regions.rds"))

gfed_data <- readRDS(file.path(input_path, "climate-fire_training_data_2001-2021_clean_varname.rds")) 

gfed_data <- left_join(gfed_data, gfed_region, by="gfed_cell_id")
gfed_data <- left_join(gfed_data, gfed_ecoregion_link, by="gfed_cell_id")

region_list <- unique(gfed_data$region)
test <- NULL 
train <- NULL 
tuning <- NULL
LOG <- F
nestedCV <- T
final_train <- F
simple_features <- T
LOOCV <- T
feature_normalize <- T

##### Load training data #####
if (LOG){
  #data_full <- readRDS(paste0(output_path, "/final_train_region_FE_log.rds"))
  data_grid <- readRDS(paste0(output_path, "/final_train_region_FE_log_1deg_grid.rds"))
} else{
  #data_full <- readRDS(paste0(output_path, "/final_train_region_FE_level.rds"))
  data_grid <- readRDS(paste0(output_path, "/final_train_region_FE_level_1deg_grid.rds"))
}

#### for model training at the grid cell level
# grid <-  expand.grid(layer1 = c(25, 50, 100, 200),
#                      layer2 = c(50, 100, 200),
#                      layer3 = c(25, 50)
#,learningrate = c(0.001, 0.01, 0.05)

## Set up resampling folds
#num_temporal_folds = length(unique(gfed_data$temporal_fold))
if (LOOCV){
  num_temporal_folds <- 21
}

##### set up the number pf hidden neurons for tuning
tmp <- expand.grid(c(10, 25, 50),
                   c(10, 25, 50)) %>% as.matrix()
# tmp <- expand.grid(c(2, 8, 16, 32),
#                    c(2, 8, 16, 32)) %>% as.matrix()
hidden_neurons <- vector(mode='list', length=nrow(tmp))
for (ii in  1:nrow(tmp)){
  hidden_neurons[[ii]] <- as.vector(tmp[ii,])
}

hyper_params <- list(
  activation = c("Rectifier"), 
  hidden = hidden_neurons,
  epochs = c(50, 200),
  # l1 = c(0, 0.00001, 0.0001), 
  # l2 = c(0, 0.00001, 0.0001),
  rate = c(0.01, 0.001)
  # rate_annealing = c(1e-8, 1e-7, 1e-6),
  # rho = c(0.9, 0.95, 0.99, 0.999),
  # epsilon = c(1e-10, 1e-8, 1e-6, 1e-4),
  # momentum_start = c(0, 0.5),
  # momentum_stable = c(0.99, 0.5, 0),
  # input_dropout_ratio = c(0, 0.1, 0.2),
  # max_w2 = c(10, 100, 1000, 3.4028235e+38)
)

search_criteria <- list(strategy = "RandomDiscrete",
                        #max_models = 100,
                        #max_runtime_secs = 900,
                        stopping_metric="RMSE",  
                        stopping_tolerance = 0.01,
                        stopping_rounds = 3,
                        seed = 100)

#### function for performing loops of nested CV
nestedcv_h2o <- function(fold = fold, data=data){
  
  if (LOOCV){
    data_train_tmp <- dplyr::filter(data, year!=fold) %>%
      mutate(row_ids = row_number()) %>%
      rename(truth = outcome_demean)
    data_test_tmp <- dplyr::filter(data, year==fold) %>% 
      rename(truth = outcome_demean)
    
    set.seed(fold)
    years <- unique(data_train_tmp$year)[sample(20)]
    # Cut into folds
    folds <- cut(seq(1, 20), breaks = 5, labels = F)
    folds <- data.frame(year=years, temporal_fold = folds)
    data_train_tmp <- data_train_tmp %>% select(-temporal_fold) %>% left_join(folds) 
      }
  else{
    data_train_tmp <- dplyr::filter(data, temporal_fold!=fold) %>%
      mutate(row_ids = row_number()) %>%
      rename(truth = outcome_demean)
    data_test_tmp <- dplyr::filter(data, temporal_fold==fold) %>% 
      rename(truth = outcome_demean)
  }
  
  data_train.hex <- as.h2o(data_train_tmp)
  data_test.hex <- as.h2o(data_test_tmp)
  
  grid_id_tmp <- paste(experiment,rrr,fold, sep="_")
  dl_grid <- h2o.grid(algorithm = "deeplearning", 
                      x = var_list,
                      y = "truth",
                      weights_column = "sd_DM",
                      grid_id = grid_id_tmp,
                      training_frame = data_train.hex,
                      #validation_frame = data_test.hex,
                      fold_column = "temporal_fold",
                      hyper_params = hyper_params,
                      regression_stop=1e-2,
                      stopping_metric="RMSE",      ## logloss is directly optimized by Deep Learning
                      stopping_tolerance=1e-2,        ## stop when validation logloss does not improve by >=1% for 2 scoring events
                      stopping_rounds=5,
                      search_criteria = search_criteria
                      )
  
  dl_grid_performance <- h2o.getGrid(grid_id_tmp, sort_by = "rmse")
  
  inner_model <- h2o.getModel(dl_grid_performance@model_ids[[1]])
  
  # Now let's evaluate the model performance on a test set
  data_test_tmp$response  <- h2o.predict(object = inner_model, newdata = data_test.hex) %>% as.vector()
  
  data_train_tmp$response  <- h2o.predict(object = inner_model, newdata = data_train.hex) %>% as.vector()
  
  # aa <- h2o.deeplearning( x = var_list,
  #                   y = "truth",
  #                   weights_column = "sd_DM",
  #                   training_frame = data_train.hex)
  
  tuning <- c(inner_model@params$input[c("hidden","epochs","rate")] %>% unlist(), year=fold)
  rm(dl_grid)
  return(list(test = data_test_tmp, train=data_train_tmp, tuning=tuning))
}


h2o.init(max_mem_size = "200g",
         nthreads = nthread)

if (nestedCV){
  for (experiment in "grid"){  ##"grid"("eco2", "eco3", "regional"
    print(paste("exp:", experiment))
    if (experiment == "regional"){
      data <- data_full$data_region
    }else if (experiment == "eco2"){
      data <- data_full$data_eco2
    }else if (experiment == "eco3"){
      data <- data_full$data_eco3
    }else if (experiment == "grid"){
      #data <- data_full$data_grid 
      data <- data_grid 
    }
    
    for (rrr in region_list[1]){
      print(rrr)
      #if (rrr == "Canada-Alaska"){next}
      data_train <- data %>% filter(region == rrr)
      
      if (rrr %in% c("Canada-Alaska","Mexico")){
        vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
                           "narr_vpd", "narr_runoff")
        if (experiment != "regional"){
          var_list <- c(vlist_climate, paste(rep(vlist_climate, each = 3),
                                             c("forest","crop", "grass"), sep="_"),
                        paste0(vlist_climate, "_temp"),
                        paste0(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff"), "_precip"),
                        paste0(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff"), "_rhum"),
                        paste0(c("narr_wspd",  "narr_vpd", "narr_runoff"), "_wspd"),
                        paste0(c("narr_vpd", "narr_runoff"), "_vpd"),
                        paste0(c("narr_runoff"), "_runoff"))
        } else {
          var_list <- c(vlist_climate,
                        paste0(vlist_climate, "_temp"),
                        paste0(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff"), "_precip"),
                        paste0(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff"), "_rhum"),
                        paste0(c("narr_wspd",  "narr_vpd", "narr_runoff"), "_wspd"),
                        paste0(c("narr_vpd", "narr_runoff"), "_vpd"),
                        paste0(c("narr_runoff"), "_runoff"))
        }
      } 
      
      if (!rrr %in% c("Canada-Alaska","Mexico")){        
        vlist_climate <- c("narr_temp", "narr_precip", "narr_rhum", "narr_wspd",        
                           "narr_vpd", "narr_runoff", "nldas_soilm")
        if (experiment != "regional"){
          var_list <- c(vlist_climate, paste(rep(vlist_climate, each = 3),
                                             c("forest","crop", "grass"), sep="_"),
                        paste0(vlist_climate, "_temp"),
                        paste0(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff", "nldas_soilm"), "_precip"),
                        paste0(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_rhum"),
                        paste0(c("narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_wspd"),
                        paste0(c("narr_vpd", "narr_runoff", "nldas_soilm"), "_vpd"),
                        paste0(c("narr_runoff", "nldas_soilm"), "_runoff"),
                        paste0(c("nldas_soilm"), "_soilm"))
        } else {
          var_list <- c(vlist_climate,
                        paste0(vlist_climate, "_temp"),
                        paste0(c("narr_precip", "narr_rhum", "narr_wspd", "narr_vpd", "narr_runoff", "nldas_soilm"), "_precip"),
                        paste0(c("narr_rhum", "narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_rhum"),
                        paste0(c("narr_wspd",  "narr_vpd", "narr_runoff", "nldas_soilm"), "_wspd"),
                        paste0(c("narr_vpd", "narr_runoff", "nldas_soilm"), "_vpd"),
                        paste0(c("narr_runoff", "nldas_soilm"), "_runoff"),
                        paste0(c("nldas_soilm"), "_soilm"))
        }
      }
      
      if (simple_features){
        var_list <- vlist_climate
      }
      
      if (feature_normalize){
        ###### compute the mean and sd of each variable ###### 
        # data_norm_mean <-  data_train %>% 
        #   summarise_at(var_list, .funs = mean, na.rm=T) %>% t() %>% as.data.frame()
        # data_norm_mean <- data_norm_mean %>% mutate(variable= row.names(data_norm_mean),stat ="mean") %>%
        #   rename(value = V1)
        # data_norm_sd <-  data_train %>% 
        #   summarise_at(var_list, .funs = sd, na.rm=T) %>% t() %>% as.data.frame()
        # data_norm_sd <- data_norm_sd %>% mutate(variable= row.names(data_norm_sd),stat ="sd") %>%
        #   rename(value = V1)
        # 
        # data_norm <- bind_rows(data_norm_mean, data_norm_sd)
        
        data_train <- data_train %>% 
          mutate_at(var_list, .fun=function(x){return((x-mean(x, na.rm=T))/sd(x, na.rm=T))})
      }
      
      #### mannually perfrom nested CV using carat and neuralnetwork
      start_time <- Sys.time()
        test_tmp <- NULL
        train_tmp <- NULL
        tuning_tmp <- NULL
        for (ii in 2001:2001){
          print(ii)
          result_tmp <- nestedcv_h2o(ii, data_train)
          test_tmp <- bind_rows(test_tmp, result_tmp$test)
          train_tmp <- bind_rows(train_tmp, result_tmp$train)
          tuning_tmp <- bind_rows(tuning_tmp, result_tmp$tuning)
             }

        end_time <- Sys.time()
        start_time - end_time
        
        test_tmp <- test_tmp %>% 
          mutate(truth = truth + outcome_mean,
                 response = response + outcome_mean)  
        
        train_tmp <- train_tmp %>% 
          mutate(truth = truth + outcome_mean,
                 response = response + outcome_mean)  
        
        if (LOG){
          test_tmp <- test_tmp %>% mutate(truth=exp(truth), response=exp(response)) 
          train_tmp <- train_tmp %>% mutate(truth=exp(truth), response=exp(response)) 
        }
        
        test_tmp <- test_tmp  %>% 
          mutate(truth = truth * area,
                 response = response * area,
                 region=rrr, experiment = experiment)
        test <- bind_rows(test, test_tmp)
        
        train_tmp <- train_tmp  %>% 
          mutate(truth = truth * area,
                 response = response * area,
                 region=rrr, experiment = experiment)
        train <- bind_rows(train, train_tmp)
        
        tuning_tmp <- tuning_tmp %>% #nn:4 xgboost:7
          mutate(region=rrr, experiment = experiment)
        tuning <- bind_rows(tuning, tuning_tmp)
      }
    }
  
  saveRDS(test, paste0(output_path, "/mlp_level_weighted_2lyr_simple_loocv_h2o_grid.rds"))
  saveRDS(train, paste0(output_path, "/mlp_level_weighted_2lyr_simple_loocv_h2o_grid_train.rds"))
  #saveRDS(tuning, paste0(output_path, "/mlp_log_weighted_2lyr_simple_loocv_h2o_tuning.rds"))
}

h2o.shutdown(prompt = F)

#####

#### Set up the hyper-params for tuning
# grid <-  expand.grid(layer1 = c(1, 2, 4, 8, 16),
#                      layer2 = c(1, 2, 4, 8, 16),
#                      layer3 = 1,
#                      decay = c(0.2,0.5,0.8))

# grid <-  expand.grid(layer1 = c(50, 100, 200),
#                      layer2 = c(50, 100, 200),
#                      layer3 = c(50, 100, 200),
#                      decay = c(0.2,0.5,0.8))

# learning.rate = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2),
# momentum = c(0.1,0.5,0.9)

# inner_model <- caret::train(truth ~ .,
#                             data = data_train_tmp[,c("truth", var_list)],
#                             method = "mlpWeightDecayML",
#                             tuneGrid = grid,
#                             metric = "RMSE",
#                             linOut = T,
#                             maxit = nnn,
#                             #weights = data_train_tmp$sd_DM,
#                             preProc = c("center", "scale"), #good idea to do this with neural nets - your error is due to non scaled data
#                             trControl = trainControl(
#                               method = "cv",
#                               index = inner_sampling,
#                               allowParallel= F,
#                               verboseIter = F))
# proc.time() - ptm
# 
# stopCluster(cl)

# start_time <- Sys.time()
# library(doParallel)
# # create the cluster for caret to use
# cl <- makeCluster(nthread)  #### change the number of clusters based on the configuration of the computing environment
# registerDoSNOW(cl)
# 
# pb <- txtProgressBar(min = 1, max = num_temporal_folds, style = 3)
# progress <- function(n) setTxtProgressBar(pb, n)
# opts <- list(progress = progress)
# 
# result_combine  <- foreach(ii=2001:2021,.options.snow = opts, #1:num_temporal_folds
#                            .packages=c("dplyr","caret","RSNNS", "neuralnet"),
#                            .combine=rbind, .multicombine=TRUE,
#                            .init=list(list())) %dopar% {
#                              nestedcv_h2o(ii, data_train)
#                            }
# test_tmp <- rbindlist(result_combine[,"test"]) %>% as.data.frame()
# train_tmp <- rbindlist(result_combine[,"train"]) %>% as.data.frame()
# tuning_tmp <- rbindlist(result_combine[,"tuning"]) %>% as.data.frame()
# 
# stopCluster(cl)
# end_time <- Sys.time()
# start_time - end_time
# 
# start_time <- Sys.time()
# library(doParallel)
# cl <- makePSOCKcluster(3)
# registerDoParallel(cl)


#nestedcv_caret(1, data_train)
# aa <- nestedcv_caret(2001, data_train)
# bb <- nestedcv_caret(2002, data_train)
# cc <- nestedcv_caret(2003, data_train)
# stopCluster(cl)
# end_time <- Sys.time()
# start_time - end_time