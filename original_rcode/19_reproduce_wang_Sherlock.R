path_project = "./"

library(dplyr)
library(purrr)
library(data.table)
library(mlr3verse)
library(xgboost)
library(future)

nthread = availableCores() - 1

# Choose experiment: reproduce, temporal_folds, or inputs
experiment = "reproduce"
exp_reproduce = (experiment == "reproduce")
exp_temporal_folds = (experiment == "temporal_folds")
exp_inputs = (experiment == "inputs")
stopifnot(exp_reproduce | exp_temporal_folds | exp_inputs)
if (!dir.exists(file.path(path_project, "output", experiment))) {
  dir.create(file.path(path_project, "output", experiment), recursive = T)
}

# Choose whether to train final model
train_final = F
replica_with_optimal_param = F
tempfold_with_optimal_param = T
# ------------------------------------------------------------------------------
# Reproduce Wang et al
# Written by: Jessica Li
# Last edited March 2023
# ------------------------------------------------------------------------------
# Load GFED grid cell IDs
gfed_ids = readRDS("/home/users/mhqiu/oak_space/mhqiu/smoke_climate/input_data/GFED4S_monthly_1997_2021_NorthAmerica_025deg_df.rds") %>% 
  distinct(gfed_cell_id, lon, lat)

# Load our training data frame
train = readRDS(file.path(path_project, "climate-fire_training_data_2001-2021_clean_varname.rds"))

## generate a split on temporal index purely randomly (10-fold)
train$date <- paste(train$year, train$month, sep="-")
date_list <- unique(train$date)
set.seed(100)
date_list <- data.frame(date=sample(date_list , length(date_list))) %>%
  mutate(temporal_fold=c(rep(1:10, 25),1,2))
train <- train %>% select(-temporal_fold) %>%
  left_join(date_list, by="date") %>% select(-date)
num_temporal_folds = length(unique(train$year)) ## for LOOCV #length(unique(train$temporal_fold))

# Load Wang et al training data frame
train_wang = readRDS(file.path(path_project, "Wang et al 2022", "Fire_emission_ML_data_2000_2020.rds")) %>% 
  select(lon = Lon, lat = Lat, year = Year, month = Month, PM2.5_emission, 
         temp, rhum, apcp, uwnd, vwnd, SPEI, FM1000, ERC, VPD, lightn_den, 
         SVD1_NCA, SVD2_NCA, SVD1_RM, SVD2_RM, SVD1_lagged2, SVD2_lagged2, 
         ET, soilm, veg_frac, LAI, fuel_load, fuel_load_nor, 
         starts_with("p_"), slope, elev, GDP, pop2)
num_folds_wang = 10
num_temporal_folds = length(unique(train_wang$year)) ## for LOOCV #length(unique(train$temporal_fold))

# Assign temporal folds to Wang et al data, as closely to our splits as possible
# temporal_folds = train %>% distinct(year, temporal_fold)
# train_wang = train_wang %>% 
#   left_join(temporal_folds %>% mutate(year = ifelse(year == 2021, 2000, year)), 
#             by = "year")
temporal_folds = train %>% distinct(year, month, temporal_fold)
train_wang = train_wang %>%
  left_join(temporal_folds, by=c("year", "month"))

# Set up task
data = train_wang
data_folds = data %>% 
  mutate(row_id = row_number()) %>% 
  select(row_id, lon, lat, year, month, temporal_fold)
task = as_task_regr(data, 
                    target = "PM2.5_emission", 
                    id = experiment)

if (exp_reproduce) {
  if (tempfold_with_optimal_param){
    # task$set_col_roles(
    #   "temporal_fold", 
    #   roles = "group"
    # )$set_col_roles(
    #   c("lon", "lat", "year", "month"), 
    #   remove_from = "feature")
    task$set_col_roles(
      "year", 
      roles = "group"
    )$set_col_roles(
      c("lon", "lat", "temporal_fold", "month"), 
      remove_from = "feature")
  }
  else{
  task$set_col_roles(
    c("lon", "lat", "year", "month", "temporal_fold"),
    remove_from = "feature")
  }
} else if (exp_temporal_folds) {
  task$set_col_roles(
    "temporal_fold", 
    roles = "group"
  )$set_col_roles(
    c("lon", "lat", "year", "month"), 
    remove_from = "feature")
} else if (exp_inputs) {
  task$set_col_roles(
    "temporal_fold", 
    roles = "group"
  )$set_col_roles(
    c("year", "SPEI", "FM1000", "ERC", "lightn_den", "SVD1_NCA", "SVD2_NCA", 
      "SVD1_RM", "SVD2_RM", "SVD1_lagged2", "SVD2_lagged2", "veg_frac", "LAI", 
      "fuel_load", "fuel_load_nor", "GDP", "pop2"), 
    remove_from = "feature")
}

# Set up learner
# learner = lrn("regr.xgboost",
#               eta = to_tune(range(c(0.1, 0.25, 0.35, 0.5, 0.8))),
#               max_depth = to_tune(range(c(3, 6, 9, 12))),
#               subsample = to_tune(range(c(0.4, 0.6, 0.8, 1))),
#               nrounds = to_tune(range(c(100, 150, 200, 250))))
learner = lrn("regr.xgboost",
              eta = 0.35,
              max_depth = 12,
              subsample = 1,
              nrounds = 150)   ####these are the best-performed hyperparameters reported in Wang et al
set_threads(learner, nthread)

# Set up resampling
if (exp_reproduce) {
  #resampling = rsmp("cv", folds = num_folds_wang)
  resampling = rsmp("cv", folds = num_temporal_folds)
} else if (exp_temporal_folds | exp_inputs) {
  outer_resampling = rsmp("cv", folds = num_temporal_folds)
  inner_resampling = rsmp("cv", folds = num_temporal_folds - 1)
  resampling = outer_resampling
}

# Set up measure
measure = msr("regr.mse")

# Set up terminator
terminator = trm("none")

# Set up tuner
tuner = tnr("grid_search")

# Conduct nested resampling
if (exp_temporal_folds | exp_inputs) {
  # Set up autotuner
  at = auto_tuner(method = tuner, 
                  learner = learner, 
                  resampling = inner_resampling, 
                  measure = measure, 
                  terminator = terminator, 
                  store_models = T)
  
  # Run nested resampling
  set.seed(364)
  rr = resample(task = task, 
                learner = at, 
                resampling = outer_resampling, 
                store_models = T)
  
  # Extract inner tuning archives
  archives = extract_inner_tuning_archives(rr)
  
  # Save raw archives
  saveRDS(archives, 
          file.path(path_project, "output", experiment, 
                    "inner_tuning_archives.rds"))
  
  # Extract inner tuning results
  results = extract_inner_tuning_results(rr)
  
  # Link iteration number to outer fold
  crosswalk_iteration_outer = 1:num_temporal_folds %>% map_dfr(function(iter) {
    inner_row_ids = archives[iteration == iter][1, resample_result][[1]]$prediction() %>% 
      as.data.table() %>% 
      pull(row_ids)
    outer_fold = data_folds %>% 
      filter(!(row_id %in% inner_row_ids)) %>% 
      pull(temporal_fold) %>% 
      unique()
    out = data.frame(iteration = iter, 
                     outer_fold = outer_fold)
    return(out)
  })
  saveRDS(crosswalk_iteration_outer, 
          file.path(path_project, "output", experiment, 
                    "crosswalk_iteration_outer_fold.rds"))
  
  # For each outer fold
  for (fold in 1:num_temporal_folds) {
    # Grab the corresponding iteration number
    iter = crosswalk_iteration_outer %>% 
      filter(outer_fold == fold) %>% 
      pull(iteration)
    
    # Save tuning history
    tuning_history = archives[iteration == iter, 
                              .(nrounds, eta, max_depth, subsample, regr.mse)]
    saveRDS(tuning_history, 
            file.path(path_project, "output", experiment, 
                      sprintf("tuning_history_fold_%s.rds", fold)))
    
    # Save optimal hyperparameters
    optimal_hyperparameters = results[iteration == iter, 
                                      .(nrounds, eta, max_depth, subsample, regr.mse)]
    saveRDS(optimal_hyperparameters, 
            file.path(path_project, "output", experiment, 
                      sprintf("optimal_hyperparameters_fold_%s.rds", fold)))
    
    # Save inner model
    inner_model = rr$learners[[iter]]$model$learner$model
    xgb.save(inner_model, 
             file.path(path_project, "output", experiment, 
                       sprintf("inner_model_fold_%s.xgb", fold)))
    
    # Save variable importance
    variable_importance = xgb.importance(model = inner_model)
    saveRDS(variable_importance, 
            file.path(path_project, "output", experiment, 
                      sprintf("variable_importance_fold_%s.rds", fold)))
    
    # Save predictions given by inner model
    predictions = data.table(
      data %>% select(lon, lat, year, month, temporal_fold, truth = PM2.5_emission), 
      response = predict(inner_model, data %>% select(task$feature_names) %>% as.matrix())
    ) %>% 
      left_join(gfed_ids, by = c("lon", "lat")) %>% 
      select(gfed_cell_id, lon, lat, year, month, temporal_fold, truth, response)
    saveRDS(predictions, 
            file.path(path_project, "output", experiment, 
                      sprintf("predictions_fold_%s.rds", fold)))
  }
}


if (replica_with_optimal_param | tempfold_with_optimal_param) {
  # Tune final model
 replica_model <- resample(task, learner, resampling)
 cv_pred <- lapply(replica_model$predictions(), 
                  function(x){as.data.table(x)}) %>% 
  rbindlist() %>%
  left_join(data_folds, by = c("row_ids" = "row_id")) %>% 
  left_join(gfed_ids, by = c("lon", "lat")) %>% 
  select(gfed_cell_id, lon, lat, year, month, truth, response)

# saveRDS(cv_pred, 
#         file.path(path_project, "output",
#                   "wang_temporal_10_cv_pred_optimal_hyperparameters.rds"))
saveRDS(cv_pred, 
        file.path(path_project, "output",
                  "wang_temporal_LOOCV_pred_optimal_hyperparameters.rds"))
  }



if (train_final) {
  # Tune final model
  instance = ti(task = task, 
                learner = learner, 
                resampling = resampling, 
                measure = measure, 
                terminator = terminator, 
                store_models = T)
  set.seed(364)
  tr = tuner$optimize(instance)
  
  # Train final model
  tuned = learner$clone()
  tuned$param_set$values = instance$result_learner_param_vals
  tuned$train(task)
  
  # Save tuning history
  tuning_history = as.data.table(instance$archive)[, .(nrounds, eta, max_depth, subsample, regr.mse)]
  saveRDS(tuning_history, 
          file.path(path_project, "output", experiment, 
                    "tuning_history_outer_model.rds"))
  
  # Save optimal hyperparameters
  optimal_hyperparameters = tr[, .(nrounds, eta, max_depth, subsample, regr.mse)]
  saveRDS(optimal_hyperparameters, 
          file.path(path_project, "output", experiment, 
                    "optimal_hyperparameters_outer_model.rds"))
  
  # Save final model
  outer_model = tuned$model
  xgb.save(outer_model, 
           file.path(path_project, "output", experiment, "outer_model.xgb"))
  
  # Save variable importance
  variable_importance = xgb.importance(model = outer_model)
  saveRDS(variable_importance, 
          file.path(path_project, "output", experiment, 
                    "variable_importance_outer_model.rds"))
  
  # Save predictions given by final model
  predictions = tuned$predict(task) %>% 
    as.data.table() %>% 
    left_join(data_folds, by = c("row_ids" = "row_id")) %>% 
    left_join(gfed_ids, by = c("lon", "lat")) %>% 
    select(gfed_cell_id, lon, lat, year, month, temporal_fold, truth, response)
  saveRDS(predictions, 
          file.path(path_project, "output", experiment, 
                    "predictions_outer_model.rds"))
}
