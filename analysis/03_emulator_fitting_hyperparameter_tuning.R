# --------------------------------------------------------
# 0. Pull data and summary plots for model stochasticity/bias etc
# --------------------------------------------------------

testna <- readRDS(file.path(here::here(), "analysis/data_derived/model_s.rds"))
testna <- testna %>% select(c("Micro.2.10","ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det", "s"))

# --------------------------------------------------------
# 1. Train model to predict hyperparameters for  xgboost
# --------------------------------------------------------

# xgboost model
library(xgboost)
library(tidyr)
library(tidyverse)
library(mlbench)
library(gbm)
library(randomForest)
library(caretEnsemble)
library(elasticnet)

# Load data
# Set up training and testing data
set.seed(123)
test <- testna %>% na.omit()
train_indices <- sample(nrow(test), nrow(test) * 0.75)
train <- test[train_indices, ]
test <- test[-train_indices, ]

# Define hyperparameter combinations to explore
param_grid <- expand.grid(eta = c(0.01, 0.1, 0.3),
                          max_depth = c(6, 8, 10),
                          subsample = c(0.75, 0.85, 0.95),
                          colsample_bytree = c(0.75, 0.85, 0.95))

# Function to train and evaluate a model with given hyperparameters
evaluate_model <- function(params, train, test) {

  # Train xgboost model with cross-validation
  xgb_cv <- xgb.cv(params = params, data = as.matrix(train[, c("Micro.2.10","ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det")]),
                   label = train$s, nfold = 10, verbose = 0, na.rm = TRUE,
                   nrounds = 200)

  # Extract best iteration
  best_iter <- which.min(xgb_cv$evaluation_log$test_rmse_mean)

  # Train final model using best iteration
  xgb_model <- xgboost::xgboost(params = xgb_params,
                                data = as.matrix(train[, c("Micro.2.10", "ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det")]),
                                label = train$s, nrounds = best_iter, verbose = 0)

  preds <- predict(xgb_model, as.matrix(test[, c("Micro.2.10", "ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det")]))
  rmse <- sqrt(mean((test$s - preds) ^ 2))
  return(rmse)
}

# Evaluate all hyperparameter combinations and store results
results <-
  param_grid %>%
  mutate(rmse = map_dbl(
    seq_len(nrow(param_grid)), ~ evaluate_model(
    list(
      objective = "reg:squarederror",
      eta = param_grid$eta[.x],
      max_depth = param_grid$max_depth[.x],
      min_child_weight = 1,
      subsample = param_grid$subsample[.x],
      colsample_bytree = param_grid$colsample_bytree[.x],
      monotone_constraints = c(-1, 1,-1,-1, 1,-1)
    ),
    train,
    test
  )), .progress = TRUE)

# Find the best set of hyperparameters
best_params <- param_grid[which.min(results$rmse), ]
saveRDS(best_params, "analysis/data_raw/xgboost_hyperparms.rds")
saveRDS(results, "analysis/data_raw/xgboost_all_hyperparms.rds")

# --------------------------------------------------------
# 2. Train model to predict hyperparameters for mars
# --------------------------------------------------------

# Define hyperparameter combinations to explore for mars model
param_grid <- expand.grid(subsample = c(0.75, 0.85, 0.95),
                          degree = seq(2:10))

# Function to train and evaluate a model with given hyperparameters
evaluate_mars_model <- function(params, train, test) {

  # Set k for k-fold cross-validation
  k <- 20

  # Create indices for k-fold cross-validation
  train_l <- nrow(train)
  folds <- replicate(k, sample(train_l, round(params$subsample * train_l), replace = FALSE))

  # Initialize variables to store results
  rmse_list <- numeric(k)
  models <- vector("list", k)

  # Perform k-fold cross-validation
  for (i in 1:k) {

    # Split the data into training and validation sets
    train_indices <- folds[, i]
    valid_indices <- setdiff(seq_len(train_l),train_indices)

    train_x <- train[train_indices, ] %>% select(-s)
    train_y <- train$s[train_indices]
    valid_x <- train[valid_indices, ] %>% select(-s)
    valid_y <- train$s[valid_indices]

    # Fit the multivariate interpolation model
    fit <- mda::mars(train_x, train_y, degree = params$degree)

    # Test the fit
    preds <- predict(fit, valid_x)
    rmse <- sqrt(mean((valid_y - preds) ^ 2))

    # Store the results
    rmse_list[i] <- rmse
    models[[i]] <- fit
  }

  # Find the best model based on the lowest RMSE
  best_index <- which.min(rmse_list)
  best_fit <- models[[best_index]]

  # Test the best model
  test_x <- as.matrix(test %>% select(-s))
  test_y <- test$s
  preds <- predict(best_fit, test_x)
  rmse <- sqrt(mean((test_y - preds) ^ 2))

  return(rmse)
}

# Evaluate all hyperparameter combinations and store results
results_mars <-
  param_grid %>%
  mutate(rmse = map_dbl(
    seq_len(nrow(param_grid)), ~ evaluate_mars_model(
      param_grid[.x,],
      train,
      test
    )), .progress = TRUE)

# Find the best set of hyperparameters
best_params <- param_grid[which.min(results_mars$rmse), ]
saveRDS(best_params, "analysis/data_raw/mars_hyperparms.rds")
saveRDS(results_mars, "analysis/data_raw/mars_all_hyperparms.rds")

# TODO: Come back and better plot hyperparams for SI
ggplot(results_mars, aes(degree, rmse, color = as.factor(subsample))) + geom_line() + theme_bw()

# --------------------------------------------------------
# 3. Train model to predict hyperparameters for brnn
# --------------------------------------------------------

# TODO:
