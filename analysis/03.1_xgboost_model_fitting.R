# --------------------------------------------------------
# 0. Pull data and summary plots for model stochasticity/bias etc
# --------------------------------------------------------

testna <- readRDS(file.path(here::here(), "analysis/data_derived/model_s.rds"))
testna <- testna %>% select(c("Micro.2.10","ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det", "s"))

# check what the mae is across stochastic reps as this will want to
# be something we aim to produce in our uncertainty from the final model
testna %>% group_by(Micro.2.10, ft, rdt.nonadherence, microscopy.use, fitness, int) %>%
  summarise(m = median(s, na.rm = TRUE), e = mean(s-m), mae = mean(abs((s-m))), n = n()) %>%
  ggplot(aes(Micro.2.10, e)) +
  geom_point() +
  scale_x_log10()

testna %>% group_by(Micro.2.10,int2) %>%
  summarise(m = median(s, na.rm = TRUE), mse = mean((s-m)^2), n = n()) %>%
  pull(mse) %>% quantile(p = 0.975, na.rm = TRUE)

# --------------------------------------------------------
# 1. Train model to predict selection coefficients
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

# Set up xgboost model
best_params <- readRDS("analysis/data_raw/xgboost_hyperparms.rds")
xgb_params <- list(objective = "reg:squarederror", # Use mean squared error as objective function
                   eta = best_params$eta, # Learning rate
                   max_depth = best_params$max_depth, # Maximum tree depth
                   min_child_weight = 1, # Minimum sum of instance weight needed in a child
                   subsample = best_params$subsample, # Subsample ratio of the training instances
                   colsample_bytree = best_params$colsample_bytree, # Subsample ratio of columns when constructing each tree
                   monotone_constraints = c(-1, 1, -1, -1, 1, -1)) # Impose monotonic constraints for known functional relationships

# Train xgboost model with cross-validation
xgb_cv <- xgb.cv(params = xgb_params, data = as.matrix(train[, c("Micro.2.10","ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det")]),
                 label = train$s, nfold = 20, verbose = 0, na.rm = TRUE,
                 nrounds = 200)

# Extract best iteration
best_iter <- which.min(xgb_cv$evaluation_log$test_rmse_mean)

# Train final model using best iteration
xgb_model <- xgboost::xgboost(params = xgb_params,
                              data = as.matrix(train[, c("Micro.2.10", "ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det")]),
                              label = train$s, nrounds = best_iter, verbose = 0)

# Evaluate performance on test set
pred <- predict(xgb_model, as.matrix(test[, c("Micro.2.10", "ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det")]))
rmse <- sqrt(mean((test$s - pred) ^ 2))
plot(pred, test$s, xlab = "XGBoost Model Predictions of Selection Coefficient", ylab = "Observed Transmission Model Selection Coefficient")
abline(0, 1, col = "red")

# Plot variable importance
importance_matrix <- xgb.importance(feature_names = c("Micro.2.10", "ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det"),
                                    model = xgb_model)
xgb.plot.importance(importance_matrix)

# Check prediction vs median values across reps
testsum <- test %>%
group_by(EIR, ft, rdt.nonadherence, microscopy.use, fitness, rdt.det, int) %>%
  summarise(s = median(s, na.rm = TRUE), Micro.2.10 = mean(Micro.2.10), n = n()) %>%
  select(ft, microscopy.use, rdt.nonadherence, fitness, rdt.det, Micro.2.10, int, n, s)

# Think the few peculiar estimates are likely stochastic variability not being captured fully with 5 reps
# But the model seems correct otherwise
pred <- predict(xgb_model, as.matrix(testsum[, c("Micro.2.10", "ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det")]))
plot(pred, testsum$s)


plots <- lapply(names(test %>% select(-s)), function(x){
  pdp::partial(xgb_model, train = train %>% select(-s), pred.var = x, grid.resolution = 20) %>%
    ggplot2::autoplot()
})

xgb_pdp <- cowplot::plot_grid(plotlist = plots)

# --------------------------------------------------------
# 2. Train model to predict standard devation of model
# --------------------------------------------------------

# first get the full data again
set.seed(123)
test <- testna %>% na.omit()

# predict absolute error using this
pred <- predict(xgb_model, as.matrix(test[, c("Micro.2.10", "ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det")]))
test$error <- abs(pred-test$s)
test <- test %>%
  group_by(EIR, ft, microscopy.use, rdt.nonadherence, fitness, rdt.det, nmf.multiplier) %>%
  summarise(Micro.2.10 = mean(Micro.2.10),
            s = median(s),
            error = sd(error, na.rm = TRUE)) %>%
  na.omit()

# now do standard train
train_indices <- sample(nrow(test), nrow(test) * 0.8)
train <- test[train_indices, ]
test <- test[-train_indices, ]

# Set up xgboost model
xgb_params <- list(objective = "reg:squarederror", # Use mean squared error as objective function
                   eta = 0.1, # Learning rate
                   max_depth = 10, # Maximum tree depth
                   min_child_weight = 1, # Minimum sum of instance weight needed in a child
                   subsample = 0.95, # Subsample ratio of the training instances
                   colsample_bytree = 0.95, # Subsample ratio of columns when constructing each tree
                   monotone_constraints = c(-1, 0, 0, 0, 0, 0)) # error increases with prevalence

# Train xgboost model with cross-validation
xgb_err_cv <- xgb.cv(params = xgb_params,
                     data = as.matrix(train[, c("Micro.2.10", "ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det")]),
                 label = train$error, nfold = 20, verbose = 0, na.rm = TRUE,
                 nrounds = 100)

# Extract best iteration
best_iter <- which.min(xgb_err_cv$evaluation_log$test_rmse_mean)

# Train final model using best iteration
xgb_err_model <- xgboost::xgboost(params = xgb_params,
                              data = as.matrix(train[, c("Micro.2.10", "ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det")]),
                              label = train$error, nrounds = best_iter, verbose = 0)

# Evaluate performance on test set
pred <- predict(xgb_err_model, as.matrix(test[, c("Micro.2.10", "ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det")]))
rmse_error <- sqrt(mean((test$error - pred) ^ 2))
plot(pred, test$error, xlab = "XGBoost Model Predictions of Error", ylab = "Observed XGBoost Error")
abline(0, 1, col = "red")

# Plot variable importance
importance_matrix <- xgb.importance(feature_names = c("Micro.2.10", "ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det"),
                                    model = xgb_err_model)
xgb.plot.importance(importance_matrix)


# --------------------------------------------------------
# 3. Create model that for given parameters produces s and uncertainty interval
# --------------------------------------------------------

hrp2_mod <- R6_hrp2_mod$new(data = testna)
hrp2_mod$add_model(xgb_model, "xbgoost")
hrp2_mod$add_model_weight(rmse, "xbgoost")
hrp2_mod$add_err_model(xgb_err_model, "xbgoost")
hrp2_mod$add_err_model_weight(rmse_error, "xbgoost")
saveRDS(hrp2_mod, file.path(here::here("analysis/data_derived/ensemble_selection_model.rds")))


final <- hrp2_mod$predict(test$Micro.2.10, test$ft, test$microscopy.use, test$rdt.nonadherence, test$fitness, test$rdt.det) %>%
  mutate(obs = test$s) %>%
  mutate(found = data.table::between(obs, smin,smax))

# if working correctly then we should recover 95% of samples.
sum(final$found)/length(final$found)

# Extra Diagnostics here
final %>%
  ggplot(aes(x = obs, y = s, ymin = smin, ymax = smax)) +
  geom_ribbon(alpha = 0.2) +
  geom_point(aes(color = found)) +
  geom_errorbar(data = . %>% filter(!found)) +
  geom_abline(slope = 1, intercept = 0, color = "red")
