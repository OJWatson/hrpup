# -------------------------------------------------------- #
# 0. Pull data and summary plots for model stochasticity/bias etc ----
# -------------------------------------------------------- #

library(tidyverse)
testna <- readRDS(file.path(here::here(), "analysis/data_derived/model_s.rds"))
testna <- testna %>% select(c("Micro.2.10","ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det", "s"))

# Load data
# Set up training and testing data
set.seed(123)
test <- testna %>% na.omit()
train_indices <- sample(nrow(test), nrow(test) * 0.75)
train <- test[train_indices, ]
test <- test[-train_indices, ]

# -------------------------------------------------------- #
# 1. Train model to predict hyperparameters for  xgboost ----
# -------------------------------------------------------- #

# xgboost model
library(xgboost)
library(tidyr)
library(tidyverse)
library(mlbench)
library(gbm)
library(randomForest)
library(caretEnsemble)
library(elasticnet)

# Define hyperparameter combinations to explore
param_grid <- expand.grid(eta = c(0.01, 0.1, 0.3),
                          max_depth = c(6, 8, 10),
                          subsample = c(0.75, 0.85, 0.95),
                          colsample_bytree = c(0.75, 0.85, 0.95))

# Function to train and evaluate a model with given hyperparameters
evaluate_model <- function(params, train, test, model = TRUE) {

  # Train xgboost model with cross-validation
  xgb_cv <- xgb.cv(params = params, data = as.matrix(train[, c("Micro.2.10","ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det")]),
                   label = train$s, nfold = 20, verbose = 0, na.rm = TRUE,
                   nrounds = 200)

  # Extract best iteration
  best_iter <- which.min(xgb_cv$evaluation_log$test_rmse_mean)

  # Train final model using best iteration
  xgb_model <- xgboost::xgboost(params = params,
                                data = as.matrix(train[, c("Micro.2.10", "ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det")]),
                                label = train$s, nrounds = best_iter, verbose = 0)

  preds <- predict(xgb_model, as.matrix(test[, c("Micro.2.10", "ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det")]))
  rmse <- sqrt(mean((test$s - preds) ^ 2))

  if(!model) {
  return(rmse)
  } else {
    return(xgb_model)
  }
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

# and get the final xgboost model
xgb_model <- evaluate_model(
  list(
    objective = "reg:squarederror",
    eta = best_params$eta[1],
    max_depth = best_params$max_depth[1],
    min_child_weight = 1,
    subsample = best_params$subsample[1],
    colsample_bytree = best_params$colsample_bytree[1],
    monotone_constraints = c(-1, 1,-1,-1, 1,-1)
  ),
  train,
  test,
  model = TRUE
)

# -------------------------------------------------------- #
# 2. Train model to predict hyperparameters for bagEarth ----
# -------------------------------------------------------- #

# Set up parallel for this task
library(doParallel)
cl <- makePSOCKcluster(8)
registerDoParallel(cl)

# set seed again
set.seed(123)
test <- testna %>% na.omit()
train_indices <- sample(nrow(test), nrow(test) * 0.75)
train <- test[train_indices, ]
test <- test[-train_indices, ]

# Define hyperparameter combinations to explore for mars model
# N.B. No way to easily enforce monotonicity here so have kept a simple bagg
train_control<- caret::trainControl(method="cv", number=20, summaryFunction=defaultSummary)
tune_grid <- expand.grid("nprune" = 100, "degree" = 8)
model <- caret::train(x = train %>% select(-s), y = train$s, method="bagEarth",
                      trControl = train_control, tuneGrid = tune_grid, keepX = FALSE,
                      B = 10)

mods <- list(model$finalModel)
names(mods) <- "brnn"
vars <- names(test %>% select(-s))
pdp_dat <- lapply(seq_along(vars), function(var){
  dat <- do.call(rbind, lapply(seq_along(mods), function(i){
    pdp::partial(mods[[i]], pred.var = as.character(vars[var]), train = train %>% select(-s), grid.resolution = 20, type = "regression") %>%
      rename(s = yhat) %>%
      mutate(model = names(mods)[i]) %>%
      setNames(c("x", "s", "model")) %>%
      mutate(var = vars[var])
  }))
})

# visual check
do.call(rbind, pdp_dat) %>% ggplot(aes(x, s)) + geom_line() + facet_wrap(~var, scales = "free")

# and get the final bagEarth model
bagEarth_model <- model$finalModel

# -------------------------------------------------------- #
# 3. Train model to predict hyperparameters for brnn ----
# -------------------------------------------------------- #

# set seed again
set.seed(123)
test <- testna %>% na.omit()
train_indices <- sample(nrow(test), nrow(test) * 0.75)
train <- test[train_indices, ]
test <- test[-train_indices, ]

# Define hyperparameter combinations to explore for brnn model
train_control <- caret::trainControl(method="cv", number=20)

# let's loop over the neurons and check for monotonicity
neurons <- 1:3
mono_checks <- lapply(1:3, function(n) {

  tune_grid <- expand.grid("neurons" = n)
  model <- caret::train(x = train %>% select(-s), y = train$s, method="brnn",
                        trControl = train_control, tuneGrid = tune_grid)

  mods <- list(model$finalModel)
  names(mods) <- "brnn"
  vars <- names(test %>% select(-s))
  pdp_dat <- lapply(seq_along(vars), function(var){
    dat <- do.call(rbind, lapply(seq_along(mods), function(i){
      pdp::partial(mods[[i]], pred.var = as.character(vars[var]), train = train %>% select(-s), grid.resolution = 20, type = "regression") %>%
        rename(s = yhat) %>%
        mutate(model = names(mods)[i]) %>%
        setNames(c("x", "s", "model")) %>%
        mutate(var = vars[var])
    }))
  })

  return(do.call(rbind, pdp_dat) %>%
           group_by(var) %>%
           summarise(mono_dec = (which.min(s) == length(s)),
                     mono_inc = (which.max(s) == length(s))))

})

# only 2 neurons possible
mono_checks

# set seed again to get the same model fit
set.seed(123)
test <- testna %>% na.omit()
train_indices <- sample(nrow(test), nrow(test) * 0.75)
train <- test[train_indices, ]
test <- test[-train_indices, ]
model <- caret::train(x = train %>% select(-s), y = train$s, method="brnn",
                                               trControl = train_control, tuneGrid = expand.grid("neurons" = 1))
model <- caret::train(x = train %>% select(-s), y = train$s, method="brnn",
                      trControl = train_control, tuneGrid = expand.grid("neurons" = 2))

# and get the final brnn model
brnn_model <- model$finalModel

# -------------------------------------------------------- #
# 4. Train model to predict hyperparameters for scam ----
# -------------------------------------------------------- #

# set seed again
set.seed(123)
test <- testna %>% na.omit()
train_indices <- sample(nrow(test), nrow(test) * 0.75)
train <- test[train_indices, ]
test <- test[-train_indices, ]

# Trialled here to identify suitable partial dependence and monotonic relationship
# with term inclusion loosely informed by our bagEarth model
fit <- scam::scam(formula = s ~
                    s(Micro.2.10, bs = "mpd", k = 4) +
                    microscopy.use +
                    fitness +
                    s(ft*Micro.2.10, k = 4) +
                    microscopy.use*Micro.2.10 +
                    rdt.nonadherence*Micro.2.10 +
                    rdt.det*Micro.2.10 +
                    rdt.det*microscopy.use +
                    Micro.2.10*ft*microscopy.use,
                  data = train)

# Visual checks on fit and performance and PDPs clear it as suitable
sqrt(mean((test$s - predict(fit, test)) ^ 2))
plot(predict(fit, test), test$s)
abline(0, 1, col = "red")

mods <- list(fit)
names(mods) <- "scam"
vars <- names(test %>% select(-s))
pdp_dat <- lapply(seq_along(vars), function(var){
  dat <- do.call(rbind, lapply(seq_along(mods), function(i){
    pdp::partial(mods[[i]], pred.var = as.character(vars[var]), train = train %>% select(-s), grid.resolution = 20, type = "regression") %>%
      rename(s = yhat) %>%
      mutate(model = names(mods)[i]) %>%
      setNames(c("x", "s", "model")) %>%
      mutate(var = vars[var])
  }))
})

# visual check
do.call(rbind, pdp_dat) %>% ggplot(aes(x, s)) + geom_line() + facet_wrap(~var, scales = "free")

# Define hyperparameter combinations to explore for scam model
param_grid <- expand.grid(gamma = c(1,2,5,10,25,50,100))

# Function to train and evaluate a model with given hyperparameters
evaluate_scam_model <- function(params, train, test, return_model = FALSE) {

  # Set k for k-fold cross-validation
  k <- 20

  # Create indices for k-fold cross-validation
  train_l <- nrow(train)
  flds <- caret::createFolds(train$s, k = k, list = TRUE, returnTrain = FALSE)

  # Initialize variables to store results
  rmse_list <- numeric(k)
  models <- vector("list", k)

  # Perform k-fold cross-validation
  for (i in 1:k) {

    # Split the data into training and validation sets
    train_indices <- unlist(flds[-i])
    valid_indices <- flds[[i]]

    train_x <- train[train_indices, ]
    valid_x <- train[valid_indices, ]
    valid_y <- train$s[valid_indices]

    # Fit the multivariate interpolation model
    gamma <- params$gamma

    fit <- scam::scam(formula = s ~
                        s(Micro.2.10, bs = "mpd", k = 4) +
                        microscopy.use +
                        fitness +
                        s(ft*Micro.2.10, k = 4) +
                        microscopy.use*Micro.2.10 +
                        rdt.nonadherence*Micro.2.10 +
                        rdt.det*Micro.2.10 +
                        rdt.det*microscopy.use +
                        Micro.2.10*ft*microscopy.use,
                      data = train_x, gamma = gamma)

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

  # return mean rmse or best model
  if(return_model == TRUE){
    return(mean(rmse_list))
  } else {
    return(best_fit)
  }
}

# Evaluate all hyperparameter combinations and store results
results_scam <-
  param_grid %>%
  mutate(rmse = furrr::future_map_dbl(
    seq_len(nrow(param_grid)), ~ evaluate_scam_model(
      param_grid[.x,],
      train,
      test
    ), .progress = TRUE))

# Find the best set of hyperparameters
best_params <- param_grid[which.min(results_scam$rmse), ]
scam_model <- evaluate_scam_model(best_params, train, test, return_model = TRUE)


# -------------------------------------------------------- #
# 5. Partial dependence plot ----
# -------------------------------------------------------- #

# Get all our models and variables to be worked on
mods <- list(xgb_model, brnn_model, scam_model, bagEarth_model)
names(mods) <- c("xgb_model", "brnn_model", "scam_model", "bagEarth_model")
vars <- names(test %>% select(-s))
names(vars) <- c("Microscopy 2-10", "Effective Treatment Seeking", "Microscopy Use",
                 "RDT Nonadherence\n", "Comparative Fitness of \npfhrp2 deleted parasites",
                 "P(RDT+ve if only infected with \npfhp2 deleted parasites)")

# Create our partial dependence
pdp_dat <- lapply(seq_along(vars), function(x){
  dat <- do.call(rbind, lapply(seq_along(mods), function(i){
    pdp::partial(mods[[i]], pred.var = as.character(vars[x]),
                 train = train %>% select(-s), grid.resolution = 20, type = "regression") %>%
      rename(s = yhat) %>%
      mutate(model = names(mods)[i]) %>%
      setNames(c("x", "s", "model"))
  }))
  dat %>% ggplot(aes(x,s,color = model)) +
    geom_line() +
    xlab(names(vars)[x]) +
    ylab("Selection Coeffiecient") +
    ggpubr::theme_pubclean(base_size = 14) +
    theme(axis.line = element_line()) +
    scale_color_viridis_d(name = "Model", end = 0.8) +
    theme(legend.key = element_rect(fill = "white"))

})

# partial dependence shows that gradient boosted tree actually not that great
ggpdp <- cowplot::plot_grid(plotlist = pdp_dat)

# Let's not use that model due to its discontinuities in pdp
mods <- list(brnn_model, scam_model, bagEarth_model)
names(mods) <- c("BRNN", "SCAM", "MARS")

# Create our partial dependence
pdp_dat <- lapply(seq_along(vars), function(x){
  dat <- do.call(rbind, lapply(seq_along(mods), function(i){
    pdp::partial(mods[[i]], pred.var = as.character(vars[x]),
                 train = train %>% select(-s), grid.resolution = 20, type = "regression") %>%
      rename(s = yhat) %>%
      mutate(model = names(mods)[i]) %>%
      setNames(c("x", "s", "model"))
  }))
  dat %>% ggplot(aes(x,s,color = model)) +
    geom_line() +
    xlab(names(vars)[x]) +
    ylab("Selection Coeffiecient") +
    ggpubr::theme_pubclean(base_size = 14) +
    theme(axis.line = element_line()) +
    scale_color_viridis_d(name = "Model", end = 0.8) +
    theme(legend.key = element_rect(fill = "white"),
          text = element_text(family = "Helvetica"))

})

ggpdp <- cowplot::plot_grid(plotlist = pdp_dat)
save_figs("partial_dependence", ggpdp, 14, 8)

# -------------------------------------------------------- #
# 6. Create hrp2_mod ----
# -------------------------------------------------------- #

# Create the hrp2 selection prediciton model
hrp2_mod <- R6_hrp2_mod$new(data = testna)

# brnn
hrp2_mod$add_model(brnn_model, "brnn")
rmse_brnn <- sqrt(mean((test$s - predict(brnn_model, test %>% select(-s))) ^ 2))
hrp2_mod$add_model_weight(rmse_brnn, "brnn")

# bagEarth
hrp2_mod$add_model(bagEarth_model, "bagEarth")
rmse_bagEarth <- sqrt(mean((test$s - predict(bagEarth_model, test %>% select(-s))) ^ 2))
hrp2_mod$add_model_weight(rmse_bagEarth, "bagEarth")

# scam
hrp2_mod$add_model(scam_model, "scam")
rmse_scam <- sqrt(mean((test$s - predict(scam_model, test %>% select(-s))) ^ 2))
hrp2_mod$add_model_weight(rmse_scam, "scam")

saveRDS(hrp2_mod, file.path(here::here("analysis/data_derived/ensemble_selection_model.rds")))

# Create a plot of model predictions plot for SI
pred_brnn <- predict(brnn_model, test %>% select(-s))
pred_bagEarth <- predict(bagEarth_model, test %>% select(-s))
pred_scam <- predict(scam_model, test %>% select(-s))
pred_ensemble <- hrp2_mod$predict_s(test)
pred_df <- rbind(
  data.frame("Prediction" = pred_brnn, "Model" = "BRNN", "Observed" = test$s),
  data.frame("Prediction" = pred_scam, "Model" = "SCAM", "Observed" = test$s),
  data.frame("Prediction" = pred_bagEarth, "Model" = "MARS", "Observed" = test$s),
  data.frame("Prediction" = pred_ensemble, "Model" = "Weighted Ensemble", "Observed" = test$s)
)

# Create comparative observation plots of OvsE
pred_gg <- pred_df %>% ggplot(aes(y = Observed, x = Prediction, color = Model)) +
  geom_point(alpha = 0.2) +
  # geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_bw() +
  facet_wrap(~Model, ncol = 4) +
  coord_equal() +
  MetBrewer::scale_color_met_d("Egypt") +
  guides(color = guide_legend(override.aes = list(alpha = 1)))
save_figs("model_prediction_performance", pred_gg, 14, 4, font_family = "Helvetica")

# Create model performance in terms of classical model performance statistics
mod_perf <- pred_df %>% group_by(Model) %>%
  summarise(RMSE = sqrt(mean((Observed - Prediction)^2)),
            MAE = mean(abs(Observed-Prediction)),
            '1-R2' = 1-summary(lm(Observed~Prediction))$r.squared) %>%
  mutate(Model =
           factor(Model, c("SCAM", "MARS", "BRNN", "Weighted Ensemble"))
  )

# Create the plot for this
mod_perf_gg <- mod_perf %>% pivot_longer(RMSE:`1-R2`) %>%
  ggplot(aes(name, value, fill = Model, group = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  MetBrewer::scale_fill_met_d("Egypt") +
  xlab("Model Statistic") +
  ylab("Value")
save_figs("model_prediction_summary_performance", mod_perf_gg, 6, 4, font_family = "Helvetica")

# -------------------------------------------------------- #
# 7. Create hrp2_mod error models ----
# -------------------------------------------------------- #

hrp2_mod <- readRDS(file.path(here::here("analysis/data_derived/ensemble_selection_model.rds")))

library(tidyverse)
testna <- readRDS(file.path(here::here(), "analysis/data_derived/model_s.rds"))
testna <- testna %>% select(c("EIR","Micro.2.10","ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det", "s", "nmf.multiplier"))

# Load data
# Set up training and testing data
set.seed(123)
test <- testna %>% na.omit()

# get out models
mods <- hrp2_mod$get_models()

# train error models
err_models <- map(seq_along(mods), function(i) {

mod <- mods[[i]]

# make error prediction
pred <-  predict(
  mod,
  test[, c("Micro.2.10", "ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det")]
  )

# work out the per simulation group error
test$error <- abs(pred-test$s)
test <- test %>%
  group_by(EIR, ft, microscopy.use, rdt.nonadherence, fitness, rdt.det, nmf.multiplier) %>%
  summarise(Micro.2.10 = mean(Micro.2.10),
            s = median(s),
            error = sd(error, na.rm = TRUE)) %>%
  na.omit() %>%
  ungroup()

train_indices <- sample(nrow(test), nrow(test) * 0.75)
train <- test[train_indices, ]
test <- test[-train_indices, ]

# now do standard train for a brnn 2 neuron model
train_control <- caret::trainControl(method="cv", number=20)
tune_grid <- expand.grid("neurons" = 2)
err_model <- caret::train(
  x = train %>% select("Micro.2.10", "ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det"),
  y = train$error, method="brnn",
  trControl = train_control, tuneGrid = tune_grid)

# Evaluate performance on test set
pred <- predict(err_model$finalModel, test[, c("Micro.2.10", "ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det")])
rmse_error <- sqrt(mean((test$error - pred) ^ 2))

return(list("err_model" = err_model, "rmse" = rmse_error))

})

# add to our model
for(i in seq_along(mods)) {
  hrp2_mod$add_err_model(err_models[[i]]$err_model, model_name = names(mods)[i])
  hrp2_mod$add_err_model_weight(err_models[[i]]$rmse, model_name = names(mods)[i])
}
saveRDS(hrp2_mod, file.path(here::here("analysis/data_derived/ensemble_selection_model.rds")))


new_mod <- R6_hrp2_mod$new(models = selection_model$get_models(),
                           err_models = selection_model$get_error_models(),
                           model_weights = selection_model$get_model_weights(),
                           err_model_weights = selection_model$get_error_model_weights(),
                           data = testna)
