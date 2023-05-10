# --------------------------------------------------------
# 0. Pull data and summary plots for model stochasticity/bias etc
# --------------------------------------------------------

testna <- readRDS(file.path(here::here(), "analysis/data_derived/model_s.rds"))
testna <- testna %>% select(c("Micro.2.10","ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det", "s"))

# --------------------------------------------------------
# 1. Train mars model  using best fitting hyperparameters
# --------------------------------------------------------

# xgboost model
library(mda)
library(tidyr)

# Load data
# Set up training and testing data
set.seed(123)
test <- testna %>% na.omit()
train_indices <- sample(nrow(test), nrow(test) * 0.75)
train <- test[train_indices, ]
test <- test[-train_indices, ]

# Function to train and return best fitting model given hyperparams
evaluate_mars_model <- function(params, train, test) {

  # Set k for k-fold cross-validation
  k <- 20

  # Create indices for k-fold cross-validation
  train_l <- nrow(train)
  folds <- replicate(k, sample(train_l, round(params$subsample * train_l), replace = FALSE))

  # Initialize variables to store results
  rmse_list <- numeric(k)
  models <- vector("list", k)
  mars_wrap <- function(x, y, degree) {
    fit <- mda::mars(x, y, degree = degree)
    fit$xin <- x
    fit$yin <- y
    class(fit) <- "mars"
    fit
  }

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
    fit <- mars_wrap(train_x, train_y, degree = params$degree)

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
  return(best_fit)
}

# Perform multivariate interpolation using
best_params <- readRDS("analysis/data_raw/mars_hyperparms.rds")
best_fit <- evaluate_mars_model(best_params, train, test)

# Test the best model
test_x <- as.matrix(test %>% select(-s))
test_y <- test$s
preds <- predict(best_fit, test_x)
rmse <- sqrt(mean((test_y - preds) ^ 2))

# Add this to our model selection object
hrp2_mod <- readRDS(file.path(here::here("analysis/data_derived/ensemble_selection_model.rds")))
hrp2_mod$add_model(best_fit, "mars")
hrp2_mod$add_model_weight(rmse, "mars")

# Save our model out again
saveRDS(hrp2_mod, file.path(here::here("analysis/data_derived/ensemble_selection_model.rds")))


# PDP plots
plots <- lapply(names(test %>% select(-s)), function(x){
  partial(best_fit, train = best_fit$xin, pred.var = x, grid.resolution = 10) %>%
    ggplot2::autoplot()
})
mars_pdp <- cowplot::plot_grid(plotlist = plots)

# turn into feature importance
earth::evimp(earth::mars.to.earth(best_fit) %>% update)
