library(torch)
library(luz)
library(dplyr)
library(ggplot2)
library(tidyr)

#' Malaria Time Series Warm-up Dataset for Torch
#'
#' A dataset class for time series modelling of malaria metrics using torch.
#' Each group (e.g. scenario or simulation replicate) is treated as a separate time series.
#' Warm-up steps are added by repeating the first observation to support autoregressive models.
#'
#' @param df A data frame containing the time series data. Must include the following numeric columns:
#' \itemize{
#'   \item \code{t}: time variable (e.g. months)
#'   \item \code{EIR}, \code{Micro.2.10}, \code{ft}, \code{microscopy.use}, \code{rdt.nonadherence}, \code{fitness}, \code{rdt.det}, \code{rep} (used for grouping)
#'   \item \code{freq} (input variable)
#'   \item \code{micro_2_10}, \code{pcr}, \code{clinical_05}, \code{clinical},
#'         \code{severe_05}, \code{severe}, \code{mortality_05}, \code{mortality_100} (output variables)
#'   \item Optionally: \code{s} if present in the normalisation ranges
#' }
#' These columns are normalised using min-max scaling based on provided or default ranges.
#'
#' @param warmup Integer specifying the number of warm-up steps to prepend by repeating the first input row (default: 24).
#' @param norm_ranges Optional named list specifying min-max normalisation ranges for each variable. If NULL, default ranges are used.
#' @param input_cols Optional vector specifying input_cols names. If NULL, default is are used.
#' @param output_cols Optional vector specifying output_cols names. If NULL, default is are used.
#' @param group_cols Optional vector specifying input_cols group_cols If NULL, default is are used.
#'
#' @return A dataset object compatible with torch dataloaders.
#' Each item returned is a list containing:
#' \describe{
#'   \item{\code{x}}{A tensor of input features with warm-up steps prepended.}
#'   \item{\code{y}}{A tensor of output features (same number of rows as \code{x}).}
#' }
#' @export
MalariaTimeSeriesWarmupDataset <- dataset(
  name = "MalariaTimeSeriesDataset",

  initialize = function(df, warmup = 24,
                        norm_ranges = list(
                          t = c(0.025, 20.074),
                          EIR = c(0, 0.5),
                          s = c(-2, 5),
                          Micro.2.10 = c(0, 1),
                          ft = c(0.01, 1),
                          microscopy.use = c(0, 1),
                          rdt.nonadherence = c(0, 1),
                          fitness = c(0.9, 1.0),
                          rdt.det = c(0, 1),
                          freq = c(0, 1),
                          micro_2_10 = c(0, 1),
                          pcr = c(0, 1),
                          clinical_05 = c(0, 4),
                          clinical = c(0, 1),
                          severe_05 = c(0, 0.04),
                          severe = c(0, 0.02),
                          mortality_05 = c(0, 0.01),
                          mortality_100 = c(0, 0.005)
                        ),
                        input_cols = c('Micro.2.10', 'ft',
                                       'microscopy.use', 'rdt.nonadherence',
                                       'fitness', 'rdt.det', 'freq'),
                        output_cols = c('micro_2_10', 'pcr', 'clinical_05', 'clinical',
                                        'severe_05', 'severe', 'mortality_05', 'mortality_100'),
                        group_cols = c('EIR', 'Micro.2.10', 'ft', 'microscopy.use',
                                       'rdt.nonadherence', 'fitness', 'rdt.det', 'rep')
                        ){

    # Data and warmup
    self$df <- df
    self$warmup <- warmup

    # Normalisation ranges (min, max) for each variable
    self$norm_ranges <- norm_ranges

    # Input cols for each variable
    self$input_cols <- input_cols

    # Input cols for each variable
    self$output_cols <- output_cols

    # Input cols for each variable
    self$group_cols <- group_cols

    # Normalisation function
    self$normalisation_fn <- function(df, norm_ranges){
      # Apply min-max normalisation to each specified column
      for (col in names(norm_ranges)) {
        range <- norm_ranges[[col]]
        df[[col]] <- (df[[col]] - range[1]) / (range[2] - range[1])
      }
      return(df)
    }
    self$df <- self$normalisation_fn(self$df, self$norm_ranges)

    # Re-Normalisation function
    self$renormalisation_fn <- function(df, norm_ranges){
      # Reverse min-max normalisation to each column
      for (col in names(norm_ranges)) {
        range <- norm_ranges[[col]]
        df[[col]] <- df[[col]] * (range[2] - range[1]) + range[1]
      }
      return(df)
    }

    # Group data by scenario and split into a list of data frames
    self$groups <- self$df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(self$group_cols))) %>%
      dplyr::group_split()

    # Warmup creation function
    self$warm_up_fn <- function(group, cols, warmup){

      # Mean from the first real time point
      w_mean <- as.numeric(group[1, cols])

      # Standard deviation from the first 12 time steps (or fewer if needed)
      n_for_stats <- min(12, nrow(group))
      stats_block <- group[1:n_for_stats, cols]
      w_sd <- apply(stats_block, 2, sd)

      # Ensure SD is non-zero and non-NA (for safe sampling)
      w_sd[is.na(w_sd) | w_sd == 0] <- 1e-6

      # Simulate warm-up inputs with N(mean, sd)
      warmup_rows <- matrix(
        rnorm(warmup * length(cols),
              mean = rep(w_mean, each = warmup),
              sd = rep(w_sd, each = warmup)),
        nrow = warmup, byrow = FALSE, dimnames = list(NULL,c(cols))
      )

      # Clamp to [0, 1] range (to match normalised inputs)
      warmup_rows <- pmin(pmax(warmup_rows, 0), 1)

      # Convert real sequence to matrices
      real <- as.matrix(group[, cols])

      # Concatenate warm-up and real inputs
      combined <- rbind(warmup_rows, real)
      return(combined)

    }

  },

  .length = function() {
    # Return number of groups (scenarios)
    length(self$groups)
  },

  .getitem = function(idx) {
    # Extract and order the current group's time series
    group <- self$groups[[idx]] %>% dplyr::arrange(t)

    # Define input and output column names
    input_cols <- self$input_cols
    output_cols <- self$output_cols

    # Get inputs and output sequences
    inputs <- self$warm_up_fn(group, input_cols, self$warmup)
    outputs <- self$warm_up_fn(group, output_cols, self$warmup)

    # Convert to torch tensors
    inputs_tensor <- torch::torch_tensor(inputs, dtype = torch::torch_float())
    outputs_tensor <- torch::torch_tensor(outputs, dtype = torch::torch_float())

    # Return input/output tensors and warm-up length
    list(
      x = inputs_tensor,
      y = outputs_tensor
    )
  }

)

#' GRU-based Neural Network Model for Malaria Time Series Prediction
#'
#' Constructs a GRU (Gated Recurrent Unit) model for time series forecasting using torch.
#' The model includes a linear projection of inputs, a GRU encoder, optional dropout,
#' layer normalisation, and a fully connected output layer. A sigmoid activation
#' is applied at the output to constrain predictions to the [0, 1] range.
#'
#' @param input_size Integer. Number of input features per timestep (default: 7).
#'        Should match the number of input columns passed in each time series step.
#' @param hidden_size Integer. Number of hidden units in each GRU layer (default: 64).
#' @param num_layers Integer. Number of stacked GRU layers (default: 2).
#' @param output_size Integer. Number of output targets predicted at each timestep (default: 8).
#'        Should match the number of output columns in the dataset.
#' @param dropout_prob Numeric. Dropout probability applied after the GRU and before the output layer (default: 0.1).
#'
#' @return An `nn_module` representing the GRU model. Can be instantiated and trained using the torch package.
#'
#' @details
#' The model includes:
#' \itemize{
#'   \item An initial \code{nn_linear} + \code{nn_relu} projection from input to hidden size.
#'   \item A GRU encoder with optional dropout between layers.
#'   \item Layer normalisation applied to the GRU output.
#'   \item Dropout applied before the output projection.
#'   \item A fully connected layer mapping hidden states to outputs.
#'   \item A final sigmoid activation constraining predictions to the [0, 1] range.
#' }
#'
#' @export
GRU_model <- function(input_size = 7,
                      hidden_size = 64,
                      num_layers = 2,
                      output_size = 8,
                      dropout_prob = 0.1) {

  MalariaGRUMixed <- nn_module(
    initialize = function(input_size = input_size,
                          hidden_size = hidden_size,
                          num_layers = num_layers,
                          output_size = output_size,
                          dropout_prob = dropout_prob) {

      # Input projection: maps raw input to hidden size
      self$input_projection <- nn_sequential(
        nn_linear(input_size, hidden_size),
        nn_relu()
      )

      # GRU layers
      self$gru <- nn_gru(
        input_size = hidden_size,
        hidden_size = hidden_size,
        dropout = ifelse(num_layers > 1, dropout_prob, 0.0),
        num_layers = num_layers,
        batch_first = TRUE
      )

      # Fully connected output layer
      self$fc <- nn_linear(hidden_size, output_size)

      # Layer normalisation and dropout
      self$ln <- nn_layer_norm(hidden_size)
      self$dropout <- nn_dropout(dropout_prob)

    },

    forward = function(x) {

      # Apply input projection
      out <- self$input_projection(x)

      # Run GRU
      out <- self$gru(out)[[1]]

      # Apply layer norm and dropout
      out <- self$ln(out)
      out <- self$dropout(out)

      # Project to output size
      out <- self$fc(out)

      # Constrain outputs to [0, 1] via sigmoid
      out <- torch_sigmoid(out)

      return(out)
    }
  )

  return(MalariaGRUMixed)
}


#' Custom Loss Function for Final Timesteps with Optional Bias Scaling
#'
#' Computes a custom loss that combines Mean Squared Error (MSE) over the final timesteps
#' of the sequence with a bias penalty on the mean prediction error per feature.
#' This is useful in forecasting tasks where a warm-up period is excluded from loss calculation,
#' and systematic bias in predictions should be penalised once the model is sufficiently accurate.
#'
#' @param pred A torch tensor of predicted values with shape \code{[batch_size, time_steps, features]}.
#' @param target A torch tensor of target values with the same shape as \code{pred}.
#' @param lambda_bias A numeric scalar controlling the maximum weight of the bias penalty term (default: 0.1).
#' @param min_mse_for_lambda A numeric threshold: bias penalty is only applied when MSE is below this value. If 0, no scaling is applied (default: 0).
#' @param final_steps An integer specifying how many timesteps at the end of the sequence to include in the loss (default: 241).
#'
#' @return A scalar tensor representing the combined loss value.
#'
#' @details
#' The total loss is computed as:
#' \deqn{
#'   \text{MSE}(Y, \hat{Y}) + \lambda \cdot \text{Bias}(Y, \hat{Y})
#' }
#' where:
#' \itemize{
#'   \item MSE is computed over the last \code{final_steps} timesteps.
#'   \item Bias is defined as the mean absolute value of the average prediction error per feature.
#' }
#' If \code{min_mse_for_lambda > 0}, the contribution of \code{lambda_bias} is scaled using a sigmoid function
#' that increases as MSE drops below the threshold.
#'
#' @export
custom_loss_final <- function(pred,
                              target,
                              lambda_bias = 0.1,
                              min_mse_for_lambda = 0,
                              final_steps = 241,
                              micro_2_10_scale = 0) {
  total_timesteps <- pred$size(2)

  # Slice final timesteps
  pred_final <- pred[, (total_timesteps - final_steps + 1):total_timesteps, ]
  target_final <- target[, (total_timesteps - final_steps + 1):total_timesteps, ]

  # Compute MSE loss
  mse_loss <- nnf_mse_loss(pred_final, target_final, reduction = "mean")

  # Compute absolute bias
  error <- pred_final - target_final
  mean_error_per_feature <- torch_mean(error, dim = c(1, 2))
  abs_bias_per_feature <- torch_abs(mean_error_per_feature)
  bias_loss <- torch_mean(abs_bias_per_feature)

  # Apply bias penalty based on MSE
  if (min_mse_for_lambda > 0) {
    # Sigmoid ramp-up as MSE drops below threshold
    scale <- torch_sigmoid(10 * (min_mse_for_lambda - mse_loss))  # sharper transition
    total_loss <- mse_loss + scale * lambda_bias * bias_loss
  } else {
    # Use full bias penalty immediately
    total_loss <- mse_loss + lambda_bias * bias_loss
  }

  # apply additional penalty for not getting micro_2_10 right:
  if (micro_2_10_scale > 0) {
    total_loss <- total_loss + (micro_2_10_scale * nnf_mse_loss(pred_final[,1:5,1], target_final[,1:5,1], reduction = "mean"))  # focus on V1
  }

  return(total_loss)
}


#' Custom Metric: Final Timesteps Mean Squared Error (MSE)
#'
#' Computes the mean squared error (MSE) over the final timesteps of a sequence
#' across all batches, useful for evaluating forecasting accuracy at the end
#' of the prediction window.
#'
#' @return A `luz_metric` object computing average MSE across batches on the final steps.
#'
#' @details
#' During training or evaluation, this metric:
#' \itemize{
#'   \item Extracts the final \code{final_steps} timesteps from each prediction/target pair.
#'   \item Computes the batch-wise MSE over those timesteps.
#'   \item Accumulates the results and returns the average over all batches seen.
#' }
#'
#' @export
custom_metric_final_mse <- function() {

  return(luz::luz_metric(
  abbrev = "Final_MSE",

  initialize = function(final_steps) {
    self$final_steps <- final_steps
  },

  update = function(pred, target) {
    total_timesteps <- pred$size(2)

    # Extract final timesteps
    pred_final <- pred[, (total_timesteps - self$final_steps + 1):total_timesteps, ]
    target_final <- target[, (total_timesteps - self$final_steps + 1):total_timesteps, ]

    # Compute MSE for current batch
    mse_batch <- nnf_mse_loss(pred_final, target_final, reduction = "mean")$item()

    # Accumulate
    self$sum_mse <- (self$sum_mse %||% 0) + mse_batch
    self$n_batches <- (self$n_batches %||% 0) + 1
  },

  compute = function() {
    self$sum_mse / self$n_batches
  }
))

}

#' Custom Metric: Final Timesteps Mean Absolute Bias
#'
#' Computes the average absolute bias (mean absolute error in the mean) across features,
#' based on the final timesteps of the sequence. Useful for measuring systematic
#' over- or under-prediction at the end of the forecast.
#'
#' @return A `luz_metric` object computing average absolute bias across batches.
#'
#' @details
#' During training or evaluation, this metric:
#' \itemize{
#'   \item Extracts the final \code{final_steps} timesteps from each prediction/target pair.
#'   \item Computes the mean signed error across time and batch for each feature.
#'   \item Computes the mean absolute value of these biases (absolute bias).
#'   \item Accumulates and averages across all batches.
#' }
#'
#' This complements the MSE by capturing **systematic** prediction error rather than overall noise.
#'
#' @export
custom_metric_final_bias <- function(){

  return(luz::luz_metric(
  abbrev = "Final_Bias",

  initialize = function(final_steps) {
    self$final_steps <- final_steps
  },

  update = function(pred, target) {
    total_timesteps <- pred$size(2)

    # Extract final timesteps
    pred_final <- pred[, (total_timesteps - self$final_steps + 1):total_timesteps, ]
    target_final <- target[, (total_timesteps - self$final_steps + 1):total_timesteps, ]

    # Compute bias: mean error per feature across batch and time
    error <- pred_final - target_final
    mean_error_per_feature <- torch_mean(error, dim = c(1, 2))
    abs_bias_per_feature <- torch_abs(mean_error_per_feature)
    bias_loss <- torch_mean(abs_bias_per_feature)$item()

    # Accumulate
    self$sum_bias <- (self$sum_bias %||% 0) + bias_loss
    self$n_batches <- (self$n_batches %||% 0) + 1
  },

  compute = function() {
    self$sum_bias / self$n_batches
  }
))

}

#' Train Malaria GRU Model Using Luz and Torch
#'
#' Trains a GRU-based malaria forecasting model using the `luz` high-level API for `torch`.
#' Supports early stopping, learning rate scheduling, checkpointing, and saving training history.
#'
#' @param model A torch `nn_module` representing the GRU model architecture (e.g. output of `GRU_model()`).
#' @param train_dl A dataloader object for the training dataset.
#' @param valid_dl A dataloader object for the validation dataset.
#' @param loss A torch loss function (default: \code{torch::nn_mse_loss()}).
#' @param metrics A list of luz-compatible metrics to evaluate (default: MSE via \code{luz_metric_mse()}).
#' @param epochs Integer. Number of training epochs (default: 100).
#' @param input_size Integer. Number of input features (should match dataset input columns; default: 7).
#' @param output_size Integer. Number of output features (should match dataset input columns; default: 8).
#' @param hidden_size Integer. Number of hidden units in each GRU layer (default: 64).
#' @param num_layers Integer. Number of GRU layers (default: 2).
#' @param dropout_prob Numeric. Dropout probability after GRU layers (default: 0.01).
#' @param max_lr Numeric. Maximum learning rate for the one-cycle scheduler (default: 0.01).
#' @param patience Integer. Number of epochs with no improvement before early stopping (default: 10).
#' @param weight_decay Numeric. Weight decay (L2 penalty) for the Adam optimiser (default: 0).
#' @param cpu Logical. If \code{TRUE}, forces training on CPU even if GPU is available (default: \code{FALSE}).
#' @param output_dir Path to directory where trained model and history will be saved (default: \code{"analysis/impact_analysis/trained_models"}).
#' @param model_name Name used for saving model checkpoint and training history (default: \code{"malaria_gru_luz"}).
#' @param verbose Logical. If \code{TRUE}, prints training progress (default: \code{TRUE}).
#'
#' @return A fitted luz model object containing training results, best model state, and training history.
#'
#' @details
#' The following callbacks are applied during training:
#' \itemize{
#'   \item \strong{Early stopping:} Stops training if validation loss does not improve for \code{patience} epochs.
#'   \item \strong{One-cycle learning rate scheduler:} Adjusts learning rate dynamically over epochs.
#'   \item \strong{Checkpointing:} Saves the best model (by validation loss) to disk as a \code{.pt} file.
#' }
#' A JSON file recording training and validation losses per epoch is also saved.
#'
#' @export
train_malaria_model <- function(model, train_dl, valid_dl,
                                loss = torch::nn_mse_loss(),
                                metrics = list(luz_metric_mse()),
                                epochs = 100,
                                input_size = 7,
                                output_size = 8,
                                hidden_size = 64,
                                num_layers = 2,
                                dropout_prob = 0.01,
                                max_lr = 0.01,
                                patience = 10,
                                weight_decay = 0,
                                cpu = FALSE,
                                output_dir = "analysis/impact_analysis/trained_models",
                                model_name = "malaria_gru_luz",
                                verbose = TRUE) {

  # Ensure output directory exists
  fs::dir_create(output_dir)

  # Train the model with early stopping, learning rate scheduler, and checkpointing
  fitted <- model %>%
    setup(
      loss = loss,
      optimizer = optim_adam,
      metrics = metrics
    ) %>%
    set_opt_hparams(weight_decay = weight_decay) %>%
    set_hparams(
      input_size = input_size,
      hidden_size = hidden_size,
      num_layers = num_layers,
      output_size = output_size,
      dropout_prob = dropout_prob
    ) %>%
    fit(
      train_dl,
      epochs = epochs,
      valid_data = valid_dl,
      verbose = verbose,
      accelerator = accelerator(cpu = cpu),
      callbacks = list(
        luz_callback_early_stopping(patience = patience),
        luz_callback_lr_scheduler(
          lr_one_cycle,
          max_lr = max_lr,
          epochs = epochs,
          steps_per_epoch = length(train_dl),
          call_on = "on_batch_end"
        ),
        luz_callback_model_checkpoint(
          path = file.path(output_dir, paste0(model_name), "/"),
          monitor = "valid_loss",
          save_best_only = FALSE
        )
      )
    )

  # Save training history as JSON
  history <- list(
    train_loss = fitted$records$metrics$train_loss,
    valid_loss = fitted$records$metrics$valid_loss,
    epochs = seq_len(length(fitted$records$metrics$train_loss)),
    best_epoch = fitted$ctx$best_epoch,
    best_val_loss = fitted$ctx$best_value
  )

  jsonlite::write_json(
    history,
    file.path(output_dir, paste0(model_name, "_history.json")),
    pretty = TRUE,
    auto_unbox = TRUE
  )

  return(fitted)
}

#' Resume Training from a Luz Checkpoint for Malaria GRU Model
#'
#' This function resumes training of a `luz` GRU model from a previously saved checkpoint.
#' It allows reconfiguration of model parameters, training duration, optimiser settings,
#' and learning rate scheduling, while preserving the training state from the checkpoint.
#'
#' @param checkpoint_path Path to the saved `.pt` checkpoint file.
#' @param train_dl Torch dataloader for training data.
#' @param valid_dl Torch dataloader for validation data.
#' @param loss Loss Function
#' @param metrics Metrics list
#' @param final_steps Integer. Number of timesteps used in the loss function (default: 241).
#' @param input_dim Number of input features.
#' @param hidden_dim Number of GRU hidden units.
#' @param num_layers Number of GRU layers.
#' @param output_dim Number of output features.
#' @param dropout_prob Dropout probability.
#' @param epochs Number of epochs to continue training.
#' @param max_lr Maximum learning rate for one-cycle scheduler.
#' @param patience Number of epochs with no improvement to trigger early stopping.
#' @param weight_decay Weight decay (L2 regularisation) for the optimiser.
#' @param cpu Logical. If TRUE, trains on CPU; otherwise on GPU if available.
#' @param model_name Optional. Used for naming the new checkpoint file.
#' @param lambda_bias Maximum weight of the bias term in the loss (default: 0.1).
#' @param min_mse_for_lambda MSE threshold before bias term ramps up (default: 1e-4).
#'
#' @return A fitted `luz` model object with training resumed from the checkpoint.
#'
#' @export
resume_malaria_training <- function(checkpoint_path,
                                    train_dl,
                                    valid_dl,
                                    loss,
                                    metrics,
                                    final_steps = 241,
                                    input_size,
                                    hidden_size,
                                    num_layers,
                                    output_size,
                                    dropout_prob,
                                    epochs,
                                    max_lr,
                                    patience,
                                    weight_decay,
                                    cpu = FALSE,
                                    model_name = "malaria_gru_resumed",
                                    lambda_bias = 0.1,
                                    min_mse_for_lambda = 1e-4) {


  # Rebuild the model architecture and training config
  model_definition <- GRU_model() %>%
    setup(
      loss = loss,
      optimizer = optim_adam,
      metrics = metrics
    ) %>%
    set_opt_hparams(weight_decay = weight_decay) %>%
    set_hparams(
      input_size = input_size,
      hidden_size = hidden_size,
      num_layers = num_layers,
      output_size = output_size,
      dropout_prob = dropout_prob
    )

  # Resume training with checkpoint callback
  model_resumed <- model_definition %>%
    fit(
      train_dl,
      epochs = epochs,
      valid_data = valid_dl,
      accelerator = accelerator(cpu = cpu),
      callbacks = list(
        luz_callback_resume_from_checkpoint(path = checkpoint_path),
        luz_callback_early_stopping(patience = patience),
        luz_callback_lr_scheduler(
          lr_one_cycle,
          max_lr = max_lr,
          epochs = epochs,
          steps_per_epoch = length(train_dl),
          call_on = "on_batch_end"
        ),
        luz_callback_model_checkpoint(
          path = file.path(
            "analysis/impact_analysis/trained_models",
            paste0(model_name, "_continued.pt")
          ),
          monitor = "valid_loss",
          save_best_only = FALSE
        )
      ),
      verbose = TRUE
    )

  return(model_resumed)
}


#' Prepare Torch Dataloaders for Malaria Time Series Dataset
#'
#' Splits a malaria time series dataset into training, validation, and test sets,
#' applies warm-up padding, and returns corresponding torch dataloaders.
#'
#' @param df A data frame containing time series input data.
#' @param warmup Integer. Number of warm-up steps to prepend (default: 24).
#' @param batch_size Integer. Batch size for all dataloaders (default: 32).
#' @param train_prop Numeric. Proportion of data to use for training (default: 0.8).
#' @param valid_prop Numeric. Proportion of data to use for validation (default: 0.1).
#' @param train_ids. Integer Vector. Indices for training data. Default = NULL,
#'  which samples based on proportions.
#' @param valid_ids Integer Vector. Indices for validation data. Default = NULL,
#'  which samples based on proprtions
#' @param test_ids Integer Vector. Indices for testing data. Default = NULL,
#'  which samples based on proprtions
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{train_dl}}{Torch dataloader for the training set.}
#'   \item{\code{valid_dl}}{Torch dataloader for the validation set.}
#'   \item{\code{test_dl}}{Torch dataloader for the test set.}
#'   \item{\code{train_ids, valid_ids, test_ids}}{Indices used for each split.}
#' }
#' @export
prepare_malaria_dataloaders <- function(df,
                                        warmup = 24,
                                        batch_size = 32,
                                        train_prop = 0.8,
                                        valid_prop = 0.1,
                                        train_ids = NULL,
                                        valid_ids = NULL,
                                        test_ids = NULL) {
  # Input checks
  stopifnot(train_prop + valid_prop <= 1)

  # Create dataset with warm-up applied
  ds <- MalariaTimeSeriesWarmupDataset(df, warmup = warmup)

  # Create ids if needed
  if(!is.null(train_ids) & !is.null(valid_ids) & !is.null(test_ids)) {} else {

  total_len <- length(ds)
  all_ids <- 1:total_len

  # Sample training IDs
  train_ids <- sample(all_ids, size = floor(train_prop * total_len))

  # Sample validation IDs from remaining
  remaining_ids <- setdiff(all_ids, train_ids)
  valid_ids <- sample(remaining_ids, size = floor(valid_prop * total_len))

  # Remaining goes to test
  test_ids <- setdiff(remaining_ids, valid_ids)

  }

  # Create subset datasets
  train_ds <- dataset_subset(ds, indices = train_ids)
  valid_ds <- dataset_subset(ds, indices = valid_ids)
  test_ds  <- dataset_subset(ds, indices = test_ids)

  # Create dataloaders
  train_dl <- dataloader(train_ds, batch_size = batch_size, shuffle = TRUE)
  valid_dl <- dataloader(valid_ds, batch_size = batch_size)
  test_dl  <- dataloader(test_ds,  batch_size = batch_size)

  return(list(
    train_dl = train_dl,
    valid_dl = valid_dl,
    test_dl  = test_dl,
    train_ids = train_ids,
    valid_ids = valid_ids,
    test_ids  = test_ids
  ))
}


#' Plot Model Predictions vs Ground Truth for a Single Test Group
#'
#' This function takes a trained model and plots its predictions against the true values
#' for a selected group (simulation replicate) from the test dataset. It automatically
#' handles input normalisation reversal, matches the group to original test data, and
#' visualises model performance across multiple malaria indicators.
#'
#' @param trained_model A fitted `luz` model object.
#' @param dataset A `dataloader` object used for testing (e.g. `test_dl`).
#' @param testna The original unnormalised test dataset (`data.frame`).
#' @param index Integer index specifying which test group to evaluate.
#' @param device String indicating whether to use `"cuda"` or `"cpu"` (default: `"cuda"`).
#'
#' @return A `ggplot` object showing prediction vs true lines for each indicator.
#'
#' @examples
#' plot_prediction(trained_model, test_dl, testna, index = 1)
#'
#' @export
plot_prediction <- function(trained_model, dataset, testna, index, device = "cuda") {

  # Extract normalisation ranges, input cols and output cols
  norm_ranges <- dataset$dataset$dataset$norm_ranges
  input_cols <- dataset$dataset$dataset$input_cols
  output_cols <- dataset$dataset$dataset$output_cols

  # Extract the indexed group from the dataset (torch tensors)
  group_data <- dataset$dataset$.getitem(index)
  replicates <- group_data$x

  # Retrieve the unnormalised metadata (for joining with original testna)
  replicate_nm <- dataset$dataset$dataset$groups[[dataset$dataset$indices[index]]]
  replicate_nm <- dataset$dataset$dataset$renormalisation_fn(replicate_nm, norm_ranges)

  # Match original true values for the same scenario (excluding freq)
  all_true <- dplyr::semi_join(
    testna, replicate_nm %>%
      dplyr::select(EIR, ft, microscopy.use, rdt.nonadherence, fitness, rdt.det),
    by = c("EIR", "ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det")
  )

  # Create average freq per timestep to use in prediction (skip normalisation as freq 0-1)
  all_true_median <- all_true %>%
    dplyr::group_by(t) %>%
    dplyr::summarise(across(all_of(input_cols), median))

  all_true_warmup <- dataset$dataset$dataset$warm_up_fn(
    all_true_median, input_cols, dataset$dataset$dataset$warmup
  )

  replicates[, ncol(replicates)] <- all_true_warmup[, "freq"]

  # Prepare input for model
  input_tensor <- replicates$unsqueeze(1)  # [1, time, features]
  device <- if (torch::cuda_is_available() && device == "cuda") torch_device("cuda") else torch_device("cpu")
  input_tensor <- input_tensor$to(device = device)

  # Predict
  trained_model$model$eval()
  preds <- trained_model$model(input_tensor)$squeeze()$detach()$cpu()
  preds <- as.data.frame(torch::as_array(preds))

  # Rename and denormalise predictions
  names(preds) <- output_cols
  preds <- dataset$dataset$dataset$renormalisation_fn(preds, norm_ranges[colnames(preds)])

  # Drop warmup period
  preds <- preds[-seq_len(dataset$dataset$dataset$warmup), ]
  preds$t <- unique(testna$t)
  preds$freq <- tail(all_true_warmup[, "freq"], -dataset$dataset$dataset$warmup)
  preds$rep <- 0

  # Long format: prediction
  pred_long <- preds %>%
    tidyr::pivot_longer(cols = -c(t, rep), names_to = "indicator", values_to = "value") %>%
    dplyr::mutate(type = "prediction")

  # Long format: true
  true_long <- all_true %>%
    dplyr::select(names(preds)) %>%
    tidyr::pivot_longer(cols = -c(t, rep), names_to = "indicator", values_to = "value") %>%
    dplyr::mutate(type = "true")

  # Combine
  final_long <- dplyr::bind_rows(pred_long, true_long)

  # Plot
  ggplot2::ggplot(final_long, ggplot2::aes(t, value, color = type, group = rep)) +
    ggplot2::geom_line(data = dplyr::filter(final_long, type == "true"), alpha = 0.5) +
    ggplot2::geom_line(data = dplyr::filter(final_long, type == "prediction")) +
    lemon::facet_rep_wrap(~stringr::str_to_title(indicator), scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::ylab("Malaria Indicator") +
    ggplot2::xlab("Year") +
    ggplot2::theme(axis.line = ggplot2::element_line(),
                   legend.position = "top") +
    ggplot2::scale_color_discrete(name = "Data Source:")
}

#' Create Comparison Data Frame of Model Predictions vs Ground Truth
#'
#' This function loops over all samples in a test dataloader, uses the trained model
#' to generate predictions, and returns a long-format data frame that includes both
#' the predicted and true values alongside their associated simulation parameters.
#'
#' @param trained_model A trained `luz` model.
#' @param dataset A torch dataloader (e.g., `test_dl`).
#' @param testna The original unnormalised test data frame used for simulation.
#'
#' @return A data frame with columns for simulation parameters, predicted and true values,
#' and the type of data (`"pred"` or `"truth"`).
#'
#' @export
create_comp_df <- function(trained_model, dataset, testna) {
  norm_ranges <- dataset$dataset$dataset$norm_ranges

  reverse_normalisation <- function(norm_values, col_name) {
    range <- norm_ranges[[col_name]]
    norm_values * (range[2] - range[1]) + range[1]
  }

  res_list <- list()

  with_no_grad({
    for (index in seq_len(dataset$dataset$.length())) {
      if (index %% 100 == 0) message(index)

      # Extract group input and move to device
      group_data <- dataset$dataset$.getitem(index)
      input_tensor <- group_data$x$unsqueeze(1)
      device <- if (cuda_is_available()) torch_device("cuda") else torch_device("cpu")
      input_tensor <- input_tensor$to(device = device)

      # produce predictions
      trained_model$model$eval()
      preds <- trained_model$model(input_tensor)$squeeze()$detach()$cpu()
      preds_df <- as.data.frame(torch::as_array(preds))
      names(preds_df) <- dataset$dataset$dataset$output_cols
      preds_df <- tail(preds_df, -dataset$dataset$dataset$warmup)

      # Extract true values
      truth_df <- as.data.frame(torch::as_array(group_data$y$detach()$cpu()))
      names(truth_df) <- names(preds_df)
      truth_df <- tail(truth_df, -dataset$dataset$dataset$warmup)

      # Get parameter values
      pars <- dataset$dataset$dataset$groups[[dataset$dataset$indices[index]]] %>%
        dplyr::select(EIR:freq)

      full <- rbind(
        cbind(pars, truth_df) %>% dplyr::mutate(type = "truth"),
        cbind(pars, preds_df) %>% dplyr::mutate(type = "pred")
      )

      res_list[[index]] <- full
    }
  })

  # Combine and reverse normalisation
  res_all <- dplyr::bind_rows(res_list)
  res_all <- dataset$dataset$dataset$renormalisation_fn(res_all, norm_ranges)

  return(res_all)
}

#' Calculate Evaluation Metrics for Predictions vs Truth
#'
#' This function calculates standard regression metrics (MAE, RMSE, ME, MAPE, SMAPE, correlation)
#' for model predictions grouped by one or more simulation parameters.
#'
#' @param df A data frame with both predicted and true values (output of `create_comp_df()`).
#' @param group_by_params A character vector of column names to group by (e.g. `c("EIR", "ft")`).
#'
#' @return A summarised data frame of metrics per group and outcome indicator.
#'
#' @export
calculate_metrics <- function(df, group_by_params = NULL) {
  df_long <- df %>%
    tidyr::pivot_longer(
      cols = c(micro_2_10, pcr, clinical_05, clinical,
               severe_05, severe, mortality_05, mortality_100),
      names_to = "name", values_to = "value"
    ) %>%
    tidyr::drop_na()

  df_wide <- df_long %>%
    tidyr::pivot_wider(
      id_cols = c(EIR, rep, s, Micro.2.10, ft, microscopy.use,
                  rdt.nonadherence, fitness, rdt.det, t, freq, name),
      names_from = type,
      values_from = value
    ) %>%
    tidyr::drop_na(pred, truth)

  compute_metrics <- function(actual, predicted) {
    mae <- mean(abs(actual - predicted))
    rmse <- sqrt(mean((actual - predicted)^2))
    me <- mean(actual - predicted)
    mape <- mean(abs((actual - predicted) / ifelse(actual == 0, NA, actual)), na.rm = TRUE) * 100
    smape <- mean(2 * abs(actual - predicted) / (abs(actual) + abs(predicted)), na.rm = TRUE) * 100
    correlation <- cor(actual, predicted, use = "complete.obs")

    tibble::tibble(MAE = mae, RMSE = rmse, ME = me,
                   MAPE = mape, SMAPE = smape, Correlation = correlation)
  }

  grouping_cols <- c(group_by_params, "name")

  df_wide %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_cols))) %>%
    dplyr::summarise(metrics = list(compute_metrics(truth, pred)), .groups = "drop") %>%
    tidyr::unnest(metrics)
}


#------------ Fitting ----------------------------


# Load and prepare data
data_path <- file.path(here::here("analysis/impact_analysis/data-derived/model_sims.rds"))
testna <- readRDS(data_path)
testna <- testna %>% dplyr::filter(!is.na(s))

# Set parameters
warmup <- 50
batch_size <- 32

# Create dataloaders and split info using the helper
loaders <- prepare_malaria_dataloaders(
  df = testna,
  warmup = warmup,
  batch_size = batch_size,
  train_prop = 0.8,
  valid_prop = 0.1
)

# Confirm CUDA is available (optional check)
torch::cuda_is_available()

# Determine final_steps based on the sequence length of one group
final_steps <- nrow(loaders$test_dl$dataset$dataset$groups[[1]])

# Define custom loss function with no bias penalty
loss_fn <- function(pred, target) {
  custom_loss_final(
    pred = pred,
    target = target,
    lambda_bias = 0, # can't seem to get the bias to move it better than just MSE
    min_mse_for_lambda = 1.5e-4,
    final_steps = final_steps,
    micro_2_10_scale = 2
  )
}

# Train the model using luz
trained_model <- train_malaria_model(
  model = GRU_model(),
  train_dl = loaders$train_dl,
  valid_dl = loaders$valid_dl,
  loss = loss_fn,
  metrics = list(
    custom_metric_final_mse()(final_steps),
    custom_metric_final_bias()(final_steps)
  ),
  input_size = 7,
  hidden_size = 64,
  num_layers = 3,
  epochs = 100,
  max_lr = 0.02,
  patience = 10,
  weight_decay = 0, # similarly seems to be difficult to get to work
  dropout_prob = 0.1,
  model_name = "malaria_gru_luz_mmse",
  verbose = TRUE
)

# function to prepare model for saving and save
save_model <- function(file, trained_model, train_ds, loaders, t_break, t_end, data_length, data_path){
  trained_model$extra <- list(
    norm_ranges = loaders$train_dl$dataset$dataset$norm_ranges,
    warmup = loaders$train_dl$dataset$dataset$warmup,
    input_cols = loaders$train_dl$dataset$dataset$input_cols,
    output_cols = loaders$train_dl$dataset$dataset$output_cols,
    group_cols = loaders$train_dl$dataset$dataset$group_cols,
    normalisation_fn = loaders$train_dl$dataset$dataset$normalisation_fn,
    renormalisation_fn = loaders$train_dl$dataset$dataset$renormalisation_fn,
    warm_up_fn = loaders$train_dl$dataset$dataset$warm_up_fn,
    data = list(
      train_ids = loaders$train_ids,
      valid_ids = loaders$valid_ids,
      test_ids = loaders$test_ids,
      data_path = data_path,
      t_break = t_break,
      t_end = t_end,
      data_length = data_length
    )
  )

  # strip expensive environments
  trained_model$extra$normalisation_fn <- with(trained_model$extra, eval(parse(text = deparse(normalisation_fn)), envir = globalenv()))
  trained_model$extra$renormalisation_fn <- with(trained_model$extra, eval(parse(text = deparse(renormalisation_fn)), envir = globalenv()))
  trained_model$extra$warm_up_fn <- with(trained_model$extra, eval(parse(text = deparse(warm_up_fn)), envir = globalenv()))

  # save
  luz_save(trained_model, file)
}
t_break <- mean(diff(unique(testna$t)))
t_end <- max(testna$t)
data_length <- length(unique(testna$t))
save_model("analysis/impact_analysis/data-derived/emulator_final.rds",
           trained_model, loaders$train_ds, loaders, t_break, t_end, data_length, data_path)

# Define custom loss function with bias penalty
# loss_fn_extra <- function(pred, target) {
#   custom_loss_final(
#     pred = pred,
#     target = target,
#     lambda_bias = 0.1,
#     min_mse_for_lambda = 1.5e-4,
#     final_steps = final_steps
#   )
# }
# trained_model_extra <- train_malaria_model(
#   model = GRU_model(),
#   train_dl = loaders$train_dl,
#   valid_dl = loaders$valid_dl,
#   loss = loss_fn_extra,
#   metrics = list(
#     custom_metric_final_mse()(final_steps),
#     custom_metric_final_bias()(final_steps)
#   ),
#   input_size = 7,
#   hidden_size = 64,
#   num_layers = 3,
#   epochs = 50,
#   max_lr = 0.02,
#   patience = 5,
#   weight_decay = 1e-7, # similarly seems to be difficult to get to work
#   dropout_prob = 0.1,
#   model_name = "malaria_gru_luz_wd_1e-7",
#   verbose = TRUE
# )
#
# save_model("analysis/impact_analysis/data-derived/emulator_final_wd_1e-7.rds",
#             trained_model, train_ds, loaders, t_break, t_end, data_length, data_path)


#------------ Resuming or loading earlier state if needed ----------------------------

checkpoint_path <- file.path("analysis/impact_analysis/trained_models", "malaria_gru_luz_mmse", "epoch-10-valid_loss-0.000.pt")
checkpoint_path <- file.path("analysis/impact_analysis/trained_models", "malaria_gru_luz_mmse_continued")
model_definition <- GRU_model() %>%
  setup(
    loss = loss_fn,
    optimizer = optim_adam,
    metrics = list(
      custom_metric_final_mse()(final_steps),
      custom_metric_final_bias()(final_steps)
    )
  ) %>%
  set_opt_hparams(weight_decay = 0) %>%
  set_hparams(
    input_size = 7,
    hidden_size = 64,
    num_layers = 3,
    output_size = 8,
    dropout_prob = 0.1
  )
model_resumed <- model_definition %>%
  fit(
    loaders$train_dl,
    epochs = 00,
    valid_data = loaders$valid_dl,
    accelerator = accelerator(cpu = FALSE),
    callbacks = list(
      luz_callback_resume_from_checkpoint(path = checkpoint_path),
      luz_callback_early_stopping(patience = 10),
      luz_callback_lr_scheduler(
        lr_one_cycle,
        max_lr = 0.005,
        epochs = 20,
        steps_per_epoch = length(loaders$train_dl),
        call_on = "on_batch_end"
      ),
      luz_callback_model_checkpoint(
        path = file.path("analysis/impact_analysis/trained_models", "malaria_gru_luz_continued"),
        monitor = "valid_loss",
        save_best_only = FALSE
      )
    ),
    verbose = TRUE
  )


#------------ Test Set Plotting  ----------------------------

trained_model <- luz::luz_load("analysis/impact_analysis/data-derived/emulator_final.rds")
trained_model$model$to(device = torch_device("cuda"))

# have a look at a span of plots
find_equally_spaced_positions <- function(vec, x) {
  stopifnot(is.numeric(vec), x >= 2)

  # Generate x equally spaced target values between min and max
  value_targets <- seq(min(vec), max(vec), length.out = x)

  # For each target value, find the index of the closest match in vec
  index_closest <- sapply(value_targets, function(target) {
    which.min(abs(vec - target))
  })

  # Return unique indices (in case of duplicate matches)
  unique(index_closest)
}
s <- unlist(lapply(loaders$test_dl$dataset$dataset$groups[loaders$test_dl$dataset$indices], function(x){x$s[1]}))
srang_ind <- find_equally_spaced_positions(s, 12)

grid_plot <- function(mod, n, test_dl) {
  s <- unlist(lapply(test_dl$dataset$dataset$groups[test_dl$dataset$indices], function(x){x$s[1]}))
  srang_ind <- tail(find_equally_spaced_positions(s, n+3), n)
  ggs <- lapply(srang_ind, function(x){plot_prediction(mod, test_dl, testna, index = x)})
  cowplot::plot_grid(plotlist = ggs, ncol = 3)
}

plot_prediction(trained_model, loaders$test_dl, testna, index = srang_ind[1])
plot_prediction(trained_model, loaders$test_dl, testna, index = srang_ind[2])
plot_prediction(trained_model, loaders$test_dl, testna, index = srang_ind[3])
plot_prediction(trained_model, loaders$test_dl, testna, index = srang_ind[4])
plot_prediction(trained_model, loaders$test_dl, testna, index = srang_ind[5])

pdf(file = "analysis/impact_analysis/plots_emulator/emulator_fits.pdf", width = 12, height = 8)
for(i in find_equally_spaced_positions(s, 100)){
  print(plot_prediction(trained_model, loaders$test_dl, testna, index = i))
}
dev.off()

#------------ Metric Evaluation ----------------------------

# 1. Generate prediction vs truth data for all test runs
comp <- create_comp_df(trained_model, loaders$test_dl, testna)

# 2. Calculate grouped metrics (e.g. by EIR and ft)
mets <- calculate_metrics(comp, c("EIR", "ft"))

# 3. Plot metrics across EIR and ft
gg <- mets %>%
  tidyr::pivot_longer(MAE:Correlation, names_to = "metric") %>%
  ggplot2::ggplot(aes(EIR, value, color = as.factor(ft))) +
  ggplot2::geom_line() +
  ggplot2::facet_wrap(name ~ metric, scales = "free", ncol = 6) +
  ggplot2::theme_minimal()

save_figs("emulator_metrics_EIR_ft", gg, width = 10, height = 14, plot_dir = "analysis/impact_analysis/plots_emulator")
