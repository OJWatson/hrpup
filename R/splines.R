#' Check if a vector is monotonically increasing
#'
#' Internal utility function that returns TRUE if the vector `x` is monotonically increasing.
#'
#' @param x Numeric vector to check.
#'
#' @return Logical TRUE if monotonically increasing, FALSE otherwise.
#' @noRd
is_increasing <- function(x) {
  all(diff(x) >= 0)
}

#' Find transition from increasing to decreasing
#'
#' Identifies the index where the numeric vector first switches from monotonic increasing to decreasing.
#'
#' @param x Numeric vector used to identify the transition.
#'
#' @return Integer index of transition point.
#' @noRd
find_transition <- function(x) {
  diffs <- diff(x)
  increasing_started <- FALSE

  for (i in seq_along(diffs)) {
    if (!increasing_started && diffs[i] > 0) {
      increasing_started <- TRUE
    } else if (increasing_started && diffs[i] <= 0) {
      return(i + 1)
    }
  }

  length(x)
}

#' Apply piecewise monotonic spline smoothing with matched frequency scaling
#'
#' Applies monotonic increasing spline smoothing to columns of `preds` until `freq_vector`
#' stops increasing, then applies monotonic decreasing smoothing scaled proportionally
#' according to the relative decrease of `freq_vector`. Specifically, each column
#' decreases back down to the value it previously had when `freq_vector` was equal to
#' its final decreased value.
#'
#' @param preds Data frame containing numeric columns to smooth.
#' @param freq_vector Numeric vector indicating the monotonic pattern.
#'
#' @return Data frame with smoothly joined monotonic constraints.
#' @noRd
apply_monotonic_spline <- function(preds, freq_vector, s) {

  # Find transition point
  transition_point <- find_transition(freq_vector)
  transition_point <- min(nrow(preds), transition_point+24)
  if(transition_point > (nrow(preds)-24)) {
    transition_point <- nrow(preds)
  }

  # First part: monotonically increasing spline
  preds_inc <- preds[1:transition_point, ] %>%
    mutate(across(everything(), ~ {
      if(s < 0.1) {
        fit_pos_guided_spline(.x, c(freq_vector[1:(transition_point-24)], rep(tail(transition_point,24)[1], 24)))
      } else {
      df <- data.frame(x = seq_along(.x), y = .x)
      my <- max(df$y)
      df$y <- df$y/my
      fit <- scam(y ~ s(x, bs = "mpi"), data = df)
      predict(fit, newdata = df, type = "response")*my
      }
    }))

  # Second part: monotonically decreasing proportional to freq_vector decrease
  if (transition_point < nrow(preds)) {
    preds_dec <- preds[(transition_point + 1):nrow(preds), ]

    # Get the corresponding final frequency after decrease
    final_freq <- tail(freq_vector, 1)

    # Find the index in the increasing phase closest to the final_freq
    if(s < 0.1) {
      freq_inc <- freq_vector[1:(transition_point-24)]
    } else {
    freq_inc <- freq_vector[1:transition_point]
    }
    idx_closest <- max(which(abs(freq_inc - final_freq) == min(abs(freq_inc - final_freq))))

    # Loop through each column explicitly to apply scaled decrease
    for (colname in colnames(preds_dec)) {
      column <- preds_dec[[colname]]

      # Starting value at transition
      start_val <- preds_inc[[colname]][transition_point]

      # Target value: matching the value at closest freq_vector index
      target_val <- preds_inc[[colname]][min(idx_closest+24,nrow(preds_inc))]

      # Fit monotonic decreasing spline to smoothly connect start_val and target_val
      df <- data.frame(x = seq_along(column), y = column)
      my <- max(df$y)
      df$y <- df$y/my

      # Monotonically decreasing spline (unscaled first)
      fit <- fit_scam_with_retries(df)
      pred_values <- predict(fit, newdata = df)*my

      # Rescale predicted values to precisely match start_val and target_val
      pred_values_rescaled <- target_val +
        (pred_values - pred_values[length(pred_values)]) /
        (pred_values[1] - pred_values[length(pred_values)]) * (start_val - target_val)

      preds_dec[[colname]] <- pred_values_rescaled
    }

    # Combine increasing and decreasing segments
    preds_smoothed <- bind_rows(preds_inc, preds_dec)
  } else {
    preds_smoothed <- preds_inc
  }

  # correct for any freq 0 locations
  if(any(freq_vector == 0)) {
    preds_smoothed <- preds_smoothed %>%
      mutate(across(everything(), ~ {
        .x[which(freq_vector == 0)] <- tail(.x[which(freq_vector == 0)],1)
        .x
      }))
  }

  # and smooth the spline joint
  index <- seq_along(preds_smoothed[[1]])
  preds_smoothed <- preds_smoothed %>%
    mutate(across(everything(), ~ {
      gam_fit <- mgcv::gam(.x ~ s(index, k = 30))
      predict(gam_fit, newdata = data.frame(index = index))
    }))

  return(preds_smoothed)
}

#' @noRd
fit_scam_with_retries <- function(df, max_retries = 3, initial_gamma = 1.4, gamma_increment = 0.3) {
  attempt <- 1
  gamma <- initial_gamma

  repeat {
    fit <- tryCatch({
      scam(y ~ s(x, bs = "micv", k = min(nrow(df)-1, 10)), data = df, gamma = gamma)
    }, error = function(e) {
      return(NULL)
    })

    if (!is.null(fit)) {
      return(fit)
    }

    if (attempt >= max_retries) {
      stop(sprintf("Failed after %d attempts. Last gamma = %.1f", attempt, gamma))
    }

    # increment gamma and retry
    gamma <- gamma + gamma_increment
    attempt <- attempt + 1
  }
}

#' Fits a guided constrained positive spline
#'
#' @param vec Numeric vector
#' @param ref_shape Guide vector
#' @param flexibility closeness to guide shape
#'
#' @noRd
fit_pos_guided_spline <- function(vec, ref_shape){
  df <- data.frame(x = seq_along(vec))

  # Calculate range scaling factor from ref_shape
  ref_range <- diff(range(ref_shape))
  scaling_factor <- min(max(ref_range, 0), 1)  # ensure it stays between 0 and 1
  scaling_factor <- ((1-exp(100*-scaling_factor)))

  # Set initial start and end points
  start_val <- mean(vec)
  end_val <- tail(vec, 1)

  reversed <- FALSE
  if (start_val > end_val) {
    reversed <- TRUE
    tmp <- start_val
    start_val <- end_val
    end_val <- tmp
  }

  # Compute scaled start point, closer to end_val based on scaling_factor
  adjusted_start <- end_val - scaling_factor * (end_val - start_val)

  # Scale ref_shape to [0,1] for guiding spline shape
  ref_scaled <- (ref_shape - min(ref_shape)) / ref_range

  # Create target curve interpolating adjusted_start to end_val guided by ref_shape
  target_curve <- adjusted_start + ref_scaled * (end_val - adjusted_start)

  # Fit monotone increasing spline to the target curve
  df$y <- target_curve

  fit <- scam(y ~ s(x, bs = "mpi", k = min(nrow(df)-1, 10)), data = df)

  fitted_curve <- predict(fit, newdata = df)

  fitted_curve <- as.numeric(fitted_curve)
  fitted_curve <- fitted_curve + (end_val-tail(fitted_curve,1))

  # correct for any freq 0 locations
  if(any(ref_shape == 0)) {

    fitted_curve[which(ref_shape == 0)] <- tail(preds_smoothed[which(ref_shape == 0)],1)

    # and smooth the spline joint
    index <- seq_along(fitted_curve)
    gam_fit <- mgcv::gam(fitted_curve ~ s(index, k = 30))
    fitted_curve <- predict(gam_fit, newdata = data.frame(index = index))

  }

  return(fitted_curve)
}

#' Fits a guided constrained negative spline
#'
#' @param vec Numeric vector
#' @param ref_shape Guide vector
#' @param flexibility closeness to guide shape
#'
#' @noRd
fit_neg_guided_spline <- function(vec, ref_shape){
  df <- data.frame(x = seq_along(vec))

  # Calculate range scaling factor from ref_shape
  ref_range <- diff(range(ref_shape))
  scaling_factor <- min(max(abs(ref_range), 0), 1)  # between 0 and 1

  # Set initial start and end points
  start_val <- mean(vec)
  end_val <- tail(vec, 1)

  reversed <- FALSE
  if (start_val < end_val) {
    reversed <- TRUE
    tmp <- start_val
    start_val <- end_val
    end_val <- tmp
  }

  # Adjusted start point moves from end towards start based on scaling_factor
  adjusted_start <- end_val + scaling_factor * (start_val - end_val)

  # Scale ref_shape to [0,1] (decreasing reference shape)
  ref_scaled <- (ref_shape - min(ref_shape)) / abs(ref_range)

  # Invert scaling for decreasing shape
  ref_scaled <- 1 - ref_scaled

  # Generate target curve
  target_curve <- adjusted_start - ref_scaled * (adjusted_start - end_val)

  df$y <- target_curve

  # Fit monotone decreasing spline (mpd)
  fit <- scam(y ~ s(x, bs = "mpd", k = min(nrow(df)-1, 10)), data = df)

  fitted_curve <- predict(fit, newdata = df)

  fitted_curve <- as.numeric(fitted_curve)
  fitted_curve <- fitted_curve + (end_val-tail(fitted_curve,1))

  return(fitted_curve)
}
