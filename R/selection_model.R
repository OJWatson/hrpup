#' R6 Class for hrp2 selection model
#'
#' Model for predicting the selective coefficient for hrp2 deletions
R6_hrp2_mod <- R6::R6Class(
  classname = "hrp2_mod",
  cloneable = FALSE,

  # PUBLIC METHODS
  public = list(

    # INITIALISATION
    #' @description
    #' Create a new hrp2 selection model using emulators.
    #' @param models Model for prediction seletion coefficient
    #' @param err_models Model for predicting sd error for predictions
    #' @param data Data used originally for model training
    #' @return A new `hrp2_mod` object.
    initialize = function(models = list(),
                          err_models  = list(),
                          model_weights = list(),
                          err_model_weights = list(),
                          data = data.frame()) {

      private$models <- models
      private$err_models <- err_models
      private$model_weights <- model_weights
      private$err_model_weights <- err_model_weights
      private$data <- data

    },

    #' @description
    #' Predict median and 95% selection coefficients and threshold times
    #' @param Micro.2.10 Microscopy prevalence in 2-10 year olds
    #' @param ft Treatment seeking rate
    #' @param microscopy.use Proportion of total testing by microscopy
    #' @param rdt.nonadherence Probability of treating negative indiviudal
    #' @param fitness Relative fitness of hrp2 deleted parasite
    #' @param rdt.det Chance of hrp2 deleted parasite yielding positive test
    #' @returns Data.frame of selection coefficients and times from 1% to 5%
    predict = function(Micro.2.10, ft, microscopy.use, rdt.nonadherence, fitness, rdt.det) {

      # get args and turn into matrix
      dat <- cbind(Micro.2.10, ft, microscopy.use, rdt.nonadherence, fitness, rdt.det)
      ret <- as.data.frame(dat)

      # predict s
      ret$s <- self$predict_s(dat)
      err <- self$predict_err(dat)
      ret$smin <- ret$s - 1.96*(err/5)
      ret$smax <- ret$s + 1.96*(err/5)

      # use to calculate times
      ret$t <- (log(0.05 / (1 - 0.05)) - log(0.01 / (1 - 0.01))) / ret$s
      ret$tmin <- (log(0.05 / (1 - 0.05)) - log(0.01 / (1 - 0.01))) / ret$smin
      ret$tmax <- (log(0.05 / (1 - 0.05)) - log(0.01 / (1 - 0.01))) / ret$smax


      return(ret)

    },

    # Predict Ensemble s
    predict_s = function(dat) {
      private$predict_internal(dat, private$models, private$model_weights)
    },

    # Predict Ensemble Error
    predict_err = function(dat) {
      private$predict_internal(dat, private$err_models, private$err_model_weights)
    },

    # Predict t for f1 and f2 given s
    predict_t = function(dat, f1, f2) {

      # catch for directionality
      if(f1 >= f2) {
        stop("f1 must be less than f2")
      }

      # create results name in data frame
      # name <- paste0("t_", f1, "_", f2)

      # create s if not available
      if(!("s" %in% names(dat))) {
        # predict s
        dat$s <- self$predict_s(dat)
        err <- self$predict_err(dat)
        dat$smin <- dat$s - 1.96*(err/5)
        dat$smax <- dat$s + 1.96*(err/5)
      }

      # create new time
      # dat[[name]] <- (log(f2 / (1 - f2)) - log(f1 / (1 - f1))) / dat$s
      # return(dat)
      ret_t <- (log(f2 / (1 - f2)) - log(f1 / (1 - f1))) / dat$s
      return(ret_t)
    },

    # Predict f2 for f1 and t given s
    predict_f2 = function(dat, f1, t) {

      # create results name in data frame
      # name <- paste0("f2_", t, "_", f1)

      # create s if not available
      if(!("s" %in% names(dat))) {
        # predict s
        dat$s <- self$predict_s(dat)
        # err <- self$predict_err(dat)
        # dat$smin <- dat$s - 1.96*err
        # dat$smax <- dat$s + 1.96*err
      }

      # create new freq
      f2 = (exp(t * dat$s) * (f1 / (1 - f1))) / (1 + exp(t * dat$s) * (f1 / (1 - f1)))
      # dat[[name]] <- f2
      return(f2)
    },


    # Add models and weights
    add_model = function(model, model_name) {
      private$models[[model_name]] <- model
    },

    add_model_weight = function(weight, model_name) {
      private$model_weights[[model_name]] <- weight
    },

    add_err_model = function(err_model, model_name) {
      private$err_models[[model_name]] <- err_model
    },

    add_err_model_weight = function(weight, model_name) {
      private$err_model_weights[[model_name]] <- weight
    },

    # GETTERS
    get_models = function() private$models,
    get_error_models = function() private$err_models,
    get_model_weights = function() private$model_weights,
    get_error_model_weights = function() private$err_model_weights,
    get_data = function() private$data,

    # SETTERS
    set_models = function(models) { private$models <- models },
    set_error_models = function(err_models) { private$err_models <- err_models },
    set_model_weights = function(model_weights) { private$model_weights <- model_weights },
    set_error_model_weights = function(err_model_weights) { private$err_model_weights <- err_model_weights },
    set_data = function(data) { private$data <- data }

  ),

  private = list(
    models = NULL,
    err_models = NULL,
    model_weights = NULL,
    err_model_weights = NULL,
    data = NULL,

    # Predict Generic
    predict_internal = function(dat, models, weights) {
      model_names <- names(models)
      model_weights <- 1 / as.numeric(weights[model_names])
      normalised_weights <- model_weights / sum(model_weights)
      dat <- dat[,c("Micro.2.10","ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det")]
      dat <- as.matrix(dat)
      predictions <- vapply(lapply(models, predict, dat), as.numeric, FUN.VALUE = numeric(nrow(dat)))
      ensemb <- apply(predictions, MARGIN = 1, weighted.mean, normalised_weights)
      return(ensemb)
    }
  )
)
