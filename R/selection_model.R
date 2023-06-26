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
    #' @param models Model for predicting selection coefficient
    #' @param err_models Model for predicting error of selection predictions
    #' @param model_weights Model for predicting selection coefficient weights
    #' @param err_model_weights Model for predicting selection error weights
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
    #' @return Data.frame of selection coefficients and times from 1% to 5%
    predict = function(Micro.2.10, ft, microscopy.use, rdt.nonadherence, fitness, rdt.det) {

      # get args and turn into matrix
      dat <- cbind(Micro.2.10, ft, microscopy.use, rdt.nonadherence, fitness, rdt.det)
      ret <- as.data.frame(dat)

      # predict s
      ret$s <- self$predict_s(dat)
      err <- self$predict_err(dat)
      # We estimated error using only 5 stochastic realisations. If we had used
      # 100 realisations, bsaed on a simple model of normally distributed variance
      # it would be 5 times smaller. Consequently, divide by 5.
      ret$smin <- ret$s - 1.96*(err/5)
      ret$smax <- ret$s + 1.96*(err/5)

      # use to calculate times
      ret$t <- (log(0.05 / (1 - 0.05)) - log(0.01 / (1 - 0.01))) / ret$s
      ret$tmin <- (log(0.05 / (1 - 0.05)) - log(0.01 / (1 - 0.01))) / ret$smin
      ret$tmax <- (log(0.05 / (1 - 0.05)) - log(0.01 / (1 - 0.01))) / ret$smax


      return(ret)

    },

    # Predict Ensemble s
    #' Predict selection coefficient for data frame of covariates
    #' @param dat Data frame of covariates
    predict_s = function(dat) {
      private$predict_internal(dat, private$models, private$model_weights)
    },

    # Predict Ensemble Error
    #' Predict error in selection coefficient estimates for covariate data frame
    #' @param dat Data frame of covariates
    predict_err = function(dat) {
      private$predict_internal(dat, private$err_models, private$err_model_weights)
    },

    # Predict t
    #' Predict t between f1 and f2 and covariate data frame
    #' @param dat Data frame of covariates
    #' @param f1 Frequency at time point 1
    #' @param f2 Frequency at time point 2
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

    # Predict f2
    #' Predict f2 given f1 and t and covariate data frame
    #' @param dat Data frame of covariates
    #' @param f1 Frequency at time point 1
    #' @param t Duration of selection
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


    #' Add models
    #' @param model Selection prediction model
    #' @param model_name Name of model
    add_model = function(model, model_name) {
      private$models[[model_name]] <- model
    },

    #' Add model weights
    #' @param weight Selection prediction model weight (RMSE)
    #' @param model_name Name of model
    add_model_weight = function(weight, model_name) {
      private$model_weights[[model_name]] <- weight
    },

    #' Add error model
    #' @param err_model Selection error prediction model
    #' @param model_name Name of model
    add_err_model = function(err_model, model_name) {
      private$err_models[[model_name]] <- err_model
    },

    #' Add error model weights
    #' @param weight Selection error prediction model weight (RMSE)
    #' @param model_name Name of model
    add_err_model_weight = function(weight, model_name) {
      private$err_model_weights[[model_name]] <- weight
    },

    # GETTERS

    #' Get all selection prediction models
    #' @return List of models
    get_models = function() private$models,

    #' Get all selection error prediction models
    #' @return List of error models
    get_error_models = function() private$err_models,

    #' Get all selection prediction model weights
    #' @return List of mode weights
    get_model_weights = function() private$model_weights,

    #' Get all selection error prediction model weights
    #' @return List of error model weights
    get_error_model_weights = function() private$err_model_weights,

    #' Get data the models were trained on
    #' @return Training data
    get_data = function() private$data,

    # SETTERS

    #' Set all selection prediction models
    #' @param models selection model list
    set_models = function(models) { private$models <- models },

    #' Set all selection error prediction models
    #' @param err_models selection error model list
    set_error_models = function(err_models) { private$err_models <- err_models },

    #' Set all selection prediction model weights
    #' @param model_weights selection model weights (RMSE) list
    set_model_weights = function(model_weights) { private$model_weights <- model_weights },

    #' Set all selection prediction model weights
    #' @param err_model_weights selection error model weights (RMSE) list
    set_error_model_weights = function(err_model_weights) { private$err_model_weights <- err_model_weights },

    #' Set data the models were trained on
    #' @param data Training data for models
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

      # get our models and their weights
      model_names <- names(models)
      model_weights <- 1 / as.numeric(weights[model_names])
      normalised_weights <- model_weights / sum(model_weights)

      # set up ur data removing NA rows
      dat <- as.data.frame(dat[,c("Micro.2.10","ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det")])
      dat_na <- na.omit(dat)

      # make ensemble prediction
      predictions <- vapply(lapply(models, predict, dat_na), as.numeric, FUN.VALUE = numeric(nrow(dat_na)))

      # Catch for when you are only requesting on one row of data
      if(!is.matrix(predictions)) {
        ensemb <- weighted.mean(predictions, normalised_weights)
      } else {
        ensemb <- apply(predictions, MARGIN = 1, weighted.mean, normalised_weights)
      }

      # return values with NAs in for missing data
      ret <- rep(NA, nrow(dat))
      ret[as.integer(which(apply(dat, 1, function(x){all(!is.na(x))})))] <- ensemb
      return(ret)
    }
  )
)
