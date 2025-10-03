## Spread Modelling
library(tidyverse)
library(brnn)
library(caret)
library(earth)
library(scam)
library(furrr)
library(future)
library(sf)
library(wpp2017)
library(torch)
library(luz)
library(quadprog)
data("UNlocations")
sf::sf_use_s2(FALSE)

#' Predict Malaria Outcomes on Real Data Using a Trained Model
#'
#' Uses a trained GRU model to predict 8 malaria indicators for each of three versions
#' of a real-world frequency input: lower confidence interval (`freq_lci`), mean (`freq`),
#' and upper confidence interval (`freq_uci`). Output is a single data frame with
#' 24 prediction columns (8 indicators × 3 frequency variants).
#'
#' @param trained_model A `luz` model object with `extra$norm_ranges` and `extra$warmup`.
#' @param real_data A data frame with one observation per row and all required input columns,
#'   including columns given by `names(freq_l)`
#' @param freq_l A list of the freq columns to use and their labels.
#'   Default: `list("freq_lci" = "lci", "freq" = "med", "freq_uci" = "uci")`
#'
#' @return A data frame of 24 predicted values, normalised back to original scale and named
#'         with `_lci`, `_mean`, and `_uci` suffixes per frequency type.
#'
#' @export
predict_from_real_data <- function(trained_model, real_data,
                                   freq_l = list("freq_lci" = "lci",
                                                 "freq_med" = "med",
                                                 "freq_uci" = "uci"),
                                   old = FALSE) {

  # Pull normalisation parameters
  norm_ranges <- trained_model$extra$norm_ranges
  warmup <- trained_model$extra$warmup
  input_cols <- trained_model$extra$input_cols
  output_cols <- trained_model$extra$output_cols

  is_monotonic_inc <- function(x) {
    all(diff(x) > 0)
  }
  is_monotonic_inc_then_flat <- function(x) {
    d <- diff(x)
    if(all(d) >= 0) {
      if(all(vapply(which(d==0), function(x){all(x > which(d > 0))}, logical(1)))) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    } else {
      return(FALSE)
    }
  }
  is_monotonic_dec <- function(x) {
    all(diff(x) < 0)
  }

  # Function to make predictions for one freq variant
  predict_one_freq <- function(freq_vector, suffix) {

    # which s range are we using
    s_var <- c("smin", "s", "smax")[match(suffix, c("lci", "med", "uci"))]
    s <- real_data[[s_var]][1]

    input_data <- real_data %>%
      dplyr::mutate(freq = freq_vector) %>%
      dplyr::select(all_of(trained_model$extra$input_cols))

    # Normalise
    norm_data <- trained_model$extra$normalisation_fn(
      input_data,  norm_ranges[colnames(input_data)]
    )

    # Build warm-up sequence
    input <- trained_model$extra$warm_up_fn(norm_data, input_cols, warmup)
    input_tensor <- torch::torch_tensor(input, dtype = torch::torch_float())$unsqueeze(1)

    # Predict
    device <- if (torch::cuda_is_available()) torch::torch_device("cuda") else torch::torch_device("cpu")
    trained_model$model$eval()
    preds <- trained_model$model(input_tensor$to(device = device))$squeeze()$cpu()$detach()
    preds <- as.data.frame(as.array(preds))

    # Drop warmup steps
    preds <- preds[-seq_len(warmup), ]

    # Name and denormalise
    colnames(preds) <- output_cols
    preds <- trained_model$extra$renormalisation_fn(
      preds, norm_ranges[colnames(preds)]
    )

    if(!old){
      # no change in freq
      if (all(freq_vector == 0)) {
        new_preds <- preds %>% mutate(across(everything(), mean))
        # no change in freq
      } else if (all(diff(freq_vector)==0)) {
        new_preds <- preds %>% mutate(across(everything(), mean))
        # monotonic decrease
      } else if (s <= 0 || is_monotonic_dec(freq_vector)) {
        new_preds <- preds %>% mutate(across(
          .cols = everything(),
          ~ fit_neg_guided_spline(.x, freq_vector)
        ))
        # monotonic increase small
      } else if (s < 0.1 && s >= 0 && is_monotonic_inc(freq_vector)) {
        new_preds <- preds %>% mutate(across(
          .cols = everything(),
          ~ fit_pos_guided_spline(.x, freq_vector)
        ))
        # monotonic increase large so just fit the first 2 years
      } else if (s >= 0.1 && (is_monotonic_inc(freq_vector) || is_monotonic_inc_then_flat(freq_vector))) {
        new_preds <- preds %>% mutate(across(everything(), ~ {
          df <- data.frame(x = seq_along(.x), y = .x)
          my <- max(df$y)
          df$y <- df$y/my
          fit <- scam(y ~ s(x, bs = "mpi"), data = df)
          newy <- predict(fit, newdata = df)*my
          oldy <- .x
          oldy[1:24] <- (newy[1:24]/max(newy[1:24]))*oldy[24]
          oldy
        }))
        # case where it is always decreasing  but not monotonic due to importation
        # } else if (s < 0 && (max(freq_vector)<0.2)) {
        #   new_preds <- preds %>% mutate(across(
        #     .cols = everything(),
        #     ~ fit_neg_guided_spline(.x, freq_vector, flexibility = 0.9)
        #   ))
        # } else if (s < 0.1 && s >= 0 && max(freq_vector) < 0.2) {
        #   new_preds <- preds %>% mutate(across(
        #     .cols = everything(),
        #     ~ fit_guided_spline(.x, freq_vector, flexibility = 0.9)
        #   ))
      } else if (s >= 0) {
        new_preds <- preds %>% apply_monotonic_spline(freq_vector, s)
      } else {
        stop("uncaught")
      }
      preds <- new_preds
    }
    colnames(preds) <- paste0(colnames(preds), "_", suffix)

    preds
  }

  # Run for each freq variant
  preds <- lapply(seq_along(freq_l), function(x){
    predict_one_freq(real_data[[names(freq_l)[x]]], freq_l[[x]])
    })
  final_preds <- do.call(cbind, preds)

  # Order correctly
  sl <- as.character(vapply(trained_model$extra$output_cols,
                            function(x){paste0(x, c("_lci", "_med", "_uci"))},
                            FUN.VALUE = character(3)))
  final_preds <- final_preds %>% select(all_of(sl))

  return(final_preds)
}

# --------------------------------------------------------------------------#
# 1. Get our needed objects to start simulating African spread of hrp2 -----
# --------------------------------------------------------------------------#

# Get the mapping object
map_obj <- readRDS("analysis/data_derived/R6_map.rds")
admin1 <- readRDS(here::here("analysis/data_derived/admin1_sf.rds"))

# Get the selection model object
mod_obj <- readRDS("analysis/data_derived/ensemble_selection_model.rds")

# Get the emulator model object
trained_model <- luz::luz_load("analysis/impact_analysis/data-derived/emulator_final.rds")
trained_model$model$to(device = torch_device("cuda"))

# subset to Africa as that is where we will do the prospective mapping
admin1$cont <- countrycode::countrycode(admin1$iso, "iso3c", "iso3n")
admin1$region <- UNlocations$area_name[match(admin1$cont, UNlocations$country_code)]
afr_map <- admin1 %>% filter(region == "Africa")

# get our non_hrp2 data
non_hrp2_use <- read.csv("analysis/data_raw/hrp2_RDT_usage.csv")
mean_afr_use <- non_hrp2_use %>% filter(iso3c %in% afr_map$iso) %>% pull(hrp2) %>% mean

# Get the scenario maps
scenario_maps <- readRDS("analysis/data_derived/scenario_maps_full.rds")

# Just noticed mal prev for DJI is wrong way round...
covars <- readRDS("analysis/data_derived/global_covariate_ranges.rds")
for(i in seq_along(scenario_maps$map_data)) {

  scenario_maps$map_data[[i]]$Micro.2.10[scenario_maps$map_data[[i]]$iso3c == "DJI"] <-
    case_when(scenario_maps$scenarios$Micro.2.10[i] == "worst" ~ covars$Micro.2.10_mean[covars$iso3c == "DJI"]*0.99,
              scenario_maps$scenarios$Micro.2.10[i] == "best" ~ covars$Micro.2.10_mean[covars$iso3c == "DJI"]*1.01,
              scenario_maps$scenarios$Micro.2.10[i] == "central" ~ covars$Micro.2.10_mean[covars$iso3c == "DJI"])
}

# update the map data with our non_hrp2_data
for(i in seq_along(scenario_maps$map_data)) {
  md <- scenario_maps$map_data[[i]]
  md <- left_join(md, non_hrp2_use %>% select(hrp2, iso3c))
  md <- md %>% mutate(hrp2 = replace_na(hrp2, mean_afr_use)) # GAB and ENQ missing volume data. Assume mean
  curr_micro <- md$microscopy.use

  # how many are not getting microscopy
  curr_non_micro <- 1-md$microscopy.use

  # what proportion of these are getting a non hrp2 and add that to the current micro
  # to now get total not getting an HRP2 RDT for diagnosis
  new_micro <- curr_micro + ((curr_non_micro) * (1-md$hrp2))
  # assume there is still always 1% of the other diagnostic option
  new_micro <- pmax(pmin(new_micro, 0.99), 0.01)

  md$microscopy.use <- new_micro
  md <- md %>% select(-hrp2)
  pred <- mod_obj$predict(Micro.2.10 = md$Micro.2.10,
                          ft = md$ft,
                          microscopy.use = md$microscopy.use,
                          rdt.nonadherence = md$rdt.nonadherence,
                          fitness = md$fitness,
                          rdt.det = md$rdt.det)
  pred <- pred %>%
    mutate(id_1 = md$id_1, .before = 1) %>%
    mutate(iso3c = md$iso3c, .before = 1)
  scenario_maps$map_data[[i]] <- pred
}

# And subset the afr_map to just those regions with selection data
map_data <- scenario_maps$map_data[[365]] %>% filter(id_1 %in% afr_map$id_1)
map_data <- filter(map_data, !is.na(s))
afr_map <- filter(afr_map, id_1 %in% map_data$id_1)

# get an adjacency matrix to figure out spread
adj_mat <- spdep::poly2nb(afr_map)
names(adj_mat) <- afr_map$id_1

# And lastly get our seed data
# Set up our initial conditions
seeds <- readxl::read_excel("analysis/data_raw/WHO MPAG map 15.4.23.xlsx",
                            col_names = TRUE, col_types = c("text", "numeric"), range = "A1:B53") %>%
  setNames(c("name","prev")) %>%
  mutate(iso = countrycode::countrycode(name, "country.name.en", "iso3c")) %>%
  mutate(continent = countrycode::countrycode(iso, "iso3c", "continent")) %>%
  filter(continent == "Africa") %>%
  filter(prev > 0) %>%
  mutate(prev = prev/100)
seeds_adm1 <- afr_map %>% sf::st_drop_geometry() %>% left_join(seeds) %>% filter(prev > 0)

# --------------------------------------------------------------------------#
# 2. Set Up RDT switch rules -----------------------------------------------------
# --------------------------------------------------------------------------#

# create per country list
switch_iso3s <- unique(afr_map$iso)
switch_info <- vector("list", length = length(switch_iso3s))
switch_iso3cs <- unique(afr_map$iso)
names(switch_info) <- switch_iso3cs

# t0 "Jan 2024"

# DJI, ETH, ERI
# assume that by now DJI, ERI, ETH ae solely using PAN, i.e microscopy.use
for (key in c("DJI", "ERI", "ETH")) {
  switch_info[[key]] <- data.frame(t = 0, microscopy.use = 0.99)
}

# others
# set to some abritraty date in the future so the only switch rules being used are those
# based on thresholds
switch_time <- as.integer((as.Date("2927-01-01") - as.Date("2024-01-01")))/365
for (key in setdiff(switch_iso3s, c("DJI", "ERI", "ETH"))) {
  switch_info[[key]] <- data.frame(t = switch_time + 0:5, microscopy.use = seq(0, 0.99, length.out = 7)[-1])
}

# --------------------------------------------------------------------------#
# 3. Run our simulation with RDT switch strategies --------------------------
# --------------------------------------------------------------------------#

# Initialise our model
spread_model <- R6_hrp2_spread$new(afr_map, mod_obj, adj_mat = adj_mat)

# Set up our initial conditions
spread_model$set_seeds(setNames(seeds_adm1$prev, seeds_adm1$id_1))

# function to run uncertainty ranges
run_lci_med_uci_sim <- function(spread_model, md, delay = 0, RDT_switch = FALSE, t_break = 0.1, t_end = 20.074, f_trig = 0.05, emulator = NULL, interp = TRUE) {

  # bring together with threshold_switch info
  switch_list <- list("iso_switch" = switch_info,
                      "threshold_switch" = data.frame(t_past = delay, microscopy.use = c(0.99), f_trig = f_trig))
  # set the data
  spread_model$set_map_data(md)

  # first simulate central
  if(RDT_switch) {
    med <- spread_model$simulate_spread_and_rdts(
      export_freq = 0.25, t_break = t_break, t_end = t_end, switch_list = switch_list, unc = "med"
    )
  } else {
    med <- spread_model$simulate_spread(export_freq = 0.25, t_break = t_break, unc = "med")
  }

  # then worst
  if(RDT_switch) {
    uci <- spread_model$simulate_spread_and_rdts(
      export_freq = 0.25, t_break = t_break, t_end = t_end, switch_list = switch_list, unc = "uci"
    )
  } else {
    uci <- spread_model$simulate_spread(export_freq = 0.25, t_break = t_break, unc = "uci")
  }

  # then best
  if(RDT_switch) {
    lci <- spread_model$simulate_spread_and_rdts(
      export_freq = 0.25, t_break = t_break, t_end = t_end, switch_list = switch_list, unc = "lci"
    )
  } else {
    lci <- spread_model$simulate_spread(export_freq = 0.25, t_break = t_break, unc = "lci")
  }

  out <- lci %>% select(id_1, t, freq_lci = freq) %>%
    left_join(med %>% select(id_1, t, freq_med = freq)) %>%
    left_join(uci %>% select(id_1, t, freq_uci = freq)) %>%
    select(id_1, t, matches("freq")) %>%
    mutate(delay = ifelse(RDT_switch, delay, -1), .after = "t")

  # if emulator was passed then add extra outputs
  if (!is.null(emulator)) {

    full <- left_join(out, md %>% select(-t), by = "id_1") %>%
      split(.$id_1) %>% lapply(head, emulator$extra$data$data_length)

    preds <- map(full, function(x){
      preds <- predict_from_real_data(emulator, x, old = !interp)
      preds <- cbind(x %>% select(names(out)), preds)
    })

    out <- do.call(rbind, preds)

  }

  return(out)

}

# our delays to explore
delays <- c(0:10)
f_trigs <- c(0.05)
param_grid <- expand.grid("delay" = delays, "f_trig" = f_trigs)

# emulator
# emulator <- NULL # dont use emulator
# t_break <- 0.1
emulator <- trained_model
t_break <- round(emulator$extra$data$t_break, 4)
t_end <- (emulator$extra$data$t_end + t_break)

# Central
central_row <- which(apply(scenario_maps$scenarios, 1, function(x){all(x == "central")}))
out_central <- list()
for(d in seq_along(param_grid$delay)) {
  type <- paste0(scales::percent(param_grid$f_trig[d]), " Threshold Strategy")
out_central[[d]] <- run_lci_med_uci_sim(spread_model, scenario_maps$map_data[[central_row]], RDT_switch = TRUE,
                                        delay = param_grid$delay[d], t_break = t_break, emulator = emulator,
                                        f_trig = param_grid$f_trig[d]) %>%
  mutate(scenario = "Central", type = type, .after = "delay")
}
# Best
best_row <- which(apply(scenario_maps$scenarios, 1, function(row) all(row == "best")))
out_best <- list()
for(d in seq_along(param_grid$delay)) {
  type <- paste0(scales::percent(param_grid$f_trig[d]), " Threshold Strategy")
out_best[[d]] <- run_lci_med_uci_sim(spread_model, scenario_maps$map_data[[best_row]], RDT_switch = TRUE,
                                     delay = param_grid$delay[d], t_break = t_break, emulator = emulator,
                                     f_trig = param_grid$f_trig[d]) %>%
  mutate(scenario = "Best", type = type, .after = "delay")
}

# Worst
worst_row <- which(apply(scenario_maps$scenarios, 1, function(row) all(row == "worst")))
out_worst <- list()
for(d in seq_along(param_grid$delay)) {
  type <- paste0(scales::percent(param_grid$f_trig[d]), " Threshold Strategy")
out_worst[[d]] <- run_lci_med_uci_sim(spread_model, scenario_maps$map_data[[worst_row]], RDT_switch = TRUE,
                                      delay = param_grid$delay[d], t_break = t_break, emulator = emulator,
                                      f_trig = param_grid$f_trig[d]) %>%
  mutate(scenario = "Worst", type = type, .after = "delay")
}

# --------------------------------------------------------------------------#
# 4. Run our simulation without RDT switch strategies -----------------------
# --------------------------------------------------------------------------#

# Here to make it comparable, let's assume that DJI, ERI and ETH all are at 100% non-hrp2 too
comp_md_func <- function(md) {
  md %>% filter(id_1 %in% afr_map$id_1) %>%
    mutate(microscopy.use_old = microscopy.use) %>%
    mutate(microscopy.use = replace(microscopy.use, iso3c %in% c("DJI", "ERI", "ETH"), 0.99)) %>%
    mod_obj$predict_dat(.) %>%
    mutate(s = if_else(iso3c %in% c("DJI", "ERI", "ETH"), pmin(s, 0), s)) %>%
    mutate(smin = if_else(iso3c %in% c("DJI", "ERI", "ETH"), pmin(smin, 0), smin)) %>%
    mutate(smax = if_else(iso3c %in% c("DJI", "ERI", "ETH"), pmin(smax, 0), smax)) %>%
    mutate(microscopy.use = microscopy.use_old) %>%
    select(-microscopy.use_old)
}

# Central
out_central_cf <- run_lci_med_uci_sim(spread_model, comp_md_func(scenario_maps$map_data[[central_row]]),
                                      t_break = t_break, emulator = emulator) %>%
  mutate(scenario = "Central", type = "No RDT Switching", .after = "delay")

# Best
out_best_cf <- run_lci_med_uci_sim(spread_model, comp_md_func(scenario_maps$map_data[[best_row]]),
                                   t_break = t_break, emulator = emulator) %>%
  mutate(scenario = "Best", type = "No RDT Switching", .after = "delay")

# Worst
out_worst_cf <- run_lci_med_uci_sim(spread_model, comp_md_func(scenario_maps$map_data[[worst_row]]),
                                    t_break = t_break, emulator = emulator) %>%
  mutate(scenario = "Worst", type = "No RDT Switching", .after = "delay")


# --------------------------------------------------------------------------#
# 5. Run our simulation fixed time RDT switch strategies -----------------------
# --------------------------------------------------------------------------#

# function to run uncertainty ranges
run_lci_med_uci_sim_fixed <- function(spread_model, md, delay = 0, RDT_switch = FALSE, t_break = 0.1, t_end = 20.074, f_trig = Inf, emulator = NULL, interp = TRUE) {


  # first set up specific switch info
  switch_iso3cs <- unique(spread_model$.__enclos_env__$private$map$iso)
  switch_info <- vector("list", length = length(switch_iso3cs))
  names(switch_info) <- switch_iso3cs

  # DJI, ETH, ERI
  # assume that by now DJI, ERI, ETH ae solely using PAN, i.e microscopy.use
  for (key in c("DJI", "ERI", "ETH")) {
    switch_info[[key]] <- data.frame(t = 0, microscopy.use = 0.99)
  }

  # others
  # set to delay years
  for (key in setdiff(switch_iso3cs, c("DJI", "ERI", "ETH"))) {
    if(delay == 0) {
    switch_info[[key]] <- data.frame(t = 0, microscopy.use = 0.99)
    } else {
      if (length(unique(md$microscopy.use[md$iso3c == key])) > 1) stop("set up issue")
      switch_info[[key]] <- data.frame(t = c(0, 0+delay), microscopy.use = c(unique(md$microscopy.use[md$iso3c == key]), 0.99))
    }
  }

  # bring together with threshold_switch info which we set to infinity for these
  switch_list <- list("iso_switch" = switch_info,
                      "threshold_switch" = data.frame(t_past = delay, microscopy.use = c(0.99), f_trig = Inf))
  # set the data
  spread_model$set_map_data(md)

  # first simulate central
  if(RDT_switch) {
    med <- spread_model$simulate_spread_and_rdts(
      export_freq = 0.25, t_break = t_break, t_end = t_end, switch_list = switch_list, unc = "med"
    )
  } else {
    med <- spread_model$simulate_spread(export_freq = 0.25, t_break = t_break, unc = "med")
  }

  # then worst
  if(RDT_switch) {
    uci <- spread_model$simulate_spread_and_rdts(
      export_freq = 0.25, t_break = t_break, t_end = t_end, switch_list = switch_list, unc = "uci"
    )
  } else {
    uci <- spread_model$simulate_spread(export_freq = 0.25, t_break = t_break, unc = "uci")
  }

  # then best
  if(RDT_switch) {
    lci <- spread_model$simulate_spread_and_rdts(
      export_freq = 0.25, t_break = t_break, t_end = t_end, switch_list = switch_list, unc = "lci"
    )
  } else {
    lci <- spread_model$simulate_spread(export_freq = 0.25, t_break = t_break, unc = "lci")
  }

  out <- lci %>% select(id_1, t, freq_lci = freq) %>%
    left_join(med %>% select(id_1, t, freq_med = freq)) %>%
    left_join(uci %>% select(id_1, t, freq_uci = freq)) %>%
    select(id_1, t, matches("freq")) %>%
    mutate(delay = ifelse(RDT_switch, delay, -1), .after = "t")

  # if emulator was passed then add extra outputs
  if (!is.null(emulator)) {

    full <- left_join(out, md %>% select(-t), by = "id_1") %>%
      split(.$id_1) %>% lapply(head, emulator$extra$data$data_length)

    preds <- map(full, function(x){
      preds <- predict_from_real_data(emulator, x, old = !interp)
      preds <- cbind(x %>% select(names(out)), preds)
    })

    out <- do.call(rbind, preds)

  }

  return(out)

}

# Central fixed delay
central_row <- which(apply(scenario_maps$scenarios, 1, function(x){all(x == "central")}))
out_central_fixed <- list()
for(d in seq_along(param_grid$delay)) {
  out_central_fixed[[d]] <- run_lci_med_uci_sim_fixed(spread_model, scenario_maps$map_data[[central_row]], RDT_switch = TRUE,
                                          delay = param_grid$delay[d], t_break = t_break, emulator = emulator,
                                          f_trig = param_grid$f_trig[d]) %>%
    mutate(scenario = "Central", type = "Immediate Switch", .after = "delay")
}
# Best fixed delay
best_row <- which(apply(scenario_maps$scenarios, 1, function(row) all(row == "best")))
out_best_fixed <- list()
for(d in seq_along(param_grid$delay)) {
  out_best_fixed[[d]] <- run_lci_med_uci_sim_fixed(spread_model, scenario_maps$map_data[[best_row]], RDT_switch = TRUE,
                                       delay = param_grid$delay[d], t_break = t_break, emulator = emulator,
                                       f_trig = param_grid$f_trig[d]) %>%
    mutate(scenario = "Best", type = "Immediate Switch", .after = "delay")
}

# Worst fixed delay
worst_row <- which(apply(scenario_maps$scenarios, 1, function(row) all(row == "worst")))
out_worst_fixed <- list()
for(d in seq_along(param_grid$delay)) {
  out_worst_fixed[[d]] <- run_lci_med_uci_sim_fixed(spread_model, scenario_maps$map_data[[worst_row]], RDT_switch = TRUE,
                                        delay = param_grid$delay[d], t_break = t_break, emulator = emulator,
                                        f_trig = param_grid$f_trig[d]) %>%
    mutate(scenario = "Worst", type = "Immediate Switch", .after = "delay")
}

# -------------------------------------------------------------------------#
# 5. Combine and create plots -----------------------
# -------------------------------------------------------------------------#

# create our full data
out_central <- do.call(rbind, out_central)
out_worst <- do.call(rbind, out_worst)
out_best <- do.call(rbind, out_best)
out_central_fixed <- do.call(rbind, out_central_fixed)
out_worst_fixed <- do.call(rbind, out_worst_fixed)
out_best_fixed <- do.call(rbind, out_best_fixed)
full_df <- do.call(rbind, list(out_central, out_worst, out_best, out_central_fixed, out_worst_fixed, out_best_fixed, out_central_cf, out_worst_cf, out_best_cf))

# loop over each variable group and correctly order low, med, high etc
for (var in c("freq", trained_model$extra$output_cols)) {
  message(var)
  lci_col <- paste0(var, "_lci")
  med_col <- paste0(var, "_med")
  uci_col <- paste0(var, "_uci")

  full_df <- full_df %>%
    mutate(across(all_of(c(lci_col, med_col, uci_col)), as.numeric)) %>%
    mutate(
      sorted = pmap(list(.data[[lci_col]], .data[[med_col]], .data[[uci_col]]), ~ sort(c(...))),
      !!lci_col := map_dbl(sorted, 1),
      !!med_col := map_dbl(sorted, 2),
      !!uci_col := map_dbl(sorted, 3)
    ) %>%
    select(-sorted)
}

# simple mean in a country
full_df %>%
  mutate(year = lubridate::dyears(t) + as.Date("2023-04-01")) %>%
  left_join(afr_map %>% sf::st_drop_geometry()) %>%
  filter(delay <= 0) %>%
  group_by(year, iso, type, scenario) %>%
  summarise(freq = mean(freq_med)) %>%
  filter(year < as.Date("2064-04-01")) %>%
  ggplot(aes(year, freq, color = scenario, linetype = type)) +
  geom_line() +
  facet_wrap(~iso, scales = "free_y")

# comparitive in one country
full_df %>%
  mutate(year = lubridate::dyears(t) + as.Date("2023-04-01")) %>%
  left_join(afr_map %>% sf::st_drop_geometry()) %>%
  filter(year < as.Date("2064-04-01")) %>%
  filter(iso == "TZA", scenario == "Worst") %>%
  ggplot(aes(year, freq_med, color = as.factor(delay), linetype = type)) +
  geom_line() +
  facet_wrap(~name_1, scales = "free_y")

# comparitive in one region
full_df %>%
  mutate(year = lubridate::dyears(t) + as.Date("2023-04-01")) %>%
  left_join(afr_map %>% sf::st_drop_geometry()) %>%
  filter(year < as.Date("2064-04-01")) %>%
  filter(iso == "UGA", scenario == "Worst") %>% filter(name_1 == "Zombo") %>%
  ggplot(aes(year, freq_med, color = as.factor(delay))) +
  geom_line() +
  facet_wrap(~type)

full_df %>% filter(id_1 == 10823084) %>% pivot_longer(starts_with(c("freq", trained_model$extra$output_cols))) %>%
  tidyr::separate(name, into = c("name", "stat"), sep = "_(?=[^_]+$)") %>% pivot_wider(names_from = stat, values_from = value) %>%
  filter(scenario == "Worst") %>%
  ggplot(aes(t, med, ymin = lci, ymax = uci, color = as.factor(delay), fill = as.factor(delay))) +
  geom_ribbon(alpha=0.4) +
  geom_line() +
  facet_wrap(type~name, scales = "free", ncol = 9)

saveRDS(full_df, "analysis/impact_analysis/data-derived/chai_sims_new.rds")

# -------------------------------------------------------------------------#
# 6. Save for partners -----------------------
# -------------------------------------------------------------------------#

# previous run with only delay equal to 5 years max
full_df <- readRDS("analysis/impact_analysis/data-derived/chai_sims.rds")
dir.create("analysis/impact_analysis/data-out")
closest_to_integer <- function(vec, integers = 0:20) {
  sapply(integers, function(i) vec[which.min(abs(vec - i))])
}

annual <- full_df %>%
  mutate(t = replace(t, t>20 & t < max(.data$t), max(.data$t))) %>%
  filter(t %in% closest_to_integer(unique(t))) %>%
  mutate(t = round(t, 1)) %>%
  filter(type %in% c("5% Threshold Strategy", "No RDT Switching"))

annual_clean <- annual %>%
  mutate(across(matches("_[a-z]"), function(x){signif(x, 4)}))

annual_clean <- annual_clean %>%
  group_by(scenario, id_1) %>%
  arrange(t) %>%
  mutate(across(matches("_[a-z]"), ~ .x - (.x[t==0] - mean(.x[t ==0])))) %>%
  arrange(id_1, delay, scenario, type, t) %>%
  ungroup()

annual_clean <- annual_clean %>% left_join(afr_map %>% sf::st_drop_geometry() %>% select(iso, id_1, name_1)) %>%
  relocate(name_1, .before = 1) %>%
  relocate(iso, .before = 1)

for(i in 0:5) {
write.csv(annual_clean %>% filter(delay %in% c(i, -1)),
          paste0("analysis/impact_analysis/data-out/longitudinal_comparison_delay_", i, ".csv"),
          row.names = FALSE)
}
