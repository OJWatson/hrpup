# ---------------------------------------------------- #
# 1. Get the world map and merge together with our covariate ranges
# ---------------------------------------------------- #

# Grab our covariate parameter ranges
covars <-  readRDS("analysis/data_derived/global_covariate_ranges.rds")

# Get our selection model
selection_model <- readRDS("analysis/data_derived/ensemble_selection_model.rds")

# MAP world map to
world_map <- malariaAtlas::getShp(ISO = na.omit(unique(prev$ISO3)), admin_level = c("admin1")) %>% sf::st_as_sf()

# make scenarios for the map
scenarios <- expand.grid(
  "prev" = c("worst", "central", "best"),
  "treat" = c("worst", "central", "best"),
  "micro" = c("worst", "central", "best"),
  "nonadherence" = c("worst", "central", "best"),
  "fitness" = c("worst", "central", "best"),
  "rdt.det"= c("worst", "central", "best")
) %>% setNames(names(formals(selection_model$predict)))

# ---------------------------------------------------- #
# 2. Estimate times to 5% for each scenario worldwide
# ---------------------------------------------------- #

# make scenario map data
map_data <- vector("list", nrow(scenarios))

# loop over and generate selection and times
for(i in seq_along(map_data)) {

  micro210 <- case_when(scenarios$Micro.2.10[i] == "central" ~ covars$Micro.2.10_mean,
                        scenarios$Micro.2.10[i] == "worst" ~ covars$Micro.2.10_low,
                        scenarios$Micro.2.10[i] == "best" ~ covars$Micro.2.10_high)

  ft <- case_when(scenarios$ft[i] == "central" ~ covars$ft_mean,
                  scenarios$ft[i] == "worst" ~ covars$ft_high,
                  scenarios$ft[i] == "best" ~ covars$ft_low)

  microscopy <- case_when(scenarios$microscopy.use[i] == "central" ~ covars$microscopy.use_mean,
                          scenarios$microscopy.use[i] == "worst" ~ covars$microscopy.use_low,
                          scenarios$microscopy.use[i] == "best" ~ covars$microscopy.use_high)

  nonadherence <- case_when(scenarios$rdt.nonadherence[i] == "central" ~ covars$rdt.nonadherence_mean,
                            scenarios$rdt.nonadherence[i] == "worst" ~ covars$rdt.nonadherence_low,
                            scenarios$rdt.nonadherence[i] == "best" ~ covars$rdt.nonadherence_high)

  fitness <- case_when(scenarios$fitness[i] == "central" ~ covars$fitness_mean,
                       scenarios$fitness[i] == "worst" ~ covars$fitness_high,
                       scenarios$fitness[i] == "best" ~ covars$fitness_low)

  rdtdet <- case_when(scenarios$rdt.det[i] == "central" ~ covars$rdt.det_mean,
                      scenarios$rdt.det[i] == "worst" ~ covars$rdt.det_low,
                      scenarios$rdt.det[i] == "best" ~ covars$rdt.det_high)


  map_data[[i]] <- selection_model$predict(micro210, ft, microscopy, nonadherence, fitness, rdtdet)
  map_data[[i]] <- map_data[[i]] %>%
    mutate(id_1 = covars$id_1, .before = 1) %>%
    mutate(iso3c = covars$iso3c, .before = 1)
}

# group together and save
scenario_maps <- list("scenarios" = scenarios,
                      "map_data" = map_data,
                      "map" = world_map)
saveRDS(scenario_maps, "analysis/data_derived/scenario_maps_full_monotonic.rds")
scenario_maps <- readRDS("analysis/data_derived/scenario_maps_full_monotonic.rds")

# now just save the t
for(i in seq_along(map_data)) {
  map_data[[i]] <- map_data[[i]] %>%
    select(id_1, t)
}
scenario_maps <- list("scenarios" = scenarios,
                      "map_data" = map_data,
                      "map" = world_map)
saveRDS(scenario_maps, "analysis/data_derived/scenario_maps_monotonic.rds")

