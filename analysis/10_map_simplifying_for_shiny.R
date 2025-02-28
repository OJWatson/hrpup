library(tidyverse)

# 1. Source objects to make mapping object
scenario_maps <- readRDS("analysis/data_derived/scenario_maps_full.rds")
who_shps <- readRDS("analysis/data_derived/who_shps.rds")

# 2. streamlined data
for(i in seq_along(scenario_maps$map_data)) {

  huh <- scenario_maps$map_data[[i]] %>% select(id_1, Micro.2.10, t) %>%
    mutate(t = replace(t, t < 0, Inf)) %>%
    mutate(t_bin = cut(t, breaks = c(0,6,12,20, Inf), include.lowest = TRUE, right = FALSE)) %>%
    select(id_1,  Micro.2.10, t_bin) %>%
    rename(hrp2_risk = t_bin)

  huh$hrp2_risk <- rev(c("Marginal", "Slight", "Moderate", "High"))[as.integer(huh$hrp2_risk)]
  huh$hrp2_risk <- factor(huh$hrp2_risk, levels = rev(c("Marginal", "Slight", "Moderate", "High")))
  scenario_maps$map_data[[i]] <- huh

}

# 3. create our smaller object for file
hrp2_map <- hrpup:::R6_hrp2_map$new(
  map = who_shps$admin1,
  map_data = scenario_maps$map_data,
  scenarios = scenario_maps$scenarios,
  map_0 = who_shps$admin0,
  lakes = who_shps$lakes,
  disputed_areas = who_shps$disputed_areas,
  disputed_borders = who_shps$disputed_borders
)

saveRDS(hrp2_map, "analysis/data_derived/R6_map.rds")
