library(tidyverse)
sf::sf_use_s2(FALSE)

# --------------------------------------------------------------------------#
# 1. Get our needed objects to start simualting African spread of hrp2 -----
# --------------------------------------------------------------------------#

# Get the mapping object
map_obj <- readRDS("analysis/data_derived/R6_map.rds")

# Get the selection model object
mod_obj <- readRDS("analysis/data_derived/ensemble_selection_model.rds")

# subset to Africa as that is where we will do the composite mapping
afr_map <- map_obj$.__enclos_env__$private$map %>% filter(region == "Africa")

# Get the scenario maps
scenario_maps <- readRDS("analysis/data_derived/scenario_maps_full_monotonic.rds")
map_data <- scenario_maps$map_data[[365]] %>% filter(id_1 %in% afr_map$id_1)
map_data <- filter(map_data, !is.na(s))

# And subset the afr_map to just those regions with selection data
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
# 2. Run our simulation -----------------------------------------------------
# --------------------------------------------------------------------------#

# Initialise our model
spread_model <- hrpup:::R6_hrp2_spread$new(afr_map, mod_obj, adj_mat = adj_mat)

# Set up our initial conditions
spread_model$set_seeds(setNames(seeds_adm1$prev, seeds_adm1$id_1))

# Set up our selection speed map data
spread_model$set_map_data(scenario_maps$map_data[[365]] %>% filter(id_1 %in% afr_map$id_1))

# Simulate our data
out <- spread_model$simulate_spread(export_freq = 0.25, t_break = 0.1)

# Now simulate spread for all our scenarios (15m to run)
out_all <- map(scenario_maps$map_data, function(x){
  spread_model$set_map_data(x %>% filter(id_1 %in% afr_map$id_1))
  central <- spread_model$simulate_spread(export_freq = 0.25, t_break = 0.1) %>%
    group_by(id_1) %>%
    summarise(t = t[freq > 0.05][1]) %>%
    left_join(x %>% filter(id_1 %in% afr_map$id_1) %>% select(id_1, s))
  spread_model$set_map_data(x %>% filter(id_1 %in% afr_map$id_1) %>% mutate(s = smin))
  worst <- spread_model$simulate_spread(export_freq = 0.25, t_break = 0.1) %>%
    group_by(id_1) %>%
    summarise(tmin = t[freq > 0.05][1]) %>%
    left_join(x %>% filter(id_1 %in% afr_map$id_1) %>% select(id_1, smin))
  spread_model$set_map_data(x %>% filter(id_1 %in% afr_map$id_1) %>% mutate(s = smax))
  best <- spread_model$simulate_spread(export_freq = 0.25, t_break = 0.1) %>%
    group_by(id_1) %>%
    summarise(tmax = t[freq > 0.05][1]) %>%
    left_join(x %>% filter(id_1 %in% afr_map$id_1) %>% select(id_1, smax))

  out_test <- left_join(central, worst) %>% left_join(best) %>%
    select(id_1, s, smin, smax, t, tmin, tmax)

  return(out_test)

}, .progress = TRUE)


# --------------------------------------------------------------------------#
# 3. Plot our output video -----------------------------------------------------
# --------------------------------------------------------------------------#

isos <- unique(countrycode::codelist$iso3c[countrycode::codelist$continent == "Africa"])

# Make a map at each t interval
pl_list <- map(unique(out$t), function(x) {
  map_obj$.__enclos_env__$private$map %>% filter(id_1 %in% out$id_1) %>%
  left_join(out %>% filter(t == x) %>% select(-t) %>% rename(t = freq)) %>%
  ggplot() +
  geom_sf(aes(fill = t), color = "grey", show.legend = TRUE, lwd = 0.1) +
  geom_sf(fill = NA, color = "black", show.legend = FALSE,
          data = map_obj$.__enclos_env__$private$map_0 %>% filter(iso %in% isos), lwd = 0.2) +
  coord_sf() +
  theme_void(base_size = 16) +
    ggtitle(as.integer(2023+x)) +
  theme(plot.caption = element_text(face = "italic"), plot.background = element_rect(fill = "white", color = "white")) +
  scale_fill_viridis_c(name = "False Negative HRP2-RDTs \namongst clinical infections \ndue to pfhrp2/3 deletions\n",
                       labels = scales::percent, limits = c(0,1)) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.key.height = unit(1, "inch"),
          legend.key.width = unit(0.5, "inch"),
          legend.text = element_text(size = 14),
          legend.title.align = 0)
})

# Plot these to the time series directory
dir.create("analysis/plots/time_series")
for(i in seq_along(pl_list)) {
  save_figs(paste0("test_",sprintf("%03d",i)), fig = pl_list[[i]],
          width = 12, height = 12, pdf_plot = FALSE, plot_dir = "analysis/plots/time_series/",
          font_family = "Helvetica")
}

# Make a movie of these
mapmate::ffmpeg(dir = "analysis/plots/time_series/",
                output_dir = "analysis/plots/time_series/",
                pattern = "test_%03d.png", output = "test_video.mp4",
                delay = 1/12, overwrite = TRUE
                )

# --------------------------------------------------------------------------#
# 4. Plot our output figure -----------------------------------------------------
# --------------------------------------------------------------------------#

# Get the starting data
who_cols <- c("#b4c6e7", "#fefd81ff", "#e9a934ff", "#c3271bff","#5e3218ff")

# coordinates for zoom in
x_coords <- c(33.24, 51.88)
y_coords <- c(1.80, 19.24)

gg_1 <- left_join(map_obj$.__enclos_env__$private$map_0, seeds) %>%
  filter(region == "Africa") %>%
  mutate(prev_break = cut(prev, c(0,1,8,15,25,100)/100, include.lowest = FALSE)) %>%
  ggplot() +
  geom_sf(fill = NA, color = "black", show.legend = FALSE,
          data = map_obj$.__enclos_env__$private$map_0 %>% filter(iso %in% isos), lwd = 0.2) +
  geom_sf(aes(fill = prev_break), color = "black", show.legend = TRUE, lwd = 0.2) +
  geom_rect(xmin = x_coords[1], ymin = y_coords[1],
            xmax = x_coords[2], ymax = y_coords[2],
            fill = NA, colour = "black", size = 2) +
  coord_sf() +
  theme_void(base_size = 16) +
  theme(plot.caption = element_text(face = "italic"), plot.background = element_rect(fill = "white", color = "white")) +
  scale_fill_manual(name = "WHO Threat Map \npfhrp2 deletion \nprevalence",
                    values = who_cols, na.value = "#f0f0f0ff",
                    labels = c("0-1%", ">1-8%",">8-15%",">15-25%",">25%", "No Data")) +
  ggtitle("2023\n") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(1, "inch"),
        legend.key.width = unit(0.5, "inch"),
        legend.text = element_text(size = 14),
        legend.title.align = 0)

# Zoom in on the starting data and start showing spread
hoa_iso <- c("ERI", "DJI", "ETH", "SOM", "KEN", "UGA", "SSD","SDN")
gg_2 <- map_obj$.__enclos_env__$private$map %>%
  filter(iso %in% hoa_iso) %>%
  left_join(out %>% filter(t_pos == 11) %>% select(-t) %>% rename(prev = freq)) %>%
  filter(region == "Africa") %>%
  mutate(prev_break = cut(prev, c(0,1,8,15,25,100)/100, include.lowest = FALSE)) %>%
  ggplot() +
  geom_sf(aes(fill = prev_break), color = "black", show.legend = TRUE, lwd = 0.1) +
  geom_sf(fill = NA, color = "black", show.legend = FALSE,
          data = map_obj$.__enclos_env__$private$map_0 %>%
            filter(iso %in% hoa_iso),
          lwd = 0.2) +
  coord_sf() +
  theme_void(base_size = 16) +
  theme(plot.caption = element_text(face = "italic"), plot.background = element_rect(fill = "white", color = "white")) +
  scale_fill_manual(name = "WHO Threat Map \npfhrp2 deletion \nprevalence",
                    values = who_cols, na.value = "#f0f0f0ff",
                    labels = c("0-1%", ">1-8%",">8-15%",">15-25%",">25%", "No Data"),
                    drop = FALSE) +
  xlim(x_coords) +
  ylim(y_coords) +
  ggtitle("2024\n") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(1, "inch"),
        legend.key.width = unit(0.5, "inch"),
        legend.text = element_text(size = 14),
        legend.title.align = 0,
        panel.border = element_rect(color = "black", linewidth = 3, fill = NA))

# Zoom in on the starting data and start showing spread
gg_3 <- map_obj$.__enclos_env__$private$map %>%
  filter(iso %in% hoa_iso) %>%
  left_join(out %>% filter(t_pos == 28) %>% select(-t) %>% rename(prev = freq)) %>%
  filter(region == "Africa") %>%
  mutate(prev_break = cut(prev, c(0,1,8,15,25,100)/100, include.lowest = FALSE)) %>%
  ggplot() +
  geom_sf(aes(fill = prev_break), color = "black", show.legend = TRUE, lwd = 0.1) +
  geom_sf(fill = NA, color = "black", show.legend = FALSE,
          data = map_obj$.__enclos_env__$private$map_0 %>%
            filter(iso %in% hoa_iso),
          lwd = 0.2) +
  coord_sf() +
  theme_void(base_size = 16) +
  theme(plot.caption = element_text(face = "italic"), plot.background = element_rect(fill = "white", color = "white")) +
  scale_fill_manual(name = "WHO Threat Map \npfhrp2 deletion \nprevalence",
                    values = who_cols, na.value = "#f0f0f0ff",
                    labels = c("0-1%", ">1-8%",">8-15%",">15-25%",">25%", "No Data"),
                    drop = FALSE) +
  xlim(x_coords) +
  ylim(y_coords) +
  ggtitle("2026\n") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(1, "inch"),
        legend.key.width = unit(0.5, "inch"),
        legend.text = element_text(size = 14),
        legend.title.align = 0,
        panel.border = element_rect(color = "black", linewidth = 3, fill = NA))

# Bring top row together
top_row_gg <- cowplot::plot_grid(
  cowplot::get_legend(gg_3 + theme(text = element_text("Helvetica"))),
  gg_1 + theme(legend.position = "none", text = element_text("Helvetica")),
  gg_2 + theme(legend.position = "none", text = element_text("Helvetica")),
  gg_3 + theme(legend.position = "none", text = element_text("Helvetica")),
  ncol = 4, rel_widths = c(0.5, 1,1,1),
  labels = c("",LETTERS[1:3]), scale = 0.9, label_size = 18)
bottom_row_gg <- cowplot::plot_grid(
  pl_list[[1]] + theme(legend.position = "none", text = element_text("Helvetica")) + ggtitle("April 2023"),
  pl_list[[as.integer(length(pl_list)/4)+1]] +
    theme(legend.position = "none", text = element_text("Helvetica")) + ggtitle("April 2033"),
  pl_list[[as.integer(length(pl_list)/2)+1]] +
    theme(legend.position = "none", text = element_text("Helvetica")) + ggtitle("April 2043"),
  cowplot::get_legend(pl_list[[1]] + theme(text = element_text("Helvetica"))),
  ncol = 4, rel_widths = c(1,1,1,0.5),
  labels = c(LETTERS[4:6], ""), label_size = 18)

spread_gg <- cowplot::plot_grid(top_row_gg, bottom_row_gg, ncol = 1, rel_heights = c(1,1.2)) +
  theme(plot.background = element_rect(fill = "white", color ="white"))
save_figs("spread_africa", spread_gg, width = 25, height = 16, pdf_plot = FALSE, font_family = "Helvetica")

# ---------------------------------------------------- #
# 5. Creating Data Outputs  ----
# ---------------------------------------------------- #

compl <- list()
for(i in seq_along(out_all)) {
  compl[[i]] <- left_join(scenario_maps$map, out_all[[i]], by = "id_1") %>%
    sf::st_drop_geometry() %>%
    select(iso, name_0, name_1, t, tmin, tmax, s, smin, smax) %>%
    mutate(Micro.2.10_scen = scenario_maps$scenarios$Micro.2.10[i],
           ft_scen = scenario_maps$scenarios$ft[i],
           microscopy.use_scen = scenario_maps$scenarios$microscopy.use[i],
           rdt.nonadherence_scen = scenario_maps$scenarios$rdt.nonadherence[i],
           fitness_scen = scenario_maps$scenarios$fitness[i],
           rdt.det_scen = scenario_maps$scenarios$rdt.det[i])
}

compl_df <- do.call(rbind, compl)
saveRDS(compl_df, "analysis/data_derived/complete_composite_times.rds")

# Times censored at 40 years
central_df <- compl_df %>% group_by(iso, name_0, name_1) %>%
  filter(if_all(ends_with("scen"), ~ . == "central")) %>%
  filter(!(iso %in% c("ATF", "UMI"))) %>%
  mutate(continent = countrycode::countrycode(iso, "iso3c", "continent"), .before = 1) %>%
  select(continent:tmax) %>%
  mutate(t_med = round(replace_na(t, Inf),2),
         t_uci = round(replace_na(tmin, Inf),2),
         t_lci = round(replace_na(tmax, Inf),2),
         scenario = "Central") %>%
  mutate(t_med = replace(t_med, t_med>40, "40+"),
         t_uci = replace(t_uci, t_uci>40, "40+"),
         t_lci = replace(t_lci, t_lci>40, "40+")) %>%
  rename(admin_0 = name_0) %>%
  rename(admin_1 = name_1) %>%
  arrange(continent, iso,  admin_1) %>%
  select(continent, iso, admin_0, admin_1, t_lci, t_med, t_uci, scenario)

worst_df <- compl_df %>% group_by(iso, name_0, name_1) %>%
  filter(if_all(ends_with("scen"), ~ . == "worst")) %>%
  filter(!(iso %in% c("ATF", "UMI"))) %>%
  mutate(continent = countrycode::countrycode(iso, "iso3c", "continent"), .before = 1) %>%
  select(continent:tmax) %>%
  mutate(t_med = round(replace_na(t, Inf),2),
         t_uci = round(replace_na(tmin, Inf),2),
         t_lci = round(replace_na(tmax, Inf),2),
         scenario = "Pessimistic") %>%
  mutate(t_med = replace(t_med, t_med>40, "40+"),
         t_uci = replace(t_uci, t_uci>40, "40+"),
         t_lci = replace(t_lci, t_lci>40, "40+")) %>%
  rename(admin_0 = name_0) %>%
  rename(admin_1 = name_1) %>%
  arrange(continent, iso,  admin_1) %>%
  select(continent, iso, admin_0, admin_1, t_lci, t_med, t_uci, scenario)

best_df <- compl_df %>% group_by(iso, name_0, name_1) %>%
  filter(if_all(ends_with("scen"), ~ . == "best")) %>%
  filter(!(iso %in% c("ATF", "UMI"))) %>%
  mutate(continent = countrycode::countrycode(iso, "iso3c", "continent"), .before = 1) %>%
  select(continent:tmax) %>%
  mutate(t_med = round(replace_na(t, Inf),2),
         t_uci = round(replace_na(tmin, Inf),2),
         t_lci = round(replace_na(tmax, Inf),2),
         scenario = "Optimistic") %>%
  mutate(t_med = replace(t_med, t_med>40, "40+"),
         t_uci = replace(t_uci, t_uci>40, "40+"),
         t_lci = replace(t_lci, t_lci>40, "40+")) %>%
  rename(admin_0 = name_0) %>%
  rename(admin_1 = name_1) %>%
  arrange(continent, iso,  admin_1) %>%
  select(continent, iso, admin_0, admin_1, t_lci, t_med, t_uci, scenario)

write.csv(central_df, "analysis/tables/central_times_composite.csv", row.names = FALSE)
write.csv(worst_df, "analysis/tables/pessimistic_times_composite.csv", row.names = FALSE)
write.csv(best_df, "analysis/tables/optimistic_times_composite.csv", row.names = FALSE)

# ---------------------------------------------------- #
# 6. Creating Key High Risk Countries ------------------
# ---------------------------------------------------- #

isos <- unique(countrycode::codelist$iso3c[countrycode::codelist$continent == "Africa"])
world <- left_join(scenario_maps$map, out_all[[365]]) %>%
  mutate(t = replace(t, t < 0, Inf)) %>%
  mutate(t = replace_na(t,Inf)) %>%
  mutate(t = replace(t, s<0, Inf)) %>%
  filter(!is.na(t)) %>%
  filter(iso %in% isos) %>%
  mutate(t_bin = cut(t, breaks = c(0,6,12,20, Inf), include.lowest = TRUE))

hrtable <- world %>% sf::st_drop_geometry() %>% group_by(name_0) %>%
  summarise(r = sum(t_bin == "[0,6]")/n()) %>%
  filter(r >= 0.5) %>%
  arrange(desc(r)) %>%
  mutate(r = scales::percent_format()(r)) %>%
  setNames(c("Country", "% Admin 1 High Risk"))

write.csv(hrtable, file = "analysis/tables/high_composite_risk_countries.csv")

# ---------------------------------------------------- #
# 7. Adding composite data to the hrp2 map object ------
# ---------------------------------------------------- #

# read the map in
hrp2_map <- readRDS("analysis/data_derived/R6_map.rds")

# create our data
cleaned_composite_map_data <- lapply(out_all, function(x){ x %>%
  mutate(t = replace(t, t < 0, Inf)) %>%
  mutate(t = replace_na(t, Inf)) %>%
  mutate(t = replace(t, s<0, Inf)) %>%
  mutate(t_bin = cut(t, breaks = c(0,6,12,20, Inf), include.lowest = TRUE, right = FALSE)) %>%
  select(id_1, t_bin) %>%
  rename(hrp2_composite_risk = t_bin) %>%
  mutate(hrp2_composite_risk = rev(c("Marginal", "Slight", "Moderate", "High"))[as.integer(hrp2_composite_risk)]) %>%
  mutate(hrp2_composite_risk = factor(hrp2_composite_risk, levels = rev(c("Marginal", "Slight", "Moderate", "High"))))
})

new_map_data <- map(seq_along(hrp2_map$.__enclos_env__$private$map_data), function(x){
  left_join(hrp2_map$.__enclos_env__$private$map_data[[x]],
            cleaned_composite_map_data[[x]], by = "id_1")
})

# set the composite data
hrp2_map$set_map_data(new_map_data)
saveRDS(hrp2_map, "analysis/data_derived/R6_map.rds")

# Simple comparison plot for ease
risk_gg1 <- hrp2_map$plot()
risk_gg2 <- hrp2_map$plot(risk = "composite")
risk_gg <- cowplot::plot_grid(risk_gg1 + theme(legend.position = "none"),
                   risk_gg2 + theme(legend.position = "none"),
                   cowplot::get_legend(risk_gg2 + theme(text = element_text("Helvetica"))),
                   labels = c("A", "B", ""),
                   rel_widths = c(1, 1, 0.3),
                   ncol = 3) +
  theme(plot.background = element_rect("white","white"))

save_figs("innate_vs_composite", risk_gg, width = 10, height = 6, pdf_plot = FALSE)
