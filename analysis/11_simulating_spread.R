library(tidyverse)
library(brnn)
library(caret)
library(earth)
library(scam)
library(furrr)
library(future)
library(wpp2017)
sf::sf_use_s2(FALSE)

# --------------------------------------------------------------------------#
# 1. Get our needed objects to start simulating African spread of hrp2 -----
# --------------------------------------------------------------------------#

# Get the maps for simulating
admin0 <- readRDS(here::here("analysis/data_derived/admin0_sf.rds"))
admin1 <- readRDS(here::here("analysis/data_derived/admin1_sf.rds"))

# and who for plotting
who_shps <- readRDS("analysis/data_derived/who_shps.rds")

# Get the selection model object
mod_obj <- readRDS("analysis/data_derived/ensemble_selection_model.rds")

# subset to Africa as that is where we will do the prospective mapping
data("UNlocations")
admin1$cont <- countrycode::countrycode(admin1$iso, "iso3c", "iso3n")
admin1$region <- UNlocations$area_name[match(admin1$cont, UNlocations$country_code)]
afr_map <- admin1 %>% filter(region == "Africa")

# get our non_hrp2 data
non_hrp2_use <- read.csv("analysis/data_raw/hrp2_RDT_usage.csv")
mean_afr_use <- non_hrp2_use %>% filter(iso3c %in% afr_map$iso) %>% pull(hrp2) %>% mean

# Get the scenario maps
scenario_maps <- readRDS("analysis/data_derived/scenario_maps_full.rds")

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

# set up parallel cluster here
n.cores <- 8
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)
future::plan(future::cluster, workers = my.cluster)


# Now simulate spread for all our scenarios (15m to run)
out_all <- furrr::future_map(scenario_maps$map_data, function(x){
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

}, .progress = TRUE, .options = furrr_options(packages = "hrpup"))
parallel::stopCluster(my.cluster)

# --------------------------------------------------------------------------#
# 3. Plot our outputs -----------------------------------------------------
# --------------------------------------------------------------------------#

# plotting args and themes
var  <- "freq"
title <- "False negative HRP2-RDTs \namongst clinical infections \ndue to pfhrp2/3 deletions\n"
scale_func <- function(...) {scale_fill_viridis_c(labels = scales::percent, limits = c(0,1), ...)}
custom_theming <- function() {
  theme(plot.background = element_rect(fill = "white", color = "white"),
        plot.title = element_text(hjust = 0.5, size = 16, family = "Helvetica"),
        legend.key.height = unit(0.4, "inch"),
        legend.key.width = unit(0.4, "inch"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.title.align = 0)
}

# add an NA for region in Morcco so we have NA for the legend
out <- out %>% rbind(data.frame(id_1 = 10314922,t = unique(out$t),freq=NA, t_pos=unique(out$t_pos)))

# Make a map at each key t interval
pl_list <- vector("list", length(unique(out$t)))
for(i in c(1, as.integer(length(unique(out$t))/4)+1, as.integer(length(unique(out$t))/2)+1)) {
  x <- unique(out$t)[i]
  pl_list[[i]] <-
    who_compliant_continuous_plot(
      out %>% filter(t == x), var = var, title = title,
      who_shps = who_shps, region = "africa",
      limits = TRUE, scale_fill_func = scale_func) +
    ggtitle(as.integer(2023+x)) +
    custom_theming()
}


# set up parallel cluster here
n.cores <- 14
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)
future::plan(future::cluster, workers = my.cluster)

# Plot all to the time series directory
dir.create("analysis/plots/time_series")
out_plots <- furrr::future_map(seq_along(unique(out$t)), function(i) {

  x <- unique(out$t)[i]
  gg_x <- who_compliant_continuous_plot(
    out %>% filter(t == x), var = var, title = title,
    who_shps = who_shps, region = "africa",
    limits = TRUE, scale_fill_func = scale_func) +
    ggtitle(as.integer(2023+x)) +
    custom_theming()

  hrpup:::save_figs(paste0("test_",sprintf("%03d",i)), fig = gg_x,
            width = 12, height = 12, pdf_plot = FALSE, plot_dir = "analysis/plots/time_series/",
            font_family = "Helvetica")

  return(1L)

}, .progress = TRUE, .options = furrr_options(packages = "hrpup"))
parallel::stopCluster(my.cluster)

# Make a movie of these
mapmate::ffmpeg(dir = "analysis/plots/time_series/",
                output_dir = "analysis/plots/time_series/",
                pattern = "test_%03d.png", output = "test_video2.mp4",
                delay = 1/12, overwrite = TRUE
                )

file.remove(grep("test_\\d", list.files("analysis/plots/time_series/", full.names = TRUE), value = TRUE))
# --------------------------------------------------------------------------#
# 4. Plot our output figure -----------------------------------------------------
# --------------------------------------------------------------------------#

## 4.a. Plot our national data -----------------------------------------------------

# Get the starting data
who_cols <- c("#b4c6e7", "#fefd81ff", "#e9a934ff", "#c3271bff","#5e3218ff")

who_scale_fill_func <- function(...){
  ggplot2::scale_fill_manual(
    values = c("#b4c6e7", "#fefd81ff", "#e9a934ff", "#c3271bff","#5e3218ff"),
    labels = (c("0-1%", ">1-8%",">8-15%",">15-25%",">25%")),...)}

nat_df <- left_join(who_shps$admin1, seeds) %>%
  select(id_1, prev) %>%
  sf::st_drop_geometry() %>%
  mutate(prev_break = cut(prev, c(0,1,8,15,25,100)/100, include.lowest = FALSE))

# national level plot
gg_1 <- who_compliant_discrete_plot(nat_df, var = "prev_break",
                                    title = "WHO Threat Map \npfhrp2 deletion \nprevalence",
                                    who_shps = who_shps, region = "africa",
                                    limits = TRUE, scale_fill_func = who_scale_fill_func) +
  ggtitle(as.integer(2023)) +
  custom_theming() +
geom_rect(xmin = x_coords[1], ymin = y_coords[1],
          xmax = x_coords[2], ymax = y_coords[2],
          fill = NA, colour = "black", size = 2)

## 4.b. Plot our initial spread -----------------------------------------------------

# countries to zoom in on
hoa_iso <- c("ERI", "DJI", "ETH", "SOM", "KEN", "UGA", "SSD","SDN")
x_coords <- c(33.24, 51.88)
y_coords <- c(1.80, 19.24)

# initial spread plot
gg_2 <- out %>% filter(t_pos == 11) %>% select(-t) %>%
  mutate(prev_break = cut(freq, c(0,1,8,15,25,100)/100, include.lowest = FALSE)) %>%
  who_compliant_discrete_plot(var = "prev_break",
                              title = "WHO Threat Map \npfhrp2 deletion \nprevalence",
                              who_shps = who_shps, region = "africa",
                              limits = TRUE, scale_fill_func = who_scale_fill_func) +
  ggtitle(as.integer(2024)) +
  custom_theming() +
  xlim(x_coords) +
  ylim(y_coords) +
  theme(panel.border = element_rect(color = "black", linewidth = 3, fill = NA)) +
  ggtitle("2024\n") +
  ggplot2::geom_sf(data = gg_2$data %>% filter(iso %in% hoa_iso), fill = NA,
                   color = "#696969ff", show.legend = FALSE, lwd = 0.05, inherit.aes = FALSE)

## 4.c. Plot our successive spread -----------------------------------------------------

gg_3 <- out %>% filter(t_pos == 29) %>% select(-t) %>%
  mutate(prev_break = cut(freq, c(0,1,8,15,25,100)/100, include.lowest = FALSE)) %>%
  who_compliant_discrete_plot(var = "prev_break",
                              title = "WHO Threat Map \npfhrp2 deletion \nprevalence",
                              who_shps = who_shps, region = "africa",
                              limits = TRUE, scale_fill_func = who_scale_fill_func) +
  ggtitle(as.integer(2024)) +
  custom_theming() +
  xlim(x_coords) +
  ylim(y_coords) +
  theme(panel.border = element_rect(color = "black", linewidth = 3, fill = NA)) +
  ggtitle("2024\n") +
  ggplot2::geom_sf(data = gg_2$data %>% filter(iso %in% hoa_iso), fill = NA,
                   color = "#696969ff", show.legend = FALSE, lwd = 0.05, inherit.aes = FALSE)

## 4.a-f. Bring all together -----------------------------------------------------

# Bring top row together
top_row_gg <- cowplot::plot_grid(
  cowplot::get_legend(gg_3 + theme(text = element_text("Helvetica"))),
  add_africa_scale(gg_1) + theme(legend.position = "none", text = element_text("Helvetica")),
  gg_2 + theme(legend.position = "none", text = element_text("Helvetica")),
  gg_3 + theme(legend.position = "none", text = element_text("Helvetica")),
  ncol = 4, rel_widths = c(0.5, 1,1,1),
  labels = c("",LETTERS[1:3]), scale = 0.9, label_size = 18)

# Bring bottom row together
bottom_row_gg <- cowplot::plot_grid(
  add_africa_scale(pl_list[[1]]) + theme(legend.position = "none", text = element_text("Helvetica")) + ggtitle("2023"),
  add_africa_scale(pl_list[[as.integer(length(pl_list)/4)+1]]) +
    theme(legend.position = "none", text = element_text("Helvetica")) + ggtitle("2033"),
  add_africa_scale(pl_list[[as.integer(length(pl_list)/2)+1]]) +
    theme(legend.position = "none", text = element_text("Helvetica")) + ggtitle("2043"),
  cowplot::get_legend(pl_list[[1]] + theme(text = element_text("Helvetica"))),
  ncol = 4, rel_widths = c(1,1,1,0.5),
  labels = c(LETTERS[4:6], ""), label_size = 18)

spread_gg <- cowplot::plot_grid(top_row_gg, bottom_row_gg, ncol = 1, rel_heights = c(1,1.3)) +
  theme(plot.background = element_rect(fill = "white", color ="white"))
save_figs("spread_africa", spread_gg, width = 22, height = 14, pdf_plot = FALSE, font_family = "Helvetica")

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
saveRDS(compl_df, "analysis/data_derived/complete_prospective_times.rds")

isos <- unique(countrycode::codelist$iso3c[countrycode::codelist$continent == "Africa"])
# Times censored at 40 years
central_df <- compl_df %>% group_by(iso, name_0, name_1) %>%
  filter(if_all(ends_with("scen"), ~ . == "central")) %>%
  filter(!(iso %in% c("ATF", "UMI"))) %>%
  mutate(continent = countrycode::countrycode(iso, "iso3c", "continent"), .before = 1) %>%
  select(continent:tmax) %>%
  mutate(t_med = round(replace(t, is.na(t) & iso %in% isos, Inf),2),
         t_uci = round(replace(tmin, is.na(t) & iso %in% isos, Inf),2),
         t_lci = round(replace(tmax, is.na(t) & iso %in% isos, Inf),2),
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
  mutate(t_med = round(replace(t, is.na(t) & iso %in% isos, Inf),2),
         t_uci = round(replace(tmin, is.na(t) & iso %in% isos, Inf),2),
         t_lci = round(replace(tmax, is.na(t) & iso %in% isos, Inf),2),
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
  mutate(t_med = round(replace(t, is.na(t) & iso %in% isos, Inf),2),
         t_uci = round(replace(tmin, is.na(t) & iso %in% isos, Inf),2),
         t_lci = round(replace(tmax, is.na(t) & iso %in% isos, Inf),2),
         scenario = "Optimistic") %>%
  mutate(t_med = replace(t_med, t_med>40, "40+"),
         t_uci = replace(t_uci, t_uci>40, "40+"),
         t_lci = replace(t_lci, t_lci>40, "40+")) %>%
  rename(admin_0 = name_0) %>%
  rename(admin_1 = name_1) %>%
  arrange(continent, iso,  admin_1) %>%
  select(continent, iso, admin_0, admin_1, t_lci, t_med, t_uci, scenario)

write.csv(central_df, "analysis/data_out/central_times_prospective.csv", row.names = FALSE)
write.csv(worst_df, "analysis/data_out/pessimistic_times_prospective.csv", row.names = FALSE)
write.csv(best_df, "analysis/data_out/optimistic_times_prospective.csv", row.names = FALSE)

# ---------------------------------------------------- #
# 6. Creating Key High Risk Countries ------------------
# ---------------------------------------------------- #

isos <- unique(countrycode::codelist$iso3c[countrycode::codelist$continent == "Africa"])
world <- left_join(scenario_maps$map, out_all[[365]]) %>%
  mutate(t = replace(t, t < 0, Inf)) %>%
  mutate(t = replace_na(t,Inf)) %>%
  #mutate(t = replace(t, s<0, Inf)) %>%
  filter(!is.na(t)) %>%
  filter(iso %in% isos) %>%
  mutate(t_bin = cut(t, breaks = c(0,6,12,20, Inf), include.lowest = TRUE, right = FALSE))

hrtable <- world %>% sf::st_drop_geometry() %>% group_by(name_0) %>%
  summarise(r = sum(t_bin == "[0,6)")/n()) %>%
  filter(r >= 0.5) %>%
  arrange(desc(r)) %>%
  mutate(r = scales::percent_format()(r)) %>%
  setNames(c("Country", "% Admin 1 High Risk"))

write.csv(hrtable, file = "analysis/tables/high_prospective_risk_countries.csv")

# ---------------------------------------------------- #
# 7. Adding prospective data to the hrp2 map object ------
# ---------------------------------------------------- #

# read the map in
hrp2_map <- readRDS("analysis/data_derived/R6_map.rds")

# create our data
cleaned_prospective_map_data <- lapply(out_all, function(x){ x %>%
  mutate(t = replace(t, t < 0, Inf)) %>%
  mutate(t = replace_na(t, Inf)) %>%
  # mutate(t = replace(t, s<0, Inf)) %>%
  mutate(t_bin = cut(t, breaks = c(0,6,12,20, Inf), include.lowest = TRUE, right = FALSE)) %>%
  select(id_1, t_bin) %>%
  rename(hrp2_prospective_risk = t_bin) %>%
  mutate(hrp2_prospective_risk = rev(c("Marginal", "Slight", "Moderate", "High"))[as.integer(hrp2_prospective_risk)]) %>%
  mutate(hrp2_prospective_risk = factor(hrp2_prospective_risk, levels = rev(c("Marginal", "Slight", "Moderate", "High"))))
})

new_map_data <- map(seq_along(hrp2_map$.__enclos_env__$private$map_data), function(x){
  left_join(hrp2_map$.__enclos_env__$private$map_data[[x]],
            cleaned_prospective_map_data[[x]], by = "id_1")
})

# set the prospective data
hrp2_map$set_map_data(new_map_data)
saveRDS(hrp2_map, "analysis/data_derived/R6_map.rds")

# Simple comparison plot for ease
risk_gg1 <- hrp2_map$plot(scale_bar_bool = FALSE) + custom_theming() +
  theme(legend.key.spacing.y = unit(0.25, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.spacing = unit(0.125, "cm"))
risk_gg2 <- hrp2_map$plot(risk = "prospective",scale_bar_bool = FALSE) + custom_theming() +
  theme(legend.key.spacing.y = unit(0.25, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.spacing = unit(0.125, "cm"))
risk_gg <- cowplot::plot_grid(add_africa_scale(risk_gg1) + theme(legend.position = "none") + ggtitle("Innate Risk"),
                              add_africa_scale(risk_gg2) + theme(legend.position = "none") + ggtitle("Prospective Risk"),
                   cowplot::get_legend(risk_gg2 + theme(text = element_text("Helvetica"))),
                   labels = c("A", "B", ""),
                   rel_widths = c(1, 1, 0.5),
                   ncol = 3) +
  theme(plot.background = element_rect("white","white"))

save_figs("innate_vs_prospective1", risk_gg, width = 12, height = 6, pdf_plot = FALSE)
