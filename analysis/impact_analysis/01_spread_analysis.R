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
data("UNlocations")
sf::sf_use_s2(FALSE)

# --------------------------------------------------------------------------#
# 1. Get our needed objects to start simulating African spread of hrp2 -----
# --------------------------------------------------------------------------#

# Get the mapping object
map_obj <- readRDS("analysis/data_derived/R6_map.rds")
admin1 <- readRDS(here::here("analysis/data_derived/admin1_sf.rds"))

# Get the selection model object
mod_obj <- readRDS("analysis/data_derived/ensemble_selection_model.rds")

# subset to Africa as that is where we will do the prospective mapping
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
# in April 2027
switch <- as.integer((as.Date("2027-01-01") - as.Date("2024-01-01")))/365
for (key in setdiff(switch_iso3s, c("DJI", "ERI", "ETH"))) {
  switch_info[[key]] <- data.frame(t = switch + 0:5, microscopy.use = seq(0, 0.99, length.out = 7)[-1])
}

# bring together with threshold_switch info
switch_list <- list("iso_switch" = switch_info,
                   "threshold_switch" = data.frame(t_past = c(0, 1, 2), microscopy.use = c(0.25, 0.5, 0.99), f_trig = 0.05))


# --------------------------------------------------------------------------#
# 3. Run our simulation with RDT switch strategies --------------------------
# --------------------------------------------------------------------------#

# Initialise our model
spread_model <- R6_hrp2_spread$new(afr_map, mod_obj, adj_mat = adj_mat)

# Set up our initial conditions
spread_model$set_seeds(setNames(seeds_adm1$prev, seeds_adm1$id_1))

# Central

# Set up our selection speed map data
spread_model$set_map_data(scenario_maps$map_data[[365]] %>% filter(id_1 %in% afr_map$id_1))

# Simulate our data
out_central <- spread_model$simulate_spread_and_rdts(export_freq = 0.25, t_break = 0.1, switch_list = switch_list) %>%
  mutate(scenario = "central")

# Best

# Set up our selection speed map data
best_row <- which(apply(scenario_maps$scenarios, 1, function(row) all(row == "best")))
spread_model$set_map_data(scenario_maps$map_data[[best_row]] %>% filter(id_1 %in% afr_map$id_1))

# Simulate our data
out_best <- spread_model$simulate_spread_and_rdts(export_freq = 0.25, t_break = 0.1, switch_list = switch_list) %>%
  mutate(scenario = "best")

# Worst

# Set up our selection speed map data
worst_row <- which(apply(scenario_maps$scenarios, 1, function(row) all(row == "worst")))
spread_model$set_map_data(scenario_maps$map_data[[worst_row]] %>% filter(id_1 %in% afr_map$id_1))

# Simulate our data
out_worst <- spread_model$simulate_spread_and_rdts(export_freq = 0.25, t_break = 0.1, switch_list = switch_list) %>%
  mutate(scenario = "worst")

# --------------------------------------------------------------------------#
# 4. Run our simulation without RDT switch strategies -----------------------
# --------------------------------------------------------------------------#

# Here to make it comparable, let's assume that DJI, ERI and ETH all are at 100% non-hrp2 too

# Central

# Set up our selection speed map data
spread_model$set_map_data(scenario_maps$map_data[[365]] %>% filter(id_1 %in% afr_map$id_1) %>%
                            mutate(microscopy.use = replace(microscopy.use, iso3c %in% c("DJI", "ERI", "ETH"), 0.99)) %>%
                            mutate(s = mod_obj$predict_s(.)) %>%
                            mutate(s = if_else(iso3c %in% c("DJI", "ERI", "ETH"), pmin(s, 0), s)))

# Simulate our data
out_central_cf <- spread_model$simulate_spread(export_freq = 0.25, t_break = 0.1) %>%
  mutate(scenario = "central")

# Best

# Set up our selection speed map data
best_row <- which(apply(scenario_maps$scenarios, 1, function(row) all(row == "best")))
spread_model$set_map_data(scenario_maps$map_data[[best_row]] %>% filter(id_1 %in% afr_map$id_1) %>%
                            mutate(microscopy.use = replace(microscopy.use, iso3c %in% c("DJI", "ERI", "ETH"), 0.99)) %>%
                            mutate(s = mod_obj$predict_s(.)) %>%
                            mutate(s = if_else(iso3c %in% c("DJI", "ERI", "ETH"), pmin(s, 0), s)))

# Simulate our data
out_best_cf <- spread_model$simulate_spread(export_freq = 0.25, t_break = 0.1) %>%
  mutate(scenario = "best")

# Worst

# Set up our selection speed map data
worst_row <- which(apply(scenario_maps$scenarios, 1, function(row) all(row == "worst")))
spread_model$set_map_data(scenario_maps$map_data[[worst_row]] %>% filter(id_1 %in% afr_map$id_1) %>%
                            mutate(microscopy.use = replace(microscopy.use, iso3c %in% c("DJI", "ERI", "ETH"), 0.99)) %>%
                            mutate(s = mod_obj$predict_s(.)) %>%
                            mutate(s = if_else(iso3c %in% c("DJI", "ERI", "ETH"), pmin(s, 0), s)))

# Simulate our data
out_worst_cf <- spread_model$simulate_spread(export_freq = 0.25, t_break = 0.1) %>%
  mutate(scenario = "worst")


# -------------------------------------------------------------------------#
# 5. Combine and create plots -----------------------
# -------------------------------------------------------------------------#

# create our full data
full_df <- rbind(
  rbind(out_central %>% mutate(type = "MedAccess RDT Strategy"),
        out_worst %>% mutate(type = "MedAccess RDT Strategy")) %>%
    rbind(out_best %>% mutate(type = "MedAccess RDT Strategy")),
  rbind(out_central_cf %>% mutate(type = "Counterfactual Strategy"),
        out_worst_cf %>% mutate(type = "Counterfactual Strategy")) %>%
    rbind(out_best_cf %>% mutate(type = "Counterfactual Strategy"))
)

# simple mean in a country

full_df %>%
  mutate(year = lubridate::year(lubridate::dyears(t) + as.Date("2023-04-01"))) %>%
  left_join(afr_map %>% sf::st_drop_geometry()) %>%
  group_by(year, iso, type, scenario) %>%
    summarise(freq = mean(freq)) %>%
  filter(year < 2034) %>%
  ggplot(aes(year, freq, color = scenario, linetype = type)) +
  geom_line() +
  facet_wrap(~iso, scales = "free_y")

saveRDS(full_df, "analysis/impact_analysis/data-derived/medaccess_sims.rds")

# ----------------------- #
# 6. Make figures and tables for Medaccess ---------
# ----------------------- #

# overall
inc_df <- readRDS("/home/oj/GoogleDrive/AcademicWork/Imperial/git/arms/analysis/data-derived/spread_mal_inc.rds")
plotting_df <- full_df %>%
  mutate(year = lubridate::year(lubridate::dyears(t) + as.Date("2023-04-01"))) %>%
  left_join(afr_map %>% sf::st_drop_geometry()) %>%
  left_join(
    inc_df %>%
       select(iso, id_1, year, scenario,starts_with("inc")) %>%
       mutate(scenario = tolower(scenario),
              year = as.integer(year)) %>%
  filter(scenario == "central") %>%
    filter(year <= 2035 & year >= 2025) %>%
  select(-scenario, -inc_low, -inc_high)
  ) %>%
  group_by(year, iso, type, scenario) %>%
  na.omit %>%
  summarise(freq = weighted.mean(freq, inc_med, na.rm = TRUE,)) %>%
  mutate(scenario = stringr::str_to_title(scenario)) %>%
  mutate(iso = paste0(iso, ": ", countrycode::countrycode(iso, "iso3c", "country.name.en")))


# all country plot
fig1 <- plotting_df %>%
  ggplot(aes(year, freq, color = scenario, linetype = type)) +
  geom_vline(xintercept = 2029, linetype = "dashed") +
  geom_line(lwd = 1) +
  lemon::facet_rep_wrap(~iso, scales = "free_y", ncol = 5) +
  theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(axis.line = element_line()) +
  scale_color_manual(values = rev(c(MetBrewer::met.brewer("Egypt", 1), "black", MetBrewer::met.brewer("Egypt", 2)[2])),
                     name = "HRP2 Selection Scenario:",
                     labels = c("Best", "Central", "Worst")) +
  scale_linetype(name = "RDT Scenario:") +
  theme(legend.position = "top", plot.background = element_rect(fill = "white")) +
  ylab("Annual False Negative Rates (%)\n") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_x_continuous(breaks = seq(2025, 2035, 2), name = "\nYear")
fig1
save_figs("all_country_comparisons", fig1, width = 14, height = 14, plot_dir = "analysis/impact_analysis/plots")

# priority country plot
fig2 <- fig1
priority_countries <- c("DJI", "ERI", "ETH", "GHA", "SSD", "SDN", "DRC", "KEN", "UGA", "ZAM", "MLI", "NGA", "SEN", "SLE", "MWI")
fig2$data <- filter(fig2$data, iso %in%
                      paste0(priority_countries, ": ", countrycode::countrycode(priority_countries, "iso3c", "country.name.en")))
fig2 <- fig2 +
  lemon::facet_rep_wrap(~iso, scales = "free_y", ncol = 3) + theme(legend.position = "right")
save_figs("priority_comparisons", fig2, width = 14, height = 11, plot_dir = "analysis/impact_analysis/plots")

# all together plot
fig3 <- full_df %>%
  mutate(year = lubridate::year(lubridate::dyears(t) + as.Date("2023-04-01"))) %>%
  left_join(afr_map %>% sf::st_drop_geometry()) %>%
  left_join(
    inc_df %>%
      select(iso, id_1, year, scenario,starts_with("inc")) %>%
      mutate(scenario = tolower(scenario),
             year = as.integer(year)) %>%
      filter(scenario == "central") %>%
      filter(year <= 2035 & year >= 2025) %>%
      select(-scenario, -inc_low, -inc_high)
  ) %>%
  group_by(year, type, scenario) %>%
  na.omit %>%
  summarise(freq = weighted.mean(freq, inc_med, na.rm = TRUE)) %>%
  mutate(scenario = stringr::str_to_title(scenario)) %>%
  ggplot(aes(year, freq, color = scenario, linetype = type)) +
  geom_vline(xintercept = 2029, linetype = "dashed") +
  geom_line(lwd = 1) +
  theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(axis.line = element_line()) +
  scale_color_manual(values = rev(c(MetBrewer::met.brewer("Egypt", 1), "black", MetBrewer::met.brewer("Egypt", 2)[2])),
                     name = "HRP2 Selection Scenario:",
                     labels = c("Best", "Central", "Worst")) +
  scale_linetype(name = "RDT Scenario:") +
  theme(legend.position = "right", plot.background = element_rect(fill = "white")) +
  ylab("Annual False Negative Rates\n") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_x_continuous(breaks = seq(2023, 2035, 2), name = "\nYear")
save_figs("africa_comparisons", fig3, width = 8, height = 6, plot_dir = "analysis/impact_analysis/plots")


# all together difference
fig4 <- full_df %>%
  mutate(year = lubridate::year(lubridate::dyears(t) + as.Date("2023-04-01"))) %>%
  left_join(afr_map %>% sf::st_drop_geometry()) %>%
  left_join(
    inc_df %>%
      select(iso, id_1, year, scenario,starts_with("inc")) %>%
      mutate(scenario = tolower(scenario),
             year = as.integer(year)) %>%
      filter(scenario == "central") %>%
      filter(year <= 2035 & year >= 2025) %>%
      select(-scenario, -inc_low, -inc_high)
  ) %>%
  group_by(year, type, scenario) %>%
  na.omit %>%
  summarise(freq = weighted.mean(freq, inc_med, na.rm = TRUE)) %>%
  mutate(scenario = stringr::str_to_title(scenario)) %>%
  pivot_wider(values_from = freq, names_from = type) %>%
  janitor::clean_names() %>%
  mutate(impact = counterfactual_strategy - med_access_rdt_strategy) %>%
  ungroup() %>%
  select(year, scenario, impact) %>%
    pivot_wider(values_from = impact, names_from = scenario) %>%
  janitor::clean_names() %>%
  pivot_longer(-year) %>%
  ggplot(aes(x = year, y = value, color = name)) +
  geom_vline(xintercept = 2029, linetype = "dashed") +
  geom_smooth(lwd = 1, se = FALSE) +
  theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(axis.line = element_line()) +
  scale_color_manual(values = rev(c(MetBrewer::met.brewer("Egypt", 1), "black", MetBrewer::met.brewer("Egypt", 2)[2])),
                     name = "HRP2 Selection Scenario:",
                     labels = c("Best", "Central", "Worst")) +
  theme(legend.position = "right", plot.background = element_rect(fill = "white")) +
  ylab("Annual Averted False Neagative Rates (%) \n(Averted = Counterfactual - MedAcces Strategy)\n") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_x_continuous(breaks = seq(2023, 2035, 2), name = "\nYear")
save_figs("africa_impact", fig4, width = 8, height = 6, plot_dir = "analysis/impact_analysis/plots")


# priority country 2029 and 2035impact
fig5 <- full_df %>%
  mutate(year = lubridate::year(lubridate::dyears(t) + as.Date("2023-04-01"))) %>%
  left_join(afr_map %>% sf::st_drop_geometry()) %>%
  left_join(
    inc_df %>%
      select(iso, id_1, year, scenario,starts_with("inc")) %>%
      mutate(scenario = tolower(scenario),
             year = as.integer(year)) %>%
      filter(scenario == "central") %>%
      filter(year <= 2035 & year >= 2025) %>%
      select(-scenario, -inc_low, -inc_high)
  ) %>%
  group_by(year, iso, type, scenario) %>%
  na.omit %>%
  summarise(freq = weighted.mean(freq, inc_med, na.rm = TRUE)) %>%
  filter(iso %in% priority_countries) %>%
  mutate(scenario = stringr::str_to_title(scenario)) %>%
  mutate(iso = paste0(iso, ": ", countrycode::countrycode(iso, "iso3c", "country.name.en"))) %>%
  janitor::clean_names() %>%
  filter(year %in% c(2029, 2035)) %>%
  ggplot(aes(x = as.factor(year), y = freq, alpha = type, group = interaction(year, type, scenario), fill = scenario)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge(width = 0.9, preserve = "single")) +
  lemon::facet_rep_wrap(~iso, scales = "free_y", ncol = 4)  +
  theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(axis.line = element_line(), plot.background = element_rect(fill = "white")) +
  scale_alpha_manual(values = c(1, 0.4), name = "RDT Scenario:") +
  scale_fill_manual(values = rev(c(MetBrewer::met.brewer("Egypt", 1), "black", MetBrewer::met.brewer("Egypt", 2)[2])),
                     name = "HRP2 Selection Scenario:",
                     labels = c("Best", "Central", "Worst")) +
  scale_y_continuous(labels = scales::percent_format(), name = "Annual False Negative Rates (%)\n") +
  xlab("Year")
fig5

save_figs("prioty_county_29_35", fig5, width = 16, height = 14, plot_dir = "analysis/impact_analysis/plots")

# -------------------------------------------------------------------------#
# 7. Example scenario plots for RDT changes -----------------------
# -------------------------------------------------------------------------#


rdt_roll_out <- data.frame("S.Times" = lubridate::year(lubridate::dyears(c(0, 1, 2, 3, 31)+1.3) + as.Date("2023-04-01")),
                           "RDT" = c(0,0.25,0.5,1,1))

complete <- full_df %>%
  mutate(year = lubridate::year(lubridate::dyears(t) + as.Date("2023-04-01"))) %>%
  left_join(afr_map %>% sf::st_drop_geometry()) %>%
  filter(iso == "KEN") %>%
  filter(id_1 == 10314353) %>%
  filter(scenario == "central") %>%
  ggplot(aes(year, freq, group = interaction(id_1,type), color = type)) +
  geom_smooth(span = 0.1, se = FALSE) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_step(aes(S.Times, RDT), data = rdt_roll_out, color = "black", inherit.aes = FALSE) +
  xlab("Year") +
  scale_x_continuous(limits = c(2023,2055), breaks = c(2023, 2030, 2040, 2050)) +
  theme_bw() +
  MetBrewer::scale_color_met_d("Egypt",name = "RDT Strategy:") +
  theme(legend.position = "top") +
  scale_y_continuous(labels = scales::percent_format(),limits = c(0,1),
                     sec.axis = sec_axis(~ . * 1, labels = scales::percent_format(),
                                         name = "Proportion of novel RDTs introduced to market")  ) +
  ylab("Percentage of False-Negative RDT Results due to pfhrp2/3 deletions")
complete

initial <- full_df %>%
  mutate(year = lubridate::year(lubridate::dyears(t) + as.Date("2023-04-01"))) %>%
  left_join(afr_map %>% sf::st_drop_geometry()) %>%
  filter(iso == "KEN") %>%
  filter(id_1 == 10314353) %>%
  filter(scenario == "central" & type == "Counterfactual Strategy") %>%
  ggplot(aes(year, freq, group = interaction(id_1,type), color = type)) +
  geom_smooth(span = 0.1, se = FALSE) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  #geom_step(aes(S.Times, RDT), data = rdt_roll_out, color = "black", inherit.aes = FALSE) +
  xlab("Year") +
  scale_x_continuous(limits = c(2023,2055), breaks = c(2023, 2030, 2040, 2050)) +
  theme_bw() +
  MetBrewer::scale_color_met_d("Egypt",name = "RDT Strategy:") +
  theme(legend.position = "top") +
  scale_y_continuous(labels = scales::percent_format(),limits = c(0,1),
                     sec.axis = sec_axis(
                       ~ . * 1, labels = scales::percent_format(),
                                         guide = guide_axis(theme = theme(axis.text = element_text(colour = "white"),
                                                                          axis.ticks = element_line(colour = "white"))),
                                         name = "")  ) +
  ylab("Percentage of False-Negative RDT Results due to pfhrp2/3 deletions")
initial

mid <- full_df %>%
  mutate(year = lubridate::year(lubridate::dyears(t) + as.Date("2023-04-01"))) %>%
  left_join(afr_map %>% sf::st_drop_geometry()) %>%
  filter(iso == "KEN") %>%
  filter(id_1 == 10314353) %>%
  filter(scenario == "central" & type == "Counterfactual Strategy") %>%
  ggplot(aes(year, freq, group = interaction(id_1,type), color = type)) +
  geom_smooth(span = 0.1, se = FALSE) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_step(aes(S.Times, RDT), data = rdt_roll_out, color = "black", inherit.aes = FALSE) +
  xlab("Year") +
  scale_x_continuous(limits = c(2023,2055), breaks = c(2023, 2030, 2040, 2050)) +
  theme_bw() +
  MetBrewer::scale_color_met_d("Egypt",name = "RDT Strategy:") +
  theme(legend.position = "top") +
  scale_y_continuous(labels = scales::percent_format(),limits = c(0,1),
                     sec.axis = sec_axis(~ . * 1, labels = scales::percent_format(),
                                         name = "Proportion of novel RDTs introduced to market")  ) +
  ylab("Percentage of False-Negative RDT Results due to pfhrp2/3 deletions")
mid

initial
mid
complete
