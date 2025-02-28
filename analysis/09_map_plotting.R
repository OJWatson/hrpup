# ---------------------------------------------------- #
# 0. Get the world map and our data  ----
# ---------------------------------------------------- #

library(tidyverse)
devtools::load_all()

# WHO world map boundaries
who_shps <- readRDS("analysis/data_derived/who_shps.rds")

# scenario data
scenario_maps <- readRDS("analysis/data_derived/scenario_maps_full.rds")

# ---------------------------------------------------- #
# 1. Plotting central time maps  ----
# ---------------------------------------------------- #

# central
scenario <- which(apply(scenario_maps$scenarios, 1, function(x){all(x == "central")}))

# cut off function
cut_off <- list(
  breaks = c(0,10,20,30,40),
  limits = c(0,40),
  color = "blue",
  cutoff = function(x){x %>% filter(t>40 | t < 0)},
  name = "Years",
  label = "40+ Years"
)

# africa t
gg_map <- who_compliant_continuous_plot(scenario_maps$map_data[[scenario]], who_shps, region = "africa",
                                        var = "t", title = "Years", limits = TRUE, cut_off = cut_off, prev_plot = TRUE, prev_var = "Micro.2.10")
save_figs("central_afr_map", add_africa_scale(gg_map), width = 10, height = 8, pdf_plot = FALSE)

# world t
gg_map2 <- who_compliant_continuous_plot(scenario_maps$map_data[[scenario]], who_shps, region = "global",
                                        var = "t", title = "Years", limits = TRUE, cut_off = cut_off, prev_plot = TRUE, prev_var = "Micro.2.10") +
  theme(legend.text = element_text(family = "Helvetica", size = 14), legend.title = element_text(family = "Helvetica", size = 16))
save_figs("central_world_map", add_global_scale(gg_map2), width = 20, height = 6, pdf_plot = FALSE)

# ---------------------------------------------------- #
# 2. Plotting central risk maps ----
# ---------------------------------------------------- #

# africa_risk first
scenario <- 365
create_risk <- function(x){
  x %>% mutate(t = replace(t, t < 0, Inf)) %>% mutate(t_bin = cut(t, breaks = c(0,6,12,20, Inf)))
}

# central
central_gg <- who_compliant_discrete_plot(create_risk(scenario_maps$map_data[[scenario]]),who_shps, region = "africa",
                                      var = "t_bin", title = "HRP2 Concern",limits = TRUE, prev_plot = TRUE, prev_var = "Micro.2.10")
save_figs("central_afr_risk", add_africa_scale(central_gg), width = 10, height = 8, pdf_plot = FALSE)

# worst
worst_gg <- who_compliant_discrete_plot(create_risk(scenario_maps$map_data[[scenario-1]]),who_shps, region = "africa",
                                        var = "t_bin", title = "HRP2 Concern",limits = TRUE, prev_plot = TRUE, prev_var = "Micro.2.10")
save_figs("worst_afr_risk", add_africa_scale(worst_gg), width = 10, height = 8, pdf_plot = FALSE)

# best
best_gg <- who_compliant_discrete_plot(create_risk(scenario_maps$map_data[[scenario+1]]),who_shps, region = "africa",
                                        var = "t_bin", title = "HRP2 Concern",limits = TRUE, prev_plot = TRUE, prev_var = "Micro.2.10")
save_figs("best_afr_risk", add_africa_scale(best_gg), width = 10, height = 8, pdf_plot = FALSE)

# world_risk
gg_map4 <- who_compliant_discrete_plot(create_risk(scenario_maps$map_data[[scenario]]),who_shps, region = "global",
                                          var = "t_bin", title = "HRP2 Concern",limits = TRUE, prev_plot = TRUE, prev_var = "Micro.2.10")
save_figs("central_world_risk1", add_global_scale(gg_map4), width = 20, height = 6, pdf_plot = FALSE)

# combine so we don't have to make by hand:
gg_map5 <- cowplot::plot_grid(
  add_africa_scale(worst_gg) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 22, family = "Helvetica")) + ggtitle("Worst Case"),
  add_africa_scale(central_gg) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 22, family = "Helvetica")) + ggtitle("Central Case"),
  add_africa_scale(best_gg) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 22, family = "Helvetica")) + ggtitle("Best Case"),
  cowplot::get_legend(central_gg + theme(legend.text = element_text(family = "Helvetica", size = 14), legend.title = element_text(family = "Helvetica", size = 16))),
  ncol = 4,
  rel_widths = c(1, 1, 1, 0.7)
) + theme(plot.background = element_rect(fill = "white", color = "white"))
save_figs("worst_central_best_africa_risk", gg_map5, width = 20, height = 6, pdf_plot = FALSE)


# ---------------------------------------------------- #
# 3. Creating Data Outputs  ----
# ---------------------------------------------------- #

compl <- list()
for(i in seq_along(scenario_maps$map_data)) {
compl[[i]] <- left_join(scenario_maps$map, scenario_maps$map_data[[i]]) %>% sf::st_drop_geometry() %>%
  select(iso, name_0, name_1, t, tmin, tmax, s, smin, smax, Micro.2.10, ft, microscopy.use, rdt.nonadherence, fitness, rdt.det) %>%
    mutate(Micro.2.10_scen = scenario_maps$scenarios$Micro.2.10[i],
           ft_scen = scenario_maps$scenarios$ft[i],
           microscopy.use_scen = scenario_maps$scenarios$microscopy.use[i],
           rdt.nonadherence_scen = scenario_maps$scenarios$rdt.nonadherence[i],
           fitness_scen = scenario_maps$scenarios$fitness[i],
           rdt.det_scen = scenario_maps$scenarios$rdt.det[i])
}

compl_df <- do.call(rbind, compl)
saveRDS(compl_df, "analysis/data_derived/complete_times.rds")

# Times censored at 40 years
central_df <- compl_df %>% group_by(iso, name_0, name_1) %>%
  filter(if_all(ends_with("scen"), ~ . == "central")) %>%
  filter(!(iso %in% c("ATF", "UMI"))) %>%
  mutate(continent = countrycode::countrycode(iso, "iso3c", "continent"), .before = 1) %>%
  na.omit %>%
  select(continent:tmax) %>%
  mutate(t_med = round(replace(t, t<0, Inf),2),
         t_uci = round(replace(tmin, tmin<0, Inf),2),
         t_lci = round(replace(tmax, tmax<0, Inf),2),
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
  na.omit %>%
  select(continent:tmax) %>%
  mutate(t_med = round(replace(t, t<0, Inf),2),
         t_uci = round(replace(tmin, tmin<0, Inf),2),
         t_lci = round(replace(tmax, tmax<0, Inf),2),
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
  na.omit %>%
  select(continent:tmax) %>%
  mutate(t_med = round(replace(t, t<0, Inf),2),
         t_uci = round(replace(tmin, tmin<0, Inf),2),
         t_lci = round(replace(tmax, tmax<0, Inf),2),
         scenario = "Optimistic") %>%
  mutate(t_med = replace(t_med, t_med>40, "40+"),
         t_uci = replace(t_uci, t_uci>40, "40+"),
         t_lci = replace(t_lci, t_lci>40, "40+")) %>%
  rename(admin_0 = name_0) %>%
  rename(admin_1 = name_1) %>%
  arrange(continent, iso,  admin_1) %>%
  select(continent, iso, admin_0, admin_1, t_lci, t_med, t_uci, scenario)

write.csv(central_df, "analysis/data_out/central_times.csv", row.names = FALSE)
write.csv(worst_df, "analysis/data_out/pessimistic_times.csv", row.names = FALSE)
write.csv(best_df, "analysis/data_out/optimistic_times.csv", row.names = FALSE)

# ---------------------------------------------------- #
# 4. Creating Key High Risk Countries
# ---------------------------------------------------- #
isos <- unique(countrycode::codelist$iso3c)#[countrycode::codelist$continent == "Africa"])
world <- left_join(scenario_maps$map, scenario_maps$map_data[[365]]) %>%
  sf::st_drop_geometry()

micro005 <- world %>% sf::st_drop_geometry() %>%
  group_by(iso3c) %>%
  summarise(nm = sum(Micro.2.10 > 0.005,na.rm=TRUE)/n()) %>%
  filter(nm > 0.5) %>%
  pull(iso3c)

world <- world %>%
  filter(!is.na(t)) %>%
  filter(iso %in% isos) %>%
  filter(iso %in% micro005) %>%
  mutate(t = replace(t, t < 0, Inf)) %>%
  mutate(t_bin = cut(t, breaks = c(0,6,12,20, Inf), right = FALSE, include.lowest = TRUE))

hrtable <- world %>% sf::st_drop_geometry() %>%
  group_by(name_0) %>%
  summarise(r = sum(t_bin == "[0,6)")/n()) %>%
  filter(r >= 0.5) %>%
  arrange(desc(r)) %>%
  mutate(r = scales::percent_format()(r)) %>%
  setNames(c("Country", "% Admin 1 High Risk"))

write.csv(hrtable, file = "analysis/tables/high_innate_risk_countries.csv", row.names = FALSE)
