# ---------------------------------------------------- #
# 0. Get the world map and our data  ----
# ---------------------------------------------------- #

library(tidyverse)

# MAP world map boundaries
available_admin <- malariaAtlas::listShp(printed = FALSE, admin_level = "admin0")
world_map_0 <- malariaAtlas::getShp(ISO = available_admin$iso, admin_level = c("admin0")) %>% sf::st_as_sf()

# scenario data
scenario_maps <- readRDS("analysis/data_derived/scenario_maps_full_monotonic.rds")

# ---------------------------------------------------- #
# 1. Plotting central time maps  ----
# ---------------------------------------------------- #

# plotting t maps -------------------------------#

# get scenario data
world_f <- function(scenario, bounds = FALSE){

world <- left_join(scenario_maps$map, scenario_maps$map_data[[scenario]]) %>%
  filter(iso %in% isos)
if(bounds){
  world <- world %>%
    filter(!is.na(t)) %>%
  mutate(t = replace(t, t>=40, NA)) %>%
  mutate(t = replace(t, t<0, NA))
}
return(world)
}

# create map
create_map_t <- function(world, world_map_0, limits = FALSE) {

gg_map <- world %>%
  #filter(Micro.2.10 > 0.0001) %>%
  ggplot() +
  #geom_sf(aes(fill = cut(ft, seq(0,1,0.2))), color = NA, show.legend = TRUE) +
  geom_sf(aes(fill = t), color = "grey", show.legend = TRUE, lwd = 0.1) +
  scale_fill_viridis_c(name = "Years", breaks = c(0,10,20,30,40), limits = c(0,40), option = "C", direction = -1, end = 0.7)


if(any((gg_map$data$Micro.2.10 < 0.0005))) {

  gg_map <- gg_map +
    ggnewscale::new_scale_fill() +
    ggpattern::geom_sf_pattern(
      pattern_fill = "grey", pattern = "stripe", fill = NA, show.legend = FALSE,
      color = NA,
      pattern_colour = NA,
      pattern_density = 0.5,
      pattern_spacing = 0.025,
      data = . %>% filter(Micro.2.10 < 0.0005), inherit.aes = FALSE) +
    scale_fill_manual(name="\nTransmission", labels="Unstable (<0.05% PfPR)", values="grey")

}

lrb <- "#4169e1"
if(any(is.na(gg_map$data$t))) {

  gg_map <- gg_map +
    ggnewscale::new_scale_fill() +
    geom_sf(aes(fill = "blue"), color = NA, show.legend = TRUE, data = . %>% filter(is.na(t))) +
    scale_fill_manual(name=" ", labels="40+ Years", values="blue")

}


xlim <- layer_scales(gg_map)$x$get_limits()
ylim <- layer_scales(gg_map)$y$get_limits()

gg_map <- gg_map +
  geom_sf(fill = NA, color = "black", show.legend = FALSE,
          data = world_map_0 %>% filter(iso %in% isos), lwd = 0.2) +
  coord_sf() +
  theme_void() +
  theme(plot.caption = element_text(face = "italic"), plot.background = element_rect(fill = "white", color = "white"))

if(limits){
  gg_map <- gg_map + xlim(xlim) + ylim(ylim)
}
gg_map
}

# central
scenario <- 365

# africa t
isos <- unique(countrycode::codelist$iso3c[countrycode::codelist$continent == "Africa"])
world <- world_f(scenario, TRUE)
gg_map <- create_map_t(world, world_map_0)
save_figs("central_afr_map", gg_map, width = 10, height = 8)

# world t
isos <- unique(countrycode::codelist$iso3c)
world <- world_f(scenario, TRUE)
gg_map2 <- create_map_t(world, world_map_0, TRUE)
save_figs("central_world_map", gg_map2, width = 30, height = 8)

# ---------------------------------------------------- #
# 2. Plotting central risk maps ----
# ---------------------------------------------------- #

# ------- Plotting Risk Maps----------------------- #

create_map_risk <- function(scenario_maps, scenario, isos, limits = FALSE)  {

# bin risk
pal <- colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=1)

world <- left_join(scenario_maps$map, scenario_maps$map_data[[scenario]]) %>%
  filter(!is.na(t)) %>%
  filter(iso %in% isos) %>%
  mutate(t = replace(t, t < 0, Inf)) %>%
  mutate(t_bin = cut(t, breaks = c(0,6,12,20, Inf)))

gg_map_risk <- world %>%
  ggplot() +
  geom_sf(aes(fill = t_bin), color = "grey", show.legend = TRUE, lwd = 0.1) +
  scale_fill_manual(name = "HRP2 Concern", values = rev(c("blue", "cyan", "yellow", "red")),
                    labels = rev(c("Marginal", "Slight", "Moderate", "High")))

if(any((gg_map_risk$data$Micro.2.10 < 0.0005))) {

  gg_map_risk <- gg_map_risk +
    ggnewscale::new_scale_fill() +
    ggpattern::geom_sf_pattern(
      pattern_fill = "grey", pattern = "stripe", fill = NA, show.legend = FALSE,
      color = NA,
      pattern_colour = NA,
      pattern_density = 0.5,
      pattern_spacing = 0.025,
      data = . %>% filter(Micro.2.10 < 0.0005), inherit.aes = FALSE) +
    scale_fill_manual(name="\nTransmission", labels="Unstable (<0.05% PfPR)", values="grey")

}

xlim <- layer_scales(gg_map_risk)$x$get_limits()
ylim <- layer_scales(gg_map_risk)$y$get_limits()

gg_map_risk <- gg_map_risk +
  geom_sf(fill = NA, color = "black", show.legend = FALSE,
          data = world_map_0 %>% filter(iso %in% isos), lwd = 0.2) +
  coord_sf() +
  theme_void() +
  theme(plot.caption = element_text(face = "italic"), plot.background = element_rect(fill = "white", color = "white"))

if(limits){
  gg_map_risk <- gg_map_risk + xlim(xlim) + ylim(ylim)
}
gg_map_risk

}

# africa_risk
scenario <- 365
isos <- unique(countrycode::codelist$iso3c[countrycode::codelist$continent == "Africa"])
gg_map3 <- create_map_risk(scenario_maps, scenario, isos, FALSE)
save_figs("central_afr_risk", gg_map3, width = 10, height = 8)

scenario <- 364
gg_map3 <- create_map_risk(scenario_maps, scenario, isos, FALSE)
save_figs("worst_afr_risk", gg_map3, width = 10, height = 8)

scenario <- 366
gg_map3 <- create_map_risk(scenario_maps, scenario, isos, FALSE)
save_figs("best_afr_risk", gg_map3, width = 10, height = 8)

# world_risk
isos <- unique(countrycode::codelist$iso3c)
gg_map4 <- create_map_risk(scenario_maps, scenario, isos, TRUE)
save_figs("central_world_risk", gg_map4, width = 30, height = 8)

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

write.csv(central_df, "analysis/tables/central_times.csv", row.names = FALSE)
write.csv(worst_df, "analysis/tables/pessimistic_times.csv", row.names = FALSE)
write.csv(best_df, "analysis/tables/optimistic_times.csv", row.names = FALSE)

# ---------------------------------------------------- #
# 4. Creating Key High Risk Countries
# ---------------------------------------------------- #

world <- left_join(scenario_maps$map, scenario_maps$map_data[[365]]) %>%
  filter(!is.na(t)) %>%
  filter(iso %in% isos) %>%
  mutate(t = replace(t, t < 0, Inf)) %>%
  mutate(t_bin = cut(t, breaks = c(0,6,12,20, Inf)))

hrtable <- world %>% sf::st_drop_geometry() %>% group_by(name_0) %>%
  summarise(r = sum(t_bin == "(0,6]")/n()) %>%
  filter(r >= 0.5) %>%
  arrange(desc(r)) %>%
  mutate(r = scales::percent_format()(r)) %>%
  setNames(c("Country", "% Admin 1 High Risk"))

write.csv(hrtable, file = "analysis/tables/high_innate_risk_countries.csv")

# -----------------------------------#
# ---------------------------------------------------- #
# 4. Creating Time Ranges  ----
# ---------------------------------------------------- #

central_time_spans <- central_df %>%
  filter(continent == "Africa") %>%
  na.omit() %>%
  filter(is.finite(t_lci) & is.finite(t_med) & is.finite(t_uci)) %>%
  mutate(uid = interaction(iso, admin_1)) %>%
  ggplot(aes(y = fct_reorder(uid,t_med), x = t_med, xmax = t_uci, xmin = t_lci)) +
  geom_errorbarh(lwd = 0.1) +
  geom_point(size = 0.1) +
  ggpubr::theme_pubclean() +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(),
        panel.grid.major.x = element_line(colour = "grey")) +
  scale_x_log10() +
  xlab("Years under central scenario for proportion of false-negative
  clinical cases to increase from 1% to 5%") +
  ylab("Arican Admin 1 Regions")
save_figs("central_time_spans", central_time_spans, 5, 5)

lvs <- levels(central_df %>% mutate(uid = interaction(iso, admin_1)) %>% mutate(uid = fct_reorder(uid,t_med)) %>% pull(uid))

all_time_spans <- central_df %>%
  rbind(worst_df) %>%
  rbind(best_df) %>%
  filter(continent == "Africa") %>%
  na.omit() %>%
  filter(is.finite(t_lci) & is.finite(t_med) & is.finite(t_uci)) %>%
  mutate(uid = interaction(iso, admin_1)) %>%
  ggplot(aes(y = fct_reorder(uid,t_med), x = t_med, xmax = t_uci, xmin = t_lci, color = scenario)) +
  geom_errorbarh(lwd = 0.1) +
  geom_point(size = 0.1) +
  ggpubr::theme_pubclean() +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(),
        panel.grid.major.x = element_line(colour = "grey")) +
  scale_x_log10() +
  xlab("Years under central scenario for proportion of false-negative
  clinical cases to increase from 1% to 5%") +
  ylab("Arican Admin 1 Regions")
save_figs("all_time_spans", all_time_spans, 5, 5)
