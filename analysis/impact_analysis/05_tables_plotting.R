## Spread Modelling
library(tidyverse)
library(sf)
library(wpp2017)
data("UNlocations")
# read in data
full_df_raw <- readRDS("analysis/impact_analysis/data-derived/chai_sims.rds")
# file.copy("/home/oj/GoogleDrive/AcademicWork/Imperial/git/arms/analysis/data-derived/spread_mal_inc.rds", "analysis/impact_analysis/data-raw/spread_mal_inc.rds")
inc_df <- readRDS("analysis/impact_analysis/data-raw/spread_mal_inc.rds")

# read in maps for iso names etc
admin1 <- readRDS(here::here("analysis/data_derived/admin1_sf.rds"))
admin1$cont <- countrycode::countrycode(admin1$iso, "iso3c", "iso3n")
admin1$region <- UNlocations$area_name[match(admin1$cont, UNlocations$country_code)]
afr_map <- admin1 %>% filter(region == "Africa")

# sort out t differences
closest_to_integer <- function(vec, integers = 0:20) {
  sapply(integers, function(i) vec[which.min(abs(vec - i))])
}

# create cleaned data
full_df <- full_df_raw %>%
  mutate(t = replace(t, t>20 & t < max(.data$t), max(.data$t))) %>%
  filter(t %in% closest_to_integer(unique(t))) %>%
  mutate(t = round(t, 1)) %>%
  filter(t %in% seq(0, 20, 1)) %>%
  group_by(scenario, id_1) %>%
  arrange(t) %>%
  mutate(across(matches("_\\w"), ~ .x - (.x[t==0] - mean(.x[t ==0])))) %>%
  mutate(across(matches("_[a-z]"), function(x){signif(x, 4)})) %>%
  arrange(delay, scenario) %>%
  ungroup() %>%
  mutate(year = lubridate::year(lubridate::dyears(t) + as.Date("2023-04-01"))) %>%
  left_join(afr_map %>% sf::st_drop_geometry() %>% select(iso, id_1)) %>%
  mutate(iso = paste0(
    iso, ": ", countrycode::countrycode(
      iso, "iso3c", "country.name.en", custom_match = c("COD" = "Democratic Republic of the Congo", "COG" == "Republic of the Congo")
    )))

# add iso
full_df <- full_df %>% mutate(
  iso = factor(iso, levels = c(
    "DJI: Djibouti", "ERI: Eritrea", "ETH: Ethiopia",
    sort(setdiff(unique(full_df$iso), c("DJI: Djibouti", "ERI: Eritrea", "ETH: Ethiopia")))
  ))
)

# priorities
priority_countries <- c("DJI", "ERI", "ETH", "GHA", "SSD", "SDN", "COD", "KEN", "UGA", "ZMB", "MLI", "NGA", "SEN", "SLE", "MWI")
priority_iso <- paste0(
  priority_countries, ": ", countrycode::countrycode(
    priority_countries, "iso3c", "country.name.en", custom_match = c("COD" = "Democratic Republic of the Congo", "COG" == "Republic of the Congo")
  ))

# population and incidence growth
inc_df <- inc_df %>%
  mutate(iso = paste0(
    iso, ": ", countrycode::countrycode(
      iso, "iso3c", "country.name.en", custom_match = c("COD" = "Democratic Republic of the Congo", "COG" == "Republic of the Congo")
    ))) %>%
  filter(iso %in% full_df$iso) %>%
  mutate(iso = factor(iso, levels = levels(full_df$iso)))

# ----------------------- #
# 1. Make hrp2 del figures and tables for CHAI Report ---------
# ----------------------- #

# overall
plotting_df <- full_df %>%
  filter(delay <= 0) %>%
  left_join(
    inc_df %>%
      select(iso, id_1, year, scenario,starts_with("inc")) %>%
      mutate(scenario = tolower(scenario),
             year = as.integer(year)) %>%
      filter(scenario == "central") %>%
      filter(year >= 2023 & year <= 2060) %>%
      select(-scenario, -inc_low, -inc_high)
  ) %>%
  group_by(year, iso, type, scenario) %>%
  na.omit %>%
  summarise(freq = weighted.mean(freq_med, inc_med, na.rm = TRUE,)) %>%
  mutate(scenario = stringr::str_to_title(scenario))

### fig1 all country plot ----

fig1 <- plotting_df %>%
  ggplot(aes(year, freq, color = scenario, linetype = type)) +
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
  scale_x_continuous(breaks = c(2023, 2030, 2040, 2050, 2060), name = "\nYear")
fig1
save_figs("all_country_comparisons", fig1, width = 14, height = 14, plot_dir = "analysis/impact_analysis/plots_chai")

### fig2 priority country plot ----

fig2 <- fig1
fig2$data <- filter(fig2$data, iso %in% priority_iso)
fig2 <- fig2 +
  lemon::facet_rep_wrap(~iso, scales = "free_y", ncol = 3) + theme(legend.position = "right")
save_figs("priority_comparisons", fig2, width = 10, height = 9, plot_dir = "analysis/impact_analysis/plots_chai")

### fig3 all together plot ----
fig3 <- full_df %>%
  filter(delay <= 0) %>%
  left_join(
    inc_df %>%
      select(iso, id_1, year, scenario,starts_with("inc")) %>%
      mutate(scenario = tolower(scenario),
             year = as.integer(year)) %>%
      filter(scenario == "central") %>%
      filter(year >= 2023) %>%
      select(-scenario, -inc_low, -inc_high)
  ) %>%
  group_by(year, type, scenario) %>%
  na.omit %>%
  summarise(freq = weighted.mean(freq_med, inc_med, na.rm = TRUE)) %>%
  mutate(scenario = stringr::str_to_title(scenario)) %>%
  ggplot(aes(year, freq, color = scenario, linetype = type)) +
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
  scale_x_continuous(breaks = c(2023, 2030, 2040, 2050, 2060), name = "\nYear")
save_figs("africa_comparisons", fig3, width = 8, height = 6, plot_dir = "analysis/impact_analysis/plots_chai")


### fig4 all together difference ----

fig4 <- full_df %>%
  filter(delay <= 0) %>%
  left_join(
    inc_df %>%
      select(iso, id_1, year, scenario,starts_with("inc")) %>%
      mutate(scenario = tolower(scenario),
             year = as.integer(year)) %>%
      filter(scenario == "central") %>%
      filter(year >= 2023) %>%
      select(-scenario, -inc_low, -inc_high)
  ) %>%
  group_by(year, type, scenario) %>%
  na.omit %>%
  summarise(freq = weighted.mean(freq_med, inc_med, na.rm = TRUE)) %>%
  mutate(scenario = stringr::str_to_title(scenario)) %>%
  pivot_wider(values_from = freq, names_from = type) %>%
  janitor::clean_names() %>%
  mutate(impact = no_rdt_switching - x5_percent_threshold_strategy) %>%
  ungroup() %>%
  select(year, scenario, impact) %>%
  pivot_wider(values_from = impact, names_from = scenario) %>%
  janitor::clean_names() %>%
  pivot_longer(-year) %>%
  ggplot(aes(x = year, y = value, color = name)) +
  geom_line(lwd = 1, se = FALSE) +
  theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(axis.line = element_line()) +
  scale_color_manual(values = rev(c(MetBrewer::met.brewer("Egypt", 1), "black", MetBrewer::met.brewer("Egypt", 2)[2])),
                     name = "HRP2 Selection Scenario:",
                     labels = c("Best", "Central", "Worst")) +
  theme(legend.position = "right", plot.background = element_rect(fill = "white")) +
  ylab("Annual Averted False Neagative Rates (%) \n(Averted = No RDT Switching - 5% Threshold)\n") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_x_continuous(breaks = c(2023, 2030, 2040, 2050, 2060), name = "\nYear")
save_figs("africa_impact", fig4, width = 8, height = 6, plot_dir = "analysis/impact_analysis/plots_chai")

### fig5 priority country 2029 and 2035impact ----

fig5 <- full_df %>%
  filter(delay <= 0) %>%
  left_join(
    inc_df %>%
      select(iso, id_1, year, scenario,starts_with("inc")) %>%
      mutate(scenario = tolower(scenario),
             year = as.integer(year)) %>%
      filter(scenario == "central") %>%
      filter(year >= 2023) %>%
      select(-scenario, -inc_low, -inc_high)
  ) %>%
  group_by(year, iso, type, scenario) %>%
  na.omit %>%
  summarise(freq = weighted.mean(freq_med, inc_med, na.rm = TRUE)) %>%
  filter(iso %in% priority_iso) %>%
  mutate(scenario = stringr::str_to_title(scenario)) %>%
  janitor::clean_names() %>%
  filter(year %in% c(2023, 2043)) %>%
  ggplot(aes(x = as.factor(year), y = freq, alpha = type, group = interaction(type, year, scenario), fill = scenario)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge(width = 0.9, preserve = "single")) +
  lemon::facet_rep_wrap(~iso, scales = "free_y", ncol = 3)  +
  theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(axis.line = element_line(), plot.background = element_rect(fill = "white")) +
  scale_alpha_manual(values = c(1, 0.4), name = "RDT Scenario:") +
  scale_fill_manual(values = rev(c(MetBrewer::met.brewer("Egypt", 1), "black", MetBrewer::met.brewer("Egypt", 2)[2])),
                    name = "HRP2 Selection Scenario:",
                    labels = c("Best", "Central", "Worst")) +
  scale_y_continuous(labels = scales::percent_format(), name = "Annual False Negative Rates (%)\n") +
  xlab("Year")
fig5

save_figs("prioty_county_23_63", fig5, width = 16, height = 16, plot_dir = "analysis/impact_analysis/plots_chai")

### fig6 impact of delay ----

fig6 <- full_df %>%
  left_join(
    inc_df %>%
      select(iso, id_1, year, scenario,starts_with("inc")) %>%
      mutate(scenario = tolower(scenario),
             year = as.integer(year)) %>%
      filter(scenario == "central") %>%
      filter(year >= 2023) %>%
      select(-scenario, -inc_low, -inc_high)
  ) %>%
  mutate(delay = replace(delay, which(delay >=0), delay[which(delay >= 0)] + 1)) %>%
  group_by(year, iso, type, scenario, delay) %>%
  na.omit %>%
  summarise(freq = weighted.mean(freq_med, inc_med, na.rm = TRUE)) %>%
  filter(iso %in% priority_iso & scenario == "Central") %>%
  mutate(scenario = stringr::str_to_title(scenario)) %>%
  janitor::clean_names() %>%
  ggplot(aes(x = year, y = freq, group = interaction(type, delay), color = as.factor(delay), linetype = type)) +
  geom_line(lwd = 1) +
  lemon::facet_rep_wrap(~iso, scales = "free_y", ncol = 3)  +
  theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(axis.line = element_line(), plot.background = element_rect(fill = "white")) +
  scale_linetype_discrete(name = "RDT Scenario:") +
  scale_color_manual(values = MetBrewer::met.brewer("Hiroshige", n = 10)[c(1, 5:10)],
                     name = "Years to \nswitch RDT:") +
  scale_x_continuous(breaks = c(2023, 2030, 2040, 2050, 2060), name = "\nYear") +
  scale_y_continuous(labels = scales::percent_format(), name = "Annual False Negative Rates (%)\n") +
  xlab("Year")
fig6
save_figs("delay_priorities", fig6, width = 11, height = 9,  plot_dir = "analysis/impact_analysis/plots_chai")

### fig7 impact of delay by impact by iso3 ----

fig7 <- full_df %>%
  filter(scenario == "Central") %>%
  left_join(
    inc_df %>%
      select(iso, id_1, year, scenario,starts_with("inc")) %>%
      mutate(scenario = tolower(scenario),
             year = as.integer(year)) %>%
      filter(scenario == "central") %>%
      filter(year >= 2023) %>%
      select(-scenario, -inc_low, -inc_high)
  ) %>%
  mutate(delay = replace(delay, which(delay >=0), delay[which(delay >= 0)] + 1)) %>%
  group_by(year, type, iso, delay) %>%
  na.omit %>%
  summarise(freq = weighted.mean(freq_med, inc_med, na.rm = TRUE)) %>%
  filter(iso %in% priority_iso) %>%
  pivot_wider(values_from = freq, names_from = type) %>%
  janitor::clean_names() %>%
  arrange(iso) %>%
  fill(no_rdt_switching, .direction = "up") %>%
  filter(delay >=0) %>%
  mutate(impact = no_rdt_switching - x5_percent_threshold_strategy) %>%
  ungroup() %>%
  select(year, iso, delay, impact) %>%
  ggplot(aes(x = year, y = impact, color = as.factor(delay))) +
  geom_line(lwd = 1) +
  lemon::facet_rep_wrap(~iso, scales = "free_y", ncol = 3)  +
  theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(axis.line = element_line()) +
  scale_color_manual(values = MetBrewer::met.brewer("Hiroshige", n = 10)[c(5:10)],
                     name = "Years to \nswitch RDT:") +
  theme(legend.position = "right", plot.background = element_rect(fill = "white")) +
  ylab("Annual Averted False Neagative Rates (%) \n(Averted = No RDT Switching - 5% Threshold)\n") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_x_continuous(breaks = c(2023, 2030, 2040, 2050, 2060), name = "\nYear")
fig7
save_figs("delay_impact_priorities", fig7, width = 11, height = 9,  plot_dir = "analysis/impact_analysis/plots_chai")

# ----------------------- #
# 2. Make malaria burden figures and tables for CHAI Report ---------
# ----------------------- #

# overall
pre_df2 <- full_df %>%
  filter(delay <= 0) %>%
  left_join(
    inc_df %>%
      select(iso, id_1, year, scenario, starts_with("pop")) %>%
      mutate(scenario = tolower(scenario),
             year = as.integer(year)) %>%
      filter(scenario == "central") %>%
      filter(year >= 2023 & year <= 2060) %>%
      select(-scenario)
  ) %>%
  na.omit %>%
  mutate(scenario = stringr::str_to_title(scenario))

iso_plotting_df2 <- pre_df2 %>%
  group_by(year, iso, type, scenario) %>%
  summarise(micro_2_10 = weighted.mean(micro_2_10_med, pop_total, na.rm = TRUE,))

all_plotting_df2 <- pre_df2 %>%
  group_by(year, type, scenario) %>%
  summarise(micro_2_10 = mean(micro_2_10_med, na.rm = TRUE))

### fig8 all country plot -----

fig8 <- iso_plotting_df2 %>%
  ggplot(aes(year, micro_2_10, color = scenario, linetype = type)) +
  geom_line(lwd = 1) +
  lemon::facet_rep_wrap(~iso, scales = "free_y", ncol = 5) +
  theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(axis.line = element_line()) +
  scale_color_manual(values = rev(c(MetBrewer::met.brewer("Egypt", 1), "black", MetBrewer::met.brewer("Egypt", 2)[2])),
                     name = "HRP2 Selection Scenario:",
                     labels = c("Best", "Central", "Worst")) +
  scale_linetype(name = "RDT Scenario:") +
  theme(legend.position = "top", plot.background = element_rect(fill = "white")) +
  ylab("Malaria Microscopy Prevalence 2-10 (%)\n") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_x_continuous(breaks = c(2023, 2030, 2040, 2050, 2060), name = "\nYear")
fig8
save_figs("micro_all_country_comparisons", fig8, width = 14, height = 14, plot_dir = "analysis/impact_analysis/plots_chai")

### fig9 priority country plot ----

fig9 <- fig8
fig9$data <- filter(fig9$data, iso %in% priority_iso)
fig9 <- fig9 +
  lemon::facet_rep_wrap(~iso, scales = "free_y", ncol = 3) + theme(legend.position = "right")
fig9
save_figs("micro_priority_comparisons", fig9, width = 10, height = 9, plot_dir = "analysis/impact_analysis/plots_chai")

### fig10 all together plot -----

fig10 <- all_plotting_df2 %>%
  ggplot(aes(year, micro_2_10, color = scenario, linetype = type)) +
  geom_line(lwd = 1) +
  theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(axis.line = element_line()) +
  scale_color_manual(values = rev(c(MetBrewer::met.brewer("Egypt", 1), "black", MetBrewer::met.brewer("Egypt", 2)[2])),
                     name = "HRP2 Selection Scenario:",
                     labels = c("Best", "Central", "Worst")) +
  scale_linetype(name = "RDT Scenario:") +
  theme(legend.position = "right", plot.background = element_rect(fill = "white")) +
  ylab("Malaria Microscopy Prevalence 2-10 (%)\n") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_x_continuous(breaks = c(2023, 2030, 2040, 2050, 2060), name = "\nYear")
save_figs("micro_africa_comparisons", fig10, width = 8, height = 6, plot_dir = "analysis/impact_analysis/plots_chai")


### fig11 all together difference -----

fig11 <- all_plotting_df2 %>%
  pivot_wider(values_from = micro_2_10, names_from = type) %>%
  janitor::clean_names() %>%
  mutate(impact = no_rdt_switching - x5_percent_threshold_strategy) %>%
  ungroup() %>%
  select(year, scenario, impact) %>%
  pivot_wider(values_from = impact, names_from = scenario) %>%
  janitor::clean_names() %>%
  pivot_longer(-year) %>%
  ggplot(aes(x = year, y = value, color = name)) +
  geom_line(lwd = 1) +
  theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(axis.line = element_line()) +
  scale_color_manual(values = rev(c(MetBrewer::met.brewer("Egypt", 1), "black", MetBrewer::met.brewer("Egypt", 2)[2])),
                     name = "HRP2 Selection Scenario:",
                     labels = c("Best", "Central", "Worst")) +
  theme(legend.position = "right", plot.background = element_rect(fill = "white")) +
  ylab("Annual Averted Difference in \nMalaria Microscopy Prevalence 2-10 (%) \n(Averted = No RDT Switching - 5% Threshold)\n") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_x_continuous(breaks = c(2023, 2030, 2040, 2050, 2060), name = "\nYear")
save_figs("micro_africa_impact", fig11, width = 8, height = 6, plot_dir = "analysis/impact_analysis/plots_chai")


### fig12 priority country 2029 and 2035impact ----

fig12 <- iso_plotting_df2 %>%
  filter(iso %in% priority_iso) %>%
  janitor::clean_names() %>%
  filter(year %in% c(2023, 2043)) %>%
  ggplot(aes(x = as.factor(year), y = micro_2_10, alpha = type, group = interaction(type, year, scenario), fill = scenario)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge(width = 0.9, preserve = "single")) +
  lemon::facet_rep_wrap(~iso, scales = "free_y", ncol = 3)  +
  theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(axis.line = element_line(), plot.background = element_rect(fill = "white")) +
  scale_alpha_manual(values = c(1, 0.4), name = "RDT Scenario:") +
  scale_fill_manual(values = rev(c(MetBrewer::met.brewer("Egypt", 1), "black", MetBrewer::met.brewer("Egypt", 2)[2])),
                    name = "HRP2 Selection Scenario:",
                    labels = c("Best", "Central", "Worst")) +
  scale_y_continuous(labels = scales::percent_format(), name = "Malaria Microscopy Prevalence 2-10 (%)\n") +
  xlab("Year")
fig12

save_figs("micro_prioty_county_23_63", fig12, width = 14, height = 10, plot_dir = "analysis/impact_analysis/plots_chai")

### fig13 impact of delay -----

fig13 <- full_df %>%
  left_join(
    inc_df %>%
      select(iso, id_1, year, scenario, starts_with("pop")) %>%
      mutate(scenario = tolower(scenario),
             year = as.integer(year)) %>%
      filter(scenario == "central") %>%
      filter(year >= 2023 & year <= 2060) %>%
      select(-scenario)
  ) %>%
  na.omit %>%
  mutate(scenario = stringr::str_to_title(scenario))  %>%
  mutate(delay = replace(delay, which(delay >=0), delay[which(delay >= 0)] + 1)) %>%
  group_by(year, iso, type, scenario, delay) %>%
  na.omit %>%
  summarise(micro_2_10 = weighted.mean(micro_2_10_med, pop_total, na.rm = TRUE)) %>%
  filter(iso %in% priority_iso & scenario == "Central") %>%
  janitor::clean_names() %>%
  ggplot(aes(x = year, y = micro_2_10, group = interaction(type, delay), color = as.factor(delay), linetype = type)) +
  geom_line(lwd = 1) +
  lemon::facet_rep_wrap(~iso, scales = "free_y", ncol = 3)  +
  theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(axis.line = element_line(), plot.background = element_rect(fill = "white")) +
  scale_linetype_discrete(name = "RDT Scenario:") +
  scale_color_manual(values = MetBrewer::met.brewer("Hiroshige", n = 10)[c(1, 5:10)],
                     name = "Years to \nswitch RDT:") +
  scale_x_continuous(breaks = c(2023, 2030, 2040, 2050, 2060), name = "\nYear") +
  scale_y_continuous(labels = scales::percent_format(), name = "Malaria Microscopy Prevalence 2-10 (%)\n") +
  xlab("Year")
fig13
save_figs("micro_delay_priorities", fig13, width = 11, height = 9,  plot_dir = "analysis/impact_analysis/plots_chai")

### fig 14 impact of delay by impact by iso3 ----

fig14 <-  full_df %>%
  left_join(
    inc_df %>%
      select(iso, id_1, year, scenario, starts_with("pop")) %>%
      mutate(scenario = tolower(scenario),
             year = as.integer(year)) %>%
      filter(scenario == "central") %>%
      filter(year >= 2023 & year <= 2060) %>%
      select(-scenario)
  ) %>%
  na.omit %>%
  mutate(scenario = stringr::str_to_title(scenario))  %>%
  mutate(delay = replace(delay, which(delay >=0), delay[which(delay >= 0)] + 1)) %>%
  group_by(year, iso, type, scenario, delay) %>%
  na.omit %>%
  summarise(micro_2_10 = signif(weighted.mean(micro_2_10_med, pop_total, na.rm = TRUE),4)) %>%
  filter(iso %in% priority_iso & scenario == "Central") %>%
  pivot_wider(values_from = micro_2_10, names_from = type) %>%
  janitor::clean_names() %>%
  arrange(iso) %>%
  fill(no_rdt_switching, .direction = "up") %>%
  filter(delay >=0) %>%
  mutate(impact = no_rdt_switching - x5_percent_threshold_strategy) %>%
  ungroup() %>%
  select(year, iso, delay, impact) %>%
  ggplot(aes(x = year, y = impact, color = as.factor(delay))) +
  geom_line(lwd = 1) +
  lemon::facet_rep_wrap(~iso, scales = "free_y", ncol = 3)  +
  theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(axis.line = element_line()) +
  scale_color_manual(values = MetBrewer::met.brewer("Hiroshige", n = 10)[c(5:10)],
                     name = "Years to \nswitch RDT:") +
  theme(legend.position = "right", plot.background = element_rect(fill = "white")) +
  ylab("Annual Averted Difference in \nMalaria Microscopy Prevalence 2-10 (%) \n(Averted = No RDT Switching - 5% Threshold)\n") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_x_continuous(breaks = c(2023, 2030, 2040, 2050, 2060), name = "\nYear")
fig14
save_figs("micro_delay_impact_priorities", fig14, width = 11, height = 9,  plot_dir = "analysis/impact_analysis/plots_chai")


# -------------------------------------------------------------------------#
# 3. Example scenario plots for RDT changes -----------------------
# -------------------------------------------------------------------------#


rdt_roll_out <- data.frame("S.Times" = lubridate::dyears(c(0, c(0)+3.3, c(0)+3.3, 36.8)) + as.Date("2023-04-01"),
                           "RDT" = c(0, 0,1, 1))

complete <- full_df %>%
  mutate(year = lubridate::dyears(t) + as.Date("2023-04-01")) %>%
  left_join(afr_map %>% sf::st_drop_geometry()) %>%
  filter(iso == "KEN" & delay <= 0 & year < as.Date("2060-01-01")) %>%
  filter(id_1 == 10314353) %>%
  filter(scenario == "Central") %>%
  ggplot(aes(year, micro_2_10, group = interaction(id_1,type), color = type)) +
  geom_line(span = 0.1, se = FALSE) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_step(aes(S.Times, RDT), data = rdt_roll_out, color = "black", inherit.aes = FALSE) +
  xlab("Year") +
  theme_bw() +
  MetBrewer::scale_color_met_d("Egypt",name = "RDT Strategy:") +
  theme(legend.position = "top") +
  scale_y_continuous(labels = scales::percent_format(),limits = c(0,1),
                     sec.axis = sec_axis(~ . * 1, labels = scales::percent_format(),
                                         name = "Proportion of novel RDTs introduced to market")  ) +
  ylab("Percentage of False-Negative RDT Results due to pfhrp2/3 deletions")
complete

initial <- full_df %>%
  mutate(year = lubridate::dyears(t) + as.Date("2023-04-01")) %>%
  left_join(afr_map %>% sf::st_drop_geometry()) %>%
  filter(iso == "KEN" & delay <= 0 & year < as.Date("2060-01-01")) %>%
  filter(id_1 == 10314353) %>%
  filter(scenario == "Central" & type == "No RDT Switching") %>%
  ggplot(aes(year, micro_2_10, group = interaction(id_1,type), color = type)) +
  geom_line(span = 0.1, se = FALSE) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  #geom_step(aes(S.Times, RDT), data = rdt_roll_out, color = "black", inherit.aes = FALSE) +
  xlab("Year") +
  theme_bw() +
  MetBrewer::scale_color_met_d("Egypt",name = "RDT Strategy:") +
  theme(legend.position = "top") +
  scale_y_continuous(labels = scales::percent_format(),limits = c(0,1),
                     sec.axis = sec_axis(~ . * 1, labels = scales::percent_format(),
                                         name = "Proportion of novel RDTs introduced to market")  ) +
  ylab("Percentage of False-Negative RDT Results due to pfhrp2/3 deletions")
initial

mid <- full_df %>%
  mutate(year = lubridate::dyears(t) + as.Date("2023-04-01")) %>%
  left_join(afr_map %>% sf::st_drop_geometry()) %>%
  filter(iso == "KEN" & delay <= 0 & year < as.Date("2060-01-01")) %>%
  filter(id_1 == 10314353) %>%
  filter(scenario == "Central" & type == "No RDT Switching") %>%
  ggplot(aes(year, micro_2_10, group = interaction(id_1,type), color = type)) +
  geom_line(span = 0.1, se = FALSE) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_step(aes(S.Times, RDT), data = rdt_roll_out, color = "black", inherit.aes = FALSE) +
  xlab("Year") +
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


# 8. example plot for one region (10314353, 10015081)

selgg <- full_df %>% filter(id_1 == 10314353) %>% select(1:5, contains("med")) %>%
  pivot_longer(contains("med")) %>%
  filter(scenario == "Central") %>%
  ggplot(aes(t, value, color = as.factor(delay))) +
  geom_line() +
  facet_grid(name~type, scales = "free_y") + theme_bw() +
  ggtitle("Setting where selection is predicted")


ids <- unique(full_df$id_1)
pdf("analysis/impact_analysis/testing_plots.pdf", width = 10, height = 8)
for(i in seq_along(ids)) {
message(i)
  gg <- full_df %>% filter(id_1 == ids[i]) %>% select(1:5, matches("freq_med|micro_2_10_med")) %>%
  pivot_longer(contains("med")) %>%
  filter(scenario == "Central") %>%
  ggplot(aes(t, value, color = as.factor(delay))) +
  geom_line() +
  facet_grid(name~type, scales = "free_y") + theme_bw() +
  ggtitle(ids[i]) +
  geom_hline(yintercept = 0.05)

  print(gg)
}
dev.off()

# Issues:


# weird one in freq falling before threshold - 604492677
# when the motnonic spline function is working on the decrease if there is very little freq decrease
# and near to fixation, then decrease in micro etc is not enough - 10313039
# immediate switch threshold weird - 10313109
# weird non monotnic post trigger - 10314195


noselgg <- full_df %>% filter(id_1 == 10015081) %>% select(1:5, contains("med")) %>%
  pivot_longer(contains("med")) %>%
  filter(scenario == "Central") %>%
  ggplot(aes(t, value, color = as.factor(delay))) +
  geom_line() +
  facet_grid(name~type, scales = "free_y") + theme_bw() +
  ggtitle("Setting where selection is not predicted/very slow")

cowplot::plot_grid(selgg, noselgg, ncol = 2)



huh[[289]] %>% filter(id_1 == 10314353 & delay == 4 & type == "1% Threshold Strategy") %>% select(1:5, contains("med")) %>%
  pivot_longer(contains("med")) %>%
  filter(scenario == "Central") %>%
  ggplot(aes(t, value, color = as.factor(delay))) +
  geom_line() +
  facet_grid(name~type, scales = "free_y") + theme_bw() +
  ggtitle("Setting where selection is not predicted/very slow")
