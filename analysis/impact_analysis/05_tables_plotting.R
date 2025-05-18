## Spread Modelling
library(tidyverse)
library(sf)
full_df <- readRDS("analysis/impact_analysis/data-derived/who_sims.rds")
# file.copy("/home/oj/GoogleDrive/AcademicWork/Imperial/git/arms/analysis/data-derived/spread_mal_inc.rds", "analysis/impact_analysis/data-raw/spread_mal_inc.rds")
inc_df <- readRDS("analysis/impact_analysis/data-raw/spread_mal_inc.rds")

full_df <- full_df %>%
  mutate(scenario = recode(scenario, "worst" = "Worst", "best" = "Best", "central" = "Central")) %>%
  mutate(t = round(t, 1)) %>%
  filter(t %in% seq(0, 20, 1))
# ----------------------- #
# 6. Make figures and tables for CHAI Report ---------
# ----------------------- #

# overall
plotting_df <- full_df %>%
  filter(delay <= 0) %>%
  mutate(year = lubridate::year(lubridate::dyears(t) + as.Date("2023-04-01"))) %>%
  left_join(afr_map %>% sf::st_drop_geometry()) %>%
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
  summarise(freq = weighted.mean(freq, inc_med, na.rm = TRUE,)) %>%
  mutate(scenario = stringr::str_to_title(scenario)) %>%
  mutate(iso = paste0(iso, ": ", countrycode::countrycode(iso, "iso3c", "country.name.en")))


# all country plot
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
save_figs("all_country_comparisons", fig1, width = 14, height = 14, plot_dir = "analysis/impact_analysis/plots_who")

# priority country plot
fig2 <- fig1
priority_countries <- c("DJI", "ERI", "ETH", "GHA", "SSD", "SDN", "DRC", "KEN", "UGA", "ZAM", "MLI", "NGA", "SEN", "SLE", "MWI")
fig2$data <- filter(fig2$data, iso %in%
                      paste0(priority_countries, ": ", countrycode::countrycode(priority_countries, "iso3c", "country.name.en")))
fig2 <- fig2 +
  lemon::facet_rep_wrap(~iso, scales = "free_y", ncol = 3) + theme(legend.position = "right")
save_figs("priority_comparisons", fig2, width = 10, height = 9, plot_dir = "analysis/impact_analysis/plots_who")

# all together plot
fig3 <- full_df %>%
  filter(delay <= 0) %>%
  mutate(year = lubridate::year(lubridate::dyears(t) + as.Date("2023-04-01"))) %>%
  left_join(afr_map %>% sf::st_drop_geometry()) %>%
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
  summarise(freq = weighted.mean(freq, inc_med, na.rm = TRUE)) %>%
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
save_figs("africa_comparisons", fig3, width = 8, height = 6, plot_dir = "analysis/impact_analysis/plots_who")


# all together difference
fig4 <- full_df %>%
  filter(delay <= 0) %>%
  mutate(year = lubridate::year(lubridate::dyears(t) + as.Date("2023-04-01"))) %>%
  left_join(afr_map %>% sf::st_drop_geometry()) %>%
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
  summarise(freq = weighted.mean(freq, inc_med, na.rm = TRUE)) %>%
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
save_figs("africa_impact", fig4, width = 8, height = 6, plot_dir = "analysis/impact_analysis/plots_who")


# priority country 2029 and 2035impact
fig5 <- full_df %>%
  filter(delay <= 0) %>%
  mutate(year = lubridate::year(lubridate::dyears(t) + as.Date("2023-04-01"))) %>%
  left_join(afr_map %>% sf::st_drop_geometry()) %>%
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
  summarise(freq = weighted.mean(freq, inc_med, na.rm = TRUE)) %>%
  filter(iso %in% priority_countries) %>%
  mutate(scenario = stringr::str_to_title(scenario)) %>%
  mutate(iso = paste0(iso, ": ", countrycode::countrycode(iso, "iso3c", "country.name.en"))) %>%
  janitor::clean_names() %>%
  filter(year %in% c(2023, 2060)) %>%
  ggplot(aes(x = as.factor(year), y = freq, alpha = type, group = interaction(type, year, scenario), fill = scenario)) +
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

save_figs("prioty_county_23_63", fig5, width = 16, height = 14, plot_dir = "analysis/impact_analysis/plots_who")

# impact of delay
fig6 <- full_df %>%
  mutate(year = lubridate::year(lubridate::dyears(t) + as.Date("2023-04-01"))) %>%
  left_join(afr_map %>% sf::st_drop_geometry()) %>%
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
  summarise(freq = weighted.mean(freq, inc_med, na.rm = TRUE)) %>%
  filter(iso %in% priority_countries & scenario == "Central") %>%
  mutate(scenario = stringr::str_to_title(scenario)) %>%
  mutate(iso = paste0(iso, ": ", countrycode::countrycode(iso, "iso3c", "country.name.en"))) %>%
  janitor::clean_names() %>%
  ggplot(aes(x = year, y = freq, group = interaction(type, delay), color = as.factor(delay), linetype = type)) +
  geom_line(lwd = 1) +
  lemon::facet_rep_wrap(~iso, scales = "free_y", ncol = 3)  +
  theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(axis.line = element_line(), plot.background = element_rect(fill = "white")) +
  scale_linetype_discrete(name = "RDT Scenario:") +
  scale_color_manual(values = MetBrewer::met.brewer("Hiroshige", n = 10)[c(1, 6:9)],
                     name = "Years to \nswitch RDT:") +
  scale_x_continuous(breaks = c(2023, 2030, 2040, 2050, 2060), name = "\nYear") +
  scale_y_continuous(labels = scales::percent_format(), name = "Annual False Negative Rates (%)\n") +
  xlab("Year")
fig6
save_figs("delay_priorities", fig6, width = 11, height = 9,  plot_dir = "analysis/impact_analysis/plots_who")

# impact of delay by impact by iso3
fig7 <- full_df %>%
  filter(scenario == "Central") %>%
  mutate(year = lubridate::year(lubridate::dyears(t) + as.Date("2023-04-01"))) %>%
  left_join(afr_map %>% sf::st_drop_geometry()) %>%
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
  summarise(freq = weighted.mean(freq, inc_med, na.rm = TRUE)) %>%
  filter(iso %in% priority_countries) %>%
  pivot_wider(values_from = freq, names_from = type) %>%
  janitor::clean_names() %>%
  arrange(iso) %>%
  fill(no_rdt_switching, .direction = "up") %>%
  filter(delay >=0) %>%
  mutate(impact = no_rdt_switching - x5_percent_threshold_strategy) %>%
  ungroup() %>%
  select(year, iso, delay, impact) %>%
  mutate(iso = paste0(iso, ": ", countrycode::countrycode(iso, "iso3c", "country.name.en"))) %>%
  ggplot(aes(x = year, y = impact, color = as.factor(delay))) +
  geom_line(lwd = 1) +
  lemon::facet_rep_wrap(~iso, scales = "free_y", ncol = 3)  +
  theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(axis.line = element_line()) +
  scale_color_manual(values = MetBrewer::met.brewer("Hiroshige", n = 10)[c(6:9)],
                     name = "Years to \nswitch RDT:") +
  theme(legend.position = "right", plot.background = element_rect(fill = "white")) +
  ylab("Annual Averted False Neagative Rates (%) \n(Averted = No RDT Switching - 5% Threshold)\n") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_x_continuous(breaks = c(2023, 2030, 2040, 2050, 2060), name = "\nYear")
fig7
save_figs("delay_impact_priorities", fig7, width = 11, height = 9,  plot_dir = "analysis/impact_analysis/plots_who")

# -------------------------------------------------------------------------#
# 7. Example scenario plots for RDT changes -----------------------
# -------------------------------------------------------------------------#


rdt_roll_out <- data.frame("S.Times" = lubridate::dyears(c(0, c(0)+3.3, c(0)+3.3, 36.8)) + as.Date("2023-04-01"),
                           "RDT" = c(0, 0,1, 1))

complete <- full_df %>%
  mutate(year = lubridate::dyears(t) + as.Date("2023-04-01")) %>%
  left_join(afr_map %>% sf::st_drop_geometry()) %>%
  filter(iso == "KEN" & delay <= 0 & year < as.Date("2060-01-01")) %>%
  filter(id_1 == 10314353) %>%
  filter(scenario == "Central") %>%
  ggplot(aes(year, freq, group = interaction(id_1,type), color = type)) +
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
  ggplot(aes(year, freq, group = interaction(id_1,type), color = type)) +
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
  ggplot(aes(year, freq, group = interaction(id_1,type), color = type)) +
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

noselgg <- full_df %>% filter(id_1 == 10015081) %>% select(1:5, contains("med")) %>%
  pivot_longer(contains("med")) %>%
  filter(scenario == "Central") %>%
  ggplot(aes(t, value, color = as.factor(delay))) +
  geom_line() +
  facet_grid(name~type, scales = "free_y") + theme_bw() +
  ggtitle("Setting where selection is not predicted/very slow")

cowplot::plot_grid(selgg, noselgg, ncol = 2)
