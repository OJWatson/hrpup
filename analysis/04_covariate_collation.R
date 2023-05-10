library(tidyverse)
# ---------------------------------------------------- #
# 1. Get the MAP 2020 Prevalence values globally
# ---------------------------------------------------- #

# use this just for matching codes for MAP data
tf <- tempfile()
download.file("https://cloud.ihme.washington.edu/s/teDKnPGcJnBjJ5F/download?path=%2F2%20-%20Data%20%5BCSV%5D%2F3%20-%20Prevalence%3A%20Plasmodium%20falciparum%20%5BGeoTIFF%5D&files=IHME_MALARIA_2000_2019_PFPR_ADMIN1_Y2020M08D31.CSV", tf)
prev <- read.csv(tf) %>% filter(year == 2019)
#
# "https://data.malariaatlas.org/global-data/json-stat/2022/v1/malaria/admin1/Pf_PR-rate_2010-2020_2-10-years_lci.json"
# "https://data.malariaatlas.org/global-data/json-stat/2022/v1/malaria/admin1/Pf_PR-rate_2010-2020_2-10-years_uci.json"
# "https://data.malariaatlas.org/global-data/json-stat/2022/v1/malaria/admin1/Pf_PR-rate_2010-2020_2-10-years_mean.json"

prev$continent <- countrycode::countrycode(prev$ISO3, "iso3c", "continent")
prev$id_1 <- prev$ADM1_Code
prev_vars <- prev %>% select(ISO3, id_1, mean, lower, upper) %>%
  rename(iso3c = ISO3,
         Micro.2.10_mean = mean,
         Micro.2.10_low = lower,
         Micro.2.10_high = upper)

# ---------------------------------------------------- #
# 2. Get the ISO level covariates
# ---------------------------------------------------- #
iso_covariates <- readRDS("analysis/data_derived/iso_covariates.rds")
iso_covariates <- rename(iso_covariates, iso3c = iso3c)

# Get the ones we need for selection
iso_covariates <-
  iso_covariates %>%
  mutate(microscopy.use_mean = micro) %>%
  select(iso3c, ft, microscopy.use_mean, non_adherence) %>%
  rename(rdt.nonadherence_mean = non_adherence) %>%
  rename(ft_mean = ft)


# ft multipliers for high and low values
upd_ft <- read.csv("analysis/data_raw/from_map/GBD2022_careseeking.csv")
ft_mults <- upd_ft %>% mutate(upp_mult = allsource.pred.high/allsource.pred,
                              low_mult = allsource.pred.low/allsource.pred) %>%
  filter(Year == 2020) %>%
  select(ISO3, upp_mult, low_mult) %>%
  group_by(ISO3) %>%
  summarise(upp_mult = mean(upp_mult),
            low_mult = mean(low_mult)) %>%
  rename(iso3c = ISO3)

# micro/adherence multipliers for high and low values
upd_test <- read.csv("analysis/data_raw/from_map/testing_rates_africa.csv")
micro_mults <- upd_test %>% mutate(upp_mult = post_upper/post_mean,
                                   low_mult = post_lower/post_mean) %>%
  filter(year == 2020) %>%
  select(ISO, upp_mult, low_mult) %>%
  group_by(ISO) %>%
  summarise(upp_mult = mean(upp_mult),
            low_mult = mean(low_mult)) %>%
  rename(iso3c = ISO)

# combine with the interval multipliers for each country to create upper and lowers
iso_covariates <- iso_covariates %>%
  left_join(ft_mults) %>%
  mutate(ft_low = ft_mean * low_mult,
         ft_high = ft_mean * upp_mult) %>%
  select(-c(upp_mult, low_mult)) %>%
  left_join(micro_mults) %>%
  mice::mice(method = "rf") %>%
  complete() %>%
  mutate(microscopy.use_low = microscopy.use_mean * low_mult,
         microscopy.use_high = microscopy.use_mean * upp_mult,
         rdt.nonadherence_low = rdt.nonadherence_mean * low_mult,
         rdt.nonadherence_high = rdt.nonadherence_mean * upp_mult) %>%
  select(-c(upp_mult, low_mult))

# ---------------------------------------------------- #
# 3. fitness and rdt.det ranges for now
# ---------------------------------------------------- #

# add in the fitness and rdt.det data range here
covars <- left_join(prev_vars, iso_covariates, by = "iso3c") %>%
  mutate(fitness_mean = 0.95, fitness_low = 0.9, fitness_high = 0.99) %>%
  mutate(rdt.det_mean = 0.225, rdt.det_low = 0.15, rdt.det_high = 0.3)

# these will then be narrowed using mcmc fitting to ERI/ETH data
saveRDS(covars, "analysis/data_derived/global_covariate_ranges.rds")

# ---------------------------------------------------
