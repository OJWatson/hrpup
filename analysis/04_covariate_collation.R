library(tidyverse)
# ---------------------------------------------------- #
# 1. Get the MAP 2020 Prevalence values globally
# ---------------------------------------------------- #

# use this just for matching codes for MAP data
tf <- tempfile()
download.file("https://cloud.ihme.washington.edu/s/teDKnPGcJnBjJ5F/download?path=%2F2%20-%20Data%20%5BCSV%5D%2F3%20-%20Prevalence%3A%20Plasmodium%20falciparum%20%5BGeoTIFF%5D&files=IHME_MALARIA_2000_2019_PFPR_ADMIN1_Y2020M08D31.CSV", tf)
prev <- read.csv(tf) %>% filter(year == 2019)

read_map_json <- function(type, prev = "Micro.2.10_mean") {

  url <- paste0(
    "https://data.malariaatlas.org/global-data/json-stat/2022/v1/malaria/admin1/Pf_PR-rate_2010-2020_2-10-years_",
    type,
    ".json")
  lci <- jsonlite::read_json(url)

  df <- data.frame("prev" = unlist(lci$value)) %>%
    cbind(expand.grid(
      "year" = unlist(lci$dimension$Year)[-1],
      "id_1" = as.character(unlist(lci$dimension$ID))[-1]
    ))

  names(df)[1] <- prev

  return(df)

}

lci <- read_map_json("lci", "Micro.2.10_low")
mean <- read_map_json("mean", "Micro.2.10_mean")
uci <- read_map_json("uci", "Micro.2.10_high")
df <- left_join(lci, mean) %>% left_join(uci) %>% mutate(id_1 = as.integer(as.character(id_1)))

prev_vars <- left_join(
  df,
  prev %>%
    select(ISO3, ADM1_Code) %>%
    rename(iso3c = ISO3, id_1 = ADM1_Code)
  ) %>%
  filter(year == 2020) %>%
  select(iso3c, id_1, Micro.2.10_mean, Micro.2.10_low, Micro.2.10_high) %>%
  arrange(iso3c)

# MAP have changed shp file since 2020 and there is
# a change in shape earlier than this in Ghana and Uganda
prev_vars <- prev_vars %>% filter(!(iso3c %in% c("GHA","UGA")))
prev_vars <- rbind(
  prev_vars,
  prev %>%
    filter(ISO3 %in% c("GHA","UGA")) %>%
    rename(iso3c = ISO3, Micro.2.10_mean = mean,
           Micro.2.10_low = lower, Micro.2.10_high = upper,
           id_1 = ADM1_Code) %>%
    mutate(year = 2020) %>%
    select(iso3c, id_1, Micro.2.10_mean, Micro.2.10_low, Micro.2.10_high)
)

# ---------------------------------------------------- #
# 2. Get the ISO level covariates
# ---------------------------------------------------- #
iso_covariates <- readRDS("analysis/data_derived/iso_covariates.rds")
iso_covariates <- rename(iso_covariates, iso3c = iso3c)

# let's not have any country assumed to fully use microscopy or RDT, there is probably still 1% of each
iso_covariates$micro <- pmax(pmin(iso_covariates$micro, 0.99), 0.01)
# same for non adherence
iso_covariates$non_adherence <- pmax(pmin(iso_covariates$non_adherence, 0.99), 0.01)

# Get the ones we need for selection
iso_covariates <-
  iso_covariates %>%
  mutate(microscopy.use_mean = micro) %>%
  select(iso3c, ft, microscopy.use_mean, non_adherence) %>%
  rename(rdt.nonadherence_mean = non_adherence) %>%
  rename(ft_mean = ft)


# ft multipliers for high and low values based on provided intervals from MAP
upd_ft <- read.csv("analysis/data_raw/from_map/careseeking.csv")
ft_mults <- upd_ft %>%
  select(ISO3, High_Mult, Low_Mult) %>%
  group_by(ISO3) %>%
  summarise(upp_mult = mean(High_Mult),
            low_mult = mean(Low_Mult)) %>%
  rename(iso3c = ISO3)

# micro/adherence multipliers for high and low values
upd_test <- read.csv("analysis/data_raw/from_map/testing.csv")
micro_mults <- upd_test %>%
  select(ISO, High_Mult, Low_Mult) %>%
  group_by(ISO) %>%
  summarise(upp_mult = mean(High_Mult),
            low_mult = mean(Low_Mult)) %>%
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
# these will then be narrowed using mcmc fitting to ERI/ETH data
covars <- left_join(prev_vars, iso_covariates, by = "iso3c") %>%
  mutate(fitness_mean = 0.9, fitness_low = 0.8, fitness_high = 0.99) %>%
  mutate(rdt.det_mean = 0.46, rdt.det_low = 0.29, rdt.det_high = 0.65)

# ---------------------------------------------------

# Last set to NA the 0 prev estimates
covars$Micro.2.10_mean[covars$Micro.2.10_mean == 0] <- NA
covars$Micro.2.10_low[covars$Micro.2.10_low == 0] <- NA
covars$Micro.2.10_high[covars$Micro.2.10_high == 0] <- NA
covars <- covars %>% filter(!is.na(iso3c))

saveRDS(covars, "analysis/data_derived/global_covariate_ranges.rds")
