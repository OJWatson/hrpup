library(tidyverse)

# -----------------------------------------------------#
# 1. Read in and format data  -----------------
# -----------------------------------------------------#

# read in our data
hrpdat <- readxl::read_xlsx("analysis/data_raw/MTM_PFHRP23_GENE_DELETIONS_20250219.xlsx", sheet = 2) %>%
  rename(HRP2_HRP3_DELETIONS_PERCENT = HRP2_HRP3_PROPORTION_DELETION,
         HRP2_HRP3_DELETIONS_TESTED_N = HRP2_HRP3_TESTED) %>%
         mutate(HRP2_HRP3_DELETIONS_POS_N = round(HRP2_HRP3_DELETIONS_PERCENT * HRP2_HRP3_DELETIONS_TESTED_N)) %>%
  filter(PATIENT_TYPE %in% c(c("Mixed symptomatic and asymptomatic", "Symptomatic")))

# Now just select data we want and filter to studies that did not just measure hrp2
hrpdat <- hrpdat %>%
  mutate(ISO_ctry = countrycode::countrycode(COUNTRY_NAME, "country.name.en", "iso3c")) %>%
  select(ISO_ctry, LATITUDE, LONGITUDE, YEAR_START, YEAR_END, PATIENT_TYPE,SURVEY_TYPE,
                             matches("HRP"), TYPE_SAMPL_ANALYZED) %>%
  mutate(across(tidyselect::everything(),.fns = str_trim)) %>%
  mutate(across(matches("TUDE"), as.numeric))

# convert NRs
conv <- names(which(unlist(lapply(hrpdat, function(x){sum(x == "NR", na.rm = TRUE)>0}))))
hrpdat <- hrpdat %>%
  mutate(across(all_of(conv), .fns = ~as.numeric(replace(.x, .x == "NR", NA))))

# convert numerics
hrpdat <- hrpdat %>%
  mutate(across(LATITUDE, LONGITUDE, PATIENT_TYPE,
                .fns = as.numeric)) %>%
  mutate(across(all_of(grep("HRP", names(hrpdat), value = TRUE)),
                .fns = as.numeric))

hrpdat <- hrpdat %>% filter(!is.na(HRP3_TESTED))

# sort out dates
my_make <- function(x){
  x[!grepl("/", x, fixed = TRUE)] <- paste0("01/", x[!grepl("/", x, fixed = TRUE)])
  x[grep("NA",x)] <- NA
  lubridate::my(x)
}

# sort out middate range for getting prevalence
hrpdat <- hrpdat %>%
  mutate(YEAR_START = my_make(YEAR_START),
         YEAR_END = my_make(YEAR_END)) %>%
  mutate(mid_date = YEAR_START + (YEAR_END - YEAR_START)/2)

# Lastly add in MAP estimated malaria prevalence for these regions
prevs <- hrpdat %>%
  split(~lubridate::year(hrpdat$mid_date)) %>%
  map(function(x){
    malariaAtlas::extractRaster(
      surface = "Plasmodium falciparum PR2 - 10 version 2020",
      df = data.frame(lat = x$LATITUDE, long = x$LONGITUDE),
      year = min(2019, lubridate::year(unique(x$mid_date)))
    )
  })

# last sort of years
hrpdat <- hrpdat %>%
  split(~lubridate::year(hrpdat$mid_date)) %>%
  list_rbind() %>%
  mutate(prev = list_rbind(prevs)$value) %>%
  mutate(prev = replace(prev, prev == -9999, NA)) %>%
  mutate(prev = replace(prev, prev <= 0, 0.00001)) # very low South American prevalence set to 1e-5

hrpdat$continent <- as.factor(countrycode::countrycode(hrpdat$ISO_ctry,"iso3c", "continent"))

hrpdat <- hrpdat %>%
  rename(HRP2_DELETION_PERCENT = HRP2_PROPORTION_DELETION,
         HRP2_DELETION_TESTED_N = HRP2_TESTED) %>%
  mutate(HRP2_DELETION_POS_N = round(HRP2_DELETION_PERCENT * HRP2_DELETION_TESTED_N)) %>%
  rename(HRP3_DELETION_PERCENT = HRP3_PROPORTION_DELETION,
         HRP3_DELETION_TESTED_N = HRP3_TESTED) %>%
  mutate(HRP3_DELETION_POS_N = round(HRP3_DELETION_PERCENT * HRP3_DELETION_TESTED_N))

hprdat <- hrpdat %>%
  select(HRP2_HRP3_DELETIONS_PERCENT, HRP2_DELETION_TESTED_N, HRP3_DELETION_TESTED_N,
         HRP2_HRP3_DELETIONS_TESTED_N, HRP2_DELETION_POS_N, HRP2_DELETION_PERCENT,
         HRP2_HRP3_DELETIONS_POS_N, prev, continent, ISO_ctry)

write.csv(hrpdat, "analysis/data_raw/WHO_hrp2_hrp3_data.csv")

