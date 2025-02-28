library(tidyverse)

sh1 <- readxl::read_excel("analysis/data_raw/PMI_GF_RDTs by brand name.xlsx", 1)
sh2 <- readxl::read_excel("analysis/data_raw/PMI_GF_RDTs by brand name.xlsx", 2, skip  = 3)

sh1 <- sh1 %>% select(Country, Quantity, `Calendar Year Product Delivered`, `Product Long Name`) %>%
  setNames(c("country", "total", "year", "product")) %>%
  mutate(product = replace(product, grep("PAN", product), "PAN")) %>%
  mutate(product = replace(product, grep("pLDH/pLDH", product), "PAN")) %>%
  mutate(product = replace(product, grep("HRP2", product), "HRP2"))


sh2 <- sh2 %>% select(`Country (simplified)`, `Total RDTs`, `PO Year`, `Description  (Source)`) %>%
  setNames(c("country", "total", "year", "product")) %>%
  mutate(product = replace(product, grep("Pan|pan|PAN", product), "PAN")) %>%
  mutate(product = replace(product, grep("pLDH/pLDH", product, fixed = TRUE), "PAN")) %>%
  mutate(product = replace(product, grep("First Response Malaria Ag. pLDH/HRP2 Combo Card Test ", product, fixed = TRUE), "PAN")) %>%
  mutate(product = replace(product, grep("LDH", product), "HRP2")) %>%
  mutate(product = replace(product, grep("Pv", product), "HRP2")) %>%
  mutate(product = replace(product, grep("Other", product), "Unknown")) %>%
  mutate(product = replace(product, which(!grepl("PAN|Unknown", product)), "HRP2"))

dat <- rbind(sh1, sh2) %>%
  filter(product != "Unknown") %>%
  group_by(country, year) %>%
  summarise(n = sum(total, na.rm = TRUE),
            hrp2 = sum(total[product == "HRP2"], na.rm = TRUE)/n)

non_hrp2_use <- dat %>% filter(year %in% 2018:2021) %>%
  mutate(country = replace(country, country == "Zanzibar", "Tanzania (United Republic)")) %>%
  mutate(iso3c = countrycode::countrycode(country, "country.name.en", "iso3c")) %>%
  group_by(iso3c) %>%
  summarise(hrp2 = weighted.mean(hrp2, n,na.rm=TRUE)) %>%
  as.data.frame() %>%
  filter(!is.na(iso3c))

write.csv(non_hrp2_use, "analysis/data_raw/hrp2_RDT_usage.csv")

