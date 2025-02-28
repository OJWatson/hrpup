# All outputs from this are saved within data_derived (MAP change their map on server frequently, so use the outputs
# if trying to replicate)

library(tidyverse)
library(wpp2017)
devtools::load_all()

# 1. MAP admin 1 world map boundaries ---------------
available_admin <- malariaAtlas::listShp(printed = FALSE, admin_level = "admin0")
map_1 <- malariaAtlas::getShp(ISO = available_admin$iso, admin_level = c("admin1")) %>% sf::st_as_sf()
map_1 <- map_1 %>% select(iso, name_0, id_0, name_1, id_1,type_1, source, country_level)

# make smaller for storage
map_small <- rmapshaper::ms_simplify(map_1)
saveRDS(map_1, "analysis/data_derived/admin1_sf.rds")

# 2. MAP lakes separate ---------------
isos <- unique(countrycode::codelist$iso3c[countrycode::codelist$continent == "Africa"])
lakes <- map_1 %>% filter(type_1 == "Water Body")
data("UNlocations")
lakes$cont <- countrycode::countrycode(lakes$iso, "iso3c", "iso3n")
lakes$region <- UNlocations$area_name[match(lakes$cont, UNlocations$country_code)]

# 3. Now to make a WHO appropriate compliant map ---------------

# first they wanted to see the map to help with GUID matching
map_who <- map_1 %>% select(name_0, name_1, id_1, type_1, source)
dir.create("analysis/data_out/shp_for_who")
sf::write_sf(map_for_who, "analysis/data_out/shp_for_who/map.shp")

# this led to the following map being shared
map_from_who <- sf::read_sf("analysis/data_raw/who_malaria_shp/WHO_forMalaria_20231201.shp")
new_map <- map_from_who %>% rename(id_1 = MalariaID) %>% rename(iso = ISO_3_CODE)
new_map$cont <- countrycode::countrycode(new_map$iso, "iso3c", "iso3n")
data("UNlocations")
new_map$region <- UNlocations$area_name[match(new_map$cont, UNlocations$country_code)]

# Changes to fix incorrect mapping from WHO
new_map <- new_map %>%
  mutate(id_1 = replace(id_1, GUID == "{D6886FA7-32DE-4141-85BA-6CC92C4606B4}",10314149)) %>%
  mutate(id_1 = replace(id_1, GUID == "{117B2F72-2F3D-4B0E-8237-84EAFC014B59}",10313882)) %>%
  mutate(id_1 = replace(id_1, GUID == "{5D25739B-8433-4036-ACE9-1DC45B76C98B}",10314646)) %>%
  mutate(id_1 = replace(id_1, GUID == "{DE4B78DA-6CF2-4168-B907-54CA006BF95F}",10316164)) %>%
  mutate(id_1 = replace(id_1, GUID == "{263F24DA-B994-4C44-AA46-3AD342DC86BC}",10313744)) %>%
  mutate(id_1 = replace(id_1, GUID == "{E3B04CAF-21F8-4277-B164-AC8F64A0A2AD}",10316356)) %>%
  mutate(id_1 = replace(id_1, GUID == "{7143C554-9A03-48C5-B51A-B2AEE32072DD}",10312917)) %>%
  mutate(id_1 = replace(id_1, GUID == "{40B93326-6640-478F-9A55-63D3B84365F3}",604100834)) %>%
  mutate(id_1 = replace(id_1, GUID == "{DA38E6BE-9E6D-4A24-94B6-EFEADBD13FFD}",604100834)) %>%
  mutate(id_1 = replace(id_1, GUID == "{867F3BE0-7A4F-43B9-93F8-DE159D54316E}",10313904)) %>%
  mutate(id_1 = replace(id_1, GUID == "{5A68645A-C112-4E72-B5FA-DED6C99DBE68}",10313523)) %>%
  mutate(id_1 = replace(id_1, GUID == "{26786DAE-445D-472E-8256-358B405321D4}",10313171)) %>%
  mutate(id_1 = replace(id_1, GUID == "{26786DAE-445D-472E-8256-358B405321D4}",10313171)) %>%
  mutate(id_1 = replace(id_1, GUID == "{D55E9F45-3E59-4811-8E6F-E320E134D57D}",10314208)) %>%
  mutate(id_1 = replace(id_1, GUID == "{3CB4667D-4886-4D7D-A1DF-A24CE2CA916F}",10315750)) %>%
  mutate(id_1 = replace(id_1, GUID == "{B0C05F5D-21ED-4248-9D1B-E9CD1405B784}",10315750))

# and assign our names
new_map$name_1 <- admi1$name_1[match(new_map$id_1, admi1$id_1)]

# fetch the disputed regions as well
disputed_areas <- sf::read_sf(
  "analysis/data_raw/who_default_shapes/Detailed_Boundary_Disputed_Areas/Detailed_Boundary_Disputed_Areas.shp"
)
disputed_areas$region <- c("Asia",             # Jammu and Kashmir
                           "Africa",           # Lake Albert
                           "Africa",           # Lake Victoria
                           "Africa",           # Lake Tanganyika
                           "Africa",           # Lake Malawi
                           "Asia",             # Aral Sea
                           NA,                 # Great Lakes of NA
                           NA,                 # Great Lakes of NA
                           NA,                 # Great Lakes of NA
                           "Latin America and the Caribbean", # Lake Titicaca
                           "Asia",             # Aksai Chin
                           "Africa",           # Abyei
                           "Africa")           # Western Sahara

disputed_borders <- sf::read_sf(
  "analysis/data_raw/who_default_shapes/Detailed_Boundary_Disputed_Borders/Detailed_Boundary_Disputed_Borders.shp"
)

disputed_borders$region <- regions <- c("Asia",              # J&K (IND Claim)
                                        "Asia",              # J&K (PAK Claim)
                                        "Asia",              # Korean DMZ
                                        "Asia",              # Gaza Strip
                                        "Asia",              # West Bank
                                        "Asia",              # J&K Line of Control
                                        "Africa",            # Bir Tawil (SDN Claim)
                                        "Africa",            # Halayib Triangle (SDN Claim)
                                        "Africa",            # Ilemi Triangle
                                        "Africa",            # Halayib Triangle (EGY Claim)
                                        "Asia",              # Aksai Chin (IND Claim)
                                        "Asia",              # Arunachal Pradesh
                                        NA,            # Kosovo
                                        "Africa",            # Sudan-South Sudan
                                        "Africa",            # Abyei (SSD Claim)
                                        "Africa",            # Sudan-South Sudan
                                        "Africa",            # Abyei (SDN Claim)
                                        "Asia",              # Jammu and Kashmir
                                        "Africa",            # Western Sahara
                                        "Africa",            # Bir Tawil (EGY Claim)
                                        "Asia",              # Aksai Chin (CHN Claim)
                                        "Africa")            # Western Sahara (coastline)

disputed_areas <- disputed_areas[!is.na(disputed_areas$region),]
disputed_borders <- disputed_borders[!is.na(disputed_borders$region),]

# Get the correct admin 0 map as well
adm0 <- sf::read_sf("analysis/data_raw/who_default_shapes/Detailed_Boundary_ADM0_2916330977523091315/GLOBAL_ADM0.shp")
adm0$cont <- countrycode::countrycode(adm0$ISO_3_CODE, "iso3c", "iso3n")
data("UNlocations")
adm0$region <- UNlocations$area_name[match(adm0$cont, UNlocations$country_code)]

# Now we have everything needed for doing WHO mapping etc
full_shp <- list(
  admin0 = adm0,
  admin1 = map_1,
  lakes = lakes,
  disputed_areas = disputed_areas,
  disputed_borders = disputed_borders
)
saveRDS(full_shp, "analysis/data_derived/who_shps.rds")

