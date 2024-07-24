library(tidyverse)

# MAP world map boundaries
available_admin <- malariaAtlas::listShp(printed = FALSE, admin_level = "admin0")
map_1 <- malariaAtlas::getShp(ISO = available_admin$iso, admin_level = c("admin1")) %>% sf::st_as_sf()
map_0 <- malariaAtlas::getShp(ISO = available_admin$iso, admin_level = c("admin0")) %>% sf::st_as_sf()

map_1 <- map_1 %>% select(iso, name_0, id_0, name_1, id_1,type_1, source, country_level)
# make smaller for storage
map_small <- rmapshaper::ms_simplify(map_1)

# save to file for easier use
saveRDS(map_0, "analysis/data_derived/admin0_sf.rds")
saveRDS(map_1, "analysis/data_derived/admin1_sf.rds")
