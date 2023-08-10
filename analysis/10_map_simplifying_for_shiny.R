# Create plotting object for the shiny
scenario_maps <- readRDS("analysis/data_derived/scenario_maps_full.rds")

new_obj <- list()

# 1. make a streamlined map object
iso_table <- (rvest::read_html("https://en.wikipedia.org/wiki/ISO_3166-1_numeric") %>% rvest::html_table())[[1]]
library(wpp2017)
data("UNlocations")

map <- scenario_maps$map
map$cont <- countrycode::countrycode(map$iso, "iso3c", "iso3n")
map$region <- UNlocations$area_name[match(map$cont, UNlocations$country_code)]
map <- map %>% select(iso, id_1, region)

new_obj$map <- map

# 2. streamlined data

for(i in seq_along(scenario_maps$map_data)) {

huh <- scenario_maps$map_data[[i]] %>% select(id_1, Micro.2.10, t) %>%
mutate(t = replace(t, t < 0, Inf)) %>%
mutate(t_bin = cut(t, breaks = c(0,6,12,20, Inf), include.lowest = TRUE, right = FALSE)) %>%
select(id_1,  Micro.2.10, t_bin) %>%
  rename(hrp2_risk = t_bin)

huh$hrp2_risk <- rev(c("Marginal", "Slight", "Moderate", "High"))[as.integer(huh$hrp2_risk)]
huh$hrp2_risk <- factor(huh$hrp2_risk, levels = rev(c("Marginal", "Slight", "Moderate", "High")))
scenario_maps$map_data[[i]] <- huh

}
new_obj$map_data <- scenario_maps$map_data

# 3. add scenarios
new_obj$scenarios <- scenario_maps$scenarios

# 4. add admin 0 map
# MAP world map boundaries
available_admin <- malariaAtlas::listShp(printed = FALSE, admin_level = "admin0")
world_map_0 <- malariaAtlas::getShp(ISO = available_admin$iso, admin_level = c("admin0")) %>% sf::st_as_sf()

world_map_0$cont <- countrycode::countrycode(world_map_0$iso, "iso3c", "iso3n")
world_map_0$region <- UNlocations$area_name[match(world_map_0$cont, UNlocations$country_code)]
world_map_0 <- world_map_0 %>% select(iso, region)
new_obj$map_0 <- world_map_0



# 5. create the object for easy plotting

# shrink the map resolution
map_small <- rmapshaper::ms_simplify(new_obj$map)

# create our smaller object for file
hrp2_map <- hrpup:::R6_hrp2_map$new(
  map = map_small,
  map_data = new_obj$map_data,
  scenarios = new_obj$scenarios,
  map_0 = new_obj$map_0
  )
saveRDS(hrp2_map, "analysis/data_derived/R6_map.rds")
