library(tidyverse)

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

# 6. Create map to share with WHO for website
map_for_who <- left_join(hrp2_map$.__enclos_env__$private$map,
          scenario_maps$map %>% select(name_0, name_1, id_1, type_1, source) %>%
            sf::st_drop_geometry())
dir.create("analysis/data_out/who_shp")
sf::write_sf(map_for_who, "analysis/data_out/who_shp/map.shp")

# 7. Alternatively we need to map these to the WHO compliant map
map_from_who <- sf::read_sf("analysis/data_out/who_malaria_shp/WHO_forMalaria.shp")
map_for_who <- sf::read_sf("analysis/data_out/shp_for_who/map.shp")
old_map <- map_for_who

new_map <- map_from_who %>% rename(id_1 = MalariaID) %>% rename(iso = ISO_3_CODE)
new_map$cont <- countrycode::countrycode(new_map$iso, "iso3c", "iso3n")
new_map$region <- UNlocations$area_name[match(new_map$cont, UNlocations$country_code)]

# Changes to fix incorrect mapping from WHO
new_map <- new_map %>%
  mutate(id_1 = replace(id_1, GUID == "{D6886FA7-32DE-4141-85BA-6CC92C4606B4}",10314149)) %>%
  mutate(id_1 = replace(id_1, GUID == "{117B2F72-2F3D-4B0E-8237-84EAFC014B59}",10313882)) %>%
  mutate(id_1 = replace(id_1, GUID == "{5D25739B-8433-4036-ACE9-1DC45B76C98B}",10314646)) %>%
  mutate(id_1 = replace(id_1, GUID == "{DE4B78DA-6CF2-4168-B907-54CA006BF95F}",10316164)) %>%
  mutate(id_1 = replace(id_1, GUID == "{263F24DA-B994-4C44-AA46-3AD342DC86BC}",10313744))


# I don't know where Djibloho, capital of GNQ is in the WHO map but it is not
# showing when I filter for iso == GNQ
# It should be id_1 = 10316356

# get the current map and data
hrp2_map_old <- readRDS("analysis/data_derived/R6_map.rds")

hrp2_map_new <- hrpup:::R6_hrp2_map$new(
  map = new_map,
  map_data = hrp2_map_old$.__enclos_env__$private$map_data,
  scenarios = hrp2_map_old$.__enclos_env__$private$scenarios,
  map_0 = hrp2_map_old$.__enclos_env__$private$map_0
)
saveRDS(hrp2_map_new, "analysis/data_derived/R6_WHO_Compliant_map.rds")





# Code for doinf the map checking that has been copied above
pl_old <- hrp2_map_old$plot(print = FALSE, region = "asia")
pl_new <- hrp2_map_new$plot(print = FALSE, region = "asia")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

plot_sub_map_comp <- function(pl_new, pl_old, iso) {

  pl_new$data <- pl_new$data[pl_new$data$iso == iso,]
  pl_new$layers[[2]] <- NULL
  pl_new$data$id_1 <- as.factor(pl_new$data$id_1)
  pl_new$layers[[1]]$mapping$fill <- as.name("id_1")
  pl_new$layers[[1]]$aes_params$linewidth <- 0.5

  pl_old$data <- pl_old$data[pl_old$data$iso == iso,]
  pl_old$layers[[2]] <- NULL
  pl_old$data$id_1 <- as.factor(pl_old$data$id_1)
  pl_old$layers[[1]]$mapping$fill <- as.name("id_1")
  pl_old$layers[[1]]$aes_params$linewidth <- 0.5

  ids <- unique(c(pl_old$data$id_1, pl_new$data$id_1))
  cols <- gg_color_hue(length(ids))
  names(cols) <- ids

  pl_old <- pl_old + scale_fill_manual(values = cols)
  pl_new <- pl_new + scale_fill_manual(values = cols)


  if(leg) {
    cowplot::plot_grid(pl_old,
                       pl_new,
                       labels = c("old", "new"))
  } else {
  cowplot::plot_grid(pl_old + theme(legend.position = "none"),
                     pl_new + theme(legend.position = "none"),
                     labels = c("old", "new"))
  }
}

plot_risk_comp <- function(pl_new, pl_old, iso) {

  pl_new$data <- pl_new$data[pl_new$data$iso == iso,]
  pl_new$layers[[2]] <- NULL
  pl_new$layers[[1]]$aes_params$linewidth <- 0.5
  pl_new$layers[[1]]$mapping$fill <- as.name("hrp2_risk")

  pl_old$data <- pl_old$data[pl_old$data$iso == iso,]
  pl_old$layers[[2]] <- NULL
  pl_old$layers[[1]]$aes_params$linewidth <- 0.5
  pl_old$layers[[1]]$mapping$fill <- as.name("hrp2_risk")

  cowplot::plot_grid(pl_old, pl_new, labels = c("old", "new"))

}

plot_full_comp <- function(pl_new, pl_old, iso) {
  print(cowplot::plot_grid(plot_sub_map_comp(pl_new, pl_old, iso),
                     plot_risk_comp(pl_new, pl_old, iso), ncol = 1))
}

plot_old_region_id <- function(pl_old, id_1) {

  pl_old$data <- pl_old$data[pl_old$data$id_1 == id_1,]
  pl_old$layers[[2]] <- NULL
  pl_old$data$id_1 <- as.factor(pl_old$data$id_1)
  pl_old$layers[[1]]$mapping$fill <- as.name("id_1")
  pl_old <- pl_old + scale_fill_discrete()
  pl_old$layers[[1]]$aes_params$linewidth <- 0.5
  print(pl_old)

}


leg <- TRUE
plot_full_comp(pl_new, pl_old, "BGD")



library(sf)
point_sf <- st_as_sf(data.frame(x = 1.583673, y = 10.838142), coords = c("x", "y"), crs = st_crs(new_map))

# Calculate distances from the point to each shape in the sf object
valid <- st_is_valid(new_map)
distances <- st_distance(point_sf, new_map[valid,] )

# Find the index of the closest shape
closest_index <- which.min(distances)

# Get the closest shape
closest_shape <- new_map[valid,][closest_index, ]
