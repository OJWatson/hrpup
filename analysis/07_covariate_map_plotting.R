# ---------------------------------------------------- #
# 0. Get the world map and our data  ----
# ---------------------------------------------------- #
library(tidyverse)

# MAP world map boundaries
available_admin <- malariaAtlas::listShp(printed = FALSE, admin_level = "admin0")
world_map_0 <- malariaAtlas::getShp(ISO = available_admin$iso, admin_level = c("admin0")) %>% sf::st_as_sf()

# covariate outputs
covars <- readRDS("analysis/data_derived/global_covariate_ranges.rds")
world_map_1 <- malariaAtlas::getShp(ISO = na.omit(unique(covars$iso3c)), admin_level = c("admin1")) %>% sf::st_as_sf()
mapped <- merge(world_map_1, covars)

# add regions to both
library(wpp2017)
data("UNlocations")
mapped$cont <- countrycode::countrycode(mapped$iso, "iso3c", "iso3n")
mapped$region <- UNlocations$area_name[match(mapped$cont, UNlocations$country_code)]
mapped <- filter(mapped, type_1 != "Water Body")

# get the admin 0 boundaries
mapped_0 <- world_map_0
mapped_0$cont <- countrycode::countrycode(mapped_0$iso, "iso3c", "iso3n")
mapped_0$region <- UNlocations$area_name[match(mapped_0$cont, UNlocations$country_code)]

# Map plotting function
create_map <- function(region, var, title, limits = FALSE) {

# filter to region we care about
if(region != "global") {
  region <- match(region, c("africa", "asia", "latam"))
  region <- c("Africa", "Asia", "Latin America and the Caribbean")[region]
  mapped <- mapped[mapped$region == region, ]
  mapped_0 <- mapped_0[mapped_0$region == region, ]
}

# generate map
gg_map <- mapped[!is.na(mapped[[var]]),]  %>%
  ggplot2::ggplot() +
  ggplot2::geom_sf(ggplot2::aes_string(fill = var), color = "grey", show.legend = TRUE, lwd = 0.1) +
  scale_fill_viridis_c(name = title, option = "C", direction = 1, end = 0.7)

xlim <- layer_scales(gg_map)$x$get_limits()
ylim <- layer_scales(gg_map)$y$get_limits()

gg_map <- gg_map +
  geom_sf(fill = NA, color = "black", show.legend = FALSE,
          data = mapped_0, lwd = 0.2) +
  coord_sf() +
  theme_void() +
  theme(plot.caption = element_text(face = "italic"), plot.background = element_rect(fill = "white", color = "white"))

if(limits){
  gg_map <- gg_map + xlim(xlim) + ylim(ylim)
}

return(gg_map + theme(text = element_text(size = 20), legend.key.height = unit(2, "cm")))

}

# Which vars are we plotting
vars <- grep("mean", names(covars), value = TRUE)
names(vars)  <- c("Microscopy \nPrevalence \n2-10 Years", "Effective \nTreatment \nSeeking", "Microscopy \nUse",
                  "RDT \nNonadherence")

# Create out maps
map_1 <- create_map(region = "global", var = vars[1], title = names(vars[1]), limits = TRUE)
map_2 <- create_map(region = "global", var = vars[2], title = names(vars[2]), limits = TRUE)
map_3 <- create_map(region = "global", var = vars[3], title = names(vars[3]), limits = TRUE)
map_4 <- create_map(region = "global", var = vars[4], title = names(vars[4]), limits = TRUE)

# Plot all four of the maps
save_figs("supp_micro210_map", map_1, width = 30, height = 8, pdf_plot = FALSE)
save_figs("supp_ft_map", map_2, width = 30, height = 8, pdf_plot = FALSE)
save_figs("supp_microuse_map", map_3, width = 30, height = 8, pdf_plot = FALSE)
save_figs("supp_rdtnonadh_map", map_4, width = 30, height = 8, pdf_plot = FALSE)


