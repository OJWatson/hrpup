#' Generate WHO-Compliant Discrete Plot
#'
#' This function creates a geographic map visualisation of a discrete variable across WHO-defined administrative regions. The map can be customised to display specific regions, highlight disputed areas and borders, and optionally include lakes or a scale bar.
#'
#' @param data A `data.frame` containing the variable of interest and corresponding administrative IDs (`id_1`) for merging with WHO shapefiles.
#' @param var A string specifying the name of the column in `data` to be used for filling regions on the map.
#' @param title A string specifying the title for the legend of the plot.
#' @param who_shps A list of WHO shapefiles that contain spatial data for administrative boundaries, lakes, and disputed regions.
#' @param region A string specifying the geographic region of interest. Options include `"global"`, `"africa"`, `"asia"`, or `"latam"`. Defaults to `"africa"`.
#' @param limits Logical. If `TRUE`, restricts the map's display to the selected region. Defaults to `FALSE`.
#' @param print_bool Logical. If `TRUE`, prints the generated map directly. Defaults to `FALSE`.
#' @param lakes_bool Logical. If `TRUE`, lakes will be displayed on the map. Defaults to `TRUE`.
#' @param disputed_areas_bool Logical. If `TRUE`, disputed areas are highlighted on the map. Defaults to `TRUE`.
#' @param disputed_borders_bool Logical. If `TRUE`, disputed borders are outlined on the map. Defaults to `TRUE`.
#' @param scale_bar_bool Logical. If `TRUE`, adds a scale bar to the map. Defaults to `FALSE`.
#' @param prev_plot Logical. If `TRUE`, overlays the new map on a previous plot. Defaults to `FALSE`.
#' @param prev_var Character. Variable name in data for malaria prevalence argument. Defaults to `Micro.2.10_mean`.
#' @param scale_fill_func Function. Function to call for discrete scale function.
#'   Defaults to `function(...){scale_fill_viridis_d(option = "C", direction = 1, end = 0.7, ...)}`.
#'
#' @return A `ggplot2` object containing the generated map. If `print_bool` is `TRUE`, the map is also displayed.
#'
#' @details
#' - Merges administrative level 1 spatial data with user-provided data for mapping.
#' - Outlines country-level administrative boundaries for visual clarity.
#' - Optional visual elements include lakes, disputed areas, disputed borders, and scale bars.
#' - Uses a discrete colour scale to represent variable values, with a legend for interpretation.
#'
#' @examples
#' \dontrun{
#' who_compliant_discrete_plot(
#'   data = example_data,
#'   var = "infection_rate",
#'   title = "Infection Rate (%)",
#'   who_shps = who_shapefiles,
#'   region = "asia",
#'   lakes_bool = TRUE,
#'   disputed_areas_bool = TRUE,
#'   disputed_borders_bool = TRUE
#' )
#' }
#'
#' @import ggplot2
#' @import ggspatial
#' @importFrom dplyr filter
#'
#' @export
who_compliant_discrete_plot = function(data,
                                       var,
                                       title,
                                       who_shps,
                                       region = "africa",
                                       limits = FALSE,
                                       print_bool = FALSE,
                                       lakes_bool = TRUE,
                                       disputed_areas_bool = TRUE,
                                       disputed_borders_bool = TRUE,
                                       scale_bar_bool = FALSE,
                                       risk = "innate",
                                       prev_plot = FALSE,
                                       prev_var = "Micro.2.10_mean",
                                       scale_fill_func = function(...){
                                         ggplot2::scale_fill_manual(
                                           values = rev(c("#6eb3deff", "#84e1c2ff", "#f8b675ff", "#fa8284ff")),
                                           labels = rev(c("No Data","Marginal", "Slight", "Moderate", "High")),...)})
{

  # checks on inputs
  stopifnot(region %in% c("global", "africa", "asia", "latam"))

  # set up line colors and constants
  adm0_brd <- "#696969ff"

  # create map data for admin regions
  mapped <- merge(who_shps$admin1, data, by = "id_1")

  # get the admin 0 boundaries
  mapped_0 <- who_shps$admin0

  # filter to region we care about
  if(region != "global") {
    region <- match(region, c("africa", "asia", "latam"))
    region <- c("Africa", "Asia", "Latin America and the Caribbean")[region]
    mapped <- mapped[mapped$region == region, ]
    mapped_0 <- mapped_0[mapped_0$region == region, ]
    lakes <- who_shps$lakes[who_shps$lakes$region == region, ]
  } else {
    lakes <- who_shps$lakes
  }

  # generate map
  gg_map <- mapped %>%
    ggplot2::ggplot()

  # get limits of relevant regions first if global
  if(region == "global") {
    # add the admin 1 mappings
    gg_map <- gg_map +
      ggplot2::geom_sf(data = mapped[!is.na(mapped[[var]]),], ggplot2::aes_string(fill = var), color = "#e6e6e6ff", show.legend = TRUE, lwd = 0.05) +
      scale_fill_func(name = title, na.value = "#e6e6e6ff")

    xlim <- layer_scales(gg_map)$x$get_limits()
    ylim <- layer_scales(gg_map)$y$get_limits()
  }

  # add the admin 0 mappings in and some simplifying themes
  gg_map <- gg_map +
    ggplot2::geom_sf(fill = "#e6e6e6ff", color = "#696969ff", show.legend = FALSE,
                     data = mapped_0, lwd = 0.2) +
    ggplot2::coord_sf() +
    ggplot2::theme_void() +
    ggplot2::theme(plot.caption = ggplot2::element_text(face = "italic"),
                   plot.background = ggplot2::element_rect(fill = "white", color = "white"))

  # add the admin 1 mappings
  gg_map <- gg_map +
    ggplot2::geom_sf(ggplot2::aes_string(fill = var), color = "#e6e6e6ff", show.legend = TRUE, lwd = 0.05) +
    scale_fill_func(name = title, na.value = "#e6e6e6ff", na.translate = FALSE)

  # add the admin 0 mappings in and some simplifying themes
  gg_map <- gg_map +
    ggplot2::geom_sf(fill = NA, color = "#696969ff", show.legend = FALSE,
                     data = mapped_0, lwd = 0.2) +
    ggplot2::coord_sf() +
    ggplot2::theme_void() +
    ggplot2::theme(plot.caption = ggplot2::element_text(face = "italic"),
                   plot.background = ggplot2::element_rect(fill = "white", color = "white"))

  # get limits of relevant regions first
  if(region != "global") {
    xlim <- layer_scales(gg_map)$x$get_limits()
    ylim <- layer_scales(gg_map)$y$get_limits()
  }

  # add lakes if needed
  if (lakes_bool) {
    gg_map <- gg_map + geom_sf(
      data = lakes, fill = "#ffffffff", color = "#ffffffff", lwd = 0
    )
  }

  # add prevalence mapping
  if(prev_plot) {
    if(any((gg_map$data[[prev_var]] < 0.0005))) {

      gg_map <- gg_map +
        ggnewscale::new_scale_fill() +
        ggpattern::geom_sf_pattern(
          pattern_fill = "grey", pattern = "stripe", fill = NA, show.legend = FALSE,
          color = NA,
          pattern_colour = NA,
          pattern_density = 0.25,
          pattern_spacing = 0.025,
          data = gg_map$data[gg_map$data[[prev_var]] < 0.0005,], inherit.aes = FALSE) +
        scale_fill_manual(name="\nTransmission", labels="Unstable (<0.05% PfPR)", values="grey")

    }
  }

  # add disputed_areas if needed
  if (disputed_areas_bool) {
    gg_map <- gg_map +
      geom_sf(data = mapped[is.na(mapped[[var]]), ],
              fill = "#e6e6e6ff", inherit.aes = FALSE,
              aes(color = "No data"), lwd = 0, show.legend = FALSE) +
      geom_sf(data = who_shps$disputed_areas %>%
                filter(NAME %in% c("Western Sahara", "Abyei", "Jammu and Kashmir")),
              fill = "#b2b2b2ff", inherit.aes = FALSE,
              aes(color = "Not applicable"), lwd = 0) +
      ggplot2::geom_sf(fill = NA, color = "#696969ff", show.legend = FALSE,
                       data = mapped_0, lwd = 0.2)+
      geom_sf(data = who_shps$disputed_areas %>% filter(grepl("Lake|Sea", NAME)),
              fill = "#ffffffff", color = "grey", lwd = 0.1, inherit.aes = FALSE) +
      guides(color = guide_legend(theme = theme(legend.title = element_text(color = "white"))))
  }

  # add disputed_areas if needed
  if (disputed_borders_bool) {
    gg_map <- gg_map +
      geom_sf(data = who_shps$disputed_borders %>%
                filter(NAME %in% c("J&K (IND Claim)", "J&K (PAK Claim)", "Korean DMZ", "Gaza Strip",
                                   "West Bank", "Bir Tawil (SDN Claim)", "Halayib Triangle (SDN Claim)",
                                   "Sudan-South Sudan")),
              color = "white", linetype = "dashed", lwd = 0.2) +
      geom_sf(data = who_shps$disputed_borders %>%
                filter(NAME %in% c("J&K Line of Control", "Ilemi Triangle", "Abyei (SSD Claim)","Abyei (SDN Claim)")),
              color = "#696969ff", linetype = "dotted", lwd = 0.2) +
      geom_sf(data = who_shps$disputed_borders %>%
                filter(NAME %in% c("Halayib Triangle (EGY Claim)", "Aksai Chin (IND Claim)",
                                   "Arunachal Pradesh", "Jammu and Kashmir", "Western Sahara",
                                   "Western Sahara (coastline)", "Bir Tawil (EGY Claim)",
                                   "Aksai Chin (CHN Claim)")),
              color = "#696969ff", linetype = "solid", lwd = 0.2)

  }

  # add scale bar if needed
  if (scale_bar_bool) {
    gg_map <- gg_map +
      ggspatial::annotation_scale(location = 'tr', width_hint = 0.15)
  }

  # set the limits
  if (limits) {
    gg_map <- gg_map + xlim(xlim) + ylim(ylim)
  }

  # last bit of styling
  gg_map <- gg_map + theme(
    plot.margin = margin(20,20,20,20),
    legend.key.spacing.y = unit(0.25, "cm"),
    legend.key.width = unit(1, "cm"),
    # Increase vertical white space between keys
    legend.key.spacing = unit(0.125, "cm")  # General spacing for the entire legend
  ) +
    guides(color = guide_legend(theme = theme(legend.title = element_text(color = "white"))))

  # print the map
  if(print_bool){
    print(gg_map)
  }
  invisible(gg_map)

}


#' Generate WHO-Compliant Continuous Plot
#'
#' This function generates a geographic map plot displaying WHO-compliant continuous cutoffs for a specified variable across administrative regions. The map can be customised to show specific geographic regions, highlight disputed areas and borders, and include lakes or scale bars.
#'
#' @param data A `data.frame` containing the variable of interest and corresponding administrative IDs (`id_1`) to merge with WHO shapefiles.
#' @param var A string specifying the column name in `data` to be plotted as the fill aesthetic.
#' @param title A string specifying the title for the plot legend.
#' @param who_shps A list of WHO shapefiles containing spatial data for administrative boundaries, lakes, and disputed regions.
#' @param region A string specifying the region of interest. Options are `"global"`, `"africa"`, `"asia"`, or `"latam"`. Defaults to `"africa"`.
#' @param limits Logical. If `TRUE`, limits the plotted data to the specified region. Defaults to `FALSE`.
#' @param print_bool Logical. If `TRUE`, prints the generated map. Defaults to `FALSE`.
#' @param lakes_bool Logical. If `TRUE`, includes lakes in the map. Defaults to `TRUE`.
#' @param disputed_areas_bool Logical. If `TRUE`, highlights disputed areas on the map. Defaults to `TRUE`.
#' @param disputed_borders_bool Logical. If `TRUE`, outlines disputed borders on the map. Defaults to `TRUE`.
#' @param scale_bar_bool Logical. If `TRUE`, adds a scale bar to the map. Defaults to `FALSE`.
#' @param prev_plot Logical. If `TRUE`, layers the new map on top of a previous plot. Defaults to `FALSE`.
#' @param prev_var Character. Variable name in data for malaria prevalence argument. Defaults to `Micro.2.10_mean`.
#' @param scale_fill_func Function. Function to call for continuous scale function.
#'   Defaults to `function(...){scale_fill_viridis_c(option = "C", direction = 1, end = 0.7, ...)}`.
#' @param cut_off List. If not null, this list used to add a layer based on cut offs for the main aesthetics. List
#'   should have named elements `c("breaks", "limits", "color", "cutoff", "name", "label")`. Defaults to `NULL`
#'
#' @return A `ggplot2` object representing the map, which can be printed or further modified. If `print_bool` is `TRUE`, the map is also printed.
#'
#' @details
#' - Administrative level 1 regions are merged with the provided data for plotting.
#' - Administrative level 0 (country-level) boundaries are overlaid for context.
#' - Options to include lakes, disputed areas, and disputed borders enhance the map’s contextual accuracy.
#' - The function supports highlighting regions where the variable exceeds specific thresholds.
#'
#' @examples
#' \dontrun{
#' who_compliant_continuous_cutoff_plot(
#'   data = my_data,
#'   var = "disease_prevalence",
#'   title = "Prevalence (%)",
#'   who_shps = my_who_shapefiles,
#'   region = "africa",
#'   lakes_bool = TRUE,
#'   disputed_areas_bool = TRUE,
#'   disputed_borders_bool = TRUE
#' )
#' }
#'
#' @import ggplot2
#' @import ggspatial
#' @importFrom dplyr filter
#' @importFrom ggnewscale new_scale_fill
#'
#' @export
who_compliant_continuous_plot = function(data,
                                         var,
                                         title,
                                         who_shps,
                                         region = "africa",
                                         limits = FALSE,
                                         print_bool = FALSE,
                                         lakes_bool = TRUE,
                                         disputed_areas_bool = TRUE,
                                         disputed_borders_bool = TRUE,
                                         scale_bar_bool = FALSE,
                                         risk = "innate",
                                         prev_plot = FALSE,
                                         prev_var = "Micro.2.10_mean",
                                         scale_fill_func = function(...){scale_fill_viridis_c(option = "C", direction = -1, end = 0.7, ...)},
                                         cut_off = NULL) {

  # checks on inputs
  stopifnot(region %in% c("global", "africa", "asia", "latam"))

  # set up line colors and constants
  adm0_brd <- "#696969ff"

  # create map data for admin regions
  mapped <- merge(who_shps$admin1, data, by = "id_1")

  # get the admin 0 boundaries
  mapped_0 <- who_shps$admin0

  # filter to region we care about
  if(region != "global") {
    region <- match(region, c("africa", "asia", "latam"))
    region <- c("Africa", "Asia", "Latin America and the Caribbean")[region]
    mapped <- mapped[mapped$region == region, ]
    mapped_0 <- mapped_0[mapped_0$region == region, ]
    lakes <- who_shps$lakes[who_shps$lakes$region == region, ]
  } else {
    lakes <- who_shps$lakes
  }

  # generate map
  gg_map <- mapped %>%
    ggplot2::ggplot()

  # get limits of relevant regions first if global
  if(region == "global") {
    # add the admin 1 mappings
    gg_map <- gg_map +
      ggplot2::geom_sf(data = mapped[!is.na(mapped[[var]]),], ggplot2::aes_string(fill = var), color = "#e6e6e6ff", show.legend = TRUE, lwd = 0.05) +
      scale_fill_func(name = title, na.value = "#e6e6e6ff")

    xlim <- layer_scales(gg_map)$x$get_limits()
    ylim <- layer_scales(gg_map)$y$get_limits()
  }

  # add the admin 0 mappings in and some simplifying themes
  gg_map <- gg_map +
    ggplot2::geom_sf(fill = "#e6e6e6ff", color = "#696969ff", show.legend = FALSE,
                     data = mapped_0, lwd = 0.2) +
    ggplot2::coord_sf() +
    ggplot2::theme_void() +
    ggplot2::theme(plot.caption = ggplot2::element_text(face = "italic"),
                   plot.background = ggplot2::element_rect(fill = "white", color = "white"))

  # add the admin 1 mappings
  gg_map <- gg_map +
    ggplot2::geom_sf(ggplot2::aes_string(fill = var), color = "#e6e6e6ff", show.legend = TRUE, lwd = 0.05) +
    scale_fill_func(name = title, na.value = "#e6e6e6ff")


  # add the admin 0 mappings in and some simplifying themes
  gg_map <- gg_map +
    ggplot2::geom_sf(fill = NA, color = "#696969ff", show.legend = FALSE,
                     data = mapped_0, lwd = 0.2) +
    ggplot2::coord_sf() +
    ggplot2::theme_void() +
    ggplot2::theme(plot.caption = ggplot2::element_text(face = "italic"),
                   plot.background = ggplot2::element_rect(fill = "white", color = "white"))

  # get limits of relevant regions first
  if(region != "global") {
    xlim <- layer_scales(gg_map)$x$get_limits()
    ylim <- layer_scales(gg_map)$y$get_limits()
  }

  # add lakes if needed
  if (lakes_bool) {
    gg_map <- gg_map + geom_sf(
      data = lakes, fill = "#ffffffff", color = "#ffffffff", lwd = 0
    )
  }

  # do our cut off mapping
  if (!is.null(cut_off)) {
    gg_map <- gg_map +
      scale_fill_func(name = "", breaks = cut_off$breaks, limits = cut_off$limits, na.value = "#e6e6e6ff") +
      ggnewscale::new_scale_fill() +
      geom_sf(aes(fill = cut_off$color), color = NA, show.legend = TRUE, data = cut_off$cutoff(gg_map$data), inherit.aes = FALSE) +
      scale_fill_manual(name = cut_off$name, labels = cut_off$label, values = cut_off$color)
  }


  # address legend order for this param combination
  if(prev_plot && is.null(cut_off)) {
    gg_map <- gg_map +
      ggnewscale::new_scale_fill() +
      scale_fill_manual(name = title, na.value = "#e6e6e6ff")
  }

  # add prevalence mapping
  if(prev_plot) {
    if(any((gg_map$data[[prev_var]] < 0.0005))) {

      gg_map <- gg_map +
        ggnewscale::new_scale_fill() +
        ggpattern::geom_sf_pattern(
          pattern_fill = "grey", pattern = "stripe", fill = NA, show.legend = FALSE,
          color = NA,
          pattern_colour = NA,
          pattern_density = 0.25,
          pattern_spacing = 0.025,
          data = gg_map$data[gg_map$data[[prev_var]] < 0.0005,], inherit.aes = FALSE) +
        scale_fill_manual(name="\nTransmission", labels="Unstable (<0.05% PfPR)", values="grey")

    }
  }


  # add disputed_areas if needed
  if (disputed_areas_bool) {
    gg_map <- gg_map +
      geom_sf(data = mapped[is.na(mapped[[var]]), ],
              fill = "#e6e6e6ff", inherit.aes = FALSE,
              aes(color = "No data"), lwd = 0) +
      geom_sf(data = who_shps$disputed_areas %>%
                filter(NAME %in% c("Western Sahara", "Abyei", "Jammu and Kashmir")),
              fill = "#b2b2b2ff", inherit.aes = FALSE,
              aes(color = "Not applicable"), lwd = 0) +
      ggplot2::geom_sf(fill = NA, color = "#696969ff", show.legend = FALSE,
                       data = mapped_0, lwd = 0.2)+
      geom_sf(data = who_shps$disputed_areas %>% filter(grepl("Lake|Sea", NAME)),
              fill = "#ffffffff", color = "grey", lwd = 0.1, inherit.aes = FALSE) +
      guides(color = guide_legend(theme = theme(legend.title = element_text(color = "white"))))
  }

  # address legend order for this param combination
  if( !prev_plot && is.null(cut_off)) {
    gg_map <- gg_map +
      guides(
        fill = guide_colorbar(order = 1),
        color = guide_legend(order = 2)
      )
  }

  # add disputed_areas if needed
  if (disputed_borders_bool) {
    gg_map <- gg_map +
      geom_sf(data = who_shps$disputed_borders %>%
                filter(NAME %in% c("J&K (IND Claim)", "J&K (PAK Claim)", "Korean DMZ", "Gaza Strip",
                                   "West Bank", "Bir Tawil (SDN Claim)", "Halayib Triangle (SDN Claim)",
                                   "Sudan-South Sudan")),
              color = "white", linetype = "dashed", lwd = 0.2) +
      geom_sf(data = who_shps$disputed_borders %>%
                filter(NAME %in% c("J&K Line of Control", "Ilemi Triangle", "Abyei (SSD Claim)","Abyei (SDN Claim)")),
              color = "#696969ff", linetype = "dotted", lwd = 0.2) +
      geom_sf(data = who_shps$disputed_borders %>%
                filter(NAME %in% c("Halayib Triangle (EGY Claim)", "Aksai Chin (IND Claim)",
                                   "Arunachal Pradesh", "Jammu and Kashmir", "Western Sahara",
                                   "Western Sahara (coastline)", "Bir Tawil (EGY Claim)",
                                   "Aksai Chin (CHN Claim)")),
              color = "#696969ff", linetype = "solid", lwd = 0.2)

  }

  # add scale bar if needed
  if (scale_bar_bool) {
    gg_map <- gg_map +
      ggspatial::annotation_scale(location = 'tr', width_hint = 0.15)
  }

  # set the limits
  if (limits) {
    gg_map <- gg_map + xlim(xlim) + ylim(ylim)
  }

  # last bit of styling
  gg_map <- gg_map + theme(
    plot.margin = margin(20,20,20,20),
    text = element_text(family = "Helvetica"),
    legend.key.spacing.y = unit(0.25, "cm"),
    legend.key.width = unit(1, "cm"),
    # Increase vertical white space between keys
    legend.key.spacing = unit(0.125, "cm")  # General spacing for the entire legend
  ) +
    guides(color = guide_legend(theme = theme(legend.title = element_text(color = "white"))))

  # print the map
  if(print_bool){
    print(gg_map)
  }
  invisible(gg_map)

}

#' Add a Custom Position Scale Bar to ggplot
#'
#' This function enables drawing a scale bar on a `ggplot` object and optionally adding an orientation arrow.
#'
#' @param lon Numeric. Longitude of the bottom-left point of the first rectangle.
#' @param lat Numeric. Latitude of the bottom-left point of the first rectangle.
#' @param distanceLon Numeric. Length of each rectangle (distance along longitude).
#' @param distanceLat Numeric. Width of each rectangle (distance along latitude).
#' @param distanceLegend Numeric. Distance between rectangles and legend texts.
#' @param dist.unit Character. Units of distance: `"km"` (default, kilometres), `"nm"` (nautical miles), or `"mi"` (statute miles).
#' @param rec.fill Character. Fill colour of the first rectangle (default: `"white"`).
#' @param rec.colour Character. Border colour of the first rectangle (default: `"black"`).
#' @param rec2.fill Character. Fill colour of the second rectangle (default: `"black"`).
#' @param rec2.colour Character. Border colour of the second rectangle (default: `"black"`).
#' @param legend.colour Character. Colour of the legend text (default: `"black"`).
#' @param legend.size Numeric. Font size of the legend text (default: `3`).
#' @param orientation Logical. If `TRUE` (default), adds an orientation arrow.
#' @param arrow.length Numeric. Length of the arrow in the same units as `dist.unit` (default: `500`).
#' @param arrow.distance Numeric. Distance between the scale bar and the bottom of the arrow (default: `300`).
#' @param arrow.North.size Numeric. Font size of the "N" label on the arrow (default: `6`).
#' @param arrow.color Character. Colour of the orientation arrow and "N" label (default: `"black"`).
#' @param box Logical. If `TRUE` (default), adds a background box for better visibility.
#' @param box.line.color Character. Colour of the box border (default: `"black"`).
#' @param box.fill.color Character. Fill colour of the box (default: `"white"`).
#' @param box.offset Numeric. Offset distance of the box from the scale bar (default: `1`).
#'
#' @return A list of `ggplot` layers that include:
#' - Rectangles for the scale bar.
#' - Legend text for the scale bar.
#'
#' @examples
#' # Add a scale bar to a ggplot map
#'
#' \dontrun{
#' ggplot(data = map_data, aes(x = lon, y = lat)) +
#'   geom_polygon() +
#'   scaleBar(lon = 0, lat = 0, distanceLon = 100, distanceLat = 10, distanceLegend = 15)
#' }
#'
scaleBar <- function(lon,
                     lat,
                     distanceLon,
                     distanceLat,
                     distanceLegend,
                     dist.unit = "km",
                     rec.fill = "white",
                     rec.colour = "black",
                     rec2.fill = "black",
                     rec2.colour = "black",
                     legend.colour = "black",
                     legend.size = 3,
                     box = TRUE,
                     box.line.color = "black",
                     box.fill.color = "white",
                     box.offset = 1) {

  # Result #
  #--------#
  # Return a list whose elements are :
  #   - rectangle : a data.frame containing the coordinates to draw the first rectangle ;
  #   - rectangle2 : a data.frame containing the coordinates to draw the second rectangle ;
  #   - legend : a data.frame containing the coordinates of the legend texts, and the texts as well.
  #
  # Arguments : #
  #-------------#
  # lon, lat : longitude and latitude of the bottom left point of the first rectangle to draw ;
  # distanceLon : length of each rectangle ;
  # distanceLat : width of each rectangle ;
  # distanceLegend : distance between rectangles and legend texts ;
  # dist.units : units of distance "km" (kilometers) (default), "nm" (nautical miles), "mi" (statute miles).
  createScaleBar <- function(lon,
                             lat,
                             distanceLon,
                             distanceLat,
                             distanceLegend,
                             dist.units = "km") {
    # First rectangle
    bottomRight <- gcDestination(
      lon = lon,
      lat = lat,
      bearing = 90,
      dist = distanceLon,
      dist.units = dist.units,
      model = "WGS84"
    )

    topLeft <- gcDestination(
      lon = lon,
      lat = lat,
      bearing = 0,
      dist = distanceLat,
      dist.units = dist.units,
      model = "WGS84"
    )
    rectangle <- cbind(
      lon = c(lon, lon, bottomRight[1, "long"], bottomRight[1, "long"], lon),
      lat = c(lat, topLeft[1, "lat"], topLeft[1, "lat"], lat, lat)
    )
    rectangle <- data.frame(rectangle, stringsAsFactors = FALSE)

    # Second rectangle t right of the first rectangle
    bottomRight2 <- gcDestination(
      lon = lon,
      lat = lat,
      bearing = 90,
      dist = distanceLon * 2,
      dist.units = dist.units,
      model = "WGS84"
    )
    rectangle2 <- cbind(
      lon = c(
        bottomRight[1, "long"],
        bottomRight[1, "long"],
        bottomRight2[1, "long"],
        bottomRight2[1, "long"],
        bottomRight[1, "long"]
      ),
      lat = c(lat, topLeft[1, "lat"], topLeft[1, "lat"], lat, lat)
    )
    rectangle2 <- data.frame(rectangle2, stringsAsFactors = FALSE)

    # Now let's deal with the text
    onTop <- gcDestination(
      lon = lon,
      lat = lat,
      bearing = 0,
      dist = distanceLegend,
      dist.units = dist.units,
      model = "WGS84"
    )
    onTop2 <- onTop3 <- onTop
    onTop2[1, "long"] <- bottomRight[1, "long"]
    onTop3[1, "long"] <- bottomRight2[1, "long"]

    legend <- rbind(onTop, onTop2, onTop3)
    legend <- data.frame(cbind(legend, text = c(0, distanceLon, distanceLon *
                                                  2)),
                         stringsAsFactors = FALSE,
                         row.names = NULL)
    return(list(
      rectangle = rectangle,
      rectangle2 = rectangle2,
      legend = legend
    ))
  }

  res <- c()
  if (box) {
    # Add a background box for better visualization on top of a base map
    topLeft <- gcDestination(
      lon = lon,
      lat = lat,
      bearing = 0,
      dist = distanceLat,
      dist.units = dist.unit,
      model = "WGS84"
    )
    bottomRight <- gcDestination(
      lon = lon,
      lat = lat,
      bearing = 90,
      dist = distanceLon * 2,
      dist.units = dist.unit,
      model = "WGS84"
    )

    boxTopLeft <- gcDestination(
      lon = topLeft[1, "long"],
      lat = topLeft[1, "lat"],
      bearing = 315,
      dist = box.offset,
      dist.units = dist.unit,
      model = "WGS84"
    )
    boxTopRight <- gcDestination(
      lon = bottomRight[1, "long"],
      lat = topLeft[1, "lat"],
      bearing = 45,
      dist = box.offset,
      dist.units = dist.unit,
      model = "WGS84"
    )
    boxBottomRight <- gcDestination(
      lon = bottomRight[1, "long"],
      lat = bottomRight[1, "lat"],
      bearing = 135,
      dist = box.offset,
      dist.units = dist.unit,
      model = "WGS84"
    )
    boxBottomLeft <- gcDestination(
      lon = topLeft[1, "long"],
      lat = bottomRight[1, "lat"],
      bearing = 225,
      dist = box.offset,
      dist.units = dist.unit,
      model = "WGS84"
    )
    bg <- cbind(
      lon = c(
        boxTopLeft[1, "long"],
        boxTopRight[1, "long"],
        boxBottomRight[1, "long"],
        boxBottomLeft[1, "long"],
        boxTopLeft[1, "long"]
      ),
      lat = c(
        boxTopLeft[1, "lat"],
        boxTopRight[1, "lat"],
        boxBottomRight[1, "lat"],
        boxBottomLeft[1, "lat"],
        boxTopLeft[1, "lat"]
      )
    )
    bgdf <- data.frame(bg, stringsAsFactors = FALSE)
    bg <- ggplot2::geom_polygon(
      data = bgdf,
      aes(x = lon, y = lat),
      fill = box.fill.color,
      colour = box.line.color,
      size = 0.2
    )
    res <- c(res, bg)
  }

  laScaleBar <- createScaleBar(
    lon = lon,
    lat = lat,
    distanceLon = distanceLon,
    distanceLat = distanceLat,
    distanceLegend = distanceLegend,
    dist.unit = dist.unit
  )
  # First rectangle
  rectangle1 <- ggplot2::geom_polygon(
    data = laScaleBar$rectangle,
    aes(x = lon, y = lat),
    fill = rec.fill,
    colour = rec.colour
  )

  # Second rectangle
  rectangle2 <- ggplot2::geom_polygon(
    data = laScaleBar$rectangle2,
    aes(x = lon, y = lat),
    fill = rec2.fill,
    colour = rec2.colour
  )

  # Legend
  dtext <- laScaleBar$legend[, "text"]
  dtext[length(dtext)] <- paste(dtext[length(dtext)], dist.unit, sep = "")
  scaleBarLegend <- ggplot2::annotate(
    "text",
    label = dtext,
    x = laScaleBar$legend[, "long"],
    y = laScaleBar$legend[, "lat"],
    size = legend.size,
    colour = legend.colour
  )

  res <- c(res, list(rectangle1, rectangle2, scaleBarLegend))

  return(res)
}



#' Find destination in geographical coordinates
#'
#' Find the destination in geographical coordinates at distance dist and for the given bearing from the starting point given by lon and lat.
#'
#' @param lon longitude (Eastings) in decimal degrees (either scalar or vector)
#' @param lat latitude (Northings) in decimal degrees (either scalar or vector)
#' @param bearing bearing from 0 to 360 degrees (either scalar or vector)
#' @param dist distance travelled (scalar)
#' @param dist.units units of distance "km" (kilometers), "nm" (nautical miles), "mi" (statute miles)
#' @param model choice of ellipsoid model ("WGS84", "GRS80", "Airy", "International", "Clarke", "GRS67")
#' @param Vincenty logical flag, default FALSE
#'
#' @details
#' The bearing argument may be a vector when lon and lat are scalar, representing a single point.
#'
#' @return
#' A matrix of decimal degree coordinates with Eastings in the first column and Northings in the second column.
#'
#' @author Sourced from maptools (archived) by Eric Archer and Roger Bivand
#'
#' @references
#' \url{http://www.movable-type.co.uk/scripts/latlong.html#ellipsoid}
#'
#' The file earlier available at \url{http://williams.best.vwh.net/avform.htm}
#'
#' \url{http://www.movable-type.co.uk/scripts/latlong-vincenty.html#direct}
#'
#' Original reference \url{https://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf}:
#'
#' Vincenty, T. 1975. Direct and inverse solutions of geodesics on the ellipsoid with application of nested equations. Survey Review 22(176):88-93
#'
#' @seealso \code{\link{gzAzimuth}}
#'
#' @examples
#' data(state)
#' res <- gcDestination(state.center$x, state.center$y, 45, 250, "km")
#' plot(state.center$x, state.center$y, asp = 1, pch = 16)
#' arrows(state.center$x, state.center$y, res[, 1], res[, 2], length = 0.05)
#'
#' llist <- vector(mode = "list", length = length(state.center$x))
#' for (i in seq(along = llist)) {
#'   llist[[i]] <- gcDestination(
#'     state.center$x[i],
#'     state.center$y[i],
#'     seq(0, 360, 5), 250, "km"
#'   )
#' }
#' plot(state.center$x, state.center$y, asp = 1, pch = 3)
#' nll <- lapply(llist, lines)
#'
gcDestination <- function(lon, lat, bearing, dist, dist.units = "km",
                          model=NULL, Vincenty=FALSE) {
  # lat, lon : lattitude and longitude in decimal degrees
  # bearing : bearing from 0 to 360 degrees
  # dist : distance travelled
  # dist.units : units of distance "km" (kilometers), "nm" (nautical
  # miles), "mi" (statute miles)
  # model : choice of ellipsoid model ("WGS84", "GRS80", "Airy",
  # "International", "Clarke", "GRS67")

  if (!is.numeric(lon)) stop("lon not numeric")
  if (!is.numeric(lat)) stop("lat not numeric")
  if (!is.numeric(bearing)) stop("bearing not numeric")
  if (!is.numeric(dist)) stop("dist not numeric")

  if (length(lon) != length(lat)) stop("lon and lat differ in length")
  if (length(bearing) > 1L && length(lon) > 1L) stop("length mismatch")
  if (length(bearing) > 1L && length(dist) > 1L) stop("length mismatch")

  as.radians <- function(degrees) degrees * pi / 180
  as.degrees <- function(radians) radians * 180 / pi
  as.bearing <- function(radians) (as.degrees(radians) + 360) %% 360

  ellipsoid <- function(model = "WGS84") {
    switch(model,
           WGS84 = c(a = 6378137, b = 6356752.3142, f = 1 / 298.257223563),
           GRS80 = c(a = 6378137, b = 6356752.3141, f = 1 / 298.257222101),
           Airy = c(a = 6377563.396, b = 6356256.909, f = 1 / 299.3249646),
           International = c(a = 6378888, b = 6356911.946, f = 1 / 297),
           Clarke = c(a = 6378249.145, b = 6356514.86955, f = 1 / 293.465),
           GRS67 = c(a = 6378160, b = 6356774.719, f = 1 / 298.25),
           c(a = NA, b = NA, f = NA)
    )}

  dist <- switch(dist.units,
                 km = dist,
                 nm = dist * 1.852,
                 mi = dist * 1.609344
  )
  lat <- as.radians(lat)
  lon <- as.radians(lon)
  bearing <- as.radians(bearing)

  if (is.null(model)) {
    # Code adapted from JavaScript by Chris Veness
    # (scripts@movable-type.co.uk) at
    # http://www.movable-type.co.uk/scripts/latlong.html#ellipsoid
    #   originally from Ed Williams' Aviation Formulary,
    # http://williams.best.vwh.net/avform.htm
    radius <- 6371
    psi <- dist / radius
    lat2 <- asin(sin(lat) * cos(psi) +  cos(lat) * sin(psi) * cos(bearing))
    lon2 <- lon + atan2(sin(bearing) * sin(psi) * cos(lat), cos(psi) -
                          sin(lat) * sin(lat2))
    if (any(is.nan(lat2)) || any(is.nan(lon2))) warning("Out of range values")
    return(cbind(long=as.degrees(lon2), lat=as.degrees(lat2)))
  }

  ellips <- ellipsoid(model)
  if (is.na(ellips["a"])) stop("no such ellipsoid model")
  if (Vincenty) {
    # Code adapted from JavaScript by Chris Veness
    # (scripts@movable-type.co.uk) at
    # http://www.movable-type.co.uk/scripts/latlong-vincenty-direct.html
    # Original reference (http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf):
    #   Vincenty, T. 1975.  Direct and inverse solutions of geodesics on
    # the ellipsoid with application of nested equations.
    #      Survey Review 22(176):88-93
    dist <- dist * 1000
    sin.alpha1 <- sin(bearing)
    cos.alpha1 <- cos(bearing)
    tan.u1 <- (1 - ellips["f"]) * tan(lat)
    cos.u1 <- 1 / sqrt(1 + (tan.u1 ^ 2))
    sin.u1 <- tan.u1 * cos.u1
    sigma1 <- atan2(tan.u1, cos.alpha1)
    sin.alpha <- cos.u1 * sin.alpha1
    cos.sq.alpha <- 1 - (sin.alpha ^ 2)
    u.sq <- cos.sq.alpha * ((ellips["a"] ^ 2) - (ellips["b"] ^ 2)) /
      (ellips["b"] ^ 2)
    cap.A <- 1 + u.sq / 16384 * (4096 + u.sq * (-768 + u.sq * (320 -
                                                                 175 * u.sq)))
    cap.B <- u.sq / 1024 * (256 + u.sq * (-128 + u.sq * (74 - 47 * u.sq)))

    sigma <- dist / (ellips["b"] * cap.A)
    sigma.p <- 2 * pi
    cos.2.sigma.m <- cos(2 * sigma1 + sigma)
    while(any(abs(sigma - sigma.p) > 1e-12)) {
      cos.2.sigma.m <- cos(2 * sigma1 + sigma)
      sin.sigma <- sin(sigma)
      cos.sigma <- cos(sigma)
      delta.sigma <- cap.B * sin.sigma * (cos.2.sigma.m + cap.B / 4 *
                                            (cos.sigma *
                                               (-1 + 2 * cos.2.sigma.m ^ 2) - cap.B / 6 * cos.2.sigma.m *
                                               (-3 + 4 * sin.sigma ^ 2) * (-3 + 4 * cos.2.sigma.m ^ 2)))
      sigma.p <- sigma
      sigma <- dist / (ellips["a"] * cap.A) + delta.sigma
    }
    tmp <- sin.u1 * sin.sigma - cos.u1 * cos.sigma * cos.alpha1
    lat2 <- atan2(sin.u1 * cos.sigma + cos.u1 * sin.sigma * cos.alpha1,
                  (1 - ellips["f"]) * sqrt(sin.alpha ^ 2 + tmp ^ 2))
    lambda <- atan2(sin.sigma * sin.alpha1, cos.u1 * cos.sigma - sin.u1 *
                      sin.sigma * cos.alpha1)
    cap.C <- ellips["f"] / 16 * cos.sq.alpha * (4 + ellips["f"] *
                                                  (ellips["f"] - 3 * cos.sq.alpha))
    cap.L <- lambda - (1 - cap.C) * ellips["f"] * sin.alpha *
      (sigma + cap.C * sin.sigma * (cos.2.sigma.m + cap.C * cos.sigma *
                                      (-1 + 2 * cos.2.sigma.m ^ 2)))
    lat2 <- as.degrees(lat2)
    lon2 <- as.degrees(lon + cap.L)
  } else {
    # Code adapted from JavaScript by Larry Bogan (larry@go.ednet.ns.ca)
    # at http://www.go.ednet.ns.ca/~larry/bsc/jslatlng.html
    e <- 0.08181922
    radius <- (ellips["a"] / 1000) * (1 - e^2) / ((1 - e^2 *
                                                     sin(lat)^2)^1.5)
    psi <- dist / radius
    phi <- pi / 2 - lat
    arc.cos <- cos(psi) * cos(phi) + sin(psi) * sin(phi) * cos(bearing)
    lat2 <- as.degrees((pi / 2) - acos(arc.cos))
    arc.sin <- sin(bearing) * sin(psi) / sin(phi)
    lon2 <- as.degrees(lon + asin(arc.sin))
  }
  return(cbind(long=lon2, lat=lat2))
}
