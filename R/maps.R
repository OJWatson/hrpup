#' @title R6 Class for hrp2 results maps
#'
#' @description A results map
#'
#' @importFrom ggplot2 ggplot aes
#' @importFrom ggnewscale new_scale_fill
#' @importFrom ggpattern geom_sf_pattern
#' @importFrom R6 R6Class
R6_hrp2_map <- R6::R6Class(
  classname = "hrp2_map",
  cloneable = FALSE,

  # PUBLIC METHODS
  public = list(

    # INITIALISATION
    #' @description
    #' Create a new hrp2 map result
    #' @param map Map sf object
    #' @param map_data Map result data for each scenario
    #' @param scenarios Scenarios data frame
    #' @param map_0 Admin 0 map sf object
    #' @param lakes Lake sf object
    #' @param disputed_areas Disputed areas sf object
    #' @param disputed_borders Disputed borders sf object
    #' @return A new `hrp2_map` object.
    initialize = function(map,
                          map_data = NULL,
                          scenarios = NULL,
                          map_0 = NULL,
                          lakes = NULL,
                          disputed_areas = NULL,
                          disputed_borders = NULL) {

      private$map <- map
      private$map_data <- map_data
      private$scenarios <- scenarios
      private$map_0 <- map_0
      private$lakes <- lakes
      private$disputed_areas <- disputed_areas
      private$disputed_borders <- disputed_borders

    },

    #' @description
    #' Plot hrp2 risk maps
    #'
    #' Plot a map for a given scenario and at a given regional scale.
    #'
    #' All scenario parameters must be one of "central", "worst", "best"
    #'
    #' @param region One of "global", "africa", "asia", "latam"
    #' @param Micro.2.10 Microscopy prevalence scenario
    #' @param ft Treatment seeking rate scenario
    #' @param microscopy.use Proportion of total testing by microscopy scenario
    #' @param rdt.nonadherence Probability of treating neg. indiviudal scenario
    #' @param fitness Relative fitness of hrp2 deleted parasite scenario
    #' @param rdt.det Chance of hrp2 deleted parasite yielding pos test scenario
    #' @param print_bool Boolean for whether to print the plot as well as return it
    #' @param lakes_bool Boolean for whether to add lakes (default = TRUE)
    #' @param disputed_areas_bool Boolean for whether to add disputed areas (default = TRUE)
    #' @param disputed_borders_bool Boolean for whether to add disputed borders (default = TRUE)
    #' @param scale_bar Position. Default = NULL, which is no scale bar
    #' @param risk What type of risk score to plot, either "innate" (default) or "prospective"
    #' @return ggplot map object silently
    plot = function(region = "africa",
                    Micro.2.10 = "central",
                    ft = "central",
                    microscopy.use = "central",
                    rdt.nonadherence = "central",
                    fitness = "central",
                    rdt.det = "central",
                    print_bool = FALSE,
                    lakes_bool = TRUE,
                    disputed_areas_bool = TRUE,
                    disputed_borders_bool = TRUE,
                    scale_bar = NULL,
                    risk = "innate",
                    prev_plot = FALSE) {

      # checks on inputs
      stopifnot(region %in% c("global", "africa", "asia", "latam"))
      stopifnot(Micro.2.10 %in% c("central", "worst", "best"))
      stopifnot(ft %in% c("central", "worst", "best"))
      stopifnot(microscopy.use %in% c("central", "worst", "best"))
      stopifnot(rdt.nonadherence %in% c("central", "worst", "best"))
      stopifnot(fitness %in% c("central", "worst", "best"))
      stopifnot(rdt.det %in% c("central", "worst", "best"))
      stopifnot(risk %in% c("innate", "prospective"))

      # set up line colors and constants
      adm0_brd <- "#696969ff"

      # find the scenario
      scenario <- which(
        private$scenarios$Micro.2.10 == Micro.2.10 &
          private$scenarios$ft == ft &
          private$scenarios$microscopy.use == microscopy.use &
          private$scenarios$rdt.nonadherence == rdt.nonadherence &
          private$scenarios$fitness == fitness &
          private$scenarios$rdt.det == rdt.det
      )

      # create map data for each region
      if(risk == "innate") {
        var <- "hrp2_risk"
        title <- "Innate Risk \nof HRP2"
      } else if (risk == "prospective") {
        var <- "hrp2_prospective_risk"
        title <- "Prospective Risk \nof HRP2"
      }
      mapped <- merge(private$map, private$map_data[[scenario]])
      #mapped <- mapped[!is.na(mapped[[var]]), ]

      # get the admin 0 boundaries
      mapped_0 <- private$map_0

      # filter to region we care about
      if(region != "global") {
        region <- match(region, c("africa", "asia", "latam"))
        region <- c("Africa", "Asia", "Latin America and the Caribbean")[region]
        mapped <- mapped[mapped$region == region, ]
        mapped_0 <- mapped_0[mapped_0$region == region, ]
        lakes <- private$lakes[private$lakes$region == region, ]
        disputed_areas <- private$disputed_areas[private$disputed_areas$region == region, ]
        disputed_borders <- private$disputed_borders[private$disputed_borders$region == region, ]

      } else {
        lakes <- private$lakes
        disputed_areas <- private$disputed_areas
        disputed_borders <- private$disputed_borders
      }

      # generate map
      gg_map_risk <- mapped %>%
        ggplot2::ggplot()

      # get limits of relevant regions first if global
      if(region == "global" && risk == "innate") {
        # add the admin 1 mappings
        gg_map_risk <- gg_map_risk +
          ggplot2::geom_sf(data = mapped[!is.na(mapped[[var]]),],
                           ggplot2::aes_string(fill = var), color = "#e6e6e6ff", show.legend = TRUE, lwd = 0.05) +
          ggplot2::scale_fill_manual(name = title, drop = FALSE,
                                     values = rev(c("#6eb3deff", "#84e1c2ff", "#f8b675ff", "#fa8284ff")),
                                     labels = rev(c("No Data","Marginal", "Slight", "Moderate", "High")),
                                     na.value = "#e6e6e6ff"
          )

        xlim <- layer_scales(gg_map_risk)$x$get_limits()
        ylim <- layer_scales(gg_map_risk)$y$get_limits()
      }

      if(region == "global" && risk == "prospective") {
        xlim <- c(-89.22461, 168.32507)
        ylim <- c(-33.74390,  38.47211)
      }

      # add the admin 0 mappings in and some simplifying themes
      gg_map_risk <- gg_map_risk +
        ggplot2::geom_sf(fill = "#e6e6e6ff", color = "#696969ff", show.legend = FALSE,
                         data = mapped_0, lwd = 0.2) +
        ggplot2::coord_sf() +
        ggplot2::theme_void() +
        ggplot2::theme(plot.caption = ggplot2::element_text(face = "italic"),
                       plot.background = ggplot2::element_rect(fill = "white", color = "white"))


      gg_map_risk <- gg_map_risk +
        ggplot2::geom_sf(ggplot2::aes_string(fill = var), color = "#e6e6e6ff", show.legend = TRUE, lwd = 0.05) +
        ggplot2::scale_fill_manual(name = title, drop = FALSE,
                                   values = rev(c("#6eb3deff", "#84e1c2ff", "#f8b675ff", "#fa8284ff")),
                                   labels = rev(c("No Data","Marginal", "Slight", "Moderate", "High")),
                                   na.value = "#e6e6e6ff"
        )

      # add prevalence mapping
      if(prev_plot) {
        if(any((gg_map_risk$data$Micro.2.10 < 0.0005))) {

          gg_map_risk <- gg_map_risk +
            ggnewscale::new_scale_fill() +
            ggpattern::geom_sf_pattern(
              pattern_fill = "grey", pattern = "stripe", fill = NA, show.legend = FALSE,
              color = NA,
              pattern_colour = NA,
              pattern_density = 0.5,
              pattern_spacing = 0.025,
              data = . %>% filter(Micro.2.10 < 0.0005), inherit.aes = FALSE) +
            ggplot2::scale_fill_manual(name="\nTransmission", labels="Unstable (<0.05% PfPR)", values="grey")

        }
      }

      # add the admin 0 mappings in and some simplifying themes
      gg_map_risk <- gg_map_risk +
        ggplot2::geom_sf(fill = NA, color = "#696969ff", show.legend = FALSE,
                         data = mapped_0, lwd = 0.2) +
        ggplot2::coord_sf() +
        ggplot2::theme_void() +
        ggplot2::theme(plot.caption = ggplot2::element_text(face = "italic"),
                       plot.background = ggplot2::element_rect(fill = "white", color = "white"))

      # get limits of relevant regions first
      if(region != "global") {
        xlim <- layer_scales(gg_map_risk)$x$get_limits()
        ylim <- layer_scales(gg_map_risk)$y$get_limits()
      }
      # add lakes if needed
      if (lakes_bool) {
        gg_map_risk <- gg_map_risk + geom_sf(
          data = lakes, fill = "#ffffffff", color = "#ffffffff", lwd = 0
        )
      }


      # add disputed_areas if needed
      if (disputed_areas_bool) {
        gg_map_risk <- gg_map_risk +
          geom_sf(data = disputed_areas %>% filter(grepl("Lake|Sea", NAME)),
                  fill = "#ffffffff", color = "grey", lwd = 0.1) +
          geom_sf(data = disputed_areas %>%
                    filter(NAME %in% c("Western Sahara", "Abyei", "Jammu and Kashmir")),
                  aes(fill = "Not applicable"),
                  color = "#b2b2b2ff", lwd = 0.2) +
          ggplot2::scale_fill_manual(name = title, drop = FALSE,
                                     values = rev(c("Not applicable"="#b2b2b2ff","Marginal"="#6eb3deff",
                                                    "Slight"="#84e1c2ff",  "Moderate"="#f8b675ff", "High"="#fa8284ff")),
                                     labels = function(breaks) {breaks[is.na(breaks)] <- "No data"; breaks},
                                     na.translate = TRUE,
                                     na.value = "#e6e6e6ff")

      }

      # add disputed_areas if needed
      if (disputed_borders_bool) {
        gg_map_risk <- gg_map_risk +
          geom_sf(data = disputed_borders %>%
                    filter(NAME %in% c("J&K (IND Claim)", "J&K (PAK Claim)", "Korean DMZ", "Gaza Strip",
                                       "West Bank", "Bir Tawil (SDN Claim)", "Halayib Triangle (SDN Claim)",
                                       "Sudan-South Sudan")),
                  color = "white", linetype = "dashed", lwd = 0.2) +
          geom_sf(data = disputed_borders %>%
                    filter(NAME %in% c("J&K Line of Control", "Ilemi Triangle", "Abyei (SSD Claim)","Abyei (SDN Claim)")),
                  color = "#696969ff", linetype = "dotted", lwd = 0.2) +
          geom_sf(data = disputed_borders %>%
                    filter(NAME %in% c("Halayib Triangle (EGY Claim)", "Aksai Chin (IND Claim)",
                                       "Arunachal Pradesh", "Jammu and Kashmir", "Western Sahara",
                                       "Western Sahara (coastline)", "Bir Tawil (EGY Claim)",
                                       "Aksai Chin (CHN Claim)")),
                  color = "#696969ff", linetype = "solid", lwd = 0.2)

      }

      # add scale bar if needed
      if (!is.null(scale_bar)) {
        gg_map_risk <- gg_map_risk +
          ggspatial::annotation_scale(location = scale_bar, width_hint = 0.15)
      }

      # set the limits
      gg_map_risk <- gg_map_risk + xlim(xlim) + ylim(ylim)

      # print the map
      if(print_bool){
        print(gg_map_risk)
      }
      invisible(gg_map_risk)

    },

    # SETTERS
    #' Set map data
    #' @param map_data Map data of selection coefficients for each scenario
    set_map_data = function(map_data) {private$map_data <- map_data},

    #' Set lakes
    #' @param lakes Map lake objects
    set_lakes = function(lakes) {private$lakes <- lakes},

    #' Set disputed areas
    #' @param lakes Map disputed_areas objects
    set_disputed_areas = function(disputed_areas) {private$disputed_areas <- disputed_areas},

    #' Set disputed_borders
    #' @param disputed_borders Map disputed_borders objects
    set_disputed_borders = function(disputed_borders) {private$disputed_borders <- disputed_borders},

    #' Set limits
    #' @param limits x and y limits for each map_0 region
    set_limits = function(limits) {private$limits <- limits},

    #' Set map data
    #' @param scenarios Scenario data frame that matches map_data
    set_scenarios = function(scenarios) {private$scenarios <- scenarios}

  ),

  private = list(
    map = NULL,
    map_data = NULL,
    scenarios = NULL,
    map_0 = NULL,
    lakes = NULL,
    disputed_areas = NULL,
    disputed_borders = NULL,
    limits = NULL,

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
    scaleBar = function(lon,
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
        bottomRight <- maptools::gcDestination(
          lon = lon,
          lat = lat,
          bearing = 90,
          dist = distanceLon,
          dist.units = dist.units,
          model = "WGS84"
        )

        topLeft <- maptools::gcDestination(
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
        bottomRight2 <- maptools::gcDestination(
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
        onTop <- maptools::gcDestination(
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
        topLeft <- maptools::gcDestination(
          lon = lon,
          lat = lat,
          bearing = 0,
          dist = distanceLat,
          dist.units = dist.unit,
          model = "WGS84"
        )
        bottomRight <- maptools::gcDestination(
          lon = lon,
          lat = lat,
          bearing = 90,
          dist = distanceLon * 2,
          dist.units = dist.unit,
          model = "WGS84"
        )

        boxTopLeft <- maptools::gcDestination(
          lon = topLeft[1, "long"],
          lat = topLeft[1, "lat"],
          bearing = 315,
          dist = box.offset,
          dist.units = dist.unit,
          model = "WGS84"
        )
        boxTopRight <- maptools::gcDestination(
          lon = bottomRight[1, "long"],
          lat = topLeft[1, "lat"],
          bearing = 45,
          dist = box.offset,
          dist.units = dist.unit,
          model = "WGS84"
        )
        boxBottomRight <- maptools::gcDestination(
          lon = bottomRight[1, "long"],
          lat = bottomRight[1, "lat"],
          bearing = 135,
          dist = box.offset,
          dist.units = dist.unit,
          model = "WGS84"
        )
        boxBottomLeft <- maptools::gcDestination(
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

  )
)

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
    bottomRight <- maptools::gcDestination(
      lon = lon,
      lat = lat,
      bearing = 90,
      dist = distanceLon,
      dist.units = dist.units,
      model = "WGS84"
    )

    topLeft <- maptools::gcDestination(
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
    bottomRight2 <- maptools::gcDestination(
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
    onTop <- maptools::gcDestination(
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
    topLeft <- maptools::gcDestination(
      lon = lon,
      lat = lat,
      bearing = 0,
      dist = distanceLat,
      dist.units = dist.unit,
      model = "WGS84"
    )
    bottomRight <- maptools::gcDestination(
      lon = lon,
      lat = lat,
      bearing = 90,
      dist = distanceLon * 2,
      dist.units = dist.unit,
      model = "WGS84"
    )

    boxTopLeft <- maptools::gcDestination(
      lon = topLeft[1, "long"],
      lat = topLeft[1, "lat"],
      bearing = 315,
      dist = box.offset,
      dist.units = dist.unit,
      model = "WGS84"
    )
    boxTopRight <- maptools::gcDestination(
      lon = bottomRight[1, "long"],
      lat = topLeft[1, "lat"],
      bearing = 45,
      dist = box.offset,
      dist.units = dist.unit,
      model = "WGS84"
    )
    boxBottomRight <- maptools::gcDestination(
      lon = bottomRight[1, "long"],
      lat = bottomRight[1, "lat"],
      bearing = 135,
      dist = box.offset,
      dist.units = dist.unit,
      model = "WGS84"
    )
    boxBottomLeft <- maptools::gcDestination(
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
