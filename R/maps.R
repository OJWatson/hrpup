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
    #' @param prev_plot Boolean for whether to overlay prevalence classification (default = FALSE)
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
    #' @param disputed_areas Map disputed_areas objects
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
    limits = NULL

  )
)

