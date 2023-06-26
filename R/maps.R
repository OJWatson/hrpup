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
    #' @return A new `hrp2_map` object.
    initialize = function(map,
                          map_data = NULL,
                          scenarios = NULL,
                          map_0 = NULL) {

      private$map <- map
      private$map_data <- map_data
      private$scenarios <- scenarios
      private$map_0 <- map_0

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
    #' @param print Boolean for whether to print the plot as well as return it
    #' @param risk What type of risk score to plot, either "innate" (default) or "prospective"
    #' @return ggplot map object silently
    plot = function(region = "africa",
                    Micro.2.10 = "central",
                    ft = "central",
                    microscopy.use = "central",
                    rdt.nonadherence = "central",
                    fitness = "central",
                    rdt.det = "central",
                    print = FALSE,
                    risk = "innate") {

      # checks on inputs
      stopifnot(region %in% c("global", "africa", "asia", "latam"))
      stopifnot(Micro.2.10 %in% c("central", "worst", "best"))
      stopifnot(ft %in% c("central", "worst", "best"))
      stopifnot(microscopy.use %in% c("central", "worst", "best"))
      stopifnot(rdt.nonadherence %in% c("central", "worst", "best"))
      stopifnot(fitness %in% c("central", "worst", "best"))
      stopifnot(rdt.det %in% c("central", "worst", "best"))
      stopifnot(risk %in% c("innate", "prospective"))

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
      } else if (risk == "prospective") {
        var <- "hrp2_prospective_risk"
      }
      mapped <- merge(private$map, private$map_data[[scenario]])
      mapped <- mapped[!is.na(mapped[[var]]), ]

      # get the admin 0 boundaries
      mapped_0 <- private$map_0

      # filter to region we care about
      if(region != "global") {
        region <- match(region, c("africa", "asia", "latam"))
        region <- c("Africa", "Asia", "Latin America and the Caribbean")[region]
        mapped <- mapped[mapped$region == region, ]
        mapped_0 <- mapped_0[mapped_0$region == region, ]
      }

      # generate map
      gg_map_risk <- mapped %>%
        ggplot2::ggplot() +
        ggplot2::geom_sf(ggplot2::aes_string(fill = var), color = "grey", show.legend = TRUE, lwd = 0.1) +
        ggplot2::scale_fill_manual(name = "HRP2 Concern", values = rev(c("blue", "cyan", "yellow", "red")),
                                   labels = rev(c("Marginal", "Slight", "Moderate", "High")),
                                   na.value = "white"
        )

      # add prevalence mapping
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

      # add the admin 0 mappings in and some simplifying themes
      gg_map_risk <- gg_map_risk +
        ggplot2::geom_sf(fill = NA, color = "black", show.legend = FALSE,
                         data = mapped_0, lwd = 0.2) +
        ggplot2::coord_sf() +
        ggplot2::theme_void() +
        ggplot2::theme(plot.caption = ggplot2::element_text(face = "italic"),
                       plot.background = ggplot2::element_rect(fill = "white", color = "white"))

      # print the map
      print(gg_map_risk)
      invisible(gg_map_risk)

    },

    # SETTERS
    #' Set map data
    #' @param map_data Map data of selection coefficients for each scenario
    set_map_data = function(map_data) {private$map_data <- map_data},

    #' Set map data
    #' @param scenarios Scenario data frame that matches map_data
    set_scenarios = function(scenarios) {private$scenarios <- scenarios}

  ),

  private = list(
    map = NULL,
    map_data = NULL,
    scenarios = NULL,
    map_0 = NULL
  )
)
