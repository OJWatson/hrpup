#' @title R6 Class for simualting spread of hrp2 on a given map
#'
#' @description A simulating model
#'
#' @importFrom R6 R6Class
#' @importFrom spdep poly2nb
#' @importFrom purrr list_rbind
R6_hrp2_spread <- R6::R6Class(
  classname = "hrp2_spread",
  cloneable = FALSE,

  # PUBLIC METHODS
  public = list(

    # INITIALISATION
    #' @description
    #' Create a new hrp2 spread model object
    #' @param map Map sf object
    #' @param hrp2_mod hrp2_mod selection model
    #' @param adj_mat Adjacency matrix for `map`. Default = NULL and will calculate internally
    initialize = function(map, hrp2_mod, adj_mat = NULL) {

      private$map <- map
      private$hrp2_mod <- hrp2_mod

      # check adj_mat has the same names
      if(is.null(adj_mat)) {
        adj_mat <- spdep::poly2nb(map)
        names(adj_mat) <- map$id_1
      } else {
        order_names(adj_mat, as.character(map$id_1))
      }
      private$adj_mat <- adj_mat
    },

    #' @description
    #' Set up simulation seeds
    #' @param seeds Named vector or list of regions seeding spread and iniial frequency
    #'   E.g. list("region_1" = 0.2, "region_2" = 0.3, "region_3" = 0.002)
    #'
    set_seeds = function(seeds){
      if(!all(names(seeds) %in% private$map$id_1)) {
        stop("All seeds names must be in spread model map")
      }
      private$seeds <- seeds
      invisible(private$seeds)
    },

    #' @description
    #' Set up map_data
    #' @param map_data Data.frame with region names ("id_1") and selection ("s")
    #'
    set_map_data = function(map_data){
      if(!all(private$map$id_1 %in% map_data$id_1)) {
        stop("map_data must include values for all spread model map regions")
      }
      private$map_data <- map_data
      invisible(private$map_data)
    },

    #' @description
    #' Simulated spread
    #' @param import_freq What frequency does importation result in. Default = 0.01
    #' @param export_freq At what frequency does exportation occur at. Default = 0.25
    #' @param t_end What year does simulation end. Default = 40
    #' @param t_break Gap between time breaks. Default = 1
    #' @param import_gap Number of years for importation to occur over. Default = 1
    #'
    simulate_spread = function(import_freq = 0.01,
                               export_freq = 0.25,
                               t_end = 40,
                               t_break = 1,
                               import_gap = 1) {

      # set up our results object
      private$set_res_list(t_end = t_end, t_break = t_break)

      # Where are we looking for deletions first
      next_pos <- which(names(private$res_list) %in% names(private$seeds))
      del_pos <- private$find_deleted_regions(t = 1, pos = next_pos)

      # Running vector of regions that have been simulated from
      simulated <- c()

      # Vector of times
      t_s <- seq(0, t_end - t_break, t_break)
      del_pos_list <- vector("list", length(t_s) + as.integer(1/t_break))
      del_pos_list[[1]] <- del_pos

      # Simulate spread process
      for(t in seq_along(t_s)) {

        # 1. Simulate Selection at deleted regions
        private$simulate_selection(t = t, del_pos = del_pos_list[[t]])

        # 2. Update which regions have been simulated
        simulated <- c(simulated, del_pos_list[[t]])

        # 3. Which regions are exporting after selections
        export_pos <- private$find_deleted_regions(t = t, del_freq = export_freq, pos = simulated)

        # 4. If a simulated region has exported then we remove it from simulated
        simulated <- setdiff(simulated, export_pos)

        # 5. Simulate Importation only if less than import gap before end
        if (float_leq(t_s[t], (t_end - import_gap))) {
        del_pos_list[[t + as.integer(1/t_break)]] <- private$simulate_importation(
          t = t, export_pos = export_pos,
          import_freq = import_freq,
          t_break = t_break, import_gap = import_gap
          )
        }
      }

      # And return our simulation
      return(purrr::list_rbind(private$res_list))

    }

  ),

  private = list(

    # Private Member Variables
    map = NULL,
    adj_mat = NULL,
    hrp2_mod = NULL,
    seeds = NULL,
    map_data = NULL,
    res_list = NULL,

    # Private Member Functions

    # Set Up
    set_res_list = function(t_end, t_break = 1){

      res <- expand.grid("id_1" = private$map$id_1, "t" = c(seq(0, t_end - t_break, t_break), t_end), "freq" = 0)
      res$t_pos <- match(res$t, c(seq(0, t_end - t_break, t_break), t_end))
      res_list <- split(res, res$id_1)

      # make it have the same ordering as our map
      res_list <- res_list[match(private$map$id_1, names(res_list))]

      # set up initial regions
      for(i in seq_along(private$seeds)){
        res_list[[names(private$seeds)[i]]]$freq[1] <- private$seeds[i]
      }

      private$res_list <- res_list

    },

    # Simulation functions
    # function to find the deleted regions
    find_deleted_regions = function(t, del_freq = NULL, pos = NULL) {

      # which positions are we finding
      if(is.null(pos)) {
        pos <- seq_along(private$res_list)
      }

      # what is our comparison criteria
      if(is.null(del_freq)) {
        del_freq_func <- function(x) {x > 0}
      } else {
        del_freq_func <- function(x) {x >= del_freq}
      }

      # Find regions if positions
      if (length(pos) > 0) {
        del <- map_lgl(private$res_list[pos], function(x){
          del_freq_func(x$freq[x$t_pos == t])
        })
        pos[which(del)]
      } else {
        integer(0L)
      }
    },

    # simulate selection
    simulate_selection = function(t, del_pos){

      if(length(del_pos) > 0) {

        # Remaining t for this time step
        t_right_pos <- which(private$res_list[[1]]$t_pos >= t)
        t_right <- private$res_list[[1]]$t[t_right_pos]
        t_forward <- t_right - t_right[1]

        # loop over the regions that need updating
        for(i in del_pos) {

          # s for our region
          s_pos <- match(
            names(private$res_list)[i],
            as.character(private$map_data$id_1)
          )
          s <- private$map_data$s[s_pos]

          # and our update positions
          private$res_list[[i]]$freq[t_right_pos] <-
            private$hrp2_mod$predict_f2(
            dat = list("s" = s),
            f1 = private$res_list[[i]]$freq[t_right_pos[1]],
            t = t_forward
            )

        }
      }
    },

    # find those that reach export freq
    simulate_importation = function(t, export_pos, import_freq, t_break = 1, import_gap = 1) {

      # imported regions
      import_pos <- unique(unlist(private$adj_mat[export_pos]))
      import_pos <- import_pos[import_pos != 0]

      # are we importing somewhere
      if(length(import_pos) > 0) {

        # Create a vector of where we are next going to be simulating
        next_del_pos <- c()

        # loop through imported regions
        for(j in import_pos) {

          # time plus one import_gap position
          tp1 <- which(private$res_list[[j]]$t_pos == t + as.integer(import_gap/t_break))

          # what is the frequency after the importation has occurred
          freq <- private$res_list[[j]]$freq[tp1]

          # if less than import then record and import
          if(freq < import_freq) {
            next_del_pos <- c(next_del_pos, j)
            private$res_list[[j]]$freq[tp1] <- import_freq
          }

        }

        # If all the regions being imported into are already higher than import_freq
        if(is.null(next_del_pos)) {
          next_del_pos <- integer(0L)
        }

      } else {

        next_del_pos <- integer(0L)

      }

      # return where we are simulating next time step
      return(next_del_pos)

    }

  )
)


# Floating less than or equals
#' @noRd
float_leq <- function(x, y, tolerance = 1e-9) {
  return(x < y | abs(x - y) <= tolerance)
}
