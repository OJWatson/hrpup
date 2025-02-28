
# to create the full data we can use the map object
hrp2_map <- readRDS("analysis/data_derived/R6_map.rds")

# this is where all the map data is for each scenario
hrp2_map$.__enclos_env__$private$map_data[[1]][1:5, ]

# with the corresponding scenario table at:
hrp2_map$.__enclos_env__$private$scenarios[1,]

# The map data id_1 column matches the map id_1 column for linking
hrp2_map$.__enclos_env__$private$map[1,]

# All the map data as one data frame can then be as:
full_dat <- list()
for(i in seq_along(hrp2_map$.__enclos_env__$private$map_data)){
  full_dat[[i]] <- hrp2_map$.__enclos_env__$private$map_data[[i]]
  full_dat[[i]] <- cbind(full_dat[[i]], hrp2_map$.__enclos_env__$private$scenarios[i,])
}

full_df <- do.call(rbind, full_dat)

# grab map with all map detail
covars <-  readRDS("analysis/data_derived/global_covariate_ranges.rds")

# MAP world map to use
world_map <- readRDS(here::here("analysis/data_derived/admin1_sf.rds"))

# write to file
full_df2 <- left_join(full_df[,-2], world_map %>% sf::st_drop_geometry() %>%
                        select(id_1, name_1, iso), by = c("id_1"))
write.csv(full_df2, "analysis/data_out/full_results.csv", row.names = FALSE)
zip::zip("analysis/data_out/full_results.csv.zip", "analysis/data_out/full_results.csv")
file.remove("analysis/data_out/full_results.csv")
