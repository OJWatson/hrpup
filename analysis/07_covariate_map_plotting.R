# ---------------------------------------------------- #
# 0. Get the world map and our data  ----
# ---------------------------------------------------- #
library(tidyverse)

# WHO world map boundaries
who_shps <- readRDS("analysis/data_derived/who_shps.rds")

# covariate outputs
covars <- readRDS("analysis/data_derived/global_covariate_ranges.rds")

# Which vars are we plotting
vars <- grep("mean", names(covars), value = TRUE)
names(vars)  <- c("Microscopy \nPrevalence \n2-10 Years", "Effective \nTreatment \nSeeking", "Microscopy \nUse",
                  "RDT \nNonadherence")

# ---------------------------------------------------- #
# 1. Create and save our maps  ----
# ---------------------------------------------------- #
map_1 <- who_compliant_continuous_plot(covars, who_shps, region = "global", var = vars[1], title = names(vars[1]), limits = TRUE)
map_2 <- who_compliant_continuous_plot(covars, who_shps, region = "global", var = vars[2], title = names(vars[2]), limits = TRUE)
map_3 <- who_compliant_continuous_plot(covars, who_shps, region = "global", var = vars[3], title = names(vars[3]), limits = TRUE)
map_4 <- who_compliant_continuous_plot(covars, who_shps, region = "global", var = vars[4], title = names(vars[4]), limits = TRUE)

# Plot all four of the maps
save_figs("supp_micro210_map", add_global_scale(map_1) + theme(plot.margin = margin(20,20,20,20)), width = 15, height = 6, pdf_plot = FALSE)
save_figs("supp_ft_map", add_global_scale(map_2) + theme(plot.margin = margin(20,20,20,20)), width = 15, height = 6, pdf_plot = FALSE)
save_figs("supp_microuse_map", add_global_scale(map_3) + theme(plot.margin = margin(20,20,20,20)), width = 15, height = 6, pdf_plot = FALSE)
save_figs("supp_rdtnonadh_map", add_global_scale(map_4) + theme(plot.margin = margin(20,20,20,20)), width = 15, height = 6, pdf_plot = FALSE)


