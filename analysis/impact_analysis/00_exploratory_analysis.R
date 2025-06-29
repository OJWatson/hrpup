library(tidyverse)
i <- 1

dat <- readRDS(file.path(
  here::here(),
  paste0(
    "analysis/data_derived/sims/true_grid_continuation_nmf_1_", i, ".rds"
  )
))


# ii <- 27 # high EIR
ii <- 2
df <- dat$res[[ii]]

# filter to after selection was implemented
df <- df %>% filter(S.Times > 365*40.05)



df %>% ggplot(aes((S.Times/365)-40, Percentange_Clin_Mono)) +
  geom_line() + geom_hline(yintercept = 0.05) +
  xlab("Years After")



library(ggplot2)
library(mgcv)

# Define a custom smoothing method for monotonic splines with adjustable smoothness
monotonic_smooth <- function(formula, data, weights = NULL, k = 10, sp = NULL, ...) {
  # Fit a monotonic GAM model with specified smoothness
  gam_model <- gam(formula, data = data,
                   family = gaussian(),
                   method = "REML",  # Optimal smoothing penalty
                   constraints = list(slope = "increasing"),
                   select = TRUE,   # Allow automatic penalty adjustment
                   sp = sp,         # Optional: User-specified smoothing parameter
                   ...)
  return(gam_model)  # Return the model
}

# Plot with adjustable smoothness
df %>%
  mutate(Years = (S.Times / 365) - 40) %>%
  filter(Years > 0) %>%
  ggplot(aes(Years, S.Incidence*100)) +
  geom_line(alpha = 0.5) +
  geom_smooth(
    method = monotonic_smooth,
    formula = y ~ s(x, bs = "cr", k = 6), # k controls smoothness (higher = smoother)
    se = FALSE,
    colour = "darkred",
    linetype = "dashed",
    lwd = 1
  ) +
  geomtextpath::geom_texthline(yintercept = 0.043, color = "darkblue", linetype = "dashed", family = "Helvetica",
                               label = "Assumed Clinical Incidence if RDT switched at 5% theshold", hjust = 0.9, lwd = 1) +
  xlab("Years After 5% Threshold is Passed") +
  ylab("Daily Clinical Incidence/1000") +
  theme_bw()

#
