# --------------------------------------------------------
# extra funcs
# --------------------------------------------------------

t_p <- function(p2, p1, s) {
  (log(p2 / (1 - p2)) - log(p1 / (1 - p1))) / s
}

dat_p <- function(dat, i, what) {
  dat$res[[i]] %>%
    ggplot(aes_string("S.Times",what)) +
    geom_line()
}

# --------------------------------------------------------
# 1. Processing raw simulation outputs
# --------------------------------------------------------

for(i in seq_along(grep("true", list.files(file.path(here::here(), "analysis/data_derived/sims/"))))) {

  dat <- readRDS(file.path(
    here::here(),
    paste0(
      "analysis/data_derived/sims/true_grid_continuation_nmf_1_", i, ".rds"
    )
  ))

res <- lapply(seq_along(dat$res), function(ii){

  if(ii %% 100 == 0) message(ii)
  df <- dat$res[[ii]]

  # filter to after selection was implemented
  df2 <- df %>% filter(S.Times < 365*40.05 & S.Times > 365*38.05)
  res <- data.frame("PCR.All" = mean(df2$S.PCR.All, na.rm = TRUE),
                    "Micro.2.10" = mean(df2$S.Micro.210, na.rm = TRUE))

  # filter to after selection was implemented
  df <- df %>% filter(S.Times > 365*40.05)

  # 1. Calculate selection
  # ---------------------------------------

  # get our x and t
  x <- df$Percentange_Clin_Mono
  t <- df$S.Times

  # remove 0 and 1 values
  rmpos <- which(x < 1 & x > 0)

  x <- x[rmpos]
  t <- t[rmpos]

  # now remove the first and last years amount (most impacted by stochasticity)
  x <- head(x, -12)
  x <- tail(x, -12)
  t <- head(t, -12)
  t <- tail(t, -12)

  # at least 2 years still remaining of data for calculating selection
  if(length(x) > 23) {
    # selection coefficient / year
    s <- lm(log(x / (1-x)) ~ t)$coefficients[2] * (365)
  } else {
    s <- NA
  }

  res$s <- s
  return(res)

})


testna <- do.call(rbind, lapply(dat$pl, as.data.frame)) %>%
  as.data.frame %>%
  cbind(do.call(rbind, res))
testna <- testna %>%
  mutate(int = interaction(ft, microscopy.use, rdt.nonadherence, fitness, rdt.det, sep = "_"))
testna$rep <- i

saveRDS(testna,
        file.path(
          here::here(),
          paste0(
            "analysis/data_derived/sims/processed_tgc_nmf_1_", i, ".rds"
          )
        ))

}

# --------------------------------------------------------
# 2. Pull together and save as one object
# --------------------------------------------------------

testna <- lapply(
  grep("processed_tgc",
       list.files(file.path(here::here(), "analysis/data_derived/sims/"), full.names = TRUE),
       value = TRUE), readRDS
) %>% do.call(rbind, .)

# check what the mse is across stochastic reps
testna$int <- interaction(testna$int, testna$EIR, sep = "_")
saveRDS(testna, file.path(here::here("analysis/data_derived/model_s.rds")))
