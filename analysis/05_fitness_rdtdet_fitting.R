# Grab our model that we will use to simulate with in our MCMC
hrp2_model <- readRDS("analysis/data_derived/ensemble_selection_model.rds")

# Grab our covariate parameter ranges
covars <-  readRDS("analysis/data_derived/global_covariate_ranges.rds")

# MAP world map boundaries
available_admin <- malariaAtlas::listShp(printed = FALSE, admin_level = "admin0")
world_map_1 <- malariaAtlas::getShp(ISO = available_admin$iso, admin_level = c("admin1")) %>% sf::st_as_sf()
world <- left_join(world_map_1, covars) %>% sf::st_drop_geometry() %>%
  rename(admin_1 = name_1)

# -----------------------------------------------------#
# Our observed data that we are trying to fit to
# -----------------------------------------------------#

# 1 year after policy change
eri_rdt_start <- as.Date("2007-01-01")
eth_rdt_start <- as.Date("2005-01-01")

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5823352/
ghindae <- data.frame(
  "x" = 21,
  "n" = 26,
  "time" = as.numeric((as.Date("2016-03-15") - eri_rdt_start)/365),
  "iso3c" = "ERI",
  "admin_1" = "Maekel"
)

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5823352/
massawa <- data.frame(
  "x" = 10,
  "n" = 24,
  "time" = as.numeric((as.Date("2016-03-15") - eri_rdt_start)/365),
  "iso3c" = "ERI",
  "admin_1" = "Semenawi Keih Bahri"
)

# https://pubmed.ncbi.nlm.nih.gov/28889944/
debub <- data.frame(
  "x" = 2,
  "n" = 9,
  "time" = as.numeric((as.Date("2013-05-15") - eri_rdt_start)/365),
  "iso3c" = "ERI",
  "admin_1" = "Debub"
)

# https://pubmed.ncbi.nlm.nih.gov/28889944/
gash_barka <- data.frame(
  "x" = 12,
  "n" = 135,
  "time" = as.numeric((as.Date("2013-05-15") - eri_rdt_start)/365),
  "iso3c" = "ERI",
  "admin_1" = "Gash Barka"
)


# https://pubmed.ncbi.nlm.nih.gov/28889944/
oromia <- data.frame(
  "x" = 50,
  "n" = 50,
  "time" = as.numeric((as.Date("2015-11-01") - eth_rdt_start)/365),
  "iso3c" = "ETH",
  "admin_1" = "Oromia"
)

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8478644/
amhara <- data.frame(
  "x" = round(0.115*1342),
  "n" = 1342,
  "time" = as.numeric((as.Date("2017-03-01") - eth_rdt_start)/365),
  "iso3c" = "ETH",
  "admin_1" = "Amhara"
)

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8478644/
tigray <- data.frame(
  "x" = round(0.149*689),
  "n" = 689,
  "time" = as.numeric((as.Date("2017-03-01") - eth_rdt_start)/365),
  "iso3c" = "ETH",
  "admin_1" = "Tigray"
)

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8478644/
gambela <- data.frame(
  "x" = round(0.011*683),
  "n" = 683,
  "time" = as.numeric((as.Date("2017-03-01") - eth_rdt_start)/365),
  "iso3c" = "ETH",
  "admin_1" = "Gambela"
)

data <- rbind(ghindae, massawa, debub, gash_barka,  amhara, tigray, gambela)
data <- left_join(data, world)
data <- data %>%
  rename(Micro.2.10 = Micro.2.10_mean,
         ft = ft_mean,
         microscopy.use = microscopy.use_mean,
         rdt.nonadherence = rdt.nonadherence_mean,
         fitness = fitness_mean,
         rdt.det = rdt.det_mean
         ) %>%
  select(x, n, time,
         Micro.2.10,
         ft,
         microscopy.use,
         rdt.nonadherence,
         fitness,
         rdt.det)

# prob of 1 and 0 will throw data
data$prob <- vapply(data$x/data$n, min, 0.999, FUN.VALUE = numeric(1))
data$prob <- vapply(data$prob, max, 0.001, FUN.VALUE = numeric(1))


# -----------------------------------------------------#


# -----------------------------------------------------#
# Likelihood function
# -----------------------------------------------------#

# define log-likelihood function
r_loglike <- function(params, data, misc) {

  # create data frame to pass to the prediction model using params
  newdat <- data$data
  newdat$fitness <- as.numeric(params["fitness"])
  newdat$rdt.det <- as.numeric(params["rdtdet"])

  # estimate new f2 using model
  f2 <- misc$hrp2_mod$predict_f2(newdat, f1 = 0.001, t = newdat$time)

  # calculate log-likelihood
  ret <- dbinom(round(newdat$n * f2), size = newdat$n, prob = newdat$prob, log = TRUE)

  # return sum of the ll
  return(sum(ret))
}

# define log-prior function
r_logprior <- function(params, misc) {

  # extract parameter values
  fitness <- as.numeric(params["fitness"])
  rdtdet <- as.numeric(params["rdtdet"])

  # calculate log-prior
  ret <- dnorm(fitness, mean = 0.95, sd = 0.02, log = TRUE) +
    dunif(rdtdet, min = 0.1, max = 0.4, log = TRUE)

  # return
  return(ret)
}


# scan suitable range
params <- expand.grid(rdtdet = seq(0.175,0.33,0.0025),
                      fitness = seq(0.93,0.97,0.0005))

params <-
  params %>%
  mutate(ll = map_dbl(
    seq_len(nrow(params)), ~ r_loglike(
      params[.x,],data_list, misc
    ), .progress = TRUE)
  ) %>%
  mutate(lp = map_dbl(
    seq_len(nrow(params)), ~ r_logprior(
      params[.x,], misc
    ), .progress = TRUE)
  ) %>%
  mutate(posterior = ll + lp)

ll_plot <- params %>%
  ggplot(aes(x = fitness, y = rdtdet, fill = -posterior)) + geom_tile() +
scale_fill_gradientn(name = "Negative Log \nLikelihood",
                     colors = c("blue", "green", "yellow", "red"),
                     values = c(0,0.2,0.4,0.6,1),
                     trans = "log", # Log transformation
                     breaks = c(20, 35, 60, 120),
                     labels = c(20, 35, 60, 120),
                     limits = c(17.5, 125)) +
  theme_minimal() +
  xlab("Relative Fitness") +
  ylab("Probability of hrp2 deleted parasite \nwith intact hrp3 and yielding +ve RDT")
ll_plot

draws <- sample(seq_len(nrow(params)), 1000, TRUE, exp(params$posterior))
params[draws,] %>%
  ggplot(aes(x = fitness, y = rdtdet)) +
  geom_density2d_filled(contour_var = "ndensity")


saveRDS(params[draws,], "analysis/data_derived/mcmc_model_fitting.rds")

# get the posterior draws from the exhaustive search
pfd <- readRDS("analysis/data_derived/mcmc_model_fitting.rds")
rdt_range <- as.numeric(round(quantile(pfd$rdtdet, c(0.025,0.5,0.975)), 2))
fitness_range <- as.numeric(round(quantile(pfd$fitness, c(0.025,0.5,0.975)), 2))

# add in the fitness and rdt.det data
covars <- covars %>%
  mutate(fitness_mean = fitness_range[2], fitness_low = fitness_range[1], fitness_high = fitness_range[3]) %>%
  mutate(rdt.det_mean = rdt_range[2], rdt.det_low = rdt_range[1], rdt.det_high = rdt_range[3])

saveRDS(covars, "analysis/data_derived/global_covariate_ranges.rds")
