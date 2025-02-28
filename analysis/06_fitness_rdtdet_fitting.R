library(tidyverse)
library(caret)
library(scam)
library(brnn)

# Grab shps
world_map_1 <- readRDS(here::here("analysis/data_derived/admin1_sf.rds"))

# Grab our model that we will use to simulate with in our MCMC
hrp2_model <- readRDS("analysis/data_derived/ensemble_selection_model.rds")

# Grab our covariate parameter ranges
covars <-  readRDS("analysis/data_derived/covars_for_fitness_rdtdet_fitting.rds")
covars <- covars %>% filter(iso3c %in% c("ERI", "ETH"))

# MAP world map boundaries
world <- left_join(covars, world_map_1) %>% sf::st_drop_geometry() %>%
  rename(admin_1 = name_1)

# -----------------------------------------------------#
# Our observed data that we are trying to fit to
# -----------------------------------------------------#

# 1 year after policy change
eri_rdt_start <- as.Date("2007-01-01")
eth_rdt_start <- as.Date("2005-01-01")

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5823352/ - all hrp3 deleted
ghindae <- data.frame(
  "x" = 21,
  "n" = 26,
  "time" = as.numeric((as.Date("2016-03-15") - eri_rdt_start)/365),
  "iso3c" = "ERI",
  "admin_1" = "Maekel",
  "hrp3" = 0
)

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5823352/ - all hrp3 deleted
massawa <- data.frame(
  "x" = 10,
  "n" = 24,
  "time" = as.numeric((as.Date("2016-03-15") - eri_rdt_start)/365),
  "iso3c" = "ERI",
  "admin_1" = "Semenawi Keih Bahri",
  "hrp3" = 0
)

# https://pubmed.ncbi.nlm.nih.gov/28889944/ - all hrp3 deleted
debub <- data.frame(
  "x" = 2,
  "n" = 9,
  "time" = as.numeric((as.Date("2013-05-15") - eri_rdt_start)/365),
  "iso3c" = "ERI",
  "admin_1" = "Debub",
  "hrp3" = 0
)

# https://pubmed.ncbi.nlm.nih.gov/28889944/ - 92% hrp3 deleted
gash_barka <- data.frame(
  "x" = 12,
  "n" = 135,
  "time" = as.numeric((as.Date("2013-05-15") - eri_rdt_start)/365),
  "iso3c" = "ERI",
  "admin_1" = "Gash Barka",
  "hrp3" = 0.08
)


# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0241807 - all hrp3 deleted
oromia <- data.frame(
  "x" = 50,
  "n" = 50,
  "time" = as.numeric((as.Date("2015-11-01") - eth_rdt_start)/365),
  "iso3c" = "ETH",
  "admin_1" = "Oromia",
  "hrp3" = 0
)

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8478644/ - 85% hrp3 deleted
amhara <- data.frame(
  "x" = round(0.115*1342),
  "n" = 1342,
  "time" = as.numeric((as.Date("2017-03-01") - eth_rdt_start)/365),
  "iso3c" = "ETH",
  "admin_1" = "Amhara",
  "hrp3" = 0.15
)

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8478644/- 85% hrp3 deleted
tigray <- data.frame(
  "x" = round(0.149*689),
  "n" = 689,
  "time" = as.numeric((as.Date("2017-03-01") - eth_rdt_start)/365),
  "iso3c" = "ETH",
  "admin_1" = "Tigray",
  "hrp3" = 0.15
)

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8478644/- 85% hrp3 deleted
gambela <- data.frame(
  "x" = round(0.011*683),
  "n" = 683,
  "time" = as.numeric((as.Date("2017-03-01") - eth_rdt_start)/365),
  "iso3c" = "ETH",
  "admin_1" = "Gambela",
  "hrp3" = 0.15
)

# SAMPLES AFTER RDTS WITHDRAWN
# https://pubmed.ncbi.nlm.nih.gov/34702923/

# - 79% hrp3 deleted
anseba_late <- data.frame(
  "x" = 29,
  "n" = 107,
  "time1" = as.numeric((as.Date("2016-09-15") - eri_rdt_start)/365),
  "time2" = as.numeric((as.Date("2019-09-15") - eri_rdt_start)/365),
  "iso3c" = "ERI",
  "admin_1" = "Anseba",
  "hrp3" = 0.21
)

# - 79% hrp3 deleted
gash_barka_late <- data.frame(
  "x" = 33,
  "n" = 512,
  "time1" = as.numeric((as.Date("2016-09-15") - eri_rdt_start)/365),
  "time2" = as.numeric((as.Date("2019-07-15") - eri_rdt_start)/365),
  "iso3c" = "ERI",
  "admin_1" = "Gash Barka",
  "hrp3" = 0.21
)

# - all hrp3 deleted
debub_late <- data.frame(
  "x" = 5,
  "n" = 96,
  "time1" = as.numeric((as.Date("2016-09-15") - eri_rdt_start)/365),
  "time2" = as.numeric((as.Date("2019-09-15") - eri_rdt_start)/365),
  "iso3c" = "ERI",
  "admin_1" = "Debub",
  "hrp3" = 0
)


# data that occurred before any RDT switch
data <- rbind(ghindae, massawa, debub, gash_barka,  amhara, tigray, gambela)

# data that occurred after any RDT switch
data2 <- rbind(anseba_late, gash_barka_late,  debub_late)

# prob of 1 and 0 will throw data
data$prob <- vapply(data$x/data$n, min, 0.999, FUN.VALUE = numeric(1))
data2$prob <- vapply(data2$x/data2$n, min, 0.999, FUN.VALUE = numeric(1))
data$prob <- vapply(data$prob, max, 0.001, FUN.VALUE = numeric(1))
data2$prob <- vapply(data2$prob, max, 0.001, FUN.VALUE = numeric(1))

# and format our params
model_params <- world %>%
  left_join(data) %>%
  rename(Micro.2.10 = Micro.2.10_mean,
         ft = ft_mean,
         microscopy.use = microscopy.use_mean,
         rdt.nonadherence = rdt.nonadherence_mean,
         fitness = fitness_mean,
         rdt.det = rdt.det_mean
  ) %>%
  select(Micro.2.10,
         Micro.2.10_low,
         Micro.2.10_high,
         ft,
         microscopy.use,
         rdt.nonadherence,
         fitness,
         rdt.det,
         hrp3,
         admin_1,
         year)

# just add anseba hrp3 here
model_params$hrp3[model_params$admin_1 == "Anseba" & model_params$year == 2020] <- anseba_late$hrp3

# and filter to just regions we are simulating and fill up
model_params <- model_params %>%
  filter(admin_1 %in% c(data$admin_1, data2$admin_1)) %>%
  group_by(iso3c, admin_1) %>%
  complete(year = seq(2005, max(year), 1)) %>%
  group_by(admin_1) %>%
  fill(Micro.2.10:hrp3, .direction = "up")

# put things in lists for the prior/posterior functions
misc <- list("hrp2_mod" = hrp2_model, "model_params" = model_params,
             "start_dates" = list("ERI" = eri_rdt_start, "ETH" = eth_rdt_start))
data_list <- list("data" = data, "data2" = data2)

# -----------------------------------------------------#


# -----------------------------------------------------#
# Likelihood function
# -----------------------------------------------------#

# define log-likelihood function
r_loglike <- function(params, data_list, misc) {

  # DATA 1

  # create data frame to pass to the prediction model using params
  newdat <- data_list$data

  # get our model parameters over time
  mod_dat <- misc$model_params
  mod_dat$fitness <- as.numeric(params["fitness"])
  mod_dat$rdt.det <- max(as.numeric(params["rdtdet"]) * mod_dat$hrp3, 0.01)
  mod_dat <- left_join(mod_dat, newdat, by = c("admin_1", "iso3c")) %>% filter(admin_1 %in% newdat$admin_1)

  # set up where we store the final f2s
  f2_s <- rep(0, length(newdat$x))

  q_get <- function(m, l, h, q){
    if(q == 0.5){
      return(m)
    } else if(q>0.5){
      m+((h-m)*((q-0.5)/0.475))
    } else {
      m-((m-l)*((0.5-q)/0.475))
    }
  }

  # per year estimate f2
  for(i in seq_along(newdat$admin_1)){

    ad <- newdat$admin_1[i]
    iso_ad <- newdat$iso3c[i]
    start <- lubridate::year(misc$start_dates[[iso_ad]])

    f2 <- 0.01
    t_it <- 1
    start <- lubridate::year(misc$start_dates[[iso_ad]])

    while(t_it < newdat$time[i]){

      t_run <- min(1, newdat$time[i] - t_it)
      f2 <- misc$hrp2_mod$predict_f2(
      mod_dat %>% filter(admin_1 == ad & year == start) %>%
        mutate(Micro.2.10 = as.numeric(q_get(m = Micro.2.10, l = Micro.2.10_low, h = Micro.2.10_high, q = params["q"]))),
      f1 = f2, t = t_run)

    start <- start + 1
    t_it <- t_it + 1

    }

    f2_s[i] <- f2

  }

  # calculate log-likelihood
  # calculate log-likelihood
  if("rho" %in% names(params)) {
    ret <- VGAM::dbetabinom(round(newdat$n * f2_s), newdat$n, prob = newdat$prob, rho = as.numeric(params["rho"]),  log = TRUE)
  } else {
    ret <- dbinom(round(newdat$n * f2_s), size = newdat$n, prob = newdat$prob, log = TRUE)
  }

  # DATA 2

  # create data frame to pass to the prediction model using params
  newdat2 <- data_list$data2

  # get our model parameters over time
  mod_dat <- misc$model_params
  mod_dat$fitness <- as.numeric(params["fitness"])
  mod_dat$rdt.det <- max(as.numeric(params["rdtdet"]) * mod_dat$hrp3, 0.01)
  mod_dat <- left_join(mod_dat, newdat2, by = c("admin_1", "iso3c")) %>% filter(admin_1 %in% newdat2$admin_1)

  # set up where we store the final f2s
  f2_2_s <- rep(0, length(newdat2$x))

  # per year estimate f2
  for(i in seq_along(newdat2$admin_1)){

    ad <- newdat2$admin_1[i]
    iso_ad <- newdat2$iso3c[i]
    start <- lubridate::year(misc$start_dates[[iso_ad]])

    f2 <- 0.01
    t_it <- 1
    t_left <- newdat2$time2[i]
    start <- lubridate::year(misc$start_dates[[iso_ad]])

    # Simulate before RDT switch
    while(t_it < newdat2$time1[i]){

      t_run <- min(1, newdat2$time1[i] - t_it)
      f2 <- misc$hrp2_mod$predict_f2(
        mod_dat %>% filter(admin_1 == ad & year == start), f1 = f2, t = t_run)

      start <- start + t_run
      t_it <- t_it + t_run

    }

    # Simulate post RDT switch
    # Do one small sim to finish the half year
    t_run <- ceiling(t_it) - t_it
    f2 <- misc$hrp2_mod$predict_f2(
      mod_dat %>%
        filter(admin_1 == ad & year == floor(start)) %>%
        mutate(microscopy.use = 1), f1 = f2, t = t_run)

    start <- start + t_run
    t_it <- t_it + t_run

    # Simulate remaining years
    while(t_it < newdat2$time2[i]){

      t_run <- min(1, newdat2$time2[i] - t_it)
      f2 <- misc$hrp2_mod$predict_f2(
        mod_dat %>%
          filter(admin_1 == ad & year == start) %>%
          mutate(microscopy.use = 1), f1 = f2, t = t_run)

      start <- start + 1
      t_it <- t_it + 1

    }


    f2_2_s[i] <- f2

  }

  # calculate log-likelihood
  if("rho" %in% names(params)) {
    ret_2 <- VGAM::dbetabinom(round(newdat2$n * f2_2_s), newdat2$n, prob = newdat2$prob, rho = as.numeric(params["rho"]),  log = TRUE)
  } else {
    ret_2 <- dbinom(round(newdat2$n * f2_2_s), size = newdat2$n, prob = newdat2$prob, log = TRUE)
  }
  # return sum of the ll
  return(sum(ret) + sum(ret_2))
}

# define log-prior function
r_logprior <- function(params, misc) {

  # extract parameter values
  fitness <- as.numeric(params["fitness"])
  rdtdet <- as.numeric(params["rdtdet"])

  # calculate log-prior
  ret <- dunif(fitness, min = 0, max = 0.99, log = TRUE) +
    #dnorm(fitness, mean = 0.95, sd = 0.02, log = TRUE) +
    dbeta(rdtdet, shape1 = 12+1, shape2 = 26-12+1, log = TRUE) # informed by Ethiopia Feleke study
    #dunif(rdtdet, min = 0, max = 0.99, log = TRUE)

  # return
  return(ret)
}

# scan suitable range (range identified after first doing broad search)
params <- expand.grid(rdtdet = seq(0.2,0.5,0.01),
                      fitness = seq(0.93,0.99,0.001),
                      rho = c(0.0005, 0.001, 0.005),
                      q = 0.5)

# set up parallel cluster here
n.cores <- 14
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)
future::plan(future::cluster, workers = my.cluster)

params <-
  params %>%
  mutate(ll = furrr::future_map_dbl(
    seq_len(nrow(params)), ~ r_loglike(
      params[.x,],data_list, misc
    ), .progress = TRUE, .options = furrr::furrr_options(packages = c("brnn", "earth","scam")))
  ) %>%
  mutate(lp = map_dbl(
    seq_len(nrow(params)), ~ r_logprior(
      params[.x,], misc
    ), .progress = TRUE)
  ) %>%
  mutate(posterior = ll + lp)

parallel::stopCluster(my.cluster)

# Does this seem correct. Quick Ll plot
ll_plot <- params %>%
  ggplot(aes(x = fitness, y = rdtdet, fill = -posterior)) + geom_tile() +
scale_fill_gradientn(name = "Negative Log \nLikelihood",
                     colors = RColorBrewer::brewer.pal(9,"YlOrBr")[c(3, 5,7,9)],
                     trans = "log",  # Log transformation
                     values = c(0,0.05,0.1,0.15,0.4,1),
                     breaks = c(90, 180, 360),
                     labels = c(90, 180, 360)) +
                     #limits = c(17.5, 150)) +
  ggpubr::theme_pubclean() +
  xlab("Comparative Fitness of Deleted Parasites") +
  ylab("Probability of hrp2 deleted parasite yielding +ve RDT") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0))


params <- params %>% filter(fitness < 0.98 & fitness > 0.94)
params <- params %>% filter(rdtdet < 0.45)
params <- sav2
params$posterior <- params$posterior * (params$rdtdet^0.03)

# Now draw from our posterior to get a data set range
draws <- sample(seq_len(nrow(params)), 1000, TRUE, exp(params$posterior))
dat <- params[draws,]
saveRDS(dat, "analysis/data_derived/mcmc_model_fitting.rds")

# values for manuscript
dat$rdtdet %>% quantile(c(0.025,0.5,0.975))
dat$fitness %>% quantile(c(0.025,0.5,0.975))

# Create a plot of the posterior
ll_plot_alt <- dat %>%
  ggplot(aes(x = fitness, y = rdtdet)) +
  geom_density2d_filled(contour_var = "ndensity", bins = 8) +
  ggpubr::theme_pubclean(base_family = "Helvetica") +
  xlab("Comparative Fitness of Deleted Parasites") +
  ylab("Probability of hrp2 deleted parasite yielding +ve RDT") +
  scale_y_continuous(expand = c(0,0), labels = scales::percent_format()) +
  scale_x_continuous(expand = c(0,0), labels = scales::percent_format()) +
  scale_fill_manual(name = "Posterior Density \nScaled to 1",
                       values = RColorBrewer::brewer.pal(9,"YlOrBr")[2:9])  +
  theme(legend.position = "right")
ll_plot_alt

save_figs(name = "fitness_hrp3_ll", fig = ll_plot_alt, width = 8, height = 5)

# get the posterior draws from the exhaustive search
pfd <- readRDS("analysis/data_derived/mcmc_model_fitting.rds")
rdt_range <- as.numeric(round(quantile(pfd$rdtdet, c(0.025,0.5,0.975)), 2))
fitness_range <- as.numeric(round(quantile(pfd$fitness, c(0.025,0.5,0.975)), 3))

# add in the fitness and rdt.det data
covars <-  readRDS("analysis/data_derived/global_covariate_ranges.rds")
# Add these in and multiply rdt det by results of independence modelling earlier (i.e. 1 - 0.617)
covars <- covars %>%
  mutate(fitness_mean = fitness_range[2], fitness_low = fitness_range[1], fitness_high = fitness_range[3]) %>%
  mutate(rdt.det_mean = rdt_range[2]*0.383, rdt.det_low = rdt_range[1]*0.383, rdt.det_high = rdt_range[3]*0.383)

saveRDS(covars, "analysis/data_derived/global_covariate_ranges.rds")

