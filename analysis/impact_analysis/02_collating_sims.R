library(tidyverse)
library(brnn)
## -------------------------- ##
## 1. Create a database from malariasimulation for equivalent simulations for incidence outputs --------
## -------------------------- ##

library(malariasimulation)

get_coeffics <- function(EIR, ft) {

# Daily simulation timesteps for 5 years
year <- 365
sim_length <- 5 * year

# With a population size of 1000
human_population <- 10000


min_ages <- c(0, 5*365)
max_upper_bound <- max(100 * 365, min_ages[length(min_ages)] + 365)
ages <- c(min_ages, max_upper_bound)
p <- malariasimulation::get_parameters(overrides = list(human_population = human_population))
p <- malariasimulation::set_epi_outputs(p, age_group = ages, incidence = ages,
                                        clinical_incidence = ages, severe_incidence = ages)

# Update parameter set with chosen drug-specific parameters (AL and DHA/PQP)
drug_params <- set_drugs(p, list(AL_params, DHA_PQP_params))

# Choose initial EIR to 10
starting_EIR <- EIR

# Set treatment program for AL (drug index = 1)
treatment_params <- set_clinical_treatment(
  parameters = drug_params,
  drug = 1,
  timesteps =  c(0), # Treatment coverage changes on day 300 and day 600
  coverages =  c(ft)) # The initial treatment coverage (0%) is the default
# and does not need to be set

# Use set_equilibrium to update the parameter set for a given initial EIR
treatment_params <- set_equilibrium(treatment_params, starting_EIR)

# Run simulation:
output <- run_simulation(sim_length, treatment_params)


output <- output %>%
  mutate(n_inc_clinical_0_36499 = n_inc_clinical_0_1824 + n_inc_clinical_1825_36499,
         n_inc_severe_0_36499 = n_inc_severe_0_1824 + n_inc_severe_1825_36499,
         n_age_0_36499 = n_age_0_1824 + n_age_1825_36499)

# get epi outputs
epi <- output |>
  postie::get_rates(scaler = 0.215)

rates <- epi |>
  dplyr::summarise(
    clinical = weighted.mean(clinical, person_days) * 365,
    severe = weighted.mean(severe, person_days) * 365,
    mortality = weighted.mean(mortality, person_days) * 365,
    .by = c("year","month", "age_lower", "age_upper")
  )

return(rates %>% tail(60) %>%
  group_by(age_lower, age_upper) %>%
  summarise(across(clinical:mortality,.fns =mean)) %>%
  mutate(across(starts_with("age"), ceiling)) %>%
  mutate(across(starts_with("age"), as.integer)) %>%
  mutate(age = paste0(age_lower, "_", age_upper)) %>%
  filter(age %in% c("0_5", "0_100")) %>%
  mutate(EIR = EIR, ft = ft))

}

eir_ft_grid <- expand.grid(EIR = magenta:::age_brackets(max_age = 200, num_age_brackets = 20)[-1], ft = seq(0.05,0.85, 0.04))

future::plan(future::multisession, workers = 14)
inc_grid <- eir_ft_grid %>% mutate(n = row_number()) %>%
  split(.$n) %>%
  furrr::future_map(function(x) {get_coeffics(x$EIR, x$ft)}, .progress = TRUE) %>%
  do.call(rbind, .)
dir.create("analysis/impact_analysis/data-derived/")
saveRDS(inc_grid, "analysis/impact_analysis/data-derived/inc_grid.rds")
inc_grid <- readRDS("analysis/impact_analysis/data-derived/inc_grid.rds")

## ALT Training

install.packages("torch")
install.packages("luz")

library(torch)
library(tabnet)
library(tidyverse)
library(tidymodels)
library(finetune) # to use tuning functions from the new finetune package
library(vip) # to plot feature importances

set.seed(777)
torch_manual_seed(777)


## -------------------------- ##
## 2. train brnn to convert from clinical incide, EIR and ft --------
## -------------------------- ##

# set seed again
set.seed(123)
test <- inc_grid %>% na.omit() %>% ungroup() %>% select(clinical, severe, mortality, age, EIR, ft)
train_indices <- sample(nrow(test), nrow(test) * 0.75)
train <- test[train_indices, ]
test <- test[-train_indices, ]

# train models
train_brnn_model <- function(train, agevar = "0_100", output = "severe", neurons = 4)  {
  train_control <- caret::trainControl(method="cv", number=20)
  model <- caret::train(x = train %>% select(-severe, -mortality) %>% filter(age == agevar) %>% select(-age) ,
                        y = (train %>% filter(age == agevar))[[output]], method="brnn",
                        trControl = train_control, tuneGrid = expand.grid("neurons" = neurons))
  return(model$finalModel)
}

model_100_severe <- train_brnn_model(train)
model_5_severe <- train_brnn_model(train, agevar = "0_5", neurons = 2)
model_100_mortality <- train_brnn_model(train, output = "mortality")
model_5_mortality <- train_brnn_model(train, agevar = "0_5",  output = "mortality", neurons = 2)

# test models
compare_model <- function(test, model, agevar = "0_100", output = "severe")  {

  plot(predict(model, test %>% select(-severe, -mortality) %>% filter(age == agevar) %>% select(-age)),
               (test %>% filter(age == agevar))[[output]])
}

compare_model(test, model_100_severe)
compare_model(test, model_5_severe, agevar = "0_5")
compare_model(test, model_100_mortality, output = "mortality")
compare_model(test, model_5_mortality, agevar = "0_5",  output = "mortality")

# put models together for use later
models <- list(model_100_severe, model_5_severe, model_100_mortality, model_5_mortality)
names(models) <- c("model_100_severe", "model_5_severe", "model_100_mortality", "model_5_mortality")
saveRDS(models, "analysis/impact_analysis/data-derived/inc_models.rds")
models <- readRDS("analysis/impact_analysis/data-derived/inc_models.rds")
# --------------------------------------------------------
# 3. Processing raw simulation outputs
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

    # filter to after selection was implemented
    df <- df %>% filter(S.Times > 365*40) %>% mutate(S.Times = round((S.Times/365)-40, 3))

    # create main outputs
    df_out <- df %>% select(
      t = S.Times,
      freq = Percentange_Clin_Mono,
      micro_2_10 = S.Micro.210,
      pcr = S.PCR.All,
      clinical_05 = S.Incidence.05,
      clinical = S.Incidence
    ) %>%
      mutate(across(starts_with("clinical"), function(x){x*365}))

    # pivot to wider
    res <- df_out %>%
      mutate(severe_05 = predict(
        models$model_5_severe, df_out %>% select(clinical) %>% mutate(EIR = dat$pl[[ii]]$EIR*365, ft = dat$pl[[ii]]$ft)
        )) %>%
      mutate(severe = predict(
        models$model_100_severe, df_out %>% select(clinical) %>% mutate(EIR = dat$pl[[ii]]$EIR*365, ft = dat$pl[[ii]]$ft)
      )) %>%
      mutate(mortality_05 = predict(
        models$model_5_mortality, df_out %>% select(clinical) %>% mutate(EIR = dat$pl[[ii]]$EIR*365, ft = dat$pl[[ii]]$ft)
      )) %>%
      mutate(mortality_100 = predict(
        models$model_100_mortality, df_out %>% select(clinical) %>% mutate(EIR = dat$pl[[ii]]$EIR*365, ft = dat$pl[[ii]]$ft)
      ))

    return(res)

  })

  # other params that will be used for training
  other_args <- readRDS(file.path(
    here::here(),
    paste0(
      "analysis/data_derived/sims/processed_tgc_nmf_1_", i, ".rds"
    )
  ))

  # bind these
  for(j in seq_along(dat$res)) {
    res[[j]] <- cbind(other_args[j,], res[[j]])
  }

  # bind the parameters
  testna <- do.call(rbind, res)
  testna$rep <- i

  # bring together and just grab the essential info
  full <- testna %>%
    select(EIR, rep, s, Micro.2.10, ft, microscopy.use, rdt.nonadherence, fitness, rdt.det, t,
           starts_with(c("freq","micro_2_10","pcr","clinical","severe","mortality")))

  dir.create("analysis/impact_analysis/data-derived/sims", showWarnings = FALSE)
  saveRDS(full,
          file.path(
            here::here(),
            paste0(
              "analysis/impact_analysis/data-derived/sims/processed_impact_nmf_1_", i, ".rds"
            )
          ))

}

# --------------------------------------------------------
# 4. Pull together and save as one object
# --------------------------------------------------------

testna <- lapply(
  grep("processed_impact_nmf",
       list.files(file.path(here::here(), "analysis/impact_analysis/data-derived/sims/"), full.names = TRUE),
       value = TRUE), readRDS
) %>% do.call(rbind, .)

testna$PCR.All <- NULL

# save the full dataset out for downstream emulator training
saveRDS(testna, file.path(here::here("analysis/impact_analysis/data-derived/model_sims.rds")))

# --------------------------------------------------------
# 5. Plot an example
# --------------------------------------------------------

semi_join(
  testna, as.data.frame(dat$pl[[77]]) %>%
    select(EIR, ft, microscopy.use, rdt.nonadherence, fitness, rdt.det),
  by = c("EIR", "ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det")
) %>%
  filter(rep == 1) %>%
  pivot_longer(freq:mortality_100) %>%
  mutate(name = stringr::str_to_title(name)) %>%
  mutate(clr = FALSE) %>%
  mutate(clr = replace(clr, name == "Freq", TRUE)) %>%
  mutate(name = replace(name, name == "Freq", "False Negative RDTs")) %>%
  ggplot(aes(t, value, color = clr)) + geom_line() + lemon::facet_rep_wrap(~name, scales = "free_y") +
  theme_minimal() + ylab("Malaria Indicator") + xlab("Year") + theme(axis.line = element_line()) +
  scale_color_discrete(name = "Monotonic Simulation Output:") +
  theme(legend.position = "top")



semi_join(
  testna, huh[1,3:7, drop = FALSE] %>% as.data.frame() %>% setNames(c("ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det")),
  by = c("ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det")
)

renorm_list <- list(
ft = c(0.01, 1),
microscopy.use = c(0, 1),
rdt.nonadherence = c(0, 1),
fitness = c(0.9, 1.0),
rdt.det = c(0, 1)
)
ch <- huh[1,3:7, drop = FALSE] %>% as.data.frame() %>% setNames(c("ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det"))

# Normalise columns
for (col in names(renorm_list)) {
  range <- renorm_list[[col]]
  ch[[col]] <-  (ch[[col]] * (range[2] - range[1])) + range[1]
}


semi_join(
  testna, ch,
  by = c("ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det")
)
