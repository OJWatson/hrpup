library(tidyverse)

# ---------------------------------------------------- #
# 1. Get access covariate dataset from MAP rasters ----
# ---------------------------------------------------- #

# rasters
rasters <- c("Walking-only travel time to healthcare map without access to motorized transport",
             "Global friction surface enumerating land-based travel speed with access to motorized transport for a nominal year 2019",
             "Global friction surface enumerating land-based travel walking-only speed without access to motorized transport for a nominal year 2019",
             "Global travel time to healthcare map with access to motorized transport")

travel_ft_map <- function(iso3c,
                          raster = "Global friction surface enumerating land-based travel speed with access to motorized transport for a nominal year 2019",
                          year = 2020) {

  # Get Admin 1 polygons from geoboundaries API
  sf_poly <- map %>% filter(iso == iso3c)

  # get countries data
  data <- cart::pull_cart(iso3c = iso3c, year = year,
                          vector_rasters = list(),
                          prevalence_rasters = list(),
                          spatial_limits_rasters = list(raster)
  )

  # extract
  extracted_data <- cart::unpack_cart(sf_poly, data)

  # population summarise by iso3c
  summary_data <- extracted_data %>%
    dplyr::select(id_1, pop, matches("layer")) %>%
    tidyr::unnest(cols = c(pop, matches("layer"))) %>%
    tidyr::drop_na() %>%
    dplyr::group_by(id_1) %>%
    dplyr::summarise(across(starts_with("layer"), function(x){weighted.mean(x, pop, na.rm = TRUE)}),
pop_total = sum(pop))


  # joing back together
  fc <- dplyr::left_join(sf_poly, summary_data, by = "id_1")

  return(fc)

}

scenario_maps <- readRDS("analysis/data_derived/scenario_maps_full.rds")
map <- scenario_maps$map

iso3cs <- unique(map$iso)
iso3cs <- iso3cs[iso3cs != "UMI"]
iso3cs <- iso3cs[iso3cs != "MEX"]
iso3cs <- iso3cs[iso3cs != "ATF"]
travel_list <- vector("list", length(iso3cs))
for(i in seq_along(travel_list)) {
  message(i)
  travel_list[[i]] <- travel_ft_map(iso3cs[i], rasters)
}
# Mexico did exist (error was misleading). add to end
travel_list[[length(travel_list)+1]] <- travel_ft_map("MEX", rasters)
saveRDS(travel_list, "analysis/data_derived/travel_list.rds")
all_travel <- do.call(rbind, travel_list)

# ---------------------------------------------------- #
# 2. Get the DHS ft values into our MAP map ----
# ---------------------------------------------------- #

# read in dhs data
shp <- sf::read_sf("analysis/data_raw/data_dhs/shps/sdr_subnational_data.shp")

# MAP world map
covars <-  readRDS("analysis/data_derived/global_covariate_ranges.rds")
#world_map <- malariaAtlas::getShp(ISO = na.omit(unique(covars$iso3c)), admin_level = c("admin1")) %>% sf::st_as_sf()
world_map <- readRDS("analysis/data_derived/admin1_sf.rds")
shp$iso3c <- countrycode::countrycode(shp$CNTRYNAMEE, origin = "country.name.en","iso3c")

# start matching it to our MAP map
iso3cs <- unique(shp$iso3c)
iso3cs <- iso3cs[iso3cs %in% world_map$iso]
map <- world_map
map$ft <- NA
ids <- map$id_1[map$iso %in% iso3cs]

sf::sf_use_s2(FALSE)

# loop through and work out ft where we have it
for(i in seq_along(ids)){

  id_i <- ids[i]
  poly <- map %>% filter(id_1 == id_i)
  centre <- sf::st_centroid(poly)
  shp_polys <- shp %>% filter(iso3c == poly$iso)
  nn <- sf::st_nearest_feature(centre, shp_polys)
  ft_i <- shp_polys[nn,]$MLFEVTCADV
  map$ft[map$id_1 == id_i] <- ft_i

}

# set NA to missings
map$ft[which(map$ft > 100)] <- NA
map$ft[which(map$ft < 0)] <- NA

# ---------------------------------------------------- #
# 3. Combine these together ----
# ---------------------------------------------------- #

new_df <- map %>% sf::st_drop_geometry() %>% select(id_1, ft) %>%
  left_join(all_travel %>% sf::st_drop_geometry() %>% select(id_1, pop_total, matches("layer")))

# Grab our covariate parameter ranges
covars <-  readRDS("analysis/data_derived/global_covariate_ranges.rds")
new_df <- covars %>% select(id_1, iso3c) %>% left_join(new_df)

# bring in a covariate dataset at the iso3c level from Economist for imputation help
preds <- readRDS("analysis/data_raw/economist_covariates.RDS")
preds <- preds %>% filter(date == 18264)
pred_df <- preds %>% select(iso3c, hospital_beds_per_thousand, life_expectancy, vdem_freedom_of_expression_score:polity_democracy_score, wdi_obs_lag:wdi_urban_pop_1m_cities_pct, percent_land_area_in_tropics:total_deaths_latest_per_100k,gdpppc_ppp_imf)

# impute missing predictors
full_df <- new_df %>% left_join(pred_df) %>% ungroup()
mic <- mice::mice(full_df %>% select(-id_1, -iso3c, -ft, -pop_total))
comp <- mice::complete(mic)
full_df[,-which(names(full_df) %in% c("id_1","ft", "iso3c", "pop_total"))] <- comp

# build brnn for this
library(caret)
set.seed(123)
train_indices <- sample(nrow(full_df), nrow(full_df) * 0.75)
train <- full_df[train_indices, ]
test <- full_df[-train_indices, ]

# model <- caret::train(x = train %>% na.omit %>% select(-ft, -iso3c, -id_1), y = train %>% na.omit %>% pull(ft),
#                       method="brnn",
#                       trControl = train_control, tuneGrid = expand.grid("neurons" =4))
#
# # check this looks okay
# plot(test %>% na.omit %>% pull(ft),
#      predict(model$finalModel, test %>% na.omit %>% select(-ft, -iso3c, -id_1)))

library(xgboost)
xgb_params <- list(objective = "reg:logistic",
                   eta = 0.05,
                   max_depth = 12, subsample = 0.85, colsample_bytree = 0.85)#, # Use mean squared error as objective function
                   # eta = best_params$eta, # Learning rate
                   # max_depth = best_params$max_depth, # Maximum tree depth
                   # min_child_weight = 1, # Minimum sum of instance weight needed in a child
                   # subsample = best_params$subsample, # Subsample ratio of the training instances
                   # colsample_bytree = best_params$colsample_bytree, # Subsample ratio of columns when constructing each tree
                   # monotone_constraints = c(-1, 1, -1, -1, 1, -1)) # Impose monotonic constraints for known functional relationships

# Train xgboost model with cross-validation
xgb_cv <- xgb.cv(params = xgb_params,
                 data = as.matrix(train %>% na.omit %>% select(-ft, -iso3c, -id_1)),
                 label = train %>% na.omit %>% mutate(ft = ft/100) %>% pull(ft),
                 nfold = 20, verbose = 0, na.rm = TRUE,
                 nrounds = 200)

# Extract best iteration
best_iter <- which.min(xgb_cv$evaluation_log$test_rmse_mean)

# Train final model using best iteration
xgb_model <- xgboost::xgboost(
  params = xgb_params,
  data = as.matrix(train %>% na.omit %>% select(-ft, -iso3c, -id_1)),
  label = train %>% na.omit %>% mutate(ft = ft/100) %>% pull(ft),
  nrounds = best_iter, verbose = 0)

# Evaluate performance on test set to check looks okay
pred <- predict(xgb_model, as.matrix(test %>% na.omit %>% select(-ft, -iso3c, -id_1)))
rmse <- sqrt(mean(((test %>% na.omit %>% mutate(ft = ft/100) %>% pull(ft)) - pred) ^ 2))
plot(pred, test %>% na.omit %>% mutate(ft = ft/100) %>% pull(ft),
     xlab = "XGBoost Model Predictions of ft", ylab = "Observed ft")
abline(0, 1, col = "red")

# use this to fill the gaps
full_df$ft[which(is.na(full_df$ft))] <- predict(
  xgb_model,
  as.matrix(full_df %>% filter(is.na(ft)) %>% select(-ft, -iso3c, -id_1))
  ) * 100

# add in the DHS informed fts
new_covars  <- covars %>% left_join(full_df %>% select(id_1, ft, pop_total) %>% rename(fs = ft))

# work out the subnational ft
new_covars <- new_covars %>%
  group_by(iso3c) %>%
  mutate(ft_high = ((ft_high[1]*sum(pop_total))/(sum(fs*pop_total))) * fs) %>%
  mutate(ft_low = ((ft_low[1]*sum(pop_total))/(sum(fs*pop_total))) * fs) %>%
  mutate(ft_mean = ((ft_mean[1]*sum(pop_total))/(sum(fs*pop_total))) * fs) %>%
  select(-fs) %>%
  select(-pop_total)

saveRDS(new_covars, "analysis/data_derived/global_covariate_ranges.rds")

# ---------------------------------------------------- #
# 4. Eth & Eri, use this to estimate historical ft and params ----
# ---------------------------------------------------- #

# This section is just creating data for the model fitting for fitness and rdtdet in script 6

# get historical country data from the MAP json file
# https://data.malariaatlas.org/map-platform-app-backend/api/v1/dnd/demand/historical
# This underpins their commodity historical treatment data
# https://data.malariaatlas.org/case-management

# get data and set up objects
dat <- jsonify::from_json("analysis/data_raw/MAP_comm_historical.json", simplify = TRUE)
isos <- c("ETH", "ERI", "UGA")
isols <- match(isos, dat$countryIso3)
isols_list <- vector("list", length(isols))

# get the iso covariates from 2020 previously collated and corrected for outiers
iso_covariates <- readRDS("analysis/data_derived/iso_covariates.rds")

for(i in seq_along(isols)) {
isols_i <- isols[i]

# convert data from json
df_i <- data.frame(year = dat$commodityDemand$year[isols_i,])
df_i$iso3c <- isos[i]

df_i$ft_publ <- dat$commodityDemand$breakdown$nCareseekingUnder5$pub[isols_i,] /
  dat$commodityDemand$breakdown$nTotalFeversUnder5[isols_i,]
df_i$ft_priv <- dat$commodityDemand$breakdown$nCareseekingUnder5$priv[isols_i,] /
  dat$commodityDemand$breakdown$nTotalFeversUnder5[isols_i,]

df_i$test_publ <- dat$commodityDemand$breakdown$nTotalTested$pub[isols_i,] /
  (dat$commodityDemand$breakdown$nCareseekingUnder5$pub[isols_i,] +
     dat$commodityDemand$breakdown$nCareseekingOver5$pub[isols_i,])
df_i$test_priv <- dat$commodityDemand$breakdown$nTotalTested$priv[isols_i,] /
  (dat$commodityDemand$breakdown$nCareseekingUnder5$priv[isols_i,] +
     dat$commodityDemand$breakdown$nCareseekingOver5$priv[isols_i,])

df_i$rdt_prop <- dat$commodityDemand$rdtDemand$pub[isols_i,]/
  dat$commodityDemand$breakdown$nTotalTested$pub[isols_i,]

# correct for RDT missingness for ETH
if(isos[i] == "ETH") {
  df_i$rdt_prop <- iso_covariates$rdt_prop[iso_covariates$iso3c == "ETH"]
}

df_i$postesttreat_publ <- dat$commodityDemand$breakdown$nTestedPosGotAm$pub[[isols_i]] /
  dat$commodityDemand$breakdown$nTestedPos$pub[[isols_i]]

df_i$postesttreat_priv <- dat$commodityDemand$breakdown$nTestedPosGotAm$priv[[isols_i]] /
  dat$commodityDemand$breakdown$nTestedPos$priv[[isols_i]]

df_i$negtesttreat_publ <- dat$commodityDemand$breakdown$nTestedNegGotAm$pub[[isols_i]] /
  dat$commodityDemand$breakdown$nTestedNeg$pub[isols_i,]

df_i$negtesttreat_priv <- dat$commodityDemand$breakdown$nTestedNegGotAm$priv[[isols_i]] /
  dat$commodityDemand$breakdown$nTestedNeg$priv[isols_i,]

df_i$untesttreat_publ <- (dat$commodityDemand$breakdown$nUntestedPosGotAm$pub[[isols_i]] +
                          dat$commodityDemand$breakdown$nUntestedNegGotAm$pub[[isols_i]]) /
  (dat$commodityDemand$breakdown$nUntestedPos$pub[[isols_i]] +
     dat$commodityDemand$breakdown$nUntestedPosO5$pub[[isols_i]] +
     dat$commodityDemand$breakdown$nUntestedNeg$pub[[isols_i]] +
     dat$commodityDemand$breakdown$nUntestedNegO5$pub[[isols_i]])

df_i$untesttreat_priv <- (dat$commodityDemand$breakdown$nUntestedPosGotAm$priv[[isols_i]] +
                          dat$commodityDemand$breakdown$nUntestedNegGotAm$priv[[isols_i]]) /
  (dat$commodityDemand$breakdown$nUntestedPos$priv[[isols_i]] +
     dat$commodityDemand$breakdown$nUntestedPosO5$priv[[isols_i]] +
     dat$commodityDemand$breakdown$nUntestedNeg$priv[[isols_i]] +
     dat$commodityDemand$breakdown$nUntestedNegO5$priv[[isols_i]])

# correct for untest outliers in ERI
if(isos[i] %in% c("ERI", "UGA")) {
  df_i$untesttreat_priv <- df_i$untesttreat_priv *
   (iso_covariates$untesttreat_priv[iso_covariates$iso3c == isos[i]] / tail(df_i$untesttreat_priv,1))
  df_i$untesttreat_publ <- df_i$untesttreat_publ *
    (iso_covariates$untesttreat_publ[iso_covariates$iso3c == isos[i]] / tail(df_i$untesttreat_publ,1))
}

# overall treatment for symptomatic cases
df_i$ft <- (df_i$ft_priv * df_i$test_priv * df_i$postesttreat_priv) +
  (df_i$ft_priv * (1-df_i$test_priv) * df_i$untesttreat_priv) +
  (df_i$ft_publ * (df_i$test_publ) * df_i$postesttreat_publ) +
  (df_i$ft_publ * (1-df_i$test_publ) * df_i$untesttreat_publ)

# use of microscopy assuming not used in private
df_i$micro <- ((1 - df_i$rdt_prop)*df_i$ft_publ)/(df_i$ft_priv + df_i$ft_publ)

# test non adherence
df_i$non_adherence <- df_i$ft_priv/(df_i$ft_priv + df_i$ft_publ) * df_i$negtesttreat_priv +
  (1- (df_i$ft_priv/(df_i$ft_priv + df_i$ft_publ))) * df_i$negtesttreat_publ

isols_list[[i]] <- df_i

}

# create covariate data.frame here
isos_df <- do.call(rbind, isols_list)

# Now add the prev
read_map_json <- function(type, prev = "Micro.2.10_mean") {

  url <- paste0(
    "https://data.malariaatlas.org/global-data/json-stat/2022/v1/malaria/admin1/Pf_PR-rate_2010-2020_2-10-years_",
    type,
    ".json")
  lci <- jsonlite::read_json(url)

  df <- data.frame("prev" = unlist(lci$value)) %>%
    cbind(expand.grid(
      "year" = unlist(lci$dimension$Year)[-1],
      "id_1" = as.character(unlist(lci$dimension$ID))[-1]
    ))

  names(df)[1] <- prev

  return(df)

}
mean <- read_map_json("mean", "Micro.2.10_mean")
low <- read_map_json("lci", "Micro.2.10_low")
high <- read_map_json("uci", "Micro.2.10_high")
df_prev <- mean %>% mutate(id_1 = as.integer(as.character(id_1))) %>%
  left_join(low %>% mutate(id_1 = as.integer(as.character(id_1)))) %>%
  left_join(high %>% mutate(id_1 = as.integer(as.character(id_1))))


  # use this just for matching codes for MAP data
  tf <- tempfile()
download.file("https://cloud.ihme.washington.edu/s/teDKnPGcJnBjJ5F/download?path=%2F2%20-%20Data%20%5BCSV%5D%2F3%20-%20Prevalence%3A%20Plasmodium%20falciparum%20%5BGeoTIFF%5D&files=IHME_MALARIA_2000_2019_PFPR_ADMIN1_Y2020M08D31.CSV", tf)
prev <- read.csv(tf)

# have to get prevalence for Uganda from this for admin 1
df_prev <- rbind(prev %>% filter(ISO3 == "UGA") %>%
  rename(Micro.2.10_mean = mean,
         Micro.2.10_low = lower,
         Micro.2.10_high = upper,
         id_1 = ADM1_Code) %>%
  select(Micro.2.10_mean, year, id_1, Micro.2.10_low, Micro.2.10_high) %>%
    arrange(id_1),
  df_prev)

prev_vars <- left_join(
  df_prev,
  prev %>% filter(year == 2019) %>% # just used for id_1 matching here
    select(ISO3, ADM1_Code) %>%
    rename(iso3c = ISO3, id_1 = ADM1_Code)
) %>%
  filter(iso3c %in% c("ERI", "ETH", "UGA")) %>%
  select(iso3c, id_1, year, Micro.2.10_mean, Micro.2.10_low, Micro.2.10_high) %>%
  arrange(iso3c)

# join with our iso level covars
spec_covars <- left_join(prev_vars %>% mutate(year = as.integer(as.character(year))) %>% filter(year >= 2010),
                         isos_df, by = c("iso3c", "year"))

# add in the DHS informed fts
new_spec_covars  <- spec_covars %>% left_join(full_df %>% select(id_1, ft, pop_total) %>% rename(fs = ft))

# work out the subnational ft
new_spec_covars <- new_spec_covars %>%
  group_by(iso3c) %>%
  mutate(ft_mean = ((ft*sum(pop_total))/(sum(fs*pop_total))) * fs) %>%
  select(iso3c, id_1, year, Micro.2.10_mean, Micro.2.10_low, Micro.2.10_high, ft_mean, micro, non_adherence) %>%
  mutate(id_1 = as.integer(id_1)) %>%
  rename(rdt.nonadherence_mean = non_adherence) %>%
  rename(microscopy.use_mean = micro) %>%
  mutate(fitness_mean = 0.9) %>%
  mutate(rdt.det_mean = 0.46)

saveRDS(new_spec_covars, "analysis/data_derived/covars_for_fitness_rdtdet_fitting.rds")
