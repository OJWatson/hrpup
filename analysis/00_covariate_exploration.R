# necessary packages
# install.packages(c("jsonify", "giscoR", "tidyverse"))

# get country data from the MAP json file
# https://data.malariaatlas.org/map-platform-app-backend/api/v1/dnd/policy/defaults
# This underpins their commodity forecasts. The default is year 2020 maintained till 2030
# https://data.malariaatlas.org/case-management
dat <- jsonify::from_json("analysis/data_raw/MAP_comm_default.json")

# convert data from json
df <- data.frame("iso3c" = dat$countryIso3)
df$ft_priv <- dat$clinicSettings$rateCareseeking$priv
df$ft_publ <- dat$clinicSettings$rateCareseeking$pub
df$test_priv <- dat$clinicSettings$rateTesting$priv
df$test_publ <- dat$clinicSettings$rateTesting$pub
df$rdt_prop <- dat$clinicSettings$rateRdt
df$postesttreat_priv <- dat$complianceSettings$rateTestedPosGotAm$priv
df$postesttreat_publ <- dat$complianceSettings$rateTestedPosGotAm$pub
df$negtesttreat_priv <- dat$complianceSettings$rateTestedNegGotAm$priv
df$negtesttreat_publ <- dat$complianceSettings$rateTestedNegGotAm$pub
df$untesttreat_priv <- dat$complianceSettings$rateUntestedGotAm$priv
df$untesttreat_publ <- dat$complianceSettings$rateUntestedGotAm$pub
df$ACT_publ <- dat$complianceSettings$rateAmAreAct$pub

# overall treatment for symptomatic cases
df$ft <- (df$ft_priv * df$test_priv * df$postesttreat_priv) +
  (df$ft_priv * (1-df$test_priv) * df$untesttreat_priv) +
  (df$ft_publ * (df$test_publ) * df$postesttreat_publ) +
  (df$ft_publ * (1-df$test_publ) * df$untesttreat_publ)

# use of microscopy assuming not used in private
df$micro <- ((1 - df$rdt_prop)*df$ft_publ)/(df$ft_priv + df$ft_publ)

# test non adherence
df$non_adherence <- df$ft_priv/(df$ft_priv + df$ft_publ) * df$negtesttreat_priv +
  (1- (df$ft_priv/(df$ft_priv + df$ft_publ))) * df$negtesttreat_publ

# huh Uganda looks high
df %>% filter(iso3c == "UGA")

# okay clearly there is something wrong with Uganda as they clearly do not always treat
# those who are not tested (clearly not accurate).
# It also looks like an outlier mathematically and also based on our literature review
df$untesttreat_priv %>% hist(breaks = 20)

# and it definitely would be classified based on IQR criteria
to_check <- names(df)[-1]
outliers <- lapply(to_check,function(x) {
  boxplot.stats(df[[x]])$out
}) %>% set_names(to_check)

# These are the full remit of relevant outliers that are not derived from other params

# The first is malaysia which given very low malaria may be correct
df %>% filter(ft_priv %in% outliers$ft_priv)

# BDI test rate in private sector will explore
df %>% filter(test_priv %in% outliers$test_priv)

# explore all untest incidences
df %>% filter(untesttreat_priv %in% outliers$untesttreat_priv)

# Set these to missing and infer using random forest
df <- df %>%
  mutate(test_priv = replace(test_priv, test_priv %in% outliers$test_priv, NA),
         untesttreat_priv = replace(untesttreat_priv, untesttreat_priv %in% outliers$untesttreat_priv, NA),
         untesttreat_publ = replace(untesttreat_publ, is.na(untesttreat_priv), NA),
         # rdt proprtion not reliably inferred with mice rf, so set to mean for Ethiopia
         # (other outliers are suitable (South America or MYS) based on literature and S. American RDT patterns)
         rdt_prop = replace(rdt_prop, rdt_prop < 0.05, mean(rdt_prop))) %>%
  mice::mice(method = "rf") %>%
  complete %>%
  mutate(untesttreat_publ = untesttreat_priv)

# overall treatment for symptomatic cases
df$ft <- (df$ft_priv * df$test_priv * df$postesttreat_priv) +
  (df$ft_priv * (1-df$test_priv) * df$untesttreat_priv) +
  (df$ft_publ * (df$test_publ) * df$postesttreat_publ) +
  (df$ft_publ * (1-df$test_publ) * df$untesttreat_publ)

# use of microscopy assuming not used much in private
df$micro <- ((1 - df$rdt_prop)*df$ft_publ)/(df$ft_priv + df$ft_publ)

# test non adherence
df$non_adherence <- df$ft_priv/(df$ft_priv + df$ft_publ) * df$negtesttreat_priv +
  (1- (df$ft_priv/(df$ft_priv + df$ft_publ))) * df$negtesttreat_publ

# save the output for the time being
saveRDS(df, "analysis/data_derived/iso_covariates.rds")
