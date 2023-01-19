# get country data
dat <- jsonify::from_json("analysis/data_raw/MAP.json")

names(dat)

dat$communitySettings

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

# suitable ranges to then be used for simulations
eir_range <- magenta:::age_brackets(num_age_brackets = 10, max_age = 130)[-1]
ft_range <- seq(0.01, 0.9, length.out = 9)
micro_range <- seq(0, max(df$micro), length.out = 4)
nonadherence_range <- seq(0, max(df$non_adherence), length.out = 4)
fitness <- c(1, 0.95, 0.90)
hrp3 <- c(0, 0.25, 0.5)
nmf <- c(0.75, 1, 1.25)

# starter param grid
param.start <- expand.grid(
  EIR = eir_range,
  ft = ft_range
) %>%
  mutate(id = seq_len(n()))
saveRDS(param.start, "analysis/data_derived/param_start.rds")

# continuation param grid
param.grid <- expand.grid(
  EIR = eir_range,
  ft = ft_range,
  microscopy.use = micro_range,
  rdt.nonadherence = nonadherence_range,
  fitness = fitness,
  rdt.det = hrp3,
  nmf = nmf
)
param.grid <- left_join(param.grid, param.start, by = c("EIR", "ft"))

dir.create("analysis/data_derived/")
saveRDS(param.grid, "analysis/data_derived/param_grid.rds")





