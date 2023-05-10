# necessary packages
# install.packages(c("jsonify", "giscoR", "tidyverse"))

# get country data from the MAP json file shared last time
dat <- jsonify::from_json("analysis/data_raw/MAP.json")

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

# overall treatment assuming that testing is near 100% for actual malaria cases
df$ft_alt <- (df$ft_priv  * df$postesttreat_priv) +
  (df$ft_publ * df$postesttreat_publ)

# use of microscopy assuming not used in private
df$micro <- ((1 - df$rdt_prop)*df$ft_publ)/(df$ft_priv + df$ft_publ)

# test non adherence
df$non_adherence <- df$ft_priv/(df$ft_priv + df$ft_publ) * df$negtesttreat_priv +
  (1- (df$ft_priv/(df$ft_priv + df$ft_publ))) * df$negtesttreat_publ


my_colors <- c("red", "blue", "green","orange", "yellow")
my_palette <- colorRampPalette(my_colors)

# quick plot of chance positive case is treated
world <- giscoR::gisco_get_countries(year = "2020", region = "Africa")
world <- left_join(world %>% mutate(iso3c = ISO3_CODE), df)



world %>%
  ggplot() +
  #geom_sf(aes(fill = cut(ft, seq(0,1,0.2))), color = NA, show.legend = TRUE) +
  geom_sf(aes(fill = ft_alt), color = NA, show.legend = TRUE) +
  # Robinson
  coord_sf(crs = "ESRI:54030") +
  theme_void() +
  labs(
    title = "  Probability of malaria clinical infection seeking treatment and being treated"
  ) +
  scale_fill_gradient(low = "#fff7ec", high = "#7f0000") +
  # scale_fill_manual(values = my_palette(5)) +
  theme(plot.caption = element_text(face = "italic"))

# huh Uganda looks high
df %>% filter(iso3c == "UGA")

# okay clearly there is something wrong with Uganda as they clearly do not always treat
# those who are not tested (clearly not accurate).
# It also looks like an outlier
df$untesttreat_priv %>% hist(breaks = 20)

# and it definitely would be classified based on IQR criteria
to_check <- names(df)[-1]
outliers <- lapply(to_check,function(x) {
  boxplot.stats(df[[x]])$out
}) %>% set_names(to_check)

# These are the full remit of relevant outliers that are not derived from other params

# first is malaysia which given very low malaria will not be a concern so will leave
df %>% filter(ft_priv %in% outliers$ft_priv)

# BDI test rate in pricate sector will explore
df %>% filter(test_priv %in% outliers$test_priv)

# explore all untest
df %>% filter(untesttreat_priv %in% outliers$untesttreat_priv)

# Set these to missing and infer using random forest
df <- df %>%
  mutate(test_priv = replace(test_priv, test_priv %in% outliers$test_priv, NA),
         untesttreat_priv = replace(untesttreat_priv, untesttreat_priv %in% outliers$untesttreat_priv, NA),
         untesttreat_publ = replace(untesttreat_publ, is.na(untesttreat_priv), NA),
         rdt_prop = replace(rdt_prop, rdt_prop < 0.05, mean(rdt_prop))) %>% # need to better work out rdt proportions for WHO edge cases
  mice::mice(method = "rf") %>%
  complete %>%
  mutate(untesttreat_publ = untesttreat_priv)

# overall treatment for symptomatic cases
df$ft <- (df$ft_priv * df$test_priv * df$postesttreat_priv) +
  (df$ft_priv * (1-df$test_priv) * df$untesttreat_priv) +
  (df$ft_publ * (df$test_publ) * df$postesttreat_publ) +
  (df$ft_publ * (1-df$test_publ) * df$untesttreat_publ)

# overall treatment assuming that testing is near 100% for actual malaria cases
df$ft_alt <- (df$ft_priv  * df$postesttreat_priv) +
  (df$ft_publ * df$postesttreat_publ)

# use of microscopy assuming not used in private
df$micro <- ((1 - df$rdt_prop)*df$ft_publ)/(df$ft_priv + df$ft_publ)

# test non adherence
df$non_adherence <- df$ft_priv/(df$ft_priv + df$ft_publ) * df$negtesttreat_priv +
  (1- (df$ft_priv/(df$ft_priv + df$ft_publ))) * df$negtesttreat_publ

# save the output for the time being
saveRDS(df, "analysis/data_derived/iso_covariates.rds")
