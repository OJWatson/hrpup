library(tidyverse)

# -----------------------------------------------------#
# 1. Read in and format data  -----------------
# -----------------------------------------------------#

# read in our data
hrpdat <- readxl::read_xlsx("analysis/data_raw/HRP2_map_13.1.23.xlsx") %>%
  rename(WHO_HRP2_HRP3_DELETIONS_PERCENT = WHO_HPR2_HRP3_DELETIONS_PERCENT,
         WHO_HRP2_HRP3_DELETIONS_TESTED_N = WHO_HPR2_HRP3_DELETIONS_TESTED_N,
         WHO_HRP2_HRP3_DELETIONS_POS_N = WHO_HPR2_HRP3_DELETIONS_POS_N)

# Fill in missing data

# fill in missing year for this one study
hrpdat$WHO_YEAR_START[is.na(hrpdat$WHO_YEAR_START) & hrpdat$ISO_ctry == "CO"] <- 1999

# Aweil state hospital location
hrpdat$WHO_LATITUDE[which(hrpdat$WHO_ADMIN2 == "Aweil City")] <- 8.7680936
hrpdat$WHO_LONGITUDE[which(hrpdat$WHO_ADMIN2 == "Aweil City")] <- 27.3288594

# Yambio
hrpdat$WHO_LATITUDE[which(hrpdat$WHO_ADMIN1 == "Yambio County")] <- 4.626894
hrpdat$WHO_LONGITUDE[which(hrpdat$WHO_ADMIN1 == "Yambio County")] <- 28.411925

# Gia Lai location average
hrpdat$WHO_LATITUDE[which(hrpdat$WHO_ADMIN2 == "Gia Lai and Dak Lak provinces")] <- 13.394084
hrpdat$WHO_LONGITUDE[which(hrpdat$WHO_ADMIN2 == "Gia Lai and Dak Lak provinces")] <- 108.706456

# Couple African locations are in the sea
hrpdat$WHO_LONGITUDE[which(hrpdat$UNIQUE_ID == "GQ-0002")] <- 8.6156703
hrpdat$WHO_LONGITUDE[which(hrpdat$UNIQUE_ID == "ER-0002")] <- 39.365273
hrpdat$WHO_LONGITUDE[which(hrpdat$UNIQUE_ID == "ET-0008")] <- 37.912280


# Now just select data we want and filter to studies that did not just measure hrp2
hrpdat <- hrpdat %>% select(ISO_ctry, WHO_LATITUDE, WHO_LONGITUDE, WHO_YEAR_START, WHO_YEAR_END, WHO_PATIENT_TYPE,
                            MEAN_AGE, WHO_RDT_TYPE, matches("WHO_HRP")) %>%
  mutate(across(tidyselect::everything(),.fns = str_trim)) %>%
  mutate(across(matches("TUDE"), as.numeric))

hrpdat <- hrpdat %>% filter(WHO_HRP2_ONLY == 0)

# convert NRs
conv <- names(which(unlist(lapply(hrpdat, function(x){sum(x == "NR", na.rm = TRUE)>0}))))
hrpdat <- hrpdat %>%
  mutate(across(all_of(conv), .fns = ~as.numeric(replace(.x, .x == "NR", NA))))

# convert numerics
hrpdat <- hrpdat %>%
  mutate(across(WHO_LATITUDE, WHO_LONGITUDE, WHO_PATIENT_TYPE,MEAN_AGE,
                .fns = as.numeric)) %>%
  mutate(across(all_of(grep("WHO_HRP", names(hrpdat), value = TRUE)),
                .fns = as.numeric))

# sort out dates
my_make <- function(x){
  x[!grepl("/", x, fixed = TRUE)] <- paste0("01/", x[!grepl("/", x, fixed = TRUE)])
  x[grep("NA",x)] <- NA
  lubridate::my(x)
}

# sort out middate range for getting prevalence
hrpdat <- hrpdat %>%
  mutate(WHO_YEAR_START = my_make(WHO_YEAR_START),
         WHO_YEAR_END = my_make(WHO_YEAR_END)) %>%
  mutate(mid_date = WHO_YEAR_START + (WHO_YEAR_END - WHO_YEAR_START)/2)

# Lastly add in MAP estimated malaria prevalence for these regions
prevs <- hrpdat %>%
  split(~lubridate::year(hrpdat$mid_date)) %>%
  map(function(x){
    malariaAtlas::extractRaster(
      surface = "Plasmodium falciparum PR2 - 10 version 2020",
      df = data.frame(lat = x$WHO_LATITUDE, long = x$WHO_LONGITUDE),
      year = min(2019, lubridate::year(unique(x$mid_date)))
    )
  })

hrpdat <- hrpdat %>%
  split(~lubridate::year(hrpdat$mid_date)) %>%
  list_rbind() %>%
  mutate(prev = list_rbind(prevs)$value) %>%
  mutate(prev = replace(prev, prev == -9999, NA)) %>%
  mutate(prev = replace(prev, prev <= 0, 0.00001)) # very low South American prevalence set to 1e-5

hrpdat$continent <- as.factor(countrycode::countrycode(hrpdat$ISO_ctry,"iso2c", "continent"))
hprdat <- hrpdat %>%
  select(WHO_HRP2_HRP3_DELETIONS_PERCENT, WHO_HRP2_DELETION_TESTED_N, WHO_HRP3_DELETION_TESTED_N,
         WHO_HRP2_HRP3_DELETIONS_TESTED_N, WHO_HRP2_DELETION_POS_N, WHO_HRP2_DELETION_PERCENT,
         WHO_HRP2_HRP3_DELETIONS_POS_N, prev, continent, ISO_ctry)

write.csv(hrpdat, "analysis/data_raw/WHO_hrp2_hrp3_data.csv")

# -----------------------------------------------------#
# 2. Estimate overall probability of hrp2 deleted parasites being pfhrp3 deleted  -----------------
# -----------------------------------------------------#

# get just the samples where all hrp2/3 testing was done
final <- hrpdat %>% filter(!is.na(WHO_HRP2_HRP3_DELETIONS_PERCENT)) %>%
  filter(WHO_HRP2_DELETION_TESTED_N == WHO_HRP3_DELETION_TESTED_N) %>%
  filter(WHO_HRP2_DELETION_TESTED_N == WHO_HRP2_HRP3_DELETIONS_TESTED_N)

# filter to just the cases with hrp2
final <- final %>% filter(WHO_HRP2_DELETION_POS_N > 0)

# also filter out where sampling was clearly not done on the same samples
# i.e. we shouldn't be able to have more samples with both deletions than
# ust one deletion assuming the deletion percentages are reported for all
# deletions, i.e. an hrp2 deletion percentage report is for any hrp2, i.e.
# with or without hrp3 and thus the double can't be more assuming this and
# assuming it was the same samples that are tested.
final <- final %>% filter(WHO_HRP2_HRP3_DELETIONS_PERCENT <= WHO_HRP2_DELETION_PERCENT)

# now fit a beta binomial model for the mean probability that hrp2 is found with hrp3 deletions
x1 <- final$WHO_HRP2_HRP3_DELETIONS_POS_N
size <- final$WHO_HRP2_DELETION_POS_N

# the model should be fit per continent as different evolution processes in S America for sure
final$continent <- as.factor(countrycode::countrycode(final$ISO_ctry,"iso2c", "continent"))

mtmp <- function(prob, theta, continent) {
  if(continent == 0) {
    -sum(emdbook::dbetabinom(
      final$WHO_HRP2_HRP3_DELETIONS_POS_N,
      prob,
      final$WHO_HRP2_DELETION_POS_N,
      theta,
      log=TRUE
    ))
  } else {
  -sum(emdbook::dbetabinom(
    final$WHO_HRP2_HRP3_DELETIONS_POS_N[as.integer(final$continent) %in% continent],
    prob,
    final$WHO_HRP2_DELETION_POS_N[as.integer(final$continent) %in% continent],
    theta,
    log=TRUE
    ))
  }
}

# summary gives 66% probability
mafrica <- bbmle::mle2(mtmp, start=list(prob=0.8,theta=9), fixed = list("continent" = 1))
bbmle::confint(mafrica)
bbmle::mle2(mtmp, start=list(prob=0.8,theta=9), fixed = list("continent" = 2))
bbmle::mle2(mtmp, start=list(prob=0.8,theta=9), fixed = list("continent" = 3))
bbmle::mle2(mtmp, start=list(prob=0.8,theta=9), fixed = list("continent" = 4))
# not enough data here
# bbmle::mle2(mtmp, start=list(prob=0.8,theta=9), fixed = list("continent" = 5))
m0 <- bbmle::mle2(mtmp, start=list(prob=0.8,theta=9), fixed = list("continent" = 0))
bbmle::confint(m0)
# -----------------------------------------------------#
# 3. Plot our overall data and results
# -----------------------------------------------------#

# add the cooccurrence onto our data
final$coocc <- final$WHO_HRP2_HRP3_DELETIONS_POS_N / final$WHO_HRP2_DELETION_POS_N
final$coocc_low <- as.numeric(Hmisc::binconf(final$WHO_HRP2_HRP3_DELETIONS_POS_N, final$WHO_HRP2_DELETION_POS_N)[,2])
final$coocc_high <- as.numeric(Hmisc::binconf(final$WHO_HRP2_HRP3_DELETIONS_POS_N, final$WHO_HRP2_DELETION_POS_N)[,3])
final$id <- seq_along(final$ISO_ctry)

# independence overall
gg1 <- ggplot(final, aes(x = coocc,y=as.factor(id), xmin = coocc_low, xmax = coocc_high, size = WHO_HRP2_DELETION_POS_N)) +
  geom_vline(xintercept = m0@coef[1], color = "red", linetype = "dashed") +
  geom_pointrange() +
  scale_size_binned(range = c(0,0.75), name = "Total *pfhrp2* deleted <br /> samples in survey:", breaks = c(0,10,50,100,150,200)) +
  labs(x = "Proportion of *pfhrp2-deleted* samples with *pfhrp3* deletions") +
  ylab("Survey") +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(axis.title.x = ggtext::element_markdown(size = 10),
        legend.title = ggtext::element_markdown(size = 10))

save_figs("hrp2_3_independence", gg1, width = 8, height = 6)

# -----------------------------------------------------#
# 4. Plot our overall and cooccurrence findings with regard prevalence for Africa
# -----------------------------------------------------#
pastellize <- function(x, p){

  hsv2rgb <- function(x){
    # convert an hsv colour to rgb
    # input:  a 3 x 1 matrix (same as output of rgb2hsv() function)
    # output: vector of length 3 with values in [0,1]

    # recover h, s, v values
    h <- x[1,1]
    s <- x[2,1]
    v <- x[3,1]

    # follow the algorithm from Wikipedia
    C <- s*v

    # in R, h takes values in [0,1] rather than [0, 360], so dividing by
    # 60 degrees is the same as multiplying by six
    hdash <- h*6
    X <- C * (1 - abs(hdash %% 2 -1))

    if (0 <= hdash & hdash <=1) RGB1 <- c(C, X, 0)
    if (1 <= hdash & hdash <=2) RGB1 <- c(X, C, 0)
    if (2 <= hdash & hdash <=3) RGB1 <- c(0, C, X)
    if (3 <= hdash & hdash <=4) RGB1 <- c(0, X, C)
    if (4 <= hdash & hdash <=5) RGB1 <- c(X, 0, C)
    if (5 <= hdash & hdash <=6) RGB1 <- c(C, 0, X)

    # the output is a vector of length 3. This is the most convenient
    # format for using as the col argument in an R plotting function
    RGB1 + (v-C)
  }

  # x is a colour
  # p is a number in [0,1]
  # p = 1 will give no pastellization

  # convert hex or letter names to rgb
  if (is.character(x)) x <- col2rgb(x)/255

  # convert vector to rgb
  if (is.numeric(x)) x <- matrix(x, nr=3)

  col <- rgb2hsv(x, maxColorValue=1)
  col[2,1] <- col[2,1]*p
  col <- hsv2rgb(col)

  # return in convenient format for plots
  rgb(col[1], col[2], col[3])
}


# independence africa
gg2 <- final %>%
  filter(continent == "Africa") %>%
  mutate(country = countrycode::countrycode(
    .data$ISO_ctry,"iso2c", "country.name",
    custom_match = c("CD" = "Democractic Republic of the Congo"))) %>%
  ggplot(aes(x = coocc,y=fct_reorder(as.factor(id), country, .desc = TRUE), group = (as.factor(country)),
             xmin = coocc_low, xmax = coocc_high, size = WHO_HRP2_DELETION_POS_N,
             color = (as.factor(country)))) +
  geom_vline(xintercept = mafrica@coef[1], color = "black", linetype = "dashed", size = 0.75) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  scale_x_continuous(labels = scales::percent) +
  scale_size_binned(range = c(0,0.85),
                    name = "Total *pfhrp2* <br /> deleted <br /> samples <br />in  survey:<br>",
                    breaks = c(0,10,50,100,150,200, 250),
                    limits = c(0, 250)) +
  scale_color_manual(name = "Country",
                     values = as.character(vapply(rev(pals::polychrome()), pastellize, character(1), 0.95))) +
  labs(x = "Percentage of *pfhrp2* deleted samples <br /> with *pfhrp3* gene deletions") +
  labs(y = "WHO Threat Map Survey of *pfhrp2*/*3* deletions") +
  theme_bw(base_family = "Helvetica", base_size = 12) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(axis.title.x = ggtext::element_markdown(size = 12),
        legend.title = ggtext::element_markdown(size = 12),
        legend.text = ggtext::element_markdown(size = 9),
        axis.title.y = ggtext::element_markdown(size = 12))

# cooccurence by Africa
gg_23_coocc <- final %>%
  filter(continent == "Africa") %>%
  ggplot(aes(x = prev, y=WHO_HRP2_HRP3_DELETIONS_POS_N/WHO_HRP2_DELETION_POS_N,
             pos=WHO_HRP2_HRP3_DELETIONS_POS_N, total=WHO_HRP2_DELETION_POS_N,
             size = WHO_HRP2_DELETION_POS_N)) +
  geom_smooth(
    method="glm",
    method.args=list(family="binomial"),
    formula = cbind(pos, total-pos) ~ x,
    fullrange = TRUE,
    show.legend = FALSE
  ) +
  geom_point() +
  scale_x_log10(limits = c(0.0001, 1), labels = scales::percent) +
  scale_size_binned(range = c(0,4),
                  name = "Total *pfhrp2* <br /> deleted <br /> samples <br />in  survey:<br>",
                  breaks = c(0,10,50,100,150,200, 250),
                  limits = c(0, 250)) +
  scale_y_continuous(labels = scales::percent)+
  labs(x = "Malaria Slide Prevalence 2-10 <br />") +
  labs(y = "Percentage of *pfhrp2* deleted samples <br /> with *pfhrp3* gene deletions") +
  theme_bw(base_family = "Helvetica", base_size = 12) +
  theme(axis.title.x = ggtext::element_markdown(size = 12),
        legend.title = ggtext::element_markdown(size = 12),
        legend.text = ggtext::element_markdown(size = 9),
        axis.title.y = ggtext::element_markdown(size = 12),
        legend.box = "horizontal")


gg_3_prev <- hrpdat %>%
  filter(continent == "Africa") %>%
  ggplot(aes(x = (prev), y=WHO_HRP3_DELETION_POS_N/WHO_HRP3_DELETION_TESTED_N,
             pos=WHO_HRP3_DELETION_POS_N, total=WHO_HRP3_DELETION_TESTED_N,
             size = WHO_HRP3_DELETION_TESTED_N)) +
  geom_smooth(
    method="glm",
    method.args=list(family="binomial"),
    formula = cbind(pos, total-pos) ~ (x),
    fullrange = TRUE,
    show.legend = FALSE
  ) +
  geom_point() +
  scale_x_log10(limits = c(0.0001, 1), labels = scales::percent) +
  scale_size_binned(range = c(0,4),
                    name = "Total <br /> samples <br /> tested for <br />*pfhrp3*: <br />",
                    breaks = c(0,10,50,100,150,200, 250),
                    limits = c(0, 250)) +
  scale_y_continuous(labels = scales::percent)+
  labs(x = "Malaria Slide Prevalence 2-10 <br />") +
  labs(y = "Percentage of samples tested for *pfhrp3* <br /> with *pfhrp3* gene deletions") +
  theme_bw(base_family = "Helvetica", base_size = 12) +
  theme(axis.title.x = ggtext::element_markdown(size = 12),
        legend.title = ggtext::element_markdown(size = 12),
        legend.text = ggtext::element_markdown(size = 9),
        axis.title.y = ggtext::element_markdown(size = 12),
        legend.box = "horizontal")

figure1 <- cowplot::plot_grid(gg2,
                   cowplot::plot_grid(gg_23_coocc, gg_3_prev, ncol = 1, labels = c("B", "C"), align = "v"),
                   labels = c("A", ""), rel_widths = c(0.55,0.45),
                   ncol = 2
)

save_figs("hrp2_3_multi_independence_africa", figure1, width = 12, height = 8)

# Significance of prevalence
glm(cbind(WHO_HRP2_HRP3_DELETIONS_POS_N, WHO_HRP2_DELETION_TESTED_N - WHO_HRP2_HRP3_DELETIONS_POS_N) ~ log10(prev),
            data = hrpdat %>% filter(continent == "Africa") %>% mutate(prev = prev), family = binomial) %>%
  broom::tidy()

glm(cbind(WHO_HRP2_HRP3_DELETIONS_POS_N, WHO_HRP2_DELETION_TESTED_N - WHO_HRP2_HRP3_DELETIONS_POS_N) ~ log10(prev),
    data = hrpdat %>% filter(continent == "Africa") %>% mutate(prev = prev), family = binomial) %>%
  broom::tidy()

# Models for each
mod1 <- glm(cbind(WHO_HRP2_HRP3_DELETIONS_POS_N, WHO_HRP2_DELETION_TESTED_N - WHO_HRP2_HRP3_DELETIONS_POS_N) ~ log10(prev),
    data = hrpdat, family = binomial)

mod2 <- glm(cbind(WHO_HRP3_DELETION_POS_N, WHO_HRP3_DELETION_TESTED_N - WHO_HRP3_DELETION_POS_N) ~ log10(prev),
    data = hrpdat, family = binomial)

star_model <- stargazer::stargazer(mod1, mod2, type = "text",
                     keep.stat = c("n","ll"),
                     report = "vcsp",
                     column.labels = c("pfhrp3- given pfhrp2-",
                                       "pfhrp3-"),
                     column.sep.width = "90pt",
                     ci=TRUE, ci.level = 0.95,
                     single.row = TRUE,
                     digits = 4,
                     dep.var.labels.include = FALSE,
                     covariate.labels = "Prevalence",
                     title = "Binomial Modelling of Impact of malaria prevalence on pfhrp2/3 status",
                     omit = "Constant", dep.var.caption = "")
write_stargazer(star_model, "analysis/tables/pfhrp3_models.csv")


# -----------------------------------------------------#
# 5. Same as above but for all samples
# -----------------------------------------------------#

# independence
gg2 <- final %>%
  filter(continent %in% c("Africa", "Americas", "Asia")) %>%
  mutate(country = countrycode::countrycode(
    .data$ISO_ctry,"iso2c", "country.name",
    custom_match = c("CD" = "Democractic Republic of the Congo"))) %>%
  ggplot(aes(x = coocc,y=fct_reorder(as.factor(id), country, .desc = TRUE), group = (as.factor(country)),
             xmin = coocc_low, xmax = coocc_high, size = WHO_HRP2_DELETION_POS_N,
             color = (as.factor(country)))) +
  geom_vline(xintercept = m0@coef[1], color = "black", linetype = "dashed", size = 0.75) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  scale_x_continuous(labels = scales::percent) +
  scale_size_binned(range = c(0,0.85),
                    #name = "Total *pfhrp2* <br /> deleted <br /> samples <br />in  survey:<br>",
                    name = "Total *pfhrp2* deleted samples in  survey:<br>",
                    breaks = c(0,10,50,100,150,200, 250),
                    limits = c(0, 250)) +
  scale_color_manual(name = "Country",
                     values = as.character(vapply(rev(pals::polychrome()), pastellize, character(1), 0.95))) +
  labs(x = "Percentage of *pfhrp2* deleted samples <br /> with *pfhrp3* gene deletions") +
  labs(y = "WHO Threat Map Survey of *pfhrp2*/*3* deletions") +
  theme_bw(base_family = "Helvetica", base_size = 12) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(axis.title.x = ggtext::element_markdown(size = 12),
        legend.title = ggtext::element_markdown(size = 12),
        legend.text = ggtext::element_markdown(size = 9),
        axis.title.y = ggtext::element_markdown(size = 12))

# cooccurence
gg_23_coocc <- final %>%
  filter(continent %in% c("Africa", "Americas", "Asia")) %>%
  ggplot(aes(x = prev, y=WHO_HRP2_HRP3_DELETIONS_POS_N/WHO_HRP2_DELETION_POS_N,
             pos=WHO_HRP2_HRP3_DELETIONS_POS_N, total=WHO_HRP2_DELETION_POS_N,
             color = continent,
             size = WHO_HRP2_DELETION_POS_N)) +
  geom_smooth(
    method="glm",
    method.args=list(family="binomial"),
    formula = cbind(pos, total-pos) ~ x,
    fullrange = TRUE,
    show.legend = FALSE
  ) +
  geom_point() +
  scale_x_log10(limits = c(0.0001, 1), labels = scales::percent) +
  scale_size_binned(range = c(0,4),
                    name = "Total *pfhrp2* <br /> deleted <br /> samples <br />in  survey:<br>",
                    breaks = c(0,10,50,100,150,200, 250),
                    limits = c(0, 250)) +
  scale_y_continuous(labels = scales::percent)+
  scale_color_viridis_d(name = "Continent", drop = TRUE) +
  labs(x = "Malaria Slide Prevalence 2-10 <br />") +
  labs(y = "Percentage of *pfhrp2* deleted samples <br /> with *pfhrp3* gene deletions") +
  theme_bw(base_family = "Helvetica", base_size = 12) +
  theme(axis.title.x = ggtext::element_markdown(size = 12),
        legend.title = ggtext::element_markdown(size = 12),
        legend.text = ggtext::element_markdown(size = 9),
        axis.title.y = ggtext::element_markdown(size = 12),
        legend.box = "horizontal")


gg_3_prev <- hrpdat %>%
  filter(continent %in% c("Africa", "Americas", "Asia")) %>%
  ggplot(aes(x = prev, y=WHO_HRP3_DELETION_POS_N/WHO_HRP3_DELETION_TESTED_N,
             pos=WHO_HRP3_DELETION_POS_N, total=WHO_HRP3_DELETION_TESTED_N,
             color = continent,
             size = WHO_HRP3_DELETION_TESTED_N)) +
  geom_smooth(
    method="glm",
    method.args=list(family="binomial"),
    formula = cbind(pos, total-pos) ~ x,
    fullrange = TRUE,
    show.legend = FALSE
  ) +
  geom_point() +
  xlim(c(0,1)) +
  scale_x_log10(limits = c(0.0001, 1), labels = scales::percent) +
  scale_size_binned(range = c(0,4),
                    name = "Total <br /> samples <br /> tested for <br />*pfhrp3*: <br />",
                    breaks = c(0,10,50,100,150,200, 250),
                    limits = c(0, 250)) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_viridis_d(name = "Continent", drop = TRUE) +
  labs(x = "Malaria Slide Prevalence 2-10 <br />") +
  labs(y = "Percentage of samples tested for *pfhrp3* <br /> with *pfhrp3* gene deletions") +
  theme_bw(base_family = "Helvetica", base_size = 12) +
  theme(axis.title.x = ggtext::element_markdown(size = 12),
        legend.title = ggtext::element_markdown(size = 12),
        legend.text = ggtext::element_markdown(size = 9),
        axis.title.y = ggtext::element_markdown(size = 12),
        legend.box = "horizontal")


# figure1 <- cowplot::plot_grid(gg2,
#                               cowplot::plot_grid(gg_23_coocc, gg_3_prev, ncol = 1, labels = c("B", "C"), align = "v"),
#                               labels = c("A", ""), rel_widths = c(0.55,0.45),
#                               ncol = 2
# )
figure1 <- cowplot::plot_grid(gg2,
                              cowplot::plot_grid(gg_23_coocc, gg_3_prev, ncol = 2, labels = c("B", "C"), align = "v"),
                              labels = c("A", ""),
                              ncol = 1, rel_heights = c(0.55,0.45)
)


save_figs("hrp2_3_multi_independence", figure1, width = 14, height = 12)

# Alternative appraoch based on 2 x 2 table

two_by_two <- final %>% filter(continent == "Africa") %>%
  select(WHO_HRP2_DELETION_TESTED_N:WHO_HRP2_HRP3_DELETIONS_PERCENT) %>%
  ungroup() %>%
  na.omit %>%
  mutate(dd = WHO_HRP2_HRP3_DELETIONS_POS_N,
            dp = WHO_HRP2_DELETION_POS_N-dd,
            pd = WHO_HRP3_DELETION_POS_N-dd,
            pp = WHO_HRP2_DELETION_TESTED_N - dd - dp - pd) %>%
  summarise(across(dd:pp, sum))

tbt_table <- two_by_two %>% as.numeric() %>% matrix(ncol = 2, nrow = 2)
chisq.test(tbt_table)

two_by_two <- two_by_two %>%
  mutate(n = sum(c(dd,dp,pd,pp))) %>%
  mutate(dd = dd/n,
         dp = dp/n,
         pd = pd/n,
         pp = pp/n)

# calculate our
p_two_d <- (two_by_two$dd+two_by_two$dp)
p_three_d <- (two_by_two$dd+two_by_two$dp)
D <- two_by_two$dd - (p_two_d*p_three_d)
D_prime <- D/min(c(p_two_d*(1-p_two_d),(p_three_d)*(1-p_three_d)))
R2 <- (D^2)/(p_two_d*(1-p_two_d)*(p_three_d)*(1-p_three_d))
