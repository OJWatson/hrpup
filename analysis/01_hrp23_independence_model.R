library(tidyverse)

# -----------------------------------------------------#
# 1. Read in formatted WHO hrp2/hrp3 data  -----------------
# -----------------------------------------------------#

hrpdat <- read.csv("analysis/data_raw/WHO_hrp2_hrp3_data.csv")

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
                                   column.sep.width = "100pt",
                                   ci=TRUE, ci.level = 0.95,
                                   single.row = TRUE,
                                   digits = 4,
                                   dep.var.labels.include = FALSE,
                                   covariate.labels = "Prevalence",
                                   title = "Binomial Modelling of Impact of malaria prevalence on pfhrp2/3 status",
                                   omit = "Constant", dep.var.caption = "")
write_stargazer(star_model, "analysis/tables/pfhrp3_models.csv", spl = "\\s{2,}|\\)\\s{1,}")


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
