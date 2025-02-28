library(tidyverse)
devtools::load_all()

# -----------------------------------------------------#
# 1. Estimate overall probability of hrp2 deleted parasites being pfhrp3 deleted  -----------------
# -----------------------------------------------------#

# read in the raw data from WHO
hrpdat <- read.csv("analysis/data_raw/WHO_hrp2_hrp3_data.csv")

# get just the samples where all hrp2/3 testing was done
final <- hrpdat %>% filter(!is.na(HRP2_HRP3_DELETIONS_PERCENT)) %>%
  filter(HRP2_DELETION_TESTED_N == HRP3_DELETION_TESTED_N) %>%
  filter(HRP2_DELETION_TESTED_N == HRP2_HRP3_DELETIONS_TESTED_N)

# filter to just the cases with hrp2
final <- final %>% filter(HRP2_DELETION_POS_N > 0)

# also filter out where sampling was clearly not done on the same samples
# i.e. we shouldn't be able to have more samples with both deletions than
# ust one deletion assuming the deletion percentages are reported for all
# deletions, i.e. an hrp2 deletion percentage report is for any hrp2, i.e.
# with or without hrp3 and thus the double can't be more assuming this and
# assuming it was the same samples that are tested.
final <- final %>% filter(HRP2_HRP3_DELETIONS_PERCENT <= HRP2_DELETION_PERCENT)
final$continent <- as.factor(final$continent)

# now fit a beta binomial model for the mean probability that hrp2 is found with hrp3 deletions
mtmp <- function(prob, theta, continent) {
  if(continent == 0) {
    -sum(emdbook::dbetabinom(
      final$HRP2_HRP3_DELETIONS_POS_N,
      prob,
      final$HRP2_DELETION_POS_N,
      theta,
      log=TRUE
    ))
  } else {
    -sum(emdbook::dbetabinom(
      final$HRP2_HRP3_DELETIONS_POS_N[as.integer(final$continent) %in% continent],
      prob,
      final$HRP2_DELETION_POS_N[as.integer(final$continent) %in% continent],
      theta,
      log=TRUE
    ))
  }
}

# summary for specific continents and CI
mafrica <- bbmle::mle2(mtmp, start=list(prob=0.6,theta=1), fixed = list("continent" = 1))
bbmle::confint(mafrica)
mafrica@coef
bbmle::mle2(mtmp, start=list(prob=0.8,theta=9), fixed = list("continent" = 2))
bbmle::mle2(mtmp, start=list(prob=0.8,theta=9), fixed = list("continent" = 3))
bbmle::mle2(mtmp, start=list(prob=0.8,theta=9), fixed = list("continent" = 4))
# not enough data here
# bbmle::mle2(mtmp, start=list(prob=0.8,theta=9), fixed = list("continent" = 5))
m0 <- bbmle::mle2(mtmp, start=list(prob=0.8,theta=9), fixed = list("continent" = 0))
bbmle::confint(m0)
m0@coef

# -----------------------------------------------------#
# 2.  Plot our overall data and results ----------
# -----------------------------------------------------#

# add the cooccurrence onto our data
final$coocc <- final$HRP2_HRP3_DELETIONS_POS_N / final$HRP2_DELETION_POS_N
final$coocc_low <- as.numeric(Hmisc::binconf(final$HRP2_HRP3_DELETIONS_POS_N, final$HRP2_DELETION_POS_N)[,2])
final$coocc_high <- as.numeric(Hmisc::binconf(final$HRP2_HRP3_DELETIONS_POS_N, final$HRP2_DELETION_POS_N)[,3])
final$id <- seq_along(final$ISO_ctry)

# independence overall
gg1 <- ggplot(final, aes(x = coocc,y=as.factor(id), xmin = coocc_low, xmax = coocc_high, size = HRP2_DELETION_POS_N)) +
  geom_vline(xintercept = m0@coef[1], color = "red", linetype = "dashed") +
  geom_pointrange() +
  scale_size_binned(range = c(0,0.75), name = "Total *pfhrp2* deleted <br /> samples in survey:", breaks = c(0,10,50,100,150,200)) +
  labs(x = "Proportion of *pfhrp2-deleted* samples with *pfhrp3* deletions") +
  ylab("Survey") +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(axis.title.x = ggtext::element_markdown(size = 10),
        legend.title = ggtext::element_markdown(size = 10))
gg1
save_figs("hrp2_3_independence", gg1, width = 8, height = 6)

# -----------------------------------------------------#
# 3. Figure 1 - Plot our overall and cooccurrence findings with regard prevalence for Africa --------
# -----------------------------------------------------#

# independence assessed just for africa for these figures

# firstly group data by country for plot a)
grouped <- final %>% group_by(ISO_ctry, continent) %>%
  summarise(
    coocc = sum(HRP2_HRP3_DELETIONS_POS_N) / sum(HRP2_DELETION_POS_N),
    coocc_low = binom.test(sum(HRP2_HRP3_DELETIONS_POS_N), sum(HRP2_DELETION_POS_N), conf.level = 0.95 )$conf.int[1],
    coocc_high = binom.test(sum(HRP2_HRP3_DELETIONS_POS_N), sum(HRP2_DELETION_POS_N), conf.level = 0.95 )$conf.int[2],
    HRP2_DELETION_POS_N = sum(HRP2_DELETION_POS_N)) %>%
  filter(continent == "Africa") %>%
  mutate(country =  countrycode::countrycode(
    .data$ISO_ctry,"iso3c", "country.name",
    custom_match = c("COD" = "Democractic Republic of the Congo")))

gg2 <- grouped %>%
  ggplot(aes(x = coocc,y=country, group = (as.factor(country)),
             xmin = coocc_low, xmax = coocc_high, size = HRP2_DELETION_POS_N,
             color = (as.factor(country)))) +
  geom_vline(xintercept = mafrica@coef[1], color = "black", linetype = "dashed", size = 0.75) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  scale_x_continuous(labels = scales::percent) +
  scale_color_manual(name = "Country",
                     values = as.character(vapply(rev(pals::polychrome()), pastellize, character(1), 0.95))) +
  scale_size_continuous(name = "Total *pfhrp2* deleted <br /> samples in  survey:<br>", range = c(0,1), breaks = c(10, 20, 50,100,300), limits = c(0,300)) +
  labs(x = "Percentage of *pfhrp2* deleted samples <br /> with *pfhrp3* gene deletions") +
  labs(y = "WHO Threat Map Survey of *pfhrp2*/*3* deletions") +
  theme_bw(base_family = "Helvetica", base_size = 12) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(axis.title.x = ggtext::element_markdown(size = 12),
        legend.title = ggtext::element_markdown(size = 12),
        legend.text = ggtext::element_markdown(size = 9),
        axis.title.y = ggtext::element_markdown(size = 12)) +
  guides(color = guide_legend(ncol = 1)) +
  scale_y_discrete(limits = rev(levels(factor(grouped$country))))
gg2

# cooccurence by Africa by prev
gg_23_coocc <- final %>%
  filter(continent == "Africa") %>%
  ggplot(aes(x = prev, y=HRP2_HRP3_DELETIONS_POS_N/HRP2_DELETION_POS_N,
             pos=HRP2_HRP3_DELETIONS_POS_N, total=HRP2_DELETION_POS_N,
             size = HRP2_DELETION_POS_N)) +
  geom_smooth(
    method="glm",
    method.args=list(family="quasibinomial"),
    formula = cbind(pos, total-pos) ~ x,
    fullrange = TRUE,
    show.legend = FALSE,
  ) +
  geom_point(aes(size = HRP2_DELETION_POS_N)) +
  scale_x_log10(limits = c(0.0001, 1), labels = scales::percent) +
  scale_size_continuous(name = "Total *pfhrp2* deleted <br /> samples in  survey:<br>", range = c(0,5), breaks = c(10, 20, 50,100,300), limits = c(0,300)) +
  scale_y_continuous(labels = scales::percent)+
  labs(x = "Malaria Slide Prevalence 2-10 <br />") +
  labs(y = "Percentage of *pfhrp2* deleted samples <br /> with *pfhrp3* gene deletions") +
  theme_bw(base_family = "Helvetica", base_size = 12) +
  theme(axis.title.x = ggtext::element_markdown(size = 12),
        legend.title = ggtext::element_markdown(size = 12),
        legend.text = ggtext::element_markdown(size = 9),
        axis.title.y = ggtext::element_markdown(size = 12),
        legend.box = "horizontal")
gg_23_coocc

# prevalence of hrp3 by prev
gg_3_prev <- hrpdat %>%
  filter(continent == "Africa") %>%
  ggplot(aes(x = (prev), y=HRP3_DELETION_POS_N/HRP3_DELETION_TESTED_N,
             pos=HRP3_DELETION_POS_N, total=HRP3_DELETION_TESTED_N,
             size = HRP3_DELETION_TESTED_N)) +
  geom_smooth(
    method="glm",
    method.args=list(family="quasibinomial"),
    formula = cbind(pos, total-pos) ~ (x),
    fullrange = TRUE,
    show.legend = FALSE
  ) +
  geom_point(aes(size = HRP3_DELETION_TESTED_N)) +
  scale_x_log10(limits = c(0.0001, 1), labels = scales::percent) +
  scale_size_continuous(range=c(0,5),
                        name = "Total samples tested <br /> for *pfhrp3* deletions: <br />",
                        breaks = c(10, 200, 500,1000,2000), limits = c(0,2000)) +
  scale_y_continuous(labels = scales::percent)+
  labs(x = "Malaria Slide Prevalence 2-10 <br />") +
  labs(y = "Percentage of samples tested for *pfhrp3* <br /> with *pfhrp3* gene deletions") +
  theme_bw(base_family = "Helvetica", base_size = 12) +
  theme(axis.title.x = ggtext::element_markdown(size = 12),
        legend.title = ggtext::element_markdown(size = 12),
        legend.text = ggtext::element_markdown(size = 9),
        axis.title.y = ggtext::element_markdown(size = 12),
        legend.box = "horizontal")
gg_3_prev

figure1 <- cowplot::plot_grid(gg2,
                              cowplot::plot_grid(gg_23_coocc, gg_3_prev, ncol = 1, labels = c("B", "C"), align = "v"),
                              labels = c("A", ""), rel_widths = c(0.55,0.45),
                              ncol = 2
)
figure1
save_figs("hrp2_3_multi_independence_africa", figure1, width = 12, height = 8)

# -----------------------------------------------------#
# 4. Statistics on findings --------
# -----------------------------------------------------#

# Significance of prevalence
# Models for each
mod1 <- glm(cbind(HRP2_HRP3_DELETIONS_POS_N, HRP2_DELETION_TESTED_N - HRP2_HRP3_DELETIONS_POS_N) ~ log10(prev),
            data = hrpdat, family = binomial)

mod2 <- glm(cbind(HRP3_DELETION_POS_N, HRP3_DELETION_TESTED_N - HRP3_DELETION_POS_N) ~ log10(prev),
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
write_stargazer(star_model, "analysis/tables/pfhrp3_models.csv", spl = "\\s{2,}|\\)\\s{1,}-")

# -----------------------------------------------------#
# 5. Supp Version of Figure 1 across continents -----
# -----------------------------------------------------#

# independence assessed for all continents for these figures

# firstly group data by country for plot a)
grouped2 <- final %>% group_by(ISO_ctry, continent) %>%
  summarise(
    coocc = sum(HRP2_HRP3_DELETIONS_POS_N) / sum(HRP2_DELETION_POS_N),
    coocc_low = binom.test(sum(HRP2_HRP3_DELETIONS_POS_N), sum(HRP2_DELETION_POS_N), conf.level = 0.95 )$conf.int[1],
    coocc_high = binom.test(sum(HRP2_HRP3_DELETIONS_POS_N), sum(HRP2_DELETION_POS_N), conf.level = 0.95 )$conf.int[2],
    HRP2_DELETION_POS_N = sum(HRP2_DELETION_POS_N)) %>%
  mutate(country =  countrycode::countrycode(
    .data$ISO_ctry,"iso3c", "country.name",
    custom_match = c("COD" = "Democractic Republic of the Congo")))

gg2 <- grouped2 %>%
  filter(continent %in% c("Africa", "Americas", "Asia")) %>%
  ggplot(aes(x = coocc,y=country, group = (as.factor(country)),
             xmin = coocc_low, xmax = coocc_high, size = HRP2_DELETION_POS_N,
             color = (as.factor(country)))) +
  geom_vline(xintercept = mafrica@coef[1], color = "black", linetype = "dashed", size = 0.75) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  scale_x_continuous(labels = scales::percent) +
  scale_color_manual(name = "Country",
                     values = as.character(vapply(rev(pals::polychrome()), pastellize, character(1), 0.95))) +
  scale_size_continuous(name = "Total *pfhrp2* deleted <br /> samples in  survey:<br>", range = c(0,1), breaks = c(10, 20, 50,100,400), limits = c(0,400)) +
  labs(x = "Percentage of *pfhrp2* deleted samples <br /> with *pfhrp3* gene deletions") +
  labs(y = "WHO Threat Map Survey of *pfhrp2*/*3* deletions") +
  theme_bw(base_family = "Helvetica", base_size = 12) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(axis.title.x = ggtext::element_markdown(size = 12),
        legend.title = ggtext::element_markdown(size = 12),
        legend.text = ggtext::element_markdown(size = 9),
        axis.title.y = ggtext::element_markdown(size = 12)) +
  guides(color = guide_legend(ncol = 2)) +
  scale_y_discrete(limits = rev(levels(factor(grouped2$country))))
gg2

# cooccurence
gg_23_coocc <- final %>%
  filter(continent %in% c("Africa", "Americas", "Asia")) %>%
  ggplot(aes(x = prev, y=HRP2_HRP3_DELETIONS_POS_N/HRP2_DELETION_POS_N,
             color = continent,
             pos=HRP2_HRP3_DELETIONS_POS_N, total=HRP2_DELETION_POS_N,
             size = HRP2_DELETION_POS_N)) +
  geom_smooth(
    method="glm",
    method.args=list(family="quasibinomial"),
    formula = cbind(pos, total-pos) ~ x,
    fullrange = TRUE,
    show.legend = FALSE,
  ) +
  geom_point(aes(size = HRP2_DELETION_POS_N)) +
  scale_x_log10(limits = c(0.0001, 1), labels = scales::percent) +
  scale_size_continuous(name = "Total *pfhrp2* deleted <br /> samples in  survey:<br>", range = c(0,5), breaks = c(10, 20, 50,100,300), limits = c(0,300)) +
  scale_color_viridis_d(name = "Continent") +
  scale_y_continuous(labels = scales::percent)+
  labs(x = "Malaria Slide Prevalence 2-10 <br />") +
  labs(y = "Percentage of *pfhrp2* deleted samples <br /> with *pfhrp3* gene deletions") +
  theme_bw(base_family = "Helvetica", base_size = 12) +
  theme(axis.title.x = ggtext::element_markdown(size = 12),
        legend.title = ggtext::element_markdown(size = 12),
        legend.text = ggtext::element_markdown(size = 9),
        axis.title.y = ggtext::element_markdown(size = 12),
        legend.box = "horizontal")
gg_23_coocc

# hrp3 prevalence
gg_3_prev <- hrpdat %>%
  filter(continent %in% c("Africa", "Americas", "Asia")) %>%
  ggplot(aes(x = (prev), y=HRP3_DELETION_POS_N/HRP3_DELETION_TESTED_N, color = continent,
             pos=HRP3_DELETION_POS_N, total=HRP3_DELETION_TESTED_N,
             size = HRP3_DELETION_TESTED_N)) +
  geom_smooth(
    method="glm",
    method.args=list(family="quasibinomial"),
    formula = cbind(pos, total-pos) ~ (x),
    fullrange = TRUE,
    show.legend = FALSE
  ) +
  geom_point(aes(size = HRP3_DELETION_TESTED_N)) +
  scale_x_log10(limits = c(0.0001, 1), labels = scales::percent) +
  scale_size_continuous(range=c(0,5),
                        name = "Total samples tested <br /> for *pfhrp3* deletions: <br />",
                        breaks = c(10, 200, 500,1000,2000), limits = c(0,2000)) +
  scale_color_viridis_d(name = "Continent") +
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
                              cowplot::plot_grid(gg_23_coocc, gg_3_prev, ncol = 2, labels = c("B", "C"), align = "v"),
                              labels = c("A", ""),
                              ncol = 1, rel_heights = c(0.55,0.45)
)
figure1

save_figs("hrp2_3_multi_independence", figure1, width = 18, height = 14)

# --------------------------------------------------- #
# 6. 2 x 2 table for cooccurrence -------
# --------------------------------------------------- #
two_by_two <- final %>% #filter(continent == "Africa") %>%
  ungroup() %>%
  na.omit %>%
  mutate(dd = HRP2_HRP3_DELETIONS_POS_N,
         dp = HRP2_DELETION_POS_N-dd,
         pd = HRP3_DELETION_POS_N-dd,
         pp = HRP2_DELETION_TESTED_N - dd - dp - pd) %>%
  summarise(across(dd:pp, sum))

tbt_table <- two_by_two %>% as.numeric() %>% matrix(ncol = 2, nrow = 2)
chisq.test(tbt_table)

two_by_two <- two_by_two %>%
  mutate(n = sum(c(dd,dp,pd,pp))) %>%
  mutate(dd = dd/n,
         dp = dp/n,
         pd = pd/n,
         pp = pp/n)

# calculate our linkage disequilbrium value
p_two_d <- (two_by_two$dd+two_by_two$dp)
p_three_d <- (two_by_two$dd+two_by_two$dp)
D <- two_by_two$dd - (p_two_d*p_three_d)
D_prime <- D/min(c(p_two_d*(1-p_two_d),(p_three_d)*(1-p_three_d)))
R2 <- (D^2)/(p_two_d*(1-p_two_d)*(p_three_d)*(1-p_three_d))
