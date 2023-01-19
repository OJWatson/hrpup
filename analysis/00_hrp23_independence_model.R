hrpdat <- readxl::read_xlsx("analysis/data_raw/HRP2_map_13.1.23.xlsx") %>%
  rename(WHO_HRP2_HRP3_DELETIONS_PERCENT = WHO_HPR2_HRP3_DELETIONS_PERCENT,
         WHO_HRP2_HRP3_DELETIONS_TESTED_N = WHO_HPR2_HRP3_DELETIONS_TESTED_N,
         WHO_HRP2_HRP3_DELETIONS_POS_N = WHO_HPR2_HRP3_DELETIONS_POS_N)
hrpdat <- hrpdat %>% select(ISO_ctry, WHO_LATITUDE, WHO_LONGITUDE, WHO_YEAR_START, WHO_YEAR_END, WHO_PATIENT_TYPE,
                 MEAN_AGE, WHO_RDT_TYPE, matches("WHO_HRP"))

hrpdat <- hrpdat %>% filter(WHO_HRP2_ONLY == 0)

# convert NRs
conv <- names(which(unlist(lapply(hrpdat, function(x){sum(x == "NR", na.rm = TRUE)>0}))))
hrpdat <- hrpdat %>%
  mutate(across(all_of(conv), .fns = ~as.numeric(replace(.x, .x == "NR", NA))))

# convert numerics
hrpdat <- hrpdat %>%
  mutate(across(WHO_LATITUDE, WHO_LONGITUDE, WHO_YEAR_START, WHO_YEAR_END, WHO_PATIENT_TYPE,MEAN_AGE,
                .fns = as.numeric)) %>%
  mutate(across(all_of(grep("WHO_HRP", names(hrpdat), value = TRUE)),
                .fns = as.numeric))

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

mtmp <- function(prob,theta) {
  -sum(emdbook::dbetabinom(final$WHO_HRP2_HRP3_DELETIONS_POS_N, prob, final$WHO_HRP2_DELETION_POS_N,theta, log=TRUE))
}

# summary gives 66% probability
m0 <- bbmle::mle2(mtmp, start=list(prob=0.6,theta=9))
summary(m0)

# plot what that looks like
final$coocc <- final$WHO_HRP2_HRP3_DELETIONS_PERCENT / final$WHO_HRP2_DELETION_PERCENT
final$coocc_low <- as.numeric(Hmisc::binconf(final$WHO_HRP2_HRP3_DELETIONS_POS_N, final$WHO_HRP2_DELETION_POS_N)[,2])
final$coocc_high <- as.numeric(Hmisc::binconf(final$WHO_HRP2_HRP3_DELETIONS_POS_N, final$WHO_HRP2_DELETION_POS_N)[,3])
final$id <- seq_along(final$ISO_ctry)

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

save_figs("hrp2_3_independence", gg1, width = 6, height = 8)
