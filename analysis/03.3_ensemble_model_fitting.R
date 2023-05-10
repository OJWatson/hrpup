# N.B. For reference, an XGBOOST model was first fit and then after having
# chosen hyperparameters for that and testing it, it was clear a more continuous
# model was required to support. multivariate additive regression spline model
# was subsequently chosen and then fit. After having fit that we selected the last
# model based on dissimilarity of the different model types as below.

# -------------------------------------------------------
# 0. Pick final model for our ensemble
# --------------------------------------------------------

# Have decided that XGB and MDA work from other scripts.
# Let's pick a 3rd model

# final model selection for ensembles
tag <- read.csv("https://topepo.github.io/caret/tag_data.csv", row.names = 1)
tag <- as.matrix(tag)

## Select only models for regression
regModels <- tag[tag[,"Regression"] == 1,]

all <- 1:nrow(regModels)
## Seed the analysis with the SVM model
start <- grep("xgbTree|bagEarth)", rownames(regModels))
pool <- all[all != start]

## Select 4 model models by maximizing the Jaccard
## dissimilarity between sets of models
nextMods <- caret::maxDissim(regModels[start,,drop = FALSE],
                      regModels[pool, ],
                      method = "Jaccard",
                      n = 4)

# what have they suggested
rownames(regModels)[c(start, nextMods)]

# let's go with a brnn for our final ensemble. First let's check it also works well

# --------------------------------------------------------
# 1. Perform brnn model fitting
# --------------------------------------------------------

library(tidyverse)
testna <- readRDS(file.path(here::here(), "analysis/data_derived/model_s.rds"))
testna <- testna %>% select(c("Micro.2.10","ft", "microscopy.use", "rdt.nonadherence", "fitness", "rdt.det", "s"))

# Split the data into training and testing sets
set.seed(123)
test <- testna %>% na.omit()
train_indices <- sample(nrow(test), nrow(test) * 0.75)
train <- test[train_indices, ]
test <- test[-train_indices, ]

# 1 is observed through pdp to enforce monotonic outputs that are needed with respect to microscopy prevalence
fit_brnn <- brnn::brnn(y = train$s, x = as.matrix(train %>% select(-s)), neurons = 1, verbose = FALSE)

# Test the fit
pred_s <- predict(fit_brnn, as.matrix(test %>% select(-s)))
mse <- mean((pred_s - test$s)^2)
rmse <- sqrt(mean((test$s - pred_s) ^ 2))
me <- mean((pred_s - test$s))
plot(pred_s, test$s)
abline(0, 1, col = "red")

# PDP plots to check
plots <- lapply(names(test %>% select(-s)), function(x){
  partial(fit_brnn, pred.var = x,  train = train %>% select(-s), grid.resolution = 20) %>%
    ggplot2::autoplot()
})
brnn_pdp <- cowplot::plot_grid(plotlist = plots)


# Add this to our model selection object
hrp2_mod <- readRDS(file.path(here::here("analysis/data_derived/ensemble_selection_model.rds")))
hrp2_mod$add_model(fit_brnn, "brnn")
hrp2_mod$add_model_weight(rmse, "brnn")

# Save our model out again
saveRDS(hrp2_mod, file.path(here::here("analysis/data_derived/ensemble_selection_model.rds")))

# --------------------------------------------------------
# 2. Partial dependence plot
# --------------------------------------------------------

mods <- hrp2_mod$get_models()
vars <- names(test %>% select(-s))
names(vars) <- c("Microscopy 2-10", "Effective Treatment Seeking", "Microscopy Use",
                 "RDT Nonadherence\n", "Comparative Fitness of \npfhrp2 deleted parasites",
                 "P(RDT+ve if only infected with \npfhp2 deleted parasites)")

pdp_dat <- lapply(seq_along(vars), function(x){
  dat <- do.call(rbind, lapply(seq_along(mods), function(i){
  partial(mods[[i]], pred.var = as.character(vars[x]), train = train %>% select(-s), grid.resolution = 20, type = "regression") %>%
      rename(s = yhat) %>%
      mutate(model = names(mods)[i]) %>%
      setNames(c("x", "s", "model"))
  }))
  dat %>% ggplot(aes(x,s,color = model)) +
    geom_line() +
    xlab(names(vars)[x]) +
    ylab("Selection Coeffiecient") +
    ggpubr::theme_pubclean(base_size = 14) +
    theme(axis.line = element_line()) +
    scale_color_viridis_d(name = "Model", end = 0.8) +
    theme(legend.key = element_rect(fill = "white"))

})

ggpdp <- cowplot::plot_grid(plotlist = pdp_dat)
save_figs("partial_dependence", ggpdp, 14, 8)


