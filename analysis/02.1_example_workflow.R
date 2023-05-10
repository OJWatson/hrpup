# Make plots to help demonstrate what selection coefficients are
# --------------------------------------------------------------- #


# 1. Get our main data objects
# ----------------------------------------#

# grab our processed data
data <- readRDS(file.path(here::here("analysis/data_derived/model_s.rds")))

# grab all our raw data
dat <- lapply(1:5, function(x) {readRDS(file.path(
  here::here(),
  paste0(
    "analysis/data_derived/sims/true_grid_continuation_nmf_1_", x, ".rds"
  )
)) })

# grab our equivalent pl df
testna <- lapply(seq_along(dat), function(x){

  testna <- do.call(rbind, lapply(dat[[x]]$pl, as.data.frame)) %>%
    as.data.frame
  testna <- testna %>%
    mutate(int = interaction(ft, microscopy.use, rdt.nonadherence, fitness, rdt.det, sep = "_"))
  testna

})

# 2. Decide on what we are plotting
# ----------------------------------------#
# what ones do we want to plot
to_extract <- data[data$EIR %in% unique(data$EIR)[c(3,4,5)] &
                     data$ft == unique(data$ft)[5] &
                     data$microscopy.use == unique(data$microscopy.use)[2] &
                     data$rdt.nonadherence == unique(data$rdt.nonadherence)[3] &
                     data$fitness == unique(data$fitness)[2] &
                     data$rdt.det == unique(data$rdt.det)[2],] %>%
  mutate(int_simple = interaction(ft, microscopy.use, rdt.nonadherence, fitness, rdt.det, sep = "_"))


make_plot <- function(to_extract, colors = c("red","green","blue"), lab = 0, size = 3.88) {

  # grab the trjectories from the raw data
  trajs <- lapply(seq_along(dat), function(x){
    dat[[x]]$res[which(testna[[x]]$int %in% to_extract$int_simple[1] & testna[[x]]$EIR %in% to_extract$EIR[1:3])]
  })
  trajs <- do.call(c, trajs)

  # and turn into a data frame for plotting
  res_all <- vector("list", length(trajs))
  for(i in seq_along(trajs)){

    # filter to after selection was implemented
    df2 <- trajs[[i]] %>% filter(S.Times > 365*38.05)
    df2 <- df2 %>% select(S.Times, S.Micro.210, Percentange_Clin_Mono)
    res_all[[i]] <- cbind(df2, to_extract[i,])

  }
  res_all <- do.call(rbind, res_all)

  binomial_smooth <- function(...) {
    geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
  }

  level_neat <- function(x) {
    lapply(split(x, nchar(x)), as.factor) %>% lapply(levels) %>% unlist
  }

  # turn this is into a data frame for easy plotting
  gg_res <- res_all %>% filter(S.Times > 365*40.05) %>%
    rename(t = S.Times, y = Percentange_Clin_Mono) %>%
    mutate(Prevalence = round(Micro.2.10/0.05)*0.05) %>%
    mutate(Prevalence = scales::percent_format(accuracy = 1)(Prevalence)) %>%
    mutate(Prevalence = factor(Prevalence, levels = level_neat(Prevalence)))

  gg_res$smooth <- glm(y ~ t*Prevalence, data = gg_res, family = "binomial") %>%
    predict(type = "response")

  gg_res$smooth_scale <- glm(y ~ t*Prevalence, data = gg_res, family = "binomial") %>%
    predict()

  p <- gg_res %>%
    ggplot(aes((t/365)-39.9, y, group = ID, color = Prevalence)) +
    geom_line(alpha = 0.5) +
    geom_line(aes(y = smooth)) +
    ggpubr::theme_pubclean() +
    theme(axis.line = element_line(), panel.grid.major.x = element_blank(),
          legend.key = element_rect(fill = "white")) +
    scale_x_continuous(expand = c(0,0), limits = c(0,20)) +
    scale_y_continuous(expand = c(0,0),  labels = scales::percent) +
    xlab("Years") +
    ylab("Percentage False-Negative HRP2-RDTs amongst \nclinical infections due to pfhrp2/3 deletions") +
    scale_color_manual(values = colors, name = "Malaria Prevalence 2-10: ")

  gg_res$label <- paste0("y==s[",as.integer(gg_res$Prevalence)+lab,"]*t+c[",as.integer(gg_res$Prevalence)+lab,"]")

  p2 <- gg_res %>%
    ggplot(aes((t/365)-39.9, log(y / (1-y)), group = ID, color = Prevalence)) +
    geom_line(alpha = 0.5) +
    geomtextpath::geom_labelline(aes(y = smooth_scale,label = label),size = size,
                                 parse = TRUE, show.legend = FALSE, label.padding = 0.15, hjust = 0.8,straight = TRUE) +
    ggpubr::theme_pubclean() +
    theme(axis.line = element_line(), panel.grid.major.x = element_blank(),
          legend.key = element_rect(fill = "white")) +
    scale_x_continuous(expand = c(0,0), limits = c(0,20)) +
    scale_y_continuous(expand = c(0,0)) +
    xlab("Years") +
    ylab("log(y / 1-y)") +
    scale_color_manual(values = colors, name = "Malaria Prevalence 2-10: ")



  example <- cowplot::plot_grid(
    p,
    grid::linesGrob(
      c(0.2, 1),
      c(2.5, 2.5),
      default.units = "inch",
      gp = grid::gpar(fill = "black"),
      arrow = grid::arrow(
        ends = "last",
        type = "closed",
        angle = 30,
        length = unit(0.5, "cm")
      )
    ),
    p2,
    ncol = 3,
    rel_widths = c(1, 0.2, 1)
  )

  vars <- names(to_extract)[2:6]
  names(vars) <- c(
    "Effective Treatment Seeking",
    "Microscopy Use",
    "RDT Nonadherence",
    "Comparative Fitness",
    "HRP3 Cross Reactivity"
  )

  cowplot::plot_grid(
    ggplot() + ggtitle(paste(
      names(vars),
      scales::percent(as.numeric(round(to_extract[1, 2:6]/0.05)*0.05)),
      sep = ": ",
      collapse = ", "
    ))  ,
    example,
    rel_heights = c(0.1, 1),
    ncol = 1) +
    theme(plot.background = element_rect(fill = "white", colour = "white"))

}

pal0 <- RColorBrewer::brewer.pal(9,"RdGy")[c(9,10,11)]
pal1 <- RColorBrewer::brewer.pal(9,"BuGn")[c(5,7,9)]
pal2 <- RColorBrewer::brewer.pal(9,"BuPu")[c(5,7,9)]
pal3 <- RColorBrewer::brewer.pal(9,"GnBu")[c(5,7,9)]
pal4 <- RColorBrewer::brewer.pal(9,"OrRd")[c(5,7,9)]
pal5 <- RColorBrewer::brewer.pal(9,"YlOrBr")[c(5,7,9)]

to_extract <- data[data$EIR %in% unique(data$EIR)[c(3,4,5)] &
                     data$ft == unique(data$ft)[5] &
                     data$microscopy.use == unique(data$microscopy.use)[2] &
                     data$rdt.nonadherence == unique(data$rdt.nonadherence)[2] &
                     data$fitness == unique(data$fitness)[2] &
                     data$rdt.det == unique(data$rdt.det)[2],] %>%
  mutate(int_simple = interaction(ft, microscopy.use, rdt.nonadherence, fitness, rdt.det, sep = "_"))
row0 <- make_plot(to_extract, pal0, size = 8)
# save a
save_figs("example_a", row0, width = 13, height = 6)
row0 <- make_plot(to_extract, pal0)

to_extract <- data[data$EIR %in% unique(data$EIR)[c(3,4,5)] &
                     data$ft == unique(data$ft)[6] &
                     data$microscopy.use == unique(data$microscopy.use)[2] &
                     data$rdt.nonadherence == unique(data$rdt.nonadherence)[2] &
                     data$fitness == unique(data$fitness)[2] &
                     data$rdt.det == unique(data$rdt.det)[2],] %>%
  mutate(int_simple = interaction(ft, microscopy.use, rdt.nonadherence, fitness, rdt.det, sep = "_"))
row1 <- make_plot(to_extract, pal1, lab = 3)


to_extract <- data[data$EIR %in% unique(data$EIR)[c(3,4,5)] &
                     data$ft == unique(data$ft)[5] &
                     data$microscopy.use == unique(data$microscopy.use)[1] &
                     data$rdt.nonadherence == unique(data$rdt.nonadherence)[2] &
                     data$fitness == unique(data$fitness)[2] &
                     data$rdt.det == unique(data$rdt.det)[2],] %>%
  mutate(int_simple = interaction(ft, microscopy.use, rdt.nonadherence, fitness, rdt.det, sep = "_"))
row2 <- make_plot(to_extract, pal2, lab = 6)

to_extract <- data[data$EIR %in% unique(data$EIR)[c(3,4,5)] &
                     data$ft == unique(data$ft)[5] &
                     data$microscopy.use == unique(data$microscopy.use)[2] &
                     data$rdt.nonadherence == unique(data$rdt.nonadherence)[1] &
                     data$fitness == unique(data$fitness)[2] &
                     data$rdt.det == unique(data$rdt.det)[2],] %>%
  mutate(int_simple = interaction(ft, microscopy.use, rdt.nonadherence, fitness, rdt.det, sep = "_"))
row3 <- make_plot(to_extract, pal3, lab = 9)

to_extract <- data[data$EIR %in% unique(data$EIR)[c(3,4,5)] &
                     data$ft == unique(data$ft)[5] &
                     data$microscopy.use == unique(data$microscopy.use)[2] &
                     data$rdt.nonadherence == unique(data$rdt.nonadherence)[2] &
                     data$fitness == unique(data$fitness)[3] &
                     data$rdt.det == unique(data$rdt.det)[2],] %>%
  mutate(int_simple = interaction(ft, microscopy.use, rdt.nonadherence, fitness, rdt.det, sep = "_"))
row4 <- make_plot(to_extract, pal4, lab = 12)


to_extract <- data[data$EIR %in% unique(data$EIR)[c(3,4,5)] &
                     data$ft == unique(data$ft)[5] &
                     data$microscopy.use == unique(data$microscopy.use)[2] &
                     data$rdt.nonadherence == unique(data$rdt.nonadherence)[2] &
                     data$fitness == unique(data$fitness)[2] &
                     data$rdt.det == unique(data$rdt.det)[3],] %>%
  mutate(int_simple = interaction(ft, microscopy.use, rdt.nonadherence, fitness, rdt.det, sep = "_"))
row5 <- make_plot(to_extract, pal5, lab = 15)


# save a
save_figs("example_a", row0, width = 13, height = 6)

# save b
save_figs("example_b",
          cowplot::plot_grid(row0,row1,row2,row3,row4,row5, ncol = 1, scale = 0.95) +
            theme(plot.background = element_rect(fill = "white", colour = "white")),
          width = 15, height = 36)
