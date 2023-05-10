`%||%` <- function(a, b) { # nolint
  if (is.null(a)) b else a
}


#' Save Figures
#'
#' @param name Name of figure
#' @param fig ggplot or similar figure object
#' @param width Width of plot in inches. Default = 6
#' @param height Height of plot in inches. Default = 6
#' @param plot_dir Plotting directory. Defaults to "analysis/plots"
#' @param pdf_plot Logical for plotting pdf too. Default = TRUE
#' @param font_family If specified, sets all font family. Default = NULL
#' @importFrom grDevices dev.off pdf
save_figs <- function(name,
                      fig,
                      width = 6,
                      height = 6,
                      plot_dir = file.path(here::here(), "analysis/plots"),
                      pdf_plot = TRUE,
                      font_family = NULL) {

  if(!is.null(font_family)) {
    fig <- fig + theme(text = element_text(family = font_family))
  }

  dir.create(plot_dir, showWarnings = FALSE)
  fig_path <- function(name) {paste0(plot_dir, "/", name)}

  cowplot::save_plot(filename = fig_path(paste0(name,".png")),
                     plot = fig,
                     base_height = height,
                     base_width = width)

  if(pdf_plot) {
  pdf(file = fig_path(paste0(name,".pdf")), width = width, height = height)
  print(fig)
  dev.off()
  }

}

#' Create path relative to root of project
#'
#' @param path Path to be appended to project root
cp_path <- function(path) {

  file.path(here::here(), path)

}


