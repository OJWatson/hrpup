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
                      font_family = "Helvetica") {

  if(!is.null(font_family)) {
    fig <- fig + ggplot2::theme(text = ggplot2::element_text(family = font_family))
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


#' Write stargazer model output to file
#'
#' @param x Stargazer model output as text
#' @param path Path for stargazer to be written to
#' @param cls Number of columns in final stargazer table. Default = NULL,
#'   which works it out based on maximum split size
#' @param spl regex string for splitting columns.
write_stargazer <- function(x, path, cls = NULL, spl = "\\s{2,}") {

  splitup <- vapply(x, strsplit, spl, FUN.VALUE = vector("list", 1))

  if(is.null(cls)) {
    cls <- max(lengths(splitup))
  }
  tbl <- do.call(rbind, splitup[lengths(splitup) == cls])
  rownames(tbl) <- NULL
  colnames(tbl) <- tbl[1,]
  utils::write.csv(tbl[-1,], path, row.names = FALSE)

}


#' Create pastel color palatte
#'
#' @param x hsv value
#' @param p Level of pastellisation. Nnumber in [0,1]
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
