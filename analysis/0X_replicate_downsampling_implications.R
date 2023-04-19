mu <- 0.5
sd <- 0.2

# find the difference between the error when estimating sd
# when using 5 vs 100 replicates

# generate 1000 samples of mae with reps 5, 100

# and span across a range of true mu and sd values
param_grid <- expand.grid(mu = seq(0,1,0.05), sd = seq(0,0.3,0.01))

diff_func <- function(params){

  mu <- params$mu
  sd <- params$sd

  top <- mean(replicate(1000, mean(abs(sd(rnorm(100, mu, sd)) - sd))))
  bottom <- mean(replicate(1000, mean(abs(sd(rnorm(5, mu, sd)) - sd))))
  diff <- bottom - top
  rel <- bottom/top
  return(data.frame(mu = mu, sd = sd, top = top, bottom = bottom, diff = diff, rel = rel))
}

param_out <- do.call(rbind, lapply(seq_len(nrow(param_grid)), function(x) {diff_func(param_grid[x,])} ))
param_out <- na.omit(param_out)

# Okay so by using only 5 reps vs 100 we over estimate the sd by 5x
param_out %>% ggplot(aes(mu,sd, fill = rel)) + geom_tile()

