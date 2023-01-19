## ------------------------------------
## 1. Setting up a cluster configuration
## ------------------------------------

# Setting Up Cluster From New

# Log in to didehpc
credentials = "C:/Users/ow813/.smbcredentials"
options(didehpc.cluster = "fi--didemrchnb",
        didehpc.username = "ow813")

# not if T is not mapped then map network drive
didehpc::didehpc_config_global(temp=didehpc::path_mapping("tmp",
                                                          "T:",
                                                          "//fi--didef3.dide.ic.ac.uk/tmp",
                                                          "T:"),
                               home=didehpc::path_mapping("OJ",
                                                          "L:",
                                                          "//fi--didenas5/malaria",
                                                          "L:"),
                               credentials=credentials,
                               cluster = "fi--didemrchnb")

# Creating a Context
context_name <- "scripts/context"

ctx <- context::context_save(
  path = context_name,
  package_sources = conan::conan_sources(
    packages = c("OJWatson/hrp2malaRia", "rrq"),
    repos = "https://mrc-ide.github.io/drat/"
  )
)

# set up a specific config for here as we need to specify the large RAM nodes
config <- didehpc::didehpc_config(use_workers = TRUE)
config$resource$parallel <- "FALSE"
config$resource$type <- "Cores"

# Configure the Queue
obj <- didehpc::queue_didehpc(ctx, config = config)

## safe submission ======
try_fail_catch <- function(expr, attempts = 3){
  r <- NULL
  attempt <- 1
  while( is.null(r) && attempt <= 3 ) {
    attempt <- attempt + 1
    try(
      r <- eval(expr)
    )
  }

}

## ------------------------------------
## 2. Steady State Runs
## ------------------------------------

pl <- readRDS("analysis/data_derived/param_start.rds")

# Submission Lists -----
paramList <- list()
for(i in seq_len(nrow(pl))){
  paramList[[i]] <- list(EIR=pl$EIR[i]/365,
                         ft=pl$ft[i],
                         strains.0=20,
                         N=100000,
                         time.step=1,
                         years=30,
                         storage=30,
                         max.age=100,
                         d.CM=67.7)
}

# submit grids
for(i in 1:5) {

  try_fail_catch(
    grp <- obj$lapply(
      X = paramList, timeout=0,
                       FUN = function(x){
                         return(hrp2_Simulation(EIR=x$EIR,strains.0 = x$strains.0, N = x$N,
                                                time.step = x$time.step, years = x$years,
                                                storage = x$storage, max.age = x$max.age,
                                                d.CM = x$d.CM, ft = x$ft,
                         ))
                       },
                       name = paste0("grid_setup_", i))
  )

}
