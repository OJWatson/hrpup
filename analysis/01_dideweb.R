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

#packages <- list(loaded = "OJWatson/hrp2malaRia")
packages <- c("OJWatson/hrp2malaRia")
ctx <- context::context_save(
  path = context_name,
  package_sources = conan::conan_sources(
    packages = packages
  )
)

# set up a specific config for here as we need to specify the large RAM nodes
config <- didehpc::didehpc_config(use_workers = TRUE)
config$resource$parallel <- "FALSE"
config$resource$type <- "Cores"

# Configure the Queue
obj <- didehpc::queue_didehpc(ctx, config = config, provision = "verylazy")

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
workers <- obj$submit_workers(200)
nmf_ranges <- unique(pl$nmf.multiplier)

for(nmf in nmf_ranges){

  # Submission Lists -----
  paramList <- list()
  plnmf <- pl %>% filter(nmf.multiplier == nmf)
  for(i in seq_len(nrow(pl))){

    paramList[[i]] <- list(EIR=plnmf$EIR[i]/365,
                           ft=plnmf$ft[i],
                           nmf.multiplier=plnmf$nmf[i],
                           include.nmf = TRUE,
                           strains.0=20,
                           N=100000,
                           time.step=1,
                           years=40,
                           storage=30,
                           max.age=100,
                           d.CM=67.7)
  }
  grp_list <- list()
  # submit grids
  for(i in 1:5) {

    try_fail_catch(
      grp_list[[i]] <- obj$lapply(
        X = paramList, timeout=0,
        FUN = function(x){
          library(hrp2malaRia)
          return(do.call(hrp2malaRia::hrp2_Simulation, x))
        },
        name = paste0("long_grid_setup_nmf_", nmf,"_",i), overwrite = TRUE)
    )

  }

}

## ------------------------------------
## 3. Continuation
## ------------------------------------

pl2 <- readRDS("analysis/data_derived/param_grid.rds")
paramList_list_final <- list()
for(nmf in seq_along(nmf_ranges)) {

  # which nmf is this
  nmf_i <- nmf
  nmf <- nmf_ranges[nmf]
  grp_start_list <- lapply(paste0("long_grid_setup_nmf_", nmf,"_",1:5), function(x){obj$task_bundle_get(x)})

  plnmf <- pl2 %>% filter(nmf.multiplier == nmf)
  paramList_list <- list()

  # Admin continuations  -----
  for(j in seq_along(grp_start_list)){

    paramList_list[[j]] <- list()
    previous_pars <- as.data.frame(do.call(rbind, lapply(grp_start_list[[j]]$X, as.data.frame)))

    for(i in seq_len(nrow(plnmf))){

      # fill our param list with the required params
      paramList_list[[j]][[i]] <- list(EIR=plnmf$EIR[i]/365,
                                       ft=plnmf$ft[i],
                                       microscopy.use=plnmf$microscopy.use[i],
                                       rdt.nonadherence=plnmf$rdt.nonadherence[i],
                                       fitness=plnmf$fitness[i],
                                       rdt.det=plnmf$rdt.det[i],
                                       nmf.multiplier=plnmf$nmf.multiplier[i],
                                       N=100000,
                                       time.step=1,
                                       years=40,
                                       storage=30,
                                       max.age=100,
                                       d.CM=67.7)

      # And the id of the previous run this relates to
      run_i <- which(previous_pars$EIR == plnmf$EIR[i]/365 &
                       previous_pars$ft == plnmf$ft[i] &
                       previous_pars$nmf.multiplier == plnmf$nmf.multiplier)
      paramList_list[[j]][[i]]$ID <- grp_start_list[[j]]$ids[run_i]
      paramList_list[[j]][[i]]$root <- context_name

      # And the adjustments to the storage and length of run time
      paramList_list[[j]][[i]]$just_storage_results <- TRUE
      paramList_list[[j]][[i]]$years <- 20
      paramList_list[[j]][[i]]$desired.freq <- 0.06
      paramList_list[[j]][[i]]$include.nmf <- TRUE

    }

  }

  paramList_list_final[[nmf_i]] <- paramList_list

}

saveRDS(paramList_list_final, "analysis/data_derived/final_param_grid.rds")

# -------------------------------------------
# Now submit our new runs to the cluster
# -------------------------------------------

paramList_list_final <- readRDS("analysis/data_derived/final_param_grid.rds")

# for now just do the list related to nmf = 1
for(nmf in seq_along(nmf_ranges)[2]){

  nmf_i <- nmf
  nmf <- nmf_ranges[nmf_i]

  paramList_list <- paramList_list_final[[nmf_i]]
  grp_list_final <- list()
  for(i in seq_along(paramList_list)[1:5]){

    try_fail_catch(
      grp_list_final[[i]] <- obj$lapply(
        X = paramList_list[[i]], timeout=0,
        FUN = function(x){
          library(hrp2malaRia)
          return(do.call(hrp2malaRia::hrp2_Simulation, x))
        },
        name = paste0("long_grid_continuation_nmf_",nmf,"_", i), overwrite = TRUE)
    )

  }

}




