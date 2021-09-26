# determine the optimal block penalty and window size 
# run like Rscript 02_optimize_block_window_ldhelmet.R simname

library("tidyverse")
library("parallel")

args <- commandArgs(trailingOnly = TRUE)
#args <- list("test")

sim_name <- args[[1]]

# list of window sizes
#window_sizes <- seq(from = 50, to = 200, by = 50)

window_sizes <- c(20, 50, 100, 1000)

# list of block penalties 
#block_penalties <- seq(from = 250, to = 2000, by = 250)

block_penalties <- c(10, 20, 50, 100, 500, 1000)

# cores 
cores_per_job <- 5
num_jobs <- 4

# function for running multiple block levels

run_ldhelmet_with_block <- function(block_size, sim_name, wind_size, cores){
  
  if(!dir.exists(paste0(sim_name, "_wind_", wind_size, "_block_", block_size))){

    ld_helm_call <- paste0("bash pipeline/run_all.bash ", sim_name, " ", wind_size, " ", block_size, " ", cores)
    system(ld_helm_call)
 } 
}

for (i in window_sizes){
  
    mclapply(block_penalties, run_ldhelmet_with_block, sim_name = sim_name, wind_size = i, cores = cores_per_job, mc.cores = num_jobs, mc.preschedule = FALSE) 
}

