# determine the optimal block penalty and window size 

library("tidyverse")

args <- commandArgs(trailingOnly = TRUE)
args <- list("test")

sim_name <- args[[1]]

# list of window sizes
window_sizes <- seq(from = 10, to = 100, by = 10)

# list of block penalties 
block_penalties <- seq(from = 10, to = 100, by = 10)

for (i in window_sizes){
  
  for (j in block_penalties){
    
    ld_helm_call <- paste0("bash run_all.bash ", sim_name, " ", i, " ", j)
    print(ld_helm_call)
    #system()
    
  }
}

