# create and launch slim simulations for range of parameter combinations
# KMS AUG 2020

library("tidyverse")
library("parallel")
library("inline")
options(scipen = 999)


# launches a single slim script with the given parameters
# dry run just prints the commands to console (for debug)
launch_slim_job <- function(model_params, out_folder, n_reps = NULL, dryrun = FALSE, array = FALSE){
  
  sim_prefix <- paste(model_params, collapse = "_")
  out_folder_base <- paste0(out_folder, sim_prefix, "/")
  
  # check if dirs exist, create them if not
  if(!dir.exists(out_folder_base)){
    dir.create(out_folder_base)
  } 
  
  # all the parameter flags to pass to the slim launcher script
  model_params$out_folder_ba <- out_folder_base
  param_flags <- c("p1si", "p2si", "mu", "rr", "rrc", "rrcd", "chrlen", "samp", "mig", "outfold")
  
  param_string <- data.frame(flag = param_flags, param = unlist(model_params)) %>%
    mutate(full_flag = paste(flag, param, sep ="=")) %>%
    pull(full_flag)
  
  param_string <- paste(param_string, collapse = " ")
  
  if (!array){
    
    if(dryrun){
      
      print(paste0("sbatch scripts/run_slim_sim.sh ", param_string))
      
    } else{
      
      system(paste0("sbatch scripts/run_slim_sim.sh ", param_string))
      
    }
    
  } else if (array){
    
    array_string = paste0("--array=0-", n_reps - 1)
    
    if(dryrun){
      
      print(paste0("sbatch ", array_string, " scripts/run_slim_sim_array.sh ", param_string))
      
    } else{
      
      system(paste0("sbatch ", array_string, " scripts/run_slim_sim_array.sh ", param_string))
      
    }
    
  }
  
  
  
}

# define the range of parameter values
#mig_rate_values <- c(seq(from = 0, to = 0.04, by = 0.01), seq(from = 0.05, to = 0.5, by = 0.05))
#recomb_rate_cha_values <- seq(from = 0.00, to = 1.0, by = 0.1)
#recomb_rate_change_dir <- c(-1, 1)

mig_rate_values <- c(0,0.01,0.05)
recomb_rate_cha_values <- c(0,0.5,1)
recomb_rate_change_dir <- c(1)

# all possible combinations of parameters
# bless u expand.grid
param_combos <- expand.grid(mig_rate_values, recomb_rate_cha_values, recomb_rate_change_dir)

# number of replicates per parameter combination
n_reps <- 1

for (i in 1:nrow(param_combos)){
#for (i in 1:3){  
  #for(j in 1:n_reps){
    
    # FIXED DROSOPHILA POP GEN PARAMETERS
    # empirically estimated population-genetic parameters for d mel
    # used for scaling the simulation parameters (see Messer 2013)
    theta <- 0.007
    mean_recomb <- 2.46 * 10^-8 # 2.46 cM/MB
    ne_estimate <- 1250000 # that's a spicy meatball
    rho_est <- 4*ne_estimate*mean_recomb
    
    # FIXED SIMULATION PARAMETERS
    p1_si <- 1000
    p2_si <- 500
    mutation_rate <- theta/(4*(p1_si)) # SCALED RATE
    recomb_rate <- rho_est/(4*(p1_si))# SCALED RATE
    chr_length <- 100000
    samp_interval <- 1000
    
    # VARIABLE
    mig_rate <- param_combos[i,][,1]
    recomb_rate_cha <- param_combos[i,][,2] 
    recomb_rate_change_dir <- param_combos[i,][,3]
    
    if(recomb_rate_cha == 1 | recomb_rate_cha == 0){
      
      recomb_rate_cha <- sprintf("%1.1f", recomb_rate_cha)
      
    }
    
    model_params <- mget(c("p1_si", "p2_si", "mutation_rate", "recomb_rate", "recomb_rate_cha", 
                           "recomb_rate_change_dir", "chr_length", "samp_interval", "mig_rate"))
    
    out_folder <- "data/slim_output/"
    
    launch_slim_job(model_params, out_folder, n_reps = n_reps, dryrun = FALSE, array = TRUE)
    Sys.sleep(1)
  
}


