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
  param_flags <- c("p1si", "p2prop", "rrc", "rrcd", "varrr", "chrlen", 
                   "samp", "mig12", "mig21", "split", "gfgen", "outfold")
  
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

####################################################
# creat parameter combos and launch jobs
####################################################

# define the range of parameter values
#mig_rate_values <- c(seq(from = 0, to = 0.04, by = 0.01), seq(from = 0.05, to = 0.5, by = 0.05))
#recomb_rate_cha_values <- seq(from = 0.00, to = 1.0, by = 0.1)
#recomb_rate_change_dir <- c(-1, 1)

# migration rate in terms of Nem
mig_rate_12_values <- c(0, 0.01, 0.1, 1, 10, 100)

mig_rate_21_values <- c(0, 0.01, 0.1, 1, 10, 100)

# the proportion of individuals that split into p2 
p2_prop_values <- c(0.5)

# the proportion of individuals that split into p2 
chr_length_values <- c(100000)

# the fold change in recomb rate (1 = +100% (1-fold) change)
recomb_rate_cha_values <- c(0, 1)

# the direction of the recombination rate change (1 = increase, -1 = decrease)
recomb_rate_change_dir <- c(1)

# variable recombination rate
variable_recomb_values <- 0

# the generation to resume gene flow after the split
# 20000 = the split generation
gene_flow_generation <- c(20000)

# all possible combinations of parameters
# bless u expand.grid
param_combos <- expand.grid(mig_rate_12_values, mig_rate_21_values, p2_prop_values, chr_length_values,
                            recomb_rate_cha_values, recomb_rate_change_dir, variable_recomb_values, gene_flow_generation)

# constrain the param combos for symmetry 
param_combos <- param_combos %>%
  filter(Var3 == 0.5) %>%
  filter(Var1 == Var2) 

# number of replicates per parameter combination
n_reps <- 100

for (i in 1:nrow(param_combos)){
#for (i in 1:3){  
  #for(j in 1:n_reps){
  
    
    # FIXED SIMULATION PARAMETERS
    p1_si <- 1000
    samp_interval <- 1
    split_generation <- 20000
    
    # VARIABLE PARAMETERS
    mig_rate_m12 <- param_combos[i,][,1]
    mig_rate_m12 <- mig_rate_m12/p1_si # convert to m from Nem
    
    mig_rate_m21 <- param_combos[i,][,2]
    mig_rate_m21 <- mig_rate_m21/p1_si 
    
    p2_prop <- param_combos[i,][,3]
    
    chr_length <- param_combos[i,][,4]
    
    recomb_rate_cha <- param_combos[i,][,5] 
    recomb_rate_change_dir <- param_combos[i,][,6]
    variable_recomb <- param_combos[i,][,7]
    
    gene_flow_generation <- param_combos[i,][,8]
    
    if(recomb_rate_cha == 1 | recomb_rate_cha == 0){
      
      recomb_rate_cha <- sprintf("%1.1f", recomb_rate_cha)
      
    }
    
    model_params <- mget(c("p1_si", "p2_prop", "recomb_rate_cha", "recomb_rate_change_dir", "variable_recomb", 
                           "chr_length", "samp_interval", "mig_rate_m12", "mig_rate_m21", "split_generation", 
                           "gene_flow_generation"))
    
    out_folder <- "data/slim_output/"
    
    
    launch_slim_job(model_params, out_folder, n_reps = n_reps, dryrun = FALSE, array = TRUE)
    Sys.sleep(1)
}
