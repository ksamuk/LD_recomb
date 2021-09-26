# convert a slim output file to a fasta file
# applies a mutation transition matrix (impossible with the VCF approach)

library("tidyverse")
library("parallel")
library("inline")
options(scipen = 999)

########################################
# SIMULATION WRAPPER(S)
########################################

# converts a list of parameters to "-d [param] = [value]" (for slim command line call)
convert_params_to_call <- function(x, model_par){
  
  # slim needs characters formatted in special quote configuration 
  if(class(unlist(model_par[x])) == "character"){
    paste0(" -d ", "\"",x, "=", "\'", unlist(model_par[x]),"\'", "\"")
  }else{
    paste0(" -d ", x, "=", unlist(model_par[x]))
  }
  
}

# generate a random string (names of simulations)
run_slim_simulation <- function(sim_num, model_params, simulation_script, out_folder){
  
  sim_prefix <- paste(model_params, collapse = "_")
  out_folder_base <- paste0(out_folder, sim_prefix, "_no_split")
  
  # check if dirs exist, create them if not
  if(!dir.exists(out_folder_base)){
    dir.create(out_folder_base)
  } 
  
  model_params$out_folder_ba <- paste0(out_folder_base, "/")
  
  # parse parameters to a system call
  param_call <- sapply(names(model_params), convert_params_to_call, model_par = model_params) 
  slim_call <- paste0("slim", paste0(param_call, collapse = ""), " ", simulation_script)
  system(slim_call)
  
  
}

# works
#run_slim_simulation(1, model_params, "slim_sim/recomb_rate_evolve.txt", out_folder_base)


# wrapper for mclapply

# wait function from 'inline'
# prevents zombie processs overflow
includes <- '#include <sys/wait.h>'
code <- 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
wait <- cfunction(body=code, includes=includes, convention='.C')

mclapply_wrapper <- function(n_reps, model_para, sim_script, output_folder, n_cores){
  
  
  for (i in 1:(n_reps / n_cores)){
  
  mclapply(1:n_cores, run_slim_simulation, 
           model_params = model_para, simulation_script = sim_script, out_folder = output_folder,
           mc.cores = n_cores, mc.preschedule = FALSE)
           
  sim_prefix <- paste(model_params, collapse = "_")
  out_folder_base <- paste0(output_folder, sim_prefix)
  system(paste0("bash functions/compress_slim_output.sh ", out_folder_base));
  
  }
  
  wait()
  
}

########################################
# RUN THE ACTUAL SIMULATIONS
########################################

# simulation parameters 

replicates <- 50
cores <- 20

theta <- 0.007
p1_si <- 1000
p2_si <- 1000
mutation_rate <- theta/(4*(p1_si))
recomb_rate <- 3 * 10^-6
recomb_rate_cha <- 0.5
recomb_rate_change_dir <- 1
chr_length <- 1000000
samp_interval <- 1000

model_params <- mget(c("p1_si", "p2_si", "mutation_rate", "recomb_rate", "recomb_rate_cha", 
                       "recomb_rate_change_dir", "chr_length", "samp_interval"))

simulation_script <- "slim_sim/recomb_rate_evolve_no_split.txt"
out_folder_base <- "data/slim_output/"

# run the actual simulation
mclapply_wrapper(replicates, model_params, "slim_sim/recomb_rate_evolve_no_split.txt", out_folder_base, cores)
