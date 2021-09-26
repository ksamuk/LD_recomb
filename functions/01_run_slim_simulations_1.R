# run slim simulations

library("tidyverse")
library("parallel")
library("inline")


########################################
# FUNCTIONS 
########################################

# (partial) function for coverting model parameters to a slim call
convert_params_to_call <- function(x, model_par){
  paste0(" -d ", x, "=", unlist(model_par[x]))
}

# main simulation control function

run_slim_simulation <- function(simulation_number, simulation_name, model_params_internal, simulation_script, max_sims = 100){
  
  
  # create directory for output (if necessary)
  simulation_root <- paste0("output/", "s_", model_params_internal$adaptive_allele_selection_coeff)
  
  if(!dir.exists(simulation_root)){
    dir.create(simulation_root)
  } 

  simulation_prefix <- paste0(simulation_root, "/", simulation_name)
  
  # check if dirs exist, create them if not
  if(!dir.exists(simulation_prefix)){
    dir.create(simulation_prefix)
    dir.create(paste0(simulation_prefix, "/slim_output"))
  } 
  
  #file name for SLiM output
  slim_file_name <- paste0(simulation_prefix, "/slim_output/", simulation_number, ".txt")

  # parse parameters to a system call
  param_call <- sapply(names(model_params_internal), convert_params_to_call, model_par = model_params_internal) 
  
  # loop until simulations succeed
  # work around for unknown seg fault :|
  inv_frq <- numeric(0)
  
  while(length(inv_frq) == 0){
    
    # build the slim call and execute it via the command line
    slim_call <- paste0("slim", paste0(param_call, collapse = ""), " ", simulation_script, " > ", slim_file_name)
    paste0("slim", paste0(param_call, collapse = ""), " ", simulation_script)
    
    # parse slim output file into a tidy data frame
    if(file.exists(slim_file_name)){
      
      system(paste0("grep -A 100 START ", slim_file_name, " > ", slim_file_name, "_tmp.txt"))
      inv_frq <- scan(paste0(slim_file_name, "_tmp.txt"), skip = 1, quiet = TRUE, what = character())
      temp <- file.remove(paste0(slim_file_name, "_tmp.txt"))
      temp <- file.remove(paste0(slim_file_name))
      
    }

    # if the simulation fails, write the call to a log file
    if(length(inv_frq) == 0){
      write.table(data.frame(call = slim_call), file = "slim_log.txt", append = TRUE)
    }
    
  }

  # write simulation output to file 
  table_out <- data.frame(simulation_name, simulation_number, model_params_internal, gen = 1:length(inv_frq), inv_frq )
  out_file_name <- paste0(simulation_prefix,"/", simulation_name, "_n=", max_sims,".csv")
  write_csv(table_out, out_file_name, append = TRUE)
  
}

# wrapper for simulation function
# allows looping when model parms are > 1 in length
# only implemented for one param (adaptive_allele_selection_coeff)

simulation_wrapper <- function(simulation_number, simulation_name, model_params, simulation_script, max_sims){
  
  if(length(model_params$adaptive_allele_selection_coeff) > 1){
    
    for (i in unlist(model_params$adaptive_allele_selection_coeff)){
      model_params_temp <- model_params
      model_params_temp$adaptive_allele_selection_coeff <- i
      run_slim_simulation(simulation_number, simulation_name, model_params_temp, simulation_script, max_sims)
    } 
    
   } else{
     
    run_slim_simulation(simulation_number, simulation_name, model_params, simulation_script, max_sims)
     
    }
}

# wrapper for mclapply

# wait function from 'inline'
# prevents zombie processs overflow
includes <- '#include <sys/wait.h>'
code <- 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
wait <- cfunction(body=code, includes=includes, convention='.C')

mclapply_wrapper <- function(n_reps, model_para, sim_name, n_cores){

  #chunks <- 1:n_reps
  #chunks <- split(chunks, ceiling(seq_along(chunks)/n_cores))
  
  #for (i in 1:length(chunks)){
    
    mclapply(1:n_reps, simulation_wrapper, simulation_name = sim_name, 
             model_params = model_para, max_sims = n_reps,
             simulation_script = sim_script, mc.cores = n_cores,
             mc.preschedule = TRUE)
    
   wait()
  #}

}

########################################
# PERFORM SIMULATIONS
########################################

# GENERAL SIMULATION PARAMETERS

sim_script <- "kirk_bart_simulation_experiment_deploy.txt"
num_replicates <- 5004
num_cores <- 20

########################################
# BASE PARAMETERS
########################################

# is the "gene flow" base model

base_parameters <- list()
base_parameters$init_inversion_freq <- 0.5
base_parameters$inversion_selection_coeff <- 0
base_parameters$inversion_underdom_penalty <- 0
base_parameters$prop_fitness_chr3 <- 0.1
base_parameters$init_freq_adpative_alleles <- 0.8
base_parameters$adaptive_allele_selection_coeff <- c(0)
base_parameters$gene_flow_exists <- 1
base_parameters$initial_proportion_migrants <- 0.5
base_parameters$migration_rate_constant <- 0.005

########################################
# CLEAR OUT OLD DATA
########################################

#system("rm -r output")
#system("rm -r slim_log.txt")
#dir.create("output")

########################################
# BASE MODELS
########################################

# GENE FLOW MODEL

model_parameters <- base_parameters
mclapply_wrapper(num_replicates, model_parameters, "gene_flow", num_cores)

# ALLOPATRY MODEL

model_parameters <- base_parameters
model_parameters$gene_flow_exists <- 0

mclapply_wrapper(num_replicates, model_parameters, "allopatry", num_cores)

########################################
# UNDERDOMINANCE
########################################

# GENE FLOW + UNDERDOMINANCE

model_parameters <- base_parameters
model_parameters$inversion_underdom_penalty <- -0.2

mclapply_wrapper(num_replicates, model_parameters, "gene_flow_underdom", num_cores)

# ALLOPATRY + UNDERDOMINANCE

model_parameters <- base_parameters
model_parameters$inversion_underdom_penalty <- -0.2
model_parameters$gene_flow_exists <- 0

mclapply_wrapper(num_replicates, model_parameters, "allopatry_underdom", num_cores)


