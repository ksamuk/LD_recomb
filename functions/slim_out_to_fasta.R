# convert a slim output file to a fasta file
# applies a mutation transition matrix (impossible with the VCF approach)

library("tidyverse")
options(scipen = 999)

slim_slugs <- list.files("data/slim_output")
slug <- slim_slugs[[1]]

########################################
# SIMULATION WRAPPER(S)
########################################

# converts a list of parameters to "-d [param] = [value]" (for slim command line call)
convert_params_to_call <- function(x, model_par){
  
  if(class(unlist(model_par[x])) == "character"){
    paste0(" -d ", "\"",x, "=", "\'", unlist(model_par[x]),"\'", "\"")
  }else{
    paste0(" -d ", x, "=", unlist(model_par[x]))
  }
  
}

# generate a random string (names of simulations)

# simulation parameters 
theta <- 0.007
p1_size <- 100
p2_size <- 100
mutation_rate <- theta/(4*(p1_size))
recomb_rate <- 3 * 10^-6
recomb_rate_change <- 0.5
recomb_rate_change_direction <- 1
chr_length <- 1000000
sample_interval <- 1000

model_params <- mget(c("p1_size", "p2_size", "mutation_rate", "recomb_rate", "recomb_rate_change", 
                       "recomb_rate_change_direction", "chr_length", "sample_interval"))

simulation_script <- "slim_sim/recomb_rate_evolve.txt"
out_folder_base <- "data/slim_output/"


run_slim_simulation <- function(model_params, simulation_script, max_sims = 10, out_folder, sim_prefix){
  
  sim_prefix <- paste(model_params, collapse = "_")
  out_folder_base <- paste0(out_folder, sim_prefix)
  
  # check if dirs exist, create them if not
  if(!dir.exists(out_folder_base)){
    dir.create(out_folder_base)
  } 
  
  # revert back to base output (slim won't accept full prefix path for some reason)
  model_params$out_folder_base <- paste0(out_folder_base, "/")
  
  # parse parameters to a system call
  param_call <- sapply(names(model_params), convert_params_to_call, model_par = model_params) 
  slim_call <- paste0("slim", paste0(param_call, collapse = ""), " ", simulation_script)
  system(slim_call)
  
}


# wrapper for mclapply

# wait function from 'inline'
# prevents zombie processs overflow
includes <- '#include <sys/wait.h>'
code <- 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
wait <- cfunction(body=code, includes=includes, convention='.C')

mclapply_wrapper <- function(n_reps, model_para, sim_script, n_cores){
  
  #chunks <- 1:n_reps
  #chunks <- split(chunks, ceiling(seq_along(chunks)/n_cores))
  
  #for (i in 1:length(chunks)){
  
  mclapply(1:n_reps, run_simulation), 
           model_params = model_para, max_sims = n_reps,
           simulation_script = sim_script, mc.cores = n_cores,
           mc.preschedule = TRUE)
  
  wait()
  #}
  
}




# initialize a blank (all A) sequnece of appropriate length
seq_len <- args[[1]]
blank_seq <- rep("A", seq_len)

###########################################################
# LIBS + INITIALIZE RAW DATA
###########################################################

library("tidyverse")
library("parallel")

# the slim output files
full_slim <- list.files("data/slim_output", pattern = "full", full.names = TRUE)
fixed_slim <- list.files("data/slim_output", pattern = "fixed", full.names = TRUE)



# - file is actually several stacked dataframe-like objects 
# - #yudothis

# read the full output lines and split based on the section headers
# we actually only care about mutations, individuals (I guess) and genomes
lines <- readLines(full_slim)

mut_line <- grep("Mutations", lines)
ind_line <- grep("Individuals", lines)
genome_line <- grep("Genomes", lines)

# for mut_df and ind_df, a function to map lines to data frames
line_to_df <- function(x){
  
  # what is the length of 1 row of x?
  x_len <- x[1] %>% strsplit(split = " ") %>% unlist %>% length
  
  # convert to df
  data.frame(col1 = x) %>% 
    separate(col1, into = paste0("col", 1:x_len), sep = " ")
  
}

# slice out the sections
mut_df <- lines[(mut_line + 1):(ind_line - 1)]
ind_df <- lines[(ind_line + 1):(genome_line - 1)]
genome_df <- lines[(genome_line + 1):(length(lines))]

##### the mutations section

# The first field is a within-file numeric identifier for the mutation, beginning at 0 and counting up (although mutations are not listed in sorted order according to this value); see below for a note on why this field exists. Second is the mutation’s id property (see section 18.8.1), a within-run unique identifier for mutations that does not change over time, and can thus be used to match up information on the same mutation within multiple output dumps made at different times. Third is the identifier of the mutation’s mutation type, such as m1. Fourth is the position of the mutation on the chromosome, as a zero-based base position. Fifth is the selection coefficient of the mutation, and sixth is its dominance coefficient (the latter being a property of the mutation type, in fact). Seventh is the identifier for the subpopulation in which the mutation originated, and eighth is the generation in which it originated. Finally, the ninth field gives the mutation’s prevalence: an integer count of the number of times that it occurs in any genome in the population.

# WARNING: careful about the difference between mut_id and mut_id_unique

mut_df <- line_to_df(mut_df)
names(mut_df) <- c("mut_id", "mut_id_unique", "type", "pos", "s", "h", "pop_orig", "gen_orig", "count_in_pop")

##### the individuals section

# This describes each individual in each subpopulation and specifies which genomes belong to it. The first field is an identifier for the individual, such as p1:i0 (which indicates the 0th individual in p1). Next is the sex of the 327 individual: H for a hermaphrodite, or if sex has been enabled, M for a male or F for a female. Following that come two genome specifiers, of the form p1:0 (indicating the 0th genome in p1).

ind_df <- line_to_df(ind_df)
names(ind_df) <- c("ind_id", "sex", "genome1", "genome2")


# select 20 (random) individuals to output to fasta (file contains ALL individuals -- only want a subset)

ind_sample <- ind_df %>%
  sample_n(10)


#### the genomes section
# ...oh boy

# filter out the individuals we want to keep

parse_genome <- function(x){
  
  x_split <- strsplit(x, split = " ") %>% unlist
  
  genome_id <- x_split[1]
  
  if(genome_id %in% ind_sample$genome1 | genome_id %in% ind_sample$genome2){
    
    # the list of mutations
    mut_id <- x_split[3:length(x_split)]
    
    # match back to an individual
    ind_id <- ind_sample$ind_id[c(which(ind_sample$genome1 == genome_id), which(ind_sample$genome2 == genome_id))]
    genome_out <- data.frame(ind_id, genome_id, mut_id)
    
    ind_id <- ind_sample$ind_id[c(which(ind_sample$genome1 == genome_id), which(ind_sample$genome2 == genome_id))]
    
    mut_df_join <- mut_df %>% 
      select(mut_id, mut_id_unique, pos, gen_orig)
    
    left_join(genome_out, mut_df_join, by = "mut_id")
    
  }
  
  
}

genome_list <- lapply(genome_df, parse_genome)
keepers <- (lapply(genome_list, length) %>% unlist) != 0
genome_df <- genome_list[keepers] %>% bind_rows
genome_df <- genome_df %>%
  arrange(ind_id, genome_id, mut_id_unique)


###########################################################
# ALL SUBSTITUTIONS (things that fixed in past gens)
###########################################################

fixed_df <- read.table(fixed_slim, skip = 2)
names(fixed_df) <- c("mut_id", "mut_id_unique", "type", "pos", "s", "h", "pop_orig", "gen_orig", "gen_fixed")

# ENGAGE THE MUTATION MATRIX!
mut_mat <- read.table("data/mut_mat.txt")
names(mut_mat) <- c("A", "C", "G", "T")

# coerce to long format
mut_mat <- data.frame(orig_nuc = c("A", "C", "G", "T"), mut_mat) %>%
  gather(key = mut_nuc, value = trans_prob, -orig_nuc) %>%
  arrange(orig_nuc)


# build mutation histories for each locus
# split, apply, combine
fixed_list <- split(fixed_df, fixed_df$pos)

# function for mutating loci based on mutation matrix
# draws from a multinomial based on the matrix
mutate_loci <- function(x){
  
  x$mut_nuc <- NA
  
  for (i in 1:nrow(x)){
    
    if(i == 1){
      original_nuc <- "A"
    } else{
      original_nuc <- x$mut_nuc[i-1]
    }
    nuc_probs <- mut_mat %>%
      filter(orig_nuc == original_nuc)
    
    draw <- rmultinom(1, size = 1, prob = nuc_probs$trans_prob) 
    x$mut_nuc[i] <- nuc_probs$mut_nuc[which(draw == 1)]
    
  }
  
  x
  
}

# apply mutations and bind to data frame
fixed_list <- lapply(fixed_list, mutate_loci)
fixed_df <- bind_rows(fixed_list)

# we actually only care about the most recent mutational state, so only retain the most recent mut
# (frequency of multi-mutants will be higher for longer sims)
fixed_df <- fixed_df %>%
  group_by(pos) %>%
  filter(gen_fixed == max(gen_fixed)) %>%
  ungroup

###########################################################
# ADD IN ALT ALLELES
###########################################################

# add in the mutation identities
genome_df$mut_id_unique <- as.numeric(genome_df$mut_id_unique)

# add in the REFERENCE alleles (based on past substitutions)
# if there was no sub, ref allele is "A" (arbitrarily)

genome_df <- left_join(genome_df, fixed_df %>% select(mut_id_unique, mut_nuc), by = "mut_id_unique")
names(genome_df)[7] <- "ref_nuc"
genome_df$ref_nuc[is.na(genome_df$ref_nuc)] <- "A"

# mutate in the alt allele for each site
# same prinicple as above
# forces all alt_alleles with same pos to be the same (as it should be)

mutate_alt_alleles <- function(alt_pos){
  
  # determine the reference allele
  ref_nuc <- genome_df %>%
    filter(pos == alt_pos) %>%
    select(ref_nuc) %>%
    .[1,]
  
  # subset the mutation matrix for the ref nuc
  nuc_probs <- mut_mat %>%
    filter(orig_nuc == ref_nuc)
  
  draw <- rmultinom(1, size = 1, prob = nuc_probs$trans_prob) 
  alt_nuc <- nuc_probs$mut_nuc[which(draw == 1)]
  
  data.frame(pos = alt_pos, ref_nuc, alt_nuc = rep(alt_nuc, nrow(genome_df %>% filter(pos == alt_pos))))
  
}
alt_nucs <- mclapply(unique(genome_df$pos), mutate_alt_alleles, mc.cores = 2)
alt_nucs <- bind_rows(alt_nucs)

alt_nucs <- alt_nucs %>%
  distinct %>%
  select(pos, alt_nuc) %>%
  arrange(pos)

genome_df <- left_join(genome_df, alt_nucs, by = c("pos")) 

# adjust positioning

genome_df$pos <- as.numeric(genome_df$pos)

# cleanup
rm(fixed_df)
rm(fixed_list)
rm(genome_list)
rm(mut_df)
rm(alt_nucs)

###########################################################
# CREATE FASTA SEQUENCES
###########################################################

# recycled from VCF script

# note that sampling has already occured, so we just output all the haplotypes in genome_df

write_fasta_out <- function(genome_id_sub, genome_df_in, file_name){
  
  hap_df_sub <- genome_df_in %>%
    filter(genome_id == genome_id_sub)
  
  hap_df_sub$pos <- hap_df_sub$pos[(hap_df_sub$pos >= 0) & (hap_df_sub$pos <= seq_len)]
  
  #apply the haplotype mask
  out_seq <- blank_seq
  out_seq[as.numeric(hap_df_sub$pos)] <- hap_df_sub$alt_nuc
  out_seq <- paste0(out_seq, collapse = "")
  
  out_header <- paste0("> ", genome_id_sub)
  
  # write to file
  cat(out_header, out_seq, file = file_name, sep = "\n",append = file.exists(file_name))
  
  
}

# total number of haplotypes
out_file_name <- paste0("data/fasta_output/", out_file_slug, ".fasta")

# process haplotypes to a fasta file
file.remove(out_file_name)
invisible(lapply(unique(genome_df$genome_id), write_fasta_out, genome_df_in = genome_df, file_name = out_file_name))

# output the ancestral allele priors (all anc alleles are A's)
# remove redundant snps (if downsampling)

genome_df <- genome_df %>%
  mutate(pos = as.numeric(pos) - 1) %>%
  arrange(pos)

anc_prior <- genome_df %>%
  select(pos) %>%
  distinct %>%
  mutate(A = 1.0, C = 0.0, G = 0.0, T = 0.0)

any(!(anc_prior$pos %in% genome_df$pos))
any(!(anc_prior$pos %in% genome_df$pos))

anc_out_name <- gsub("\\.fasta", "_anc_prior.txt", out_file_name)

write.table(anc_prior, file  = anc_out_name, col.names = FALSE, row.names = FALSE, quote = FALSE)
