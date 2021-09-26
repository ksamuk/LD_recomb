# compute ld from VCFs 

library("vcfR")
library("tidyverse")
library("parallel")
library("SNPRelate")

# function for dividing phased genotypes in two haplotypes
split_column <- function(x, vcf_df_source){
  
  lab <- names(vcf_df_source[,-1])[x]
  x <- vcf_df_source[,-1][,x]
  
  col_names <- c(paste0(lab, "_hap1"), paste0(lab, "_hap2"))
  
  data.frame(x) %>%
    separate(x, into = col_names, sep = "\\|")
  
}

# summarize LD for one window
summarize_ld <- function(x){
  
  # remove naan bread and nas
  x <- x[!is.nan(x)]
  x <- x[!is.na(x)]
  
  #
  mean_ld <- mean(x^2)
  mean_ld_mag <- mean(x)
  sd_ld <- sd(x^2)
  n_snps <- length(x)
  
  data.frame(mean_ld, mean_ld_mag, sd_ld, n_snps)
  
}

# master function for computing LD
compute_ld <- function(vcf_file, sim_param, window_size = 100, downsample = 50){
  
  message(vcf_file)
  message("Processing VCF...")
  # read a vcf
  vcf_df <- read.vcfR(vcf_file)
  
  # get the pos field
  vcf_df <- data.frame(pos = as.numeric(vcf_df@fix[,2]), data.frame(vcf_df@gt)[,-1]) %>%
    arrange(pos)
  
  # if no valid vcf, skip  
  if (length(vcf_df) == 1){
  
	return(data.frame(x = "failed!"))
  }
  
  # downsample vcf (if specified)
  if(is.numeric(downsample)){
    
    name_sample <- sample(names(vcf_df)[-1], size = downsample)
    vcf_df <- vcf_df[,c("pos", name_sample)]
  }
  
  
  # split out haplotypes
  haps <- lapply(1:length(vcf_df[,-1]), split_column, vcf_df_source = vcf_df) %>%
    bind_cols()
  
  haps <- data.frame(pos = vcf_df$pos, haps) %>%
    filter(!duplicated(pos))
  
  message("Computing LD...")  
  
  # create a temporary gds file
  geno_slug <- strsplit(vcf_file, split = "/")[[1]][4:5]
  geno_slug[2] <- geno_slug[2] %>% gsub(".vcf.gz", "", .)
  geno_slug <- paste0(geno_slug[1], "_", geno_slug[2])
  snpgdsCreateGeno(paste0("data/", geno_slug, ".gds"), genmat = data.matrix(haps)[,-1], 
                   sample.id = names(haps)[-1], snp.id = haps$pos, snp.position = haps$pos)
  
  geno_file <- snpgdsOpen(paste0("data/", geno_slug, ".gds"))
  ld_mat <- snpgdsLDMat(geno_file, slide = window_size)
  snpgdsClose(geno_file)
  file.remove(paste0("data/", geno_slug, ".gds"))
  
  ld_df <- data.frame(ld_mat$LD)
  names(ld_df) <- paste0("snp_", ld_mat$snp.id)
  
  ld_list <- lapply(ld_df, summarize_ld)
  ld_df <- data.frame(sim_param, sim_run = geno_slug, snp_id = ld_mat$snp.id, bind_rows(ld_list))
  
  return(ld_df)
  
}

vcf_files <- list.files("slim_sim/bak/recomb_out1602306361870", pattern = ".vcf", full.names = TRUE)
compute_ld(vcf_files[1])

###########################
# DO THE THING
###########################

sim_params <- list.files("data/slim_output", full.names = TRUE)

for (i in sim_params){
  
  message(i)
  sub_sims <- list.files(i, full.names = TRUE)
  
  for (j in sub_sims){
    
    param_slug <- i %>% gsub("data/slim_output/", "", .)
    run_slug <- j %>% gsub("data/slim_output/", "", .) %>% gsub(paste0(param_slug, "/"), "", .)
    
    if(!file.exists(paste0("data/sim_ld/", param_slug, "/", run_slug, "_ld.txt.gz"))){
      
      message(j)
      vcf_files <- list.files(j, full.names = TRUE, pattern = "vcf.gz")
      
      #sim_ld_df <- lapply(vcf_files, compute_ld, sim_param = i, window_size = 1000, downsample = 200)
      sim_ld_df <- mclapply(vcf_files, compute_ld, sim_param = i, window_size = 1000, downsample = 100, mc.cores = 10)
      
      # clean list of failed runs
      names_list <- lapply(sim_ld_df, function(x)length(names(x)))
	  sim_ld_df <- sim_ld_df[!(unlist(names_list) == 1)]
      
      # format to a tidy data frame
      sim_ld_df <- bind_rows(sim_ld_df)
      
      sim_ld_df <- sim_ld_df %>%
        separate(sim_run, sep = "_", into = c("run", "pop","gen"))
      
      # output file
      dir.create("data/sim_ld")
      
      
      dir.create(paste0("data/sim_ld/", param_slug))
      
      write.table(sim_ld_df, gzfile(paste0("data/sim_ld/", param_slug, "/", unique(sim_ld_df$run), "_ld.txt.gz")), row.names = FALSE)
      
    }

  }
  
}


