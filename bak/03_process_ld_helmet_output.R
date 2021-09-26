# process ld helmet output
# runlike Rscript 03_process_ld_helmet_output.R simname

# read in files

library("tidyverse")
library("cowplot")

#args <- commandArgs(trailingOnly = TRUE)
args <- list("test")

# the rho scaling factor (rNe4)
rho_simulated <- (3e-8)*500*4

options(scipen = 999)

# function for processing ld_helmet files to plots

process_sub_sim <- function(sub_sim, sim_name, real_rho = NA, rates = NA, ends = NA){
  
  # convert ld helmet data to df
  ldh_file_name <- paste0("data/ldhelmet_output/", sim_name, "/", sub_sim, "/", sub_sim, ".txt")
  ld_df <- read.table(ldh_file_name, skip = 3)
  names(ld_df) <- c("left_snp", "right_snp", "mean_rho", "p0.025", "p0.500", "p0.975")
  
  #### VARIABLE RATES
  if (length(rates) > 1){
    
    # build a df of the actualy (simulated rates)
    recomb_sim_df <- data.frame(left_snp = ends[1:11], right_snp = ends[seq(from = 2, to = 12, by = 1)], mean_rho = rates) %>%
      rowwise %>%
      mutate(mid_snp = mean(c(left_snp, right_snp))) %>%
      ungroup
    
    # interpolate + assess fit
    
    # assign the windows
    ld_df$left_wind <- NA
    ld_df$right_wind <- NA
    ld_df$real_rho <- NA
    
    for (i in 1:length(ends)){
      
      ld_df$right_wind[is.na(ld_df$right_wind) & ld_df$left_snp <= ends[i]] <- ends[i]
      ld_df$real_rho[is.na(ld_df$left_wind) & ld_df$right_snp <= ends[i]] <- rates[i-1]
      ld_df$left_wind[is.na(ld_df$left_wind) & ld_df$right_snp <= ends[i]] <- ends[i-1]

    }
    
    p <- ld_df %>%
      rowwise %>%
      mutate(mid_snp = mean(c(left_snp, right_snp))) %>%
      ungroup %>%
      ggplot(aes(x = mid_snp, y = mean_rho))+
      geom_line()+
      geom_smooth(se = FALSE)+
      geom_line(aes(x = mid_snp, y= real_rho), color = "red", size = 1)
    
    p2 <- ld_df %>%
      rowwise %>%
      mutate(mid_snp = mean(c(left_snp, right_snp))) %>%
      ungroup %>%
      filter(mean_rho < max(real_rho)) %>%
      ggplot(aes(x = mid_snp, y = mean_rho))+
      geom_line()+
      geom_smooth(se = FALSE)+
      geom_line(aes(x = mid_snp, y= real_rho), color = "red", size = 1)
    
   
    resid <- data.frame(sub_sim = sub_sim, resid = sum(abs(ld_df$mean_rho - real_rho) / length(ld_df$mean_rho)))
    
  } else{
    
    #### VARIABLE RATES
    #assess fit of each window
    
    
    p <- ld_df %>%
      rowwise %>%
      mutate(mid_snp = mean(c(left_snp, right_snp))) %>%
      ungroup %>%
      ggplot(aes(x = mid_snp, y = mean_rho))+
      geom_line()+
      geom_smooth(se = FALSE)+
      geom_hline(yintercept = real_rho, color = "red", linetype = 2)
    
    p2 <- ld_df %>%
      rowwise %>%
      mutate(mid_snp = mean(c(left_snp, right_snp))) %>%
      ungroup %>%
      filter(mean_rho < max(real_rho)) %>%
      ggplot(aes(x = mid_snp, y = mean_rho))+
      geom_line()+
      geom_smooth(se = FALSE)+
      geom_line(aes(x = mid_snp, y= real_rho), color = "red", size = 1)
    
    resid <- data.frame(sub_sim = sub_sim, resid = sum(abs(ld_df$mean_rho - real_rho) / length(ld_df$mean_rho)))
    
    
  }

  ggsave(p, filename = paste0("data/block_optim_output/",sim_name,"/plots/",sub_sim, ".pdf"), height = 5, width = 5, device = "pdf")
  return(list(ld_df, resid))
}

rho_simulated <- (3e-7)*500*4

# rate variation
rates <- c(10, 8, 6, 4, 2, 1, 2, 4, 6, 8, 10) * rho_simulated
ends <-  c(0, 909090 ,1818180, 2727270, 3636360, 4545450, 5454540, 6363630, 7272720, 8181810, 9090900, 10000000)
# get list of ldhelmet output files
sim_name <- list.files("data/ldhelmet_output")[3] 
sub_sims <- list.files(paste0("data/ldhelmet_output/", sim_name))
dir.create(paste0("data/block_optim_output/",sim_name))
dir.create(paste0("data/block_optim_output/",sim_name,"/plots/"))
sub_sim_resid <- lapply(sub_sims, process_sub_sim, sim_name = sim_name, real_rho = rho_simulated, rates = rates, ends = ends)


# one rate
# get list of ldhelmet output files
sim_name <- list.files("data/ldhelmet_output")[2] 
sub_sims <- list.files(paste0("data/ldhelmet_output/", sim_name))
dir.create(paste0("data/block_optim_output/",sim_name))
dir.create(paste0("data/block_optim_output/",sim_name,"/plots/"))
sub_sim_resid <- lapply(sub_sims, process_sub_sim, sim_name = sim_name, real_rho = rho_simulated)

