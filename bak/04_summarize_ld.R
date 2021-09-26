# analyze LD data

library("tidyverse")
library("parallel")

ld_files <- list.files("data/sim_ld", full.names = TRUE, recursive = TRUE)
ld_files <- ld_files[-c(grep("summary", ld_files))]

dir.create("data/sim_ld/summary")

summarize_ld <- function(x, overwrite = FALSE){
  
  # determine the summary file slug
  file_slug <- gsub(".*/", "", x) %>%
    gsub("_.*", "", .) %>%
    paste0("data/sim_ld/summary/",., "_", "ld_summary.txt")
  
  # proceed only if file exists
  if(!file.exists(file_slug) & !overwrite){
    
    message(x)
    
    # read in the ld info
    ld_df <- read.table(x, h = T)
    
    # catch cases where 'pop' wasn't properly specified
    if(is.na(ld_df$gen[1])){
      
      ld_df$gen <- ld_df$pop
      ld_df$pop <- "p1"
      
    }
    
    if(length(ld_df) > 0){
      
      # summarise the ld calculations
      ld_df <- ld_df %>%
        mutate(gen = as.numeric(gen)) %>%
        group_by(sim_param, run, pop, gen) %>%
        filter(n_snps > 10) %>%
        summarise(mean_ld = mean(mean_ld, na.rm = T), 
                  mean_ld_mag = mean(mean_ld_mag, na.rm = T), 
                  sd_ld = mean(sd_ld, na.rm = T))
      
      
      # old method
      #print(ld_df)
      #return(ld_df)
      
      # output summary to file
      write.table(ld_df, file_slug, row.names = FALSE)
      
      
    } else{
      print("FAILED")
      return(data.frame(paste0(x, " failed")))
      
    }
    
  }
}

ld_list <- mclapply(ld_files, summarize_ld, mc.cores = 5)
#ld_list <- lapply(ld_files, summarize_ld)

#if (length(grep("failed", ld_list) >=1)){

#  ld_list <- ld_list[!grepl("failed", ld_list)]

#}

ld_list <- lapply(list.files("data/sim_ld/summary", full.names = T), read.table, h = T)

# simple validation the df list
names_list <- lapply(ld_list, function(x)length(names(x)))
ld_list <- ld_list[unlist(names_list) == 7]

ld_df <- bind_rows(ld_list)

write.table(ld_df, "data/sim_ld/ld_summary_out.txt", row.names = FALSE)
