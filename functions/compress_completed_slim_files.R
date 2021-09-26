library("tidyverse")
library("parallel")

sub_params <- list.files("data/slim_output", full.names = TRUE)
sub_sims <- lapply(sub_params, list.files, full.names = TRUE) %>% unlist

x = sub_sims[100]

zip_file <- function(x){

system(paste0("gzip ", x))

}

determine_if_compression_is_needed <- function(x){
  
  # find all the files in the simulation folder
  all_files <- list.files(x, full.names = TRUE, pattern = "\\.vcf$|\\.txt$")
  
  if(length(all_files) > 0){
  
  file_ages <- lapply(all_files, function(x)file.info(x)$ctime) %>% do.call("c", .)
  old_files <- (difftime(Sys.time(), file_ages, units = "hours") %>% as.numeric) > 4
  zip_files <- all_files[old_files]
  
  if(length(zip_files) > 0){
  
      print(paste0("compressing ", length(zip_files), " files in ", x, "..."))
      #print(zip_files)
      mclapply(zip_files, zip_file, mc.cores = 10)
      
  } else{
      print(paste0("nothing to compress in ", x))
  }

  }
}

# check all folders
lapply(sub_sims, determine_if_compression_is_needed)
