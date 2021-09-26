# purge old ld files

# the current list of completed simulations
sub_params <- list.files("data/slim_output", full.names = TRUE)
sub_sims <- lapply(sub_params, list.files) %>% unlist

# the list of processed ld files
ld_file_params <- grep(list.files("data/sim_ld", full.names = TRUE), pattern = "summary", inv = T, value = T)
sub_ld_files <- lapply(ld_file_params, list.files, full.names = TRUE) %>% unlist

# 
sub_ld_names <- sub_ld_files %>% gsub(".*/", "", .) %>% gsub("_ld.txt.gz", "", .)