library("tidyverse")

sim_params <- list.files("data/slim_output", full.names = TRUE)

sub_sims <- lapply(sim_params, list.files, full.names = TRUE)

vcf_files <- lapply(sub_sims, list.files, full.names = TRUE)

vcf_list <- lapply(unlist(vcf_files), function(x) strsplit(x, split = "/") %>% unlist %>% data.frame %>% t %>% data.frame)

vcf_df <- bind_rows(vcf_list)

names(vcf_df) <- c("data_folder", "slim_folder", "sim_params", "seed", "file")

vcf_df <- vcf_df %>%
  arrange(sim_params, seed, file)

write.table(vcf_df, "simulation_report.txt", row.names = FALSE)