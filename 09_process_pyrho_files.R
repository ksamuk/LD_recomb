library("tidyverse")
library("parallel")
#library("extrafont")
#library("plotly")
#library("officer")
#library("aod")
#library("rvg")
#source("functions/theme_presentation.R")
#loadfonts(device = "win", quiet = TRUE)

rmap_files <- list.files("optimize", full.names = TRUE, recursive = TRUE, 
                         pattern = "*.rmap")

#rmap_files <- c(rmap_files[grepl("\\/p[12]+_19998", rmap_files)],
#                rmap_files[grepl("\\/p[12]+_20000", rmap_files)],
#                rmap_files[grepl("\\/p[12]+_20100", rmap_files)],
#                rmap_files[grepl("\\/p[12]+_21000", rmap_files)],
#                rmap_files[grepl("\\/p[12]+_25000", rmap_files)],
#                rmap_files[grepl("\\/p[12]+_30000", rmap_files)],
#                rmap_files[grepl("\\/p[12]+_31000", rmap_files)],
#                rmap_files[grepl("\\/p[12]+_50000", rmap_files)])

#rmap_file <- "data/optimize/1000_0.1_0.0_-1_0_100000_1000_0_0_20000_20000/1615085801584/p1_19998_no_split_recomb.rmap"

cat("files read, processing...")

tidy_rmap <- function(rmap_file){
  
  cat(rmap_file)
  
  rep <- rmap_file %>% gsub("/p.*", "", .) %>% gsub(".*/", "", .)
  sim <- rmap_file %>% gsub("/[0-9]*/p.*", "", .) %>% gsub(".*/", "", .)
  ne_control <- ifelse(grepl("no_split", rmap_file), "no_split", "split_aware")
  gen <- rmap_file %>% gsub(".*/p[1-2]*_", "", .) %>% gsub("_.*", "", .)
  pop <- rmap_file %>% gsub(".*/p", "", .) %>% gsub("_.*", "", .)
  
  outfile <- paste0("data/recomb_maps/", sim, "_", rep, "_", gen, "_", pop, "_", ne_control, "_recomb_map_1000bp.txt")
  
  if(!file.exists(outfile)){
    
    rmap_df <- read.table(rmap_file)
    names(rmap_df) <- c("pos1", "pos2", "recomb_est")
    rmap_df$rep <- rep
    rmap_df$sim <- sim
    rmap_df$pop <- pop
    rmap_df$gen <- gen
    rmap_df$ne_control <- ne_control
    
    window_size <- 1000
    
    rmap_df <- rmap_df %>% 
      mutate(window_pos1 = floor(pos1 / window_size) * window_size) %>%
      mutate(window_pos2 = floor(pos1 / window_size) * 
               window_size + window_size) %>%
      group_by(sim, rep, ne_control, gen, pop, window_pos1) %>%
      summarise(recomb_est_avg = mean(recomb_est, na.rm = TRUE)) %>%
      ungroup
    
    write.table(rmap_df, outfile, quote = FALSE, row.names = FALSE)

    
  }
  
}

mclapply(rmap_files, tidy_rmap, mc.cores = 16)

rmap_processed_files <- list.files("data/recomb_maps", full.names = TRUE)
rmap_out <- lapply(rmap_processed_files, read.table, header = TRUE)
rmap_out <- bind_rows(rmap_out)

write.table(rmap_out, "data/recomb_maps_1000bp.txt", quote = FALSE, row.names = FALSE)

