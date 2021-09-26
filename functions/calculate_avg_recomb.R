library("tidyverse")

recomb_dat <- lapply(list.files("data/comeron_stdpopsim_maps", full.names = TRUE, pattern = "[LR]+.txt"), read.table, h = T) %>%
  bind_rows

x <- recomb_files[1]

process_recomb_file <- function(x){
  
  recomb_dat <-  read.table(x, h = T) 
  start_rows <- recomb_dat[seq(from = 1, to = nrow(recomb_dat), by = 2),]
  end_rows <- recomb_dat[seq(from = 2, to = nrow(recomb_dat), by = 2),]
  
  if (nrow(end_rows) < nrow(start_rows)){
    
    blank_row <- data.frame(Chromosome = NA, Position.bp. = NA, Rate.cM.Mb. = NA)
    end_rows <- bind_rows(end_rows, blank_row)
    
  }
  
  recomb_dat <- bind_cols(start_rows, end_rows)
  names(recomb_dat) <- c("chr", "pos1", "rr1", "chr2", "pos2", "rr2")
  
  recomb_dat
  
}

recomb_files <- list.files("data/comeron_stdpopsim_maps", full.names = TRUE, pattern = "[LR]+.txt")
recomb_dat <- lapply(recomb_files, process_recomb_file) %>%
  bind_rows

recomb_dat %>%
  group_by(chr) %>%
  summarise(avg_recomb = mean(c(rr1, rr2), na.rm = TRUE))

# per chromosome
# chr   avg_recomb
# <fct>      <dbl>
# 1 chr2L       2.40
# 2 chr2R       2.66
# 3 chr3L       1.79
# 4 chr3R       1.97

recomb_dat %>%
  summarise(avg_recomb = mean(c(rr1, rr2), na.rm = TRUE))

# genome wide
#avg_recomb
# 1   2.176645

(2.176645/100)/1000000
# [1] 2.176645e-08

