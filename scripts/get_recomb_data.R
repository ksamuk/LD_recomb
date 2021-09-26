library("tidyverse")

start_pos <- c(0, seq(from = 5000, to = 1000000, by  = 5000)) + 9000000
end_pos <- start_pos + 4999

pos_dat <- data.frame(start_pos, end_pos)

pos_dat %>%
  mutate(pos_string = paste0("2L:", start_pos, "..", end_pos)) %>%
  select(pos_string) %>%
  write.table(file = "meta/recomb_windows.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)

## obtained recombination data by uploading the text file to:
# https://petrov.stanford.edu/cgi-bin/recombination-rates_updateR5.pl
# Fiston-Lavier AS and Petrov DA. Drosophila melanogaster Recombination Rate Calculator: http://petrov.stanford.edu/cgi-bin/recombination-rates_updateR5.pl

recomb_dat <- read.csv("meta/comeron_recomb.csv", h = T)

# rescale for input into slim
recomb_dat$locus <- (c(0, seq(from = 5000, to = 1000000, by  = 5000))+1)[-201]

# output 
write.table(recomb_dat, file = "meta/recomb_windows_slim_5k.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)


