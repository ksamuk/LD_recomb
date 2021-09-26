library("tidyverse")
library("extrafont")
library("plotly")
library("officer")
library("aod")
library("rvg")
source("functions/theme_presentation.R")
loadfonts(device = "win", quiet = TRUE)

# read in windowed pyrho output
rmap_df <- read_delim("data/recomb_maps_1000bp.txt", delim = " ") %>%
  filter(gen != "gen", pop != "pop", window_pos1 != "window_pos1", recomb_est_avg != "recomb_est_avg", rep != "rep", sim != "sim") %>%
  mutate(rep = paste0("r", rep)) %>%
  type_convert()

# reduce to GW averages
rmap_df <- rmap_df %>%
  group_by(sim, rep, ne_control, pop, gen) %>%
  summarise(recomb_est = mean(recomb_est_avg, na.rm = TRUE))

# the model parameters
model_params <- c("p1_si", "p2_prop", "recomb_rate_cha", "recomb_rate_change_dir", "variable_recomb", 
                  "chr_length", "samp_interval", "mig_rate_m12", "mig_rate_m21", "split_generation", 
                  "gene_flow_generation")

rmap_df <- rmap_df %>%
  separate(sim, into = model_params, sep = "_") %>%
  select(-variable_recomb, -split_generation, -chr_length, -samp_interval, -p1_si) 

write_rds(rmap_df, path = "data/pyrho_processed_df.rds")
