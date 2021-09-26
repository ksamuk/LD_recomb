##########################################
# Statistical analyses and plots 
# of pyrho recombation rate estimates
# from SLiM simulations
# 
# Code for FIGURE 4 and associated
# statistical tests
#
# Kieran Samuk
# April 2021 
# ksamuk@gmail.com
##########################################\

##########################################
# Initialize libraries and read in 
# raw pyrho estimates
##########################################

#extrafont::loadfonts(device="win")
library("tidyverse")
library("patchwork")
library("grid")
library("ggthemes")

# the model parameters
model_params <- c("p1_si", "p2_prop", "recomb_rate_cha", "recomb_rate_change_dir", "variable_recomb", 
                  "chr_length", "samp_interval", "mig_rate_m12", "mig_rate_m21", "split_generation", 
                  "gene_flow_generation")

recomb_df <- read.table("data/recomb_maps_1000bp.txt.gz", h = T, stringsAsFactors = TRUE) %>%
  filter(gen != "gen", pop != "pop", window_pos1 != "window_pos1", recomb_est_avg != "recomb_est_avg", rep != "rep", sim != "sim") %>%
  mutate(rep = factor(paste0("r", rep))) %>%
  select(-ne_control) %>%
  type_convert() %>%
  distinct 

#recomb_df <- recomb_df %>%
#  separate(sim, into = model_params, sep = "_") %>%
#  filter(gene_flow_generation %in% c(20000, 21000)) %>%
#  select(-variable_recomb, -split_generation, -chr_length, -samp_interval, -p1_si) 

recomb_sub <- recomb_df %>%
  filter(sim %in% c("1000_0.5_0.0_1_0_100000_1_0.01_0.01_20000_20000", "1000_0.5_0.0_1_0_100000_1_0.01_0.01_20000_21000")) %>%
  filter(window_pos1 <= 95000) %>%
  separate(sim, into = model_params, sep = "_") %>%
  type_convert() %>%
  mutate(gen_since_split = gen - as.numeric(gene_flow_generation)) %>%
  filter(gen_since_split >= 0, gen_since_split <= 30) %>%
  mutate(rep_pop = paste0(rep, "_", pop))

emp_rho <- 2.23e-08


recomb_sub %>%
  filter(pop == 1) %>%
  filter(gen_since_split %in% c(1, 5, 10, 20)) %>%
  ggplot(aes(x = window_pos1, y = recomb_est_avg, group = rep_pop, color = pop))+
  geom_step(alpha = 0.25)+
  geom_hline(yintercept = emp_rho) +
  facet_grid(gene_flow_generation~gen_since_split)+
  theme(legend.position = "none")

recomb_sub %>%
  filter(gen_since_split %in% c(1, 10, 20)) %>%
  ggplot(aes(x = as.factor(window_pos1), y = recomb_est_avg))+
  geom_boxplot()+
  facet_grid(gene_flow_generation~gen_since_split)+
  theme(legend.position = "none")+
  ylim(c(1.5e-08, 3e-08))

recomb_sd <- recomb_df %>% 
  group_by(sim, rep, gen, pop) %>%
  summarize(mean_recomb = mean(recomb_est_avg, na.rm = TRUE), stdev_recomb = sd(recomb_est_avg, na.rm = TRUE))

recomb_sd <- ungroup(recomb_sd) %>%
  separate(sim, into = model_params, sep = "_") 

recomb_sd %>%
  filter(pop == 1) %>%
  filter(gene_flow_generation == 20000 | gene_flow_generation == 21000) %>% 
  mutate(gen_since_split = gen - as.numeric(gene_flow_generation)) %>%
  filter(gen_since_split >= 0, gen_since_split <= 30) %>%
  mutate(gf_gen_label = ifelse(gene_flow_generation == 20000, "Continuous Migration", "Secondary Contact")) %>%
  ggplot(aes(x = mig_rate_m12, y = stdev_recomb, color = gf_gen_label))+
  stat_summary(fun = "mean")+
  facet_grid(gf_gen_label~mig_rate_m12)
  
