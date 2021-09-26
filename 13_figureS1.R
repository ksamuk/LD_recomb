##########################################
# Statistical analyses and plots 
# of pyrho recombation rate estimates
# from SLiM simulations
# 
# Code for FIGURE S1 
#
# Kieran Samuk
# April 2021 
# ksamuk@gmail.com
##########################################\

extrafont::loadfonts(device="win")
library("tidyverse")
library("ggrastr")
library("grid")
library("ggthemes")
#remotes::install_git("https://git.rud.is/hrbrmstr/hrbrthemes.git")
#remotes::install_git("https://github.com/brooke-watson/bplots.git")
library("hrbrthemes")
library("bplots")

lapply(list.files("functions", pattern = "theme", full.names = TRUE), 
       source, echo = FALSE)

fst_df <- read_delim("data/fst_comparisons.txt.gz", col_names = TRUE, delim = " ")
names(fst_df) <- c("sim", "rep", "file1", "file2", "fst")

model_params <- c("p1_si", "p2_prop", "recomb_rate_cha", "recomb_rate_change_dir", "variable_recomb", 
                  "chr_length", "samp_interval", "mig_rate_m12", "mig_rate_m21", "split_generation", 
                  "gene_flow_generation")

fst_df <- fst_df %>%
  mutate(gen = gsub(".*/p[12]_", "", file1) %>% gsub("_.*", "", .)) %>%
  separate(sim, into = model_params, sep = "_")

figureS1 <- fst_df %>%
  #filter(gen < 20040) %>%
  filter(gene_flow_generation == 20000 | gene_flow_generation == 21000) %>%
  mutate(gen_since_split = as.numeric(gen) - as.numeric(gene_flow_generation)) %>%
  filter(gen_since_split > 0, gen_since_split <=30) %>%
  mutate(gf_gen_label = ifelse(gene_flow_generation == 20000, "Continuous Migration", "Secondary Contact")) %>%
  sample_frac(1) %>%
  ggplot(aes(x = gen_since_split, color = mig_rate_m12, group = mig_rate_m12, y = fst)) +
  geom_point_rast(size = 0.75, alpha = 0.15, shape = 20)+
  #geom_point(size = 0.5, alpha = 0.1)+
  geom_smooth(color = "black", se = FALSE)+
  labs(x = "Generations Since Migration Resumed (thousands)", y = "Weir and Cockerham's FST (pop1 vs pop2)")+
  facet_grid(gf_gen_label~mig_rate_m12, scales = "free_y")+
  theme_publication()+
  theme(legend.position = "none")+
  theme(panel.grid.minor = element_line(colour = "grey75", size = 0.25),
        panel.grid.major = element_line(colour = "grey75", size = 0.25))+
  scale_color_viridis_d()

ggsave("figures/FigureS1.pdf", figureS1, device = "pdf", height = 8, width = 8)
