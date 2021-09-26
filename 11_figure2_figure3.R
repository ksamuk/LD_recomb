##########################################
# Statistical analyses and plots 
# of pyrho recombation rate estimates
# from SLiM simulations
# 
# Code for FIGURES 2 & 3 and associated
# statistical tests
#
# Kieran Samuk
# April 2021 
# ksamuk@gmail.com
##########################################

##########################################
# Initialize libraries and read in 
# summarized pyrho estimates
##########################################

#extrafont::loadfonts(device="win")
library("tidyverse")
library("patchwork")
library("grid")
library("ggthemes")

lapply(list.files("functions", pattern = "theme", full.names = TRUE), 
       source, echo = FALSE)

# summarized pyrho estimates
rmap_df <- read_rds("data/pyrho_processed_df.rds") %>%
  ungroup

# the true simulated value of rho
emp_rho <- 2.23e-08

# the expected empirical values of rho
exp_rho_df_1 <- data.frame(pop = rep(1, 3), 
                           recomb_rate_cha = c(0,0.5,1), 
                           exp_rho = emp_rho, linecol = "1")
exp_rho_df_2 <- data.frame(pop = rep(2 ,3), 
                           recomb_rate_cha = c(0,0.5,1), 
                           exp_rho = c(emp_rho, emp_rho * 1.5, emp_rho * 2), 
                           linecol = "2")

exp_rho_df <- bind_rows(exp_rho_df_1, exp_rho_df_2) %>%
  distinct 

# data frame of the expected (average) recombination rates with *no migration*
allo_exp <- rmap_df %>%
  filter(mig_rate_m12 == 0, mig_rate_m21 == 0) %>% 
  group_by(ne_control, gene_flow_generation, recomb_rate_change_dir, 
           recomb_rate_cha, p2_prop, pop, gen) %>%
  summarise(allo_exp_rho = mean(recomb_est, na.rm = TRUE)) %>%
  type_convert() %>%
  ungroup

# in the absence of a change in recombination, 
# average the expected value for pop1 and 2

allo_exp <- allo_exp %>%
  filter(recomb_rate_cha == 0) %>%
  group_by(ne_control, gene_flow_generation, recomb_rate_change_dir, 
           recomb_rate_cha, p2_prop, gen) %>%
  summarise(allo_exp_rho2 = mean(allo_exp_rho, na.rm = TRUE)) %>%
  left_join(allo_exp, .) %>%
  mutate(allo_exp_rho = ifelse(is.na(allo_exp_rho2), 
                               allo_exp_rho, allo_exp_rho2)) %>%
  select(-allo_exp_rho2)

# join in empirical and allopatric expected values
rmap_df_summary <- left_join(type_convert(rmap_df), exp_rho_df)
rmap_df_summary <- left_join(rmap_df_summary, allo_exp)

# ungroup and appaned string to rep to force factorization
rmap_df_summary <- rmap_df_summary %>%
  ungroup %>%
  mutate(rep = paste0("r",rep)) %>%
  type_convert

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

##########################################
# Figure 2: recomb estimates without 
# a change in recombination rate
##########################################

# figure 2a: simulation results

fig2_sub <- rmap_df_summary %>%
  filter(recomb_rate_cha == 0) %>%
  filter(mig_rate_m12 == mig_rate_m21) %>%
  mutate(mig_rate_m12 = as.numeric(mig_rate_m12) * 1000) %>%
  mutate(gen_since_split = gen - gene_flow_generation) %>%
  filter(gen_since_split >= 0, gen_since_split <= 30) %>%
  filter(p2_prop == 0.5) %>%
  filter(gene_flow_generation == 20000 | gene_flow_generation == 21000) %>% 
  mutate(gf_gen_label = ifelse(gene_flow_generation == 20000, "Continuous Migration", "Secondary Contact")) %>%
  mutate(rep_pop = paste0(rep, pop))

fig2_allo_exp <- fig2_sub %>%
  select(-rep_pop, -rep, -recomb_est) %>%
  distinct

fig2a <- fig2_sub %>%
  ggplot(aes(x = gen_since_split, y = recomb_est, group = rep_pop, color = as.factor(pop)))+
  geom_line(alpha = 0.1, size = 0.5)+
  stat_smooth(data = fig2_allo_exp, aes(x = gen_since_split, y = allo_exp_rho, group = 1), linetype = 1, size = 1, se = FALSE, color = "black", span = 1.0, n = 10) +
  facet_grid(gf_gen_label~mig_rate_m12)+
  labs(x = "Generations since divergence (thousands)", y = "Inferred Recomb. Rate (cM/Mb)")+
  scale_color_brewer(palette = "Set1")+
  theme_publication()+
  theme(legend.position = "none")+
  scale_y_continuous(label=scientific_10)

# figure 2b: deviation from allopatry

fig2b <- rmap_df_summary %>%
  filter(recomb_rate_cha == 0) %>%
  filter(mig_rate_m12 == mig_rate_m21) %>%
  mutate(mig_rate_m12 = as.numeric(mig_rate_m12) * 1000) %>%
  mutate(gen_since_split = gen - gene_flow_generation) %>%
  filter(gen_since_split >= 1, gen_since_split <= 30) %>%
  filter(p2_prop == 0.5) %>%
  filter(gene_flow_generation == 20000 | gene_flow_generation == 21000) %>%
  mutate(gf_gen_label = ifelse(gene_flow_generation == 20000, "Continuous Migration", "Secondary Contact")) %>%
  mutate(rep_pop = paste0(rep, pop)) %>%
  select(gf_gen_label, rep, gen, recomb_rate_cha, mig_rate_m12, gen_since_split, recomb_est, allo_exp_rho) %>%
  mutate(recomb_ratio = recomb_est/allo_exp_rho) %>%
  group_by(gf_gen_label, mig_rate_m12) %>%
  summarise(mean_diff = mean(recomb_ratio, na.rm = TRUE), sd_diff = sd(recomb_ratio, na.rm = TRUE)) %>%
  mutate(sd_min = mean_diff - sd_diff) %>%
  mutate(sd_max = mean_diff + sd_diff) %>%
  ggplot(aes(x = as.factor(mig_rate_m12), y = mean_diff, 
             ymin = sd_min, ymax = sd_max, group = 1))+
  geom_point()+
  geom_errorbar(width = 0)+
  geom_line()+
  geom_hline(yintercept = 1, linetype = 2, color = "grey50")+
  facet_grid(gf_gen_label~.)+
  xlab("Migration Rate (Nem)")+
  ylab("Inferred / Allopatric Recomb Rate.")+
  scale_color_brewer(palette = "Set1")+
  theme_publication()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# figure 2c: magnitude of difference

pop1 <- rmap_df_summary %>%
  filter(recomb_rate_cha == 0) %>%
  filter(p2_prop == 0.5) %>%
  filter(pop == 1) %>%
  filter(mig_rate_m12 == mig_rate_m21) %>%
  mutate(mig_rate_m12 = as.numeric(mig_rate_m12) * 1000) %>%
  mutate(gen_since_split = gen - gene_flow_generation) %>%
  filter(gene_flow_generation == 20000 | gene_flow_generation == 21000) %>%
  mutate(gf_gen_label = ifelse(gene_flow_generation == 20000, "Continuous Migration", "Secondary Contact")) %>%
  rename(recomb_est_pop1 = recomb_est) %>%
  select(-allo_exp_rho, -exp_rho, -linecol, -pop)

pop2 <- rmap_df_summary %>%
  filter(recomb_rate_cha == 0) %>%
  filter(p2_prop == 0.5) %>%
  filter(pop == 2) %>%
  filter(mig_rate_m12 == mig_rate_m21) %>%
  mutate(mig_rate_m12 = as.numeric(mig_rate_m12) * 1000) %>%
  mutate(gen_since_split = gen - gene_flow_generation) %>%
  filter(gene_flow_generation == 20000 | gene_flow_generation == 21000) %>%
  mutate(gf_gen_label = ifelse(gene_flow_generation == 20000, "Continuous Migration", "Secondary Contact")) %>%
  rename(recomb_est_pop2 = recomb_est) %>%
  select(-allo_exp_rho, -exp_rho, -linecol, -pop)

left_join(pop1, pop2) %>%
  mutate(pop_diff = recomb_est_pop2/recomb_est_pop1) %>%
  select(gf_gen_label, rep, mig_rate_m12, gen, pop_diff) %>%
  group_by(gf_gen_label, mig_rate_m12) %>%
  summarise(mean_diff = mean(pop_diff, na.rm = TRUE), sd_diff = sd(pop_diff, na.rm = TRUE)) %>%
  mutate(sd_min = mean_diff - sd_diff) %>%
  mutate(sd_max = mean_diff + sd_diff) %>%
  ungroup %>%
  mutate(group_id = 1) %>%
  ggplot(aes(x = as.factor(mig_rate_m12), y = mean_diff, 
             ymin = sd_min, ymax = sd_max, group = 1))+
  geom_hline(yintercept = 2, color = "grey50", linetype = 2)+
  geom_point()+
  geom_errorbar(width = 0)+
  geom_line()+
  facet_grid(gf_gen_label~.)+
  xlab("Migration Rate (Nem)")+
  ylab("Inferred Recomb. Difference (cM/Mb)") +
  theme_publication()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

fig2c <- left_join(pop1, pop2) %>%
  mutate(pop_diff = recomb_est_pop2/recomb_est_pop1) %>%
  select(gf_gen_label, rep, mig_rate_m12, gen, pop_diff) %>%
  group_by(gf_gen_label, mig_rate_m12) %>%
  summarise(mean_diff = mean(pop_diff, na.rm = TRUE), sd_diff = sd(pop_diff, na.rm = TRUE)) %>%
  mutate(sd_min = mean_diff - sd_diff) %>%
  mutate(sd_max = mean_diff + sd_diff) %>%
  ungroup %>%
  mutate(group_id = 1) %>%
  ggplot(aes(x = as.factor(mig_rate_m12), y = mean_diff, 
             ymin = sd_min, ymax = sd_max, group = 1))+
  geom_hline(yintercept = 1, color = "grey50", linetype = 2)+
  geom_point()+
  geom_errorbar(width = 0)+
  geom_line()+
  facet_grid(gf_gen_label~.)+
  xlab("Migration Rate (Nem)")+
  ylab("Pop 2 / Pop 1 Recomb. Rate") +
  theme_publication()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

figure2 <- fig2a / (fig2b | fig2c)
ggsave("figures/Figure2.pdf", figure2, device = "pdf", height = 8, width = 8, useDingbats=FALSE)

##########################################
# Figure 3: recomb estimates *with*
# a change in recombination rate
##########################################

# figure 3a: simulation results

fig3_sub <- rmap_df_summary %>%
  filter(recomb_rate_cha == 1) %>%
  filter(mig_rate_m12 == mig_rate_m21) %>%
  mutate(mig_rate_m12 = as.numeric(mig_rate_m12) * 1000) %>%
  mutate(gen_since_split = gen - gene_flow_generation) %>%
  filter(gen_since_split >= 0, gen_since_split <= 30) %>%
  filter(p2_prop == 0.5) %>%
  filter(gene_flow_generation == 20000 | gene_flow_generation == 21000) %>% 
  mutate(gf_gen_label = ifelse(gene_flow_generation == 20000, "Continuous Migration", "Secondary Contact")) %>%
  mutate(rep_pop = paste0(rep, pop)) %>%
  filter(!is.na(pop))

fig3_allo_exp <- fig3_sub %>%
  select(-rep_pop, -rep, -recomb_est) %>%
  distinct

fig3a <- fig3_sub %>%
  ggplot(aes(x = gen_since_split, y = recomb_est, group = rep_pop, color = as.factor(pop)))+
  geom_line(alpha = 0.1, size = 0.5)+
  stat_smooth(data = fig3_allo_exp, aes(x = gen_since_split, y = allo_exp_rho, group = pop, color = as.factor(pop)), linetype = 1, size = 1, se = FALSE, span = 1.0, n = 10) +
  facet_grid(gf_gen_label~mig_rate_m12)+
  labs(x = "Generations since divergence (thousands)", y = "Inferred Recomb. Rate (cM/Mb)")+
  scale_color_brewer(palette = "Set1")+
  theme_publication()+
  theme(legend.position = "none")+
  scale_y_continuous(label=scientific_10)

# figure 3b: deviation from p1 allopatry

# subset the allopatric expectations for the first population only
# (for 3b below)
p1_exp <- allo_exp %>%
  filter(pop == 1) %>%
  select(-pop, -ne_control) %>%
  rename(allo_exp_rho_p1 = allo_exp_rho)

fig3b <- rmap_df_summary %>%
  filter(recomb_rate_cha == 1) %>%
  filter(mig_rate_m12 == mig_rate_m21) %>%
  mutate(mig_rate_m12 = as.numeric(mig_rate_m12) * 1000) %>%
  mutate(gen_since_split = gen - gene_flow_generation) %>%
  filter(p2_prop == 0.5) %>%
  left_join(p1_exp) %>%
  filter(gene_flow_generation == 20000 | gene_flow_generation == 21000) %>%
  mutate(gf_gen_label = ifelse(gene_flow_generation == 20000, "Continuous Migration", "Secondary Contact")) %>%
  mutate(rep_pop = paste0(rep, pop)) %>%
  mutate(recomb_ratio = recomb_est/allo_exp_rho_p1) %>%
  group_by(gf_gen_label, mig_rate_m12, pop) %>%
  summarise(mean_diff = mean(recomb_ratio, na.rm = TRUE), sd_diff = sd(recomb_ratio, na.rm = TRUE)) %>%
  mutate(sd_min = mean_diff - sd_diff) %>%
  mutate(sd_max = mean_diff + sd_diff) %>%
  ggplot(aes(x = as.factor(mig_rate_m12), y = mean_diff, 
             ymin = sd_min, ymax = sd_max, group = pop, color = as.factor(pop)))+
  geom_point()+
  geom_errorbar(width = 0)+
  geom_line()+
  geom_hline(yintercept = 1, linetype = 2, color = "#E41A1C")+
  geom_hline(yintercept = 2, linetype = 2, color = "#377EB8")+
  facet_grid(gf_gen_label~.)+
  xlab("Migration Rate (Nem)")+
  ylab("Inferred / (Pop1 Allopatric) Recomb Rate.")+
  scale_color_brewer(palette = "Set1")+
  theme_publication()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# figure 3c: magnitude of difference

pop1 <- rmap_df_summary %>%
  filter(recomb_rate_cha == 1) %>%
  filter(p2_prop == 0.5) %>%
  filter(pop == 1) %>%
  filter(mig_rate_m12 == mig_rate_m21) %>%
  mutate(mig_rate_m12 = as.numeric(mig_rate_m12) * 1000) %>%
  mutate(gen_since_split = gen - gene_flow_generation) %>%
  filter(gene_flow_generation == 20000 | gene_flow_generation == 21000) %>%
  mutate(gf_gen_label = ifelse(gene_flow_generation == 20000, "Continuous Migration", "Secondary Contact")) %>%
  rename(recomb_est_pop1 = recomb_est) %>%
  select(-allo_exp_rho, -exp_rho, -linecol, -pop)

pop2 <- rmap_df_summary %>%
  filter(recomb_rate_cha == 1) %>%
  filter(p2_prop == 0.5) %>%
  filter(pop == 2) %>%
  filter(mig_rate_m12 == mig_rate_m21) %>%
  mutate(mig_rate_m12 = as.numeric(mig_rate_m12) * 1000) %>%
  mutate(gen_since_split = gen - gene_flow_generation) %>%
  filter(gene_flow_generation == 20000 | gene_flow_generation == 21000) %>%
  mutate(gf_gen_label = ifelse(gene_flow_generation == 20000, "Continuous Migration", "Secondary Contact")) %>%
  rename(recomb_est_pop2 = recomb_est) %>%
  select(-allo_exp_rho, -exp_rho, -linecol, -pop)

fig3c <- left_join(pop1, pop2) %>%
  mutate(pop_diff = recomb_est_pop2/recomb_est_pop1) %>%
  select(gf_gen_label, rep, mig_rate_m12, gen, pop_diff) %>%
  group_by(gf_gen_label, mig_rate_m12) %>%
  summarise(mean_diff = mean(pop_diff, na.rm = TRUE), sd_diff = sd(pop_diff, na.rm = TRUE)) %>%
  mutate(sd_min = mean_diff - sd_diff) %>%
  mutate(sd_max = mean_diff + sd_diff) %>%
  ungroup %>%
  mutate(group_id = 1) %>%
  ggplot(aes(x = as.factor(mig_rate_m12), y = mean_diff, 
             ymin = sd_min, ymax = sd_max, group = 1))+
  geom_hline(yintercept = 2, color = "grey50", linetype = 2)+
  geom_point()+
  geom_errorbar(width = 0)+
  geom_line()+
  facet_grid(gf_gen_label~.)+
  xlab("Migration Rate (Nem)")+
  ylab("Inferred Recomb. Difference (cM/Mb)") +
  theme_publication()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


figure3 <- fig3a / (fig3b | fig3c)
ggsave("figures/Figure3.pdf", figure3, device = "pdf", height = 8, width = 8, useDingbats=FALSE)

# create slides of figures
read_pptx() %>% 
  add_slide(layout = "Title and Content", master = "Office Theme") %>% 
  ph_with_vg_at(code = print(fig1), left = 0, top = 0, width = 8, height = 6) %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>% 
  ph_with_vg_at(code = print(fig2), left = 0, top = 0, width = 8, height = 6) %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>% 
  ph_with_vg_at(code = print(fig3), left = 0, top = 0, width = 8, height = 6) %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>% 
  ph_with_vg_at(code = print(fig4), left = 0, top = 0, width = 8, height = 6) %>%
  print(target = "lab_meeting.pptx")



