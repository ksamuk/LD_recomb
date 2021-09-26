library("tidyverse")
library("extrafont")
library("plotly")
library("officer")
library("aod")
library("rvg")
source("functions/theme_presentation.R")
loadfonts(device = "win", quiet = TRUE)

rmap_df <- read_rds("data/pyrho_processed_df.rds")

# re-calculate the scaled empirical estimate of recomb 
emp_rho <- 8.4e-09
#dmel_ne <- 1720600
#pi_size <- 1000
#scaling_factor <- dmel_ne/p1_size
#emp_recomb_per_bp <- 0.5*(1-(1-2*(dmel_recomb)))^scaling_factor

exp_rho_df_1 <- data.frame(pop = rep(1, 3), recomb_rate_change = c(0,0.5,1), exp_rho = emp_rho*0.5, linecol = "1")
exp_rho_df_2 <- data.frame(pop = rep(2 ,3), recomb_rate_change = c(0,0.5,1), exp_rho = c(emp_rho*0.5, emp_rho * 0.75, emp_rho), linecol = "2")
exp_rho_df <-bind_rows(exp_rho_df_1, exp_rho_df_2) %>%
  mutate(recomb_rate_change = paste0(recomb_rate_change + 1, "x" )) %>%
  mutate(recomb_rate_change = factor(recomb_rate_change, levels = c("1x", "1.5x", "2x")))

rmap_df_summary <- left_join(rmap_df, exp_rho_df)

rmap_df_summary <- rmap_df_summary %>%
  ungroup %>%
  type_convert

# rrate change / migration rate
rmap_df_summary %>%
  filter(recomb_rate_change_dir == 1) %>%
  filter(mig_rate_m12 == mig_rate_m21) %>%
  filter(recomb_rate_cha == 0.5) %>%
  filter(p2_prop == 0.5) %>%
  #filter(pop == 2) %>%
  mutate(rep_pop = paste0(rep, pop)) %>%
  ggplot(aes(x = as.factor(gen), y = recomb_est, group = rep_pop, color = as.factor(pop)))+
  geom_line()+
  #geom_jitter(width = 0.1)+
  facet_grid(mig_rate_m12~gene_flow_generation)+
  ylab("Recombination Rate Estimate") +
  xlab("Generation")+
  geom_hline(yintercept = emp_rho, linetype = 2, color = "red", size = 1, alpha = 0.2) +
  geom_hline(yintercept = emp_rho*1.5, linetype = 2, color = "blue", size = 1, alpha = 0.2) +
  scale_color_brewer(palette = "Set1")+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# rrate change / migration rate
rmap_df_summary %>%
  filter(recomb_rate_change_dir == 1) %>%
  filter(mig_rate_m12 == mig_rate_m21) %>%
  #filter(recomb_rate_cha == 0.5) %>%
  filter(gene_flow_generation == 20000) %>%
  #filter(p2_prop == 0.5) %>%
  filter(pop == 1) %>%
  mutate(rep_pop = paste0(rep, pop)) %>%
  ggplot(aes(x = as.factor(gen), y = recomb_est, group = rep_pop, color = recomb_rate_cha))+
  geom_line()+
  #geom_jitter(width = 0.1)+
  facet_grid(mig_rate_m12~p2_prop)+
  ylab("Recombination Rate Estimate") +
  xlab("Generation")+
  geom_hline(yintercept = emp_rho, linetype = 2, color = "blue", size = 1, alpha = 0.2) +
  scale_color_viridis_c()+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#

# no gene flow 
# fig 1

fig1 <-  rmap_df_summary %>%
  filter(mig_rate_m12 == 0, mig_rate_m21 == 0) %>%
  filter(recomb_rate_cha == 0.0 ) %>%
  filter(recomb_rate_change_dir == 1) %>%
  filter(gene_flow_generation %in% c(20000)) %>%
  mutate(rep_pop = paste0(rep, pop)) %>%
  ggplot(aes(x = as.factor(gen), y = recomb_est, group = rep_pop, color = factor(pop)))+
  geom_line()+
  #facet_wrap(~recomb_rate_change)+
  geom_hline(yintercept = emp_rho*0.5, linetype = 2, color = "firebrick", size = 1, alpha = 0.75) +
  geom_hline(yintercept = emp_rho*0.5, linetype = 2, color = "cornflowerblue", size = 1, alpha = 0.75) +
  ylab("Recombination Rate Estimate") +
  xlab("Generation")+
  scale_color_brewer(palette = "Set1")+
  scale_y_continuous(breaks = scales::pretty_breaks())+
  labs(color = "Population")+
  labs(title = "No Migration, No Change in Recombination", subtitle = "Gene Flow t=20000")


fig2 <- rmap_df_summary %>%
  filter(migration_rate == 0) %>%
  mutate(migration_rate = paste0("Nm = ", migration_rate*500)) %>%
  #filter(recomb_rate_change == 0.0 ) %>%
  filter(recomb_rate_change_dir == 1) %>%
  filter(migration_generation %in% c(20000)) %>%
  mutate(rep_pop = paste0(rep, pop)) %>%
  mutate(recomb_rate_change = paste0(recomb_rate_change + 1, "x" )) %>%
  mutate(recomb_rate_change = factor(recomb_rate_change, levels = c("1x", "1.5x", "2x"))) %>%
  ggplot(aes(x = as.factor(gen), y = recomb_est_avg, group = rep_pop, color = pop))+
  geom_line(size = 0.75)+
  facet_wrap(~recomb_rate_change)+
  ylab("Recombination Rate Estimate") +
  xlab("Generation")+
  geom_hline(aes(yintercept =exp_rho, color = linecol), linetype = 2, size = 0.5, alpha = 0.85, data= exp_rho_df) +
  scale_color_brewer(palette = "Set1")+
  theme_presentation(base_size = 22, title_size = 24, legend_size = 18, axis_text_size = 16, 
                     panel_col = "black", bg_col = "black", grid_col = "grey25", x_margin = 0) +
  scale_y_continuous(breaks = scales::pretty_breaks())+
  labs(color = "Population")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, margin = margin(0,0,0,0)))+
  labs(title = "No Migration, Recombination Increase", subtitle = "Gene Flow t=20000")

fig3 <- rmap_df_summary %>%
  filter(migration_rate %in% c(0,2e-05)) %>%
  mutate(migration_rate = paste0("Nm = ", migration_rate*500)) %>%
  mutate(recomb_rate_change = paste0(recomb_rate_change + 1, "x" )) %>%
  mutate(recomb_rate_change = factor(recomb_rate_change, levels = c("1x", "1.5x", "2x"))) %>%
  filter(recomb_rate_change_dir == 1) %>%
  filter(migration_generation %in% c(20000)) %>%
  mutate(rep_pop = paste0(rep, pop)) %>%
  ggplot(aes(x = as.factor(gen), y = recomb_est_avg, group = rep_pop, color = pop))+
  geom_line(size = 0.5)+
  facet_grid(migration_rate~recomb_rate_change)+
  ylab("Recombination Rate Estimate") +
  xlab("Generation")+
  geom_hline(aes(yintercept =exp_rho, color = linecol), linetype = 2, size = 0.33, alpha = 1.0, data= exp_rho_df) +
  scale_color_brewer(palette = "Set1")+
  theme_presentation(base_size = 22, title_size = 24, legend_size = 18, axis_text_size = 16, 
                     panel_col = "black", bg_col = "black", grid_col = "grey25", x_margin = 0) +
  scale_y_continuous(breaks = scales::pretty_breaks())+
  labs(color = "Population")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, margin = margin(0,0,0,0)))+
  labs(title = "Migration, Recombination Increase", subtitle = "Gene Flow t=20000")

fig4 <- rmap_df_summary %>%
  filter(migration_rate %in% c(0,2e-05, 2e-04, 0.2)) %>%
  mutate(migration_rate = paste0("Nm = ", migration_rate*500)) %>%
  mutate(recomb_rate_change = paste0(recomb_rate_change + 1, "x" )) %>%
  mutate(recomb_rate_change = factor(recomb_rate_change, levels = c("1x", "1.5x", "2x"))) %>%
  filter(recomb_rate_change_dir == 1) %>%
  filter(migration_generation %in% c(20000)) %>%
  mutate(rep_pop = paste0(rep, pop)) %>%
  ggplot(aes(x = as.factor(gen), y = recomb_est_avg, group = rep_pop, color = pop))+
  geom_line(size = 0.5)+
  facet_grid(migration_rate~recomb_rate_change)+
  ylab("Recombination Rate Estimate") +
  xlab("Generation")+
  geom_hline(aes(yintercept =exp_rho, color = linecol), linetype = 2, size = 0.33, alpha = 1.0, data= exp_rho_df) +
  scale_color_brewer(palette = "Set1")+
  theme_presentation(base_size = 22, title_size = 24, legend_size = 18, axis_text_size = 16, 
                     panel_col = "black", bg_col = "black", grid_col = "grey25", x_margin = 0) +
  scale_y_continuous(breaks = scales::pretty_breaks())+
  labs(color = "Population")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, margin = margin(0,0,0,0)))+
  labs(title = "Migration, Recombination Increase", subtitle = "Gene Flow t=20000")

fig5 <- rmap_df_summary %>%
  filter(migration_rate %in% c(0,2e-05, 2e-04, 0.2)) %>%
  mutate(migration_rate = paste0("Nm = ", migration_rate*500)) %>%
  mutate(recomb_rate_change = paste0(recomb_rate_change + 1, "x" )) %>%
  mutate(recomb_rate_change = factor(recomb_rate_change, levels = c("1x", "1.5x", "2x"))) %>%  filter(recomb_rate_change_dir == 1) %>%
  filter(migration_generation %in% c(30000)) %>%
  mutate(rep_pop = paste0(rep, pop)) %>%
  ggplot(aes(x = as.factor(gen), y = recomb_est_avg, group = rep_pop, color = pop))+
  geom_line(size = 0.5)+
  facet_grid(migration_rate~recomb_rate_change)+
  ylab("Recombination Rate Estimate") +
  xlab("Generation")+
  geom_hline(aes(yintercept =exp_rho, color = linecol), linetype = 2, size = 0.33, alpha = 1.0, data= exp_rho_df) +
  scale_color_brewer(palette = "Set1")+
  theme_presentation(base_size = 22, title_size = 24, legend_size = 18, axis_text_size = 16, 
                     panel_col = "black", bg_col = "black", grid_col = "grey25", x_margin = 0) +
  scale_y_continuous(breaks = scales::pretty_breaks())+
  labs(color = "Population")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, margin = margin(0,0,0,0)))+
  labs(title = "Migration, Recombination Increase", subtitle = "Gene Flow t=30000")

read_pptx() %>% 
  add_slide(layout = "Title and Content", master = "Office Theme") %>% 
  ph_with_vg_at(code = print(fig1), left = 0, top = 0, width = 13.33, height = 7.5) %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>% 
  ph_with_vg_at(code = print(fig2), left = 0, top = 0, width = 13.33, height = 7.5) %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>% 
  ph_with_vg_at(code = print(fig3), left = 0, top = 0, width = 13.33, height = 7.5) %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>% 
  ph_with_vg_at(code = print(fig4), left = 0, top = 0, width = 13.33, height = 7.5) %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>% 
  ph_with_vg_at(code = print(fig5), left = 0, top = 0, width = 13.33, height = 7.5) %>%
  print(target = "lab_meeting.pptx")


ggsave(fig1, filename = "fig1.pdf", height = 8, width = 9, device = "pdf", useDingbats= FALSE)


rmap_df_summary %>%
  #filter(migration_rate == 0.00002) %>%
  filter(recomb_rate_change == 0.5 ) %>%
  filter(recomb_rate_change_dir == 1) %>%
  filter(migration_generation %in% c(20000)) %>%
  mutate(rep_pop = paste0(rep, pop)) %>%
  ggplot(aes(x = as.factor(gen), y = recomb_est_avg, group = rep_pop, color = pop))+
  geom_line()+
  #geom_jitter(width = 0.1)+
  facet_grid(migration_rate~recomb_rate_change)+
  ylab("Recombination Rate Estimate") +
  xlab("Generation")+
  geom_hline(yintercept = emp_rho, linetype = 2, color = "blue", size = 1, alpha = 0.2) +
  scale_color_brewer(palette = "Set1")+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))








ggsave(tmo, filename = "out.pdf", height = 8, width = 9, device = "pdf", useDingbats= FALSE)

