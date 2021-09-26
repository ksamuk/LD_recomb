

library("tidyverse")
library("plotly")
# plot summarized simulation results

options(scipen = 999)

ld_df <- read.table("data/sim_ld/ld_summary_out.txt", h = T)

ld_df <- ld_df %>%
  mutate(run_pop = paste0(pop , "_", run)) %>%
  mutate(param_pop = paste0(sim_param, "_", pop)) 

#######

ld_df_sumry <- ld_df %>%
  group_by(sim_param, pop, param_pop, gen) %>%
  summarise(sim_avg_ld = mean(mean_ld)) %>%
  filter(grepl("1000_1000_0.00000175_0.00003075_0.5_1_100000_1000", param_pop)) %>%
  filter(grepl("1000_1000_0.00000175_0.00003075_0.5_1_100000_1000", param_pop)) %>%
  #filter(grepl("0.00003075|500", param_pop)) %>%
  filter(!grepl("no_split", param_pop)) %>%
  filter(!grepl("gene_flow", param_pop)) %>%
  filter(gen %% 1000 == 0)

# whole trace
tmp <- ld_df %>%
  #filter(grepl("8455", run)) %>%
  #filter(grepl("1000_1000_0.00000175_0.00003075_0.5_-1_100000_1000", param_pop)) %>%
  filter(grepl("1000_1000_0.00000175_0.00003075_0.5_1_100000_1000", param_pop)) %>%
  #filter(grepl("0.00003075|500", param_pop)) %>%
  filter(!grepl("no_split", param_pop)) %>%
  filter(!grepl("gene_flow", param_pop)) %>%
  filter(gen %% 1000 == 0) %>%
  #filter(param_pop == "data/slim_output/1000_1000_0.00000175_0.000003_0.5_1_1000000_1000_p1") %>%
    ggplot(aes(x = gen, y = mean_ld, color = pop)) +
    geom_line(aes(group = run_pop), alpha = 0.2)+
    #geom_line(data = ld_df_sumry, aes(x = gen, y = sim_avg_ld), size = 2)+
    #geom_smooth(se = F)+
    #geom_smooth(method= "lm", formula= (y ~ exp(x)))+
    #facet_wrap(~sim_param, scale = "free_y")+
    facet_grid(sim_param~run, scale = "free_y")+
    theme(legend.position = "none")+
    scale_color_brewer(type = "qual", palette = "Set1")+
    xlab("Generation")+
    ylab("Mean LD")

ggplotly(tmp)

ld_df_sumry <- ld_df %>%
  group_by(sim_param, pop, param_pop, gen) %>%
  summarise(sim_avg_ld = mean(mean_ld)) %>%
  filter(grepl("0.00003075|500", param_pop)) %>%
  filter(!grepl("no_split", param_pop))
    
ld_df %>%
  filter(grepl("0.00003075|500", param_pop)) %>%
  filter(!grepl("no_split", param_pop)) %>%
  #filter(param_pop == "data/slim_output/1000_1000_0.00000175_0.000003_0.5_1_1000000_1000_p1") %>%
  ggplot(aes(x = gen, y = mean_ld, color = pop)) +
  geom_line(aes(group = run_pop), alpha = 0.2)+
  geom_line(data = ld_df_sumry, aes(x = gen, y = sim_avg_ld), size = 2)+
  #geom_smooth(se = F)+
  #geom_smooth(method= "lm", formula= (y ~ exp(x)))+
  facet_wrap(~sim_param)+
  theme(legend.position = "none")+
  scale_color_brewer(type = "qual", palette = "Set1")+
  xlim(c(20000,20100))

# GENE FLOW

# whole trace

ld_df_sumry <- ld_df %>%
  group_by(sim_param, pop, param_pop, gen) %>%
  summarise(sim_avg_ld = mean(mean_ld)) %>%
  filter(grepl("1000_1000_0.00000175_0.00003075", param_pop)) %>%
  filter(!grepl("1000_1000_0.00000175_0.00003075_0.5_1_100000_1000_no_split_p1", param_pop)) %>%
  filter(gen %% 1000 == 0)


ld_df %>%
  filter(grepl("1000_1000_0.00000175_0.00003075", param_pop)) %>%
  filter(!grepl("1000_1000_0.00000175_0.00003075_0.5_1_100000_1000_no_split_p1", param_pop)) %>%
  filter(gen %% 1000 == 0) %>%
  #filter(grepl("gene_flow", param_pop)) %>%
  #filter(!grepl("no_split", param_pop)) %>%
  #filter(param_pop == "data/slim_output/1000_1000_0.00000175_0.000003_0.5_1_1000000_1000_p1") %>%
  ggplot(aes(x = gen, y = mean_ld, color = pop)) +
  geom_line(aes(group = run_pop), alpha = 0.2)+
  geom_line(data = ld_df_sumry, aes(x = gen, y = sim_avg_ld), size = 2)+
  #geom_smooth(se = F)+
  #geom_smooth(method= "lm", formula= (y ~ exp(x)))+
  facet_wrap(~sim_param)+
  theme(legend.position = "none")+
  scale_color_brewer(type = "qual", palette = "Set1")+
  xlab("Generation")+
  ylab(expression("Mean linkage disequilibrium ("*r^{2}*")"))


# just change point

ld_df_sumry <- ld_df %>%
  group_by(sim_param, pop, param_pop, gen) %>%
  summarise(sim_avg_ld = mean(mean_ld)) %>%
  filter(grepl("1000_1000_0.00000175_0.00003075", param_pop)) %>%
  filter(!grepl("1000_1000_0.00000175_0.00003075_0.5_1_100000_1000_no_split_p1", param_pop))


ld_df %>%
  filter(grepl("1000_1000_0.00000175_0.00003075", param_pop)) %>%
  filter(!grepl("1000_1000_0.00000175_0.00003075_0.5_1_100000_1000_no_split_p1", param_pop)) %>%
  #filter(grepl("gene_flow", param_pop)) %>%
  #filter(!grepl("no_split", param_pop)) %>%
  #filter(param_pop == "data/slim_output/1000_1000_0.00000175_0.000003_0.5_1_1000000_1000_p1") %>%
  ggplot(aes(x = gen, y = mean_ld, color = pop)) +
  geom_line(aes(group = run_pop), alpha = 0.2)+
  geom_line(data = ld_df_sumry, aes(x = gen, y = sim_avg_ld), size = 2)+
  #geom_smooth(se = F)+
  #geom_smooth(method= "lm", formula= (y ~ exp(x)))+
  facet_wrap(~sim_param)+
  theme(legend.position = "none")+
  scale_color_brewer(type = "qual", palette = "Set1")+
  xlim(c(19995,20200))+
  xlab("Generation")+
  ylab(expression("Mean linkage disequilibrium ("*r^{2}*")"))

### CONTRACTION

ld_df_sumry <- ld_df %>%
  group_by(sim_param, pop, param_pop, gen) %>%
  summarise(sim_avg_ld = mean(mean_ld)) %>%
  filter(gen %% 1000 == 0) %>%
  filter(grepl("1000_500", param_pop))

ld_df %>%
  filter(grepl("1000_500", param_pop)) %>%
  filter(gen %% 1000 == 0) %>%
  ggplot(aes(x = gen, y = mean_ld, color = pop)) +
  geom_line(aes(group = run_pop), alpha = 0.2)+
  geom_line(data = ld_df_sumry, aes(x = gen, y = sim_avg_ld), size = 2)+
  facet_wrap(~sim_param)+
  theme(legend.position = "none")+
  scale_color_brewer(type = "qual", palette = "Set1")+
  #xlim(c(19990,21000))+
  xlab("Generation")+
  ylab(expression("Mean linkage disequilibrium ("*r^{2}*")"))
