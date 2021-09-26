library("tidyverse")

options(scipen = 999)

sim_df <- read.table("simulation_report.txt", h = T)

head(sim_df)

sim_df_sum <- sim_df %>%
  mutate(gen = as.numeric(gsub("^p._", "", file) %>% gsub("[^0-9]", "", .))) %>%
  group_by(sim_params, seed) %>%
  summarise(max_gen = as.numeric(max(gen, na.rm = TRUE)))

sim_df_sum$seed %>% unique

#View(sim_df_sum)

sim_df_sum %>%
  mutate(seed = reorder(seed, max_gen)) %>%
  ggplot(aes(x = seed, y = max_gen))+
  geom_bar(stat = "identity")+
  facet_wrap(~sim_params, scales = "free_x")