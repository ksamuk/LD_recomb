

library("tidyverse")
# plot summarized simulation results

options(scipen = 999)

ld_df <- read.table("data/sim_ld/ld_summary_out.txt", h = T)

ld_df %>%
  mutate(run_pop = paste0(pop , "_", run)) %>%
  mutate(param_pop = paste0(sim_param, "_", pop)) %>%
 # filter(param_pop == "data/slim_output/1000_1000_0.00000175_0.0000015_0.5_1_1000000_1000_p1") %>%
  #filter(param_pop == "data/slim_output/1000_1000_0.00000175_0.000003_0.5_1_1000000_1000_p1") %>%
    ggplot(aes(x = gen, y = mean_ld, color = param_pop)) +
    geom_line(aes(group = run_pop), alpha = 0.1)+
    geom_smooth(se = F)+
    facet_wrap(~param_pop)+
    theme(legend.position = "none")

runs <- ld_df$run %>% unique
runs == 1509379654005
