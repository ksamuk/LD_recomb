##########################################
# Initialize libraries and read in 
# summarized pyrho estimates
##########################################

library("tidyverse")
library("lme4")
library("car")
library("lmerTest")
library("visreg")
library("sjPlot")
library("report")

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

dat_sub <- rmap_df_summary %>%
  filter(pop == 1) %>%
  filter(gene_flow_generation == 20000 | gene_flow_generation == 21000) %>%
  filter(ne_control == "split_aware") %>%
  filter(recomb_rate_cha == 0)

##########################################
# statistical comparisons of results
##########################################

# PROSE:
# First, in the model of continuous gene flow, when NeM >= 1, 
# we observed a systematic increase (overestimate) in estimated rates
# of recombination in both populations (Figure 2A, 2B, top row, Nem = 1-100). 
# This increase was statistically significant

# STATS:

dat_sub <- rmap_df_summary %>%
  filter(pop == 1) %>%
  filter(gene_flow_generation == 20000) %>%
  filter(ne_control == "split_aware") %>%
  filter(recomb_rate_cha == 0)

dat_sub <- dat_sub %>%
  mutate(recomb_est = scale(recomb_est)) %>%
  mutate(mig_rate_m12 = scale(mig_rate_m12)) %>%
  mutate(gene_flow_generation = as.factor(gene_flow_generation))

mod1 <- lmer(data = dat_sub, recomb_est ~ mig_rate_m12 + (1|gen) + (1|rep)) 
Anova(mod1, type = "III")

tab_model(mod1, show.df = TRUE, show.fstat = TRUE, df.method = "satterthwaite", file = "figures/TableS1.html")

report(mod1)

# REPORT:
#We fitted a linear mixed model (estimated using REML and nloptwrap optimizer) to predict recomb_est with mig_rate_m12 (formula: recomb_est ~ mig_rate_m12). The model included gen and rep as random effects (formula: list(~1 | gen, ~1 | rep)). The model's total explanatory power is substantial (conditional R2 = 0.54) and the part related to the fixed effects alone (marginal R2) is of 0.41. The model's intercept, corresponding to mig_rate_m12 = 0, is at -0.05 (95% CI [-0.16, 0.05], t(19495) = -1.00, p = 0.316). Within this model:
#- The effect of mig_rate_m12 is statistically significant and positive (beta = 0.65, 95% CI [0.63, 0.67], t(19495) = 71.34, p < .001; Std. beta = 0.65, 95% CI [0.63, 0.67])
#Standardized parameters were obtained by fitting the model on a standardized version of the dataset. 95% Confidence Intervals (CIs) and p-values were computed using the Wald approximation.


################################################################################

# PROSE:
# In contrast to the continuous gene flow case, under a model of secondary 
# contact, there was a marked systematic decrease (underestimate) of recombination 
# rates, which also became visible when Nem >=1 (Figure 2A, 2B, 
# bottom row, Nem = 1-100). This decrease was statistically significant

# STATS:

dat_sub <- rmap_df_summary %>%
  filter(pop == 1) %>%
  filter(gene_flow_generation == 21000) %>%
  filter(ne_control == "split_aware") %>%
  filter(recomb_rate_cha == 0)

dat_sub <- dat_sub %>%
  mutate(recomb_est = scale(recomb_est)) %>%
  mutate(mig_rate_m12 = scale(mig_rate_m12)) %>%
  mutate(gene_flow_generation = as.factor(gene_flow_generation))

mod1 <- lmer(data = dat_sub, recomb_est ~ mig_rate_m12 + (1|gen) + (1|rep)) 
Anova(mod1, type = "III")
report(mod1)



tab_model(mod1, show.df = TRUE, show.fstat = TRUE, df.method = "satterthwaite", file = "figures/TableS2.html")

# REPORT:
# Response: recomb_est
#Chisq Df Pr(>Chisq)    
#(Intercept)      0  1     0.9995    
#mig_rate_m12  1512  1     <2e-16 ***

#We fitted a linear mixed model (estimated using REML and nloptwrap optimizer) to predict recomb_est with mig_rate_m12 (formula: recomb_est ~ mig_rate_m12). The model included gen and rep as random effects (formula: list(~1 | gen, ~1 | rep)). The model's total explanatory power is substantial (conditional R2 = 0.53) and the part related to the fixed effects alone (marginal R2) is of 0.27. The model's intercept, corresponding to mig_rate_m12 = 0, is at -3.55e-05 (95% CI [-0.12, 0.12], t(31846) = -5.94e-04, p > .999). Within this model:
# - The effect of mig_rate_m12 is statistically significant and negative (beta = -0.52, 95% CI [-0.54, -0.49], t(31846) = -38.88, p < .001; Std. beta = -0.52, 95% CI [-0.54, -0.49])
# Standardized parameters were obtained by fitting the model on a standardized version of the dataset. 95% Confidence Intervals (CIs) and p-values were computed using the Wald approximation.

################################################################################


# PROSE
# This decrease was accompanied by a clear increase in the variance of 
# recombination rate estimates, especially for Nem = 1 and Nem = 10 (Fig 2A, bottom row)

dat_sub <- rmap_df_summary %>%
  filter(pop == 1) %>%
  filter(gene_flow_generation == 21000) %>%
  filter(ne_control == "split_aware") %>%
  filter(recomb_rate_cha == 0)

dat_sub <- dat_sub %>%
  mutate(recomb_est = scale(recomb_est)) %>%
  mutate(mig_rate_m12_scale = scale(mig_rate_m12)) %>%
  mutate(gene_flow_generation = as.factor(gene_flow_generation))


nem_lt_1 <- dat_sub %>%
  filter(mig_rate_m12 < 1e-04) %>%
  pull(recomb_est)

nem_1_10 <- dat_sub %>%
  filter(mig_rate_m12 == 1e-04 | mig_rate_m12 == 1e-03) %>%
  pull(recomb_est)

var.test(nem_lt_1, nem_1_10, alternative = "two.sided")

# REPORT:
# F test to compare two variances
# 
# data:  nem_lt_1 and nem_1_10
# F = 0.20863, num df = 10429, denom df = 13860, p-value < 2.2e-16
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   0.2012793 0.2162766
# sample estimates:
#   ratio of variances 
# 0.2086336 


################################################################################

# PROSE
# When recombination rates diverged between populations, we also observed the 
# two forms of bias described above (Figure 3). The estimates from the 
# continuous gene flow scenario exhibited a statistically significant 
# increase () whereas estimates from the secondary contact model exhibited a 
# statistically significant decrease ().

dat_sub <- rmap_df_summary %>%
  filter(pop == 1) %>%
  filter(gene_flow_generation == 20000) %>%
  filter(ne_control == "split_aware") %>%
  filter(recomb_rate_cha == 1)

dat_sub <- dat_sub %>%
  mutate(recomb_est = scale(recomb_est)) %>%
  mutate(mig_rate_m12_scale = scale(mig_rate_m12)) %>%
  mutate(gene_flow_generation = as.factor(gene_flow_generation))


mod1 <- lmer(data = dat_sub, recomb_est ~ mig_rate_m12_scale + (1|gen) + (1|rep)) 

dat_sub <- rmap_df_summary %>%
  filter(pop == 1) %>%
  filter(gene_flow_generation == 21000) %>%
  filter(ne_control == "split_aware") %>%
  filter(recomb_rate_cha == 1)

dat_sub <- dat_sub %>%
  mutate(recomb_est = scale(recomb_est)) %>%
  mutate(mig_rate_m12_scale = scale(mig_rate_m12)) %>%
  mutate(gene_flow_generation = as.factor(gene_flow_generation))


mod2 <- lmer(data = dat_sub, recomb_est ~ mig_rate_m12_scale + (1|gen) + (1|rep)) 

Anova(mod1, type = "III")
as.report_text(report(mod1), summary = TRUE)

Anova(mod2, type = "III")
as.report_text(report(mod2), summary = TRUE)


# continuous gene flow

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: recomb_est
# Chisq Df Pr(>Chisq)    
# (Intercept)    42.485  1  7.121e-11 ***
# mig_rate_m12 8936.441  1  < 2.2e-16 ***
#   
# We fitted a linear mixed model to predict recomb_est with mig_rate_m12_scale. 
# The model included gen and rep as random effects. The model's total explanatory power is substantial 
# (conditional R2 = 0.56) and the part related to the fixed effects alone (marginal R2) is of 0.43. 
# The model's intercept is at -0.03 (95% CI [-0.14, 0.08]). Within this model:
  
# - The effect of mig_rate_m12_scale is statistically significant and positive 
# (beta = 0.66, 95% CI [0.65, 0.67], t(19495) = 94.53, p < .001, Std. beta = 0.66)> 


# secondary contact
# 
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: recomb_est
# Chisq Df Pr(>Chisq)    
# (Intercept)    2.1733  1     0.1404    
# mig_rate_m12 538.9673  1     <2e-16 ***
#   ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# We fitted a linear mixed model to predict recomb_est with mig_rate_m12. 

# The model included gen and rep as random effects. 
# The model's total explanatory power is substantial (conditional R2 = 0.43) 
# and the part related to the fixed effects alone (marginal R2) is of 0.06. 
# The model's intercept is at 5.33e-13 (95% CI [-0.18, 0.18]). Within this model:
  
# - The effect of mig_rate_m12 is statistically significant and negative 
# (beta = -0.25, 95% CI [-0.27, -0.22], t(34505) = -23.22, p < .001, Std. beta = -0.25)

################################################################################

# In keeping with the models with no recombination divergence, 
# starting at Nem ~ 1,  migration-associated LD resulted in the systematic 
# underestimation and increase in variance for estimated recombination rates 
 # within both p1 and p2 (Figure 3B, Secondary Contact)

# bias

# STATS:

dat_sub <- rmap_df_summary %>%
  filter(pop == 1) %>%
  filter(gene_flow_generation == 21000) %>%
  filter(ne_control == "split_aware") %>%
  filter(recomb_rate_cha == 1)

dat_sub <- dat_sub %>%
  mutate(recomb_est = scale(recomb_est)) %>%
  mutate(mig_rate_m12 = scale(mig_rate_m12)) %>%
  mutate(gene_flow_generation = as.factor(gene_flow_generation))

mod1 <- lmer(data = dat_sub, recomb_est ~ mig_rate_m12 + (1|gen) + (1|rep)) 
Anova(mod1, type = "III")
report(mod1)

# REPORT:
# Response: recomb_est
#Chisq Df Pr(>Chisq)    
#(Intercept)      0  1     0.9995    
#mig_rate_m12  1512  1     <2e-16 ***

#We fitted a linear mixed model (estimated using REML and nloptwrap optimizer) to predict recomb_est with mig_rate_m12 (formula: recomb_est ~ mig_rate_m12). The model included gen and rep as random effects (formula: list(~1 | gen, ~1 | rep)). The model's total explanatory power is substantial (conditional R2 = 0.53) and the part related to the fixed effects alone (marginal R2) is of 0.27. The model's intercept, corresponding to mig_rate_m12 = 0, is at -3.55e-05 (95% CI [-0.12, 0.12], t(31846) = -5.94e-04, p > .999). Within this model:
# - The effect of mig_rate_m12 is statistically significant and negative (beta = -0.52, 95% CI [-0.54, -0.49], t(31846) = -38.88, p < .001; Std. beta = -0.52, 95% CI [-0.54, -0.49])
# Standardized parameters were obtained by fitting the model on a standardized version of the dataset. 95% Confidence Intervals (CIs) and p-values were computed using the Wald approximation.



# variance increase

dat_sub <- rmap_df_summary %>%
  filter(pop == 1) %>%
  filter(gene_flow_generation == 21000) %>%
  filter(ne_control == "split_aware") %>%
  filter(recomb_rate_cha == 1)

dat_sub <- dat_sub %>%
  mutate(recomb_est = scale(recomb_est)) %>%
  mutate(mig_rate_m12_scale = scale(mig_rate_m12)) %>%
  mutate(gene_flow_generation = as.factor(gene_flow_generation))


nem_lt_1 <- dat_sub %>%
  filter(mig_rate_m12 < 1e-04) %>%
  pull(recomb_est)

nem_1_10 <- dat_sub %>%
  filter(mig_rate_m12 == 1e-04 | mig_rate_m12 == 1e-03) %>%
  pull(recomb_est)

var.test(nem_lt_1, nem_1_10, alternative = "two.sided")

# REPORT:
# F test to compare two variances
# 
# data:  nem_lt_1 and nem_1_10
# F = 0.20863, num df = 10429, denom df = 13860, p-value < 2.2e-16
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   0.2012793 0.2162766
# sample estimates:
#   ratio of variances 
# 0.2086336 



