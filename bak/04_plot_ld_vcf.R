# plot LD for a slim sim

library("tidyverse")

# the actual (scaled) recomb rates (real values are rates*1e-8)
rates <-  c(10, 8, 6, 4, 2, 1, 2, 4, 6, 8, 10)
ends <-  c(0, 909090 ,1818180, 2727270, 3636360, 4545450, 5454540, 6363630, 7272720, 8181810, 9090900, 10000000)

ld_dat <- read.table("data/vcftools_output/ld_windows.txt", h = T)

ld_dat <- ld_dat %>%
  select(POS1, POS2, R.2) 
  
names(ld_dat) <- c("pos1", "pos2", "r2")

ld_dat <- ld_dat %>%
  group_by(pos1) %>%
  summarise(mean_r2 = mean(r2, na.rm = TRUE))



ld_dat_slope <- ld_dat %>%
  filter(complete.cases(ld_dat)) %>%
  group_by(pos1) %>%
  do(slope_r2 = lm(.$r2 ~ .$pos2, data = .))

ld_dat_coeffs <- 

# assign the windows
ld_dat$window1 <- NA
ld_dat$window2 <- NA

for (i in 1:length(ends)){
  
  ld_dat$window2[is.na(ld_dat$window2) & ld_dat$pos1 <= ends[i]] <- ends[i]
  ld_dat$window1[is.na(ld_dat$window1) & ld_dat$pos1 <= ends[i]] <- ends[i-1]
  
}

ld_dat %>%
  rowwise %>%
  mutate(wind_mind = mean(window1, window2)) %>%
  ungroup %>%
  ggplot(aes(x = pos1, y = mean_r2)) +
  geom_point()+
  geom_smooth()

