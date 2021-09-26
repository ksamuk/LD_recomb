library("tidyverse")



sample_with_replacement <- function(index, len = 1000){
  
  x <- 1:len
  df <- data.frame(num = sample(x, len, replace = TRUE)) 
  df$num %>% unique %>% length
  
}

len_list <- lapply(1:10000, sample_with_replacement)
len_list %>% unlist %>% hist


df %>%
  group_by(num) %>%
  count %>%
  ggplot(aes(x = n))+
  geom_histogram()