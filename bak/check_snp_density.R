# check snp density

library("vcfR")
library("dplyr")
library("tidyr")

# read in the VCF (dropping unnecessary fields)

args <- list(10000000, "data/slim_output/output_vcf.vcf", "test")
slim_vcf <- read.vcfR(args[[2]])

slim_vcf <- data.frame(pos = as.numeric(slim_vcf@fix[,2]), data.frame(slim_vcf@gt)[,-1]) %>%
  arrange(pos)

slim_vcf %>%
  ggplot(aes(x = pos))+
  geom_histogram()

