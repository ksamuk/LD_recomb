# convert a VCF created by SLiM 2.0 to a fasta file
# KS Aug 2017

library("vcfR")
library("dplyr")
library("tidyr")

# parse args
args <- commandArgs(trailingOnly = TRUE)

if(is.null(args) | length(args) == 0){

	args <- list(10000000, "data/slim_output/output_vcf.vcf", "test")
}

# initialize a blank (all A) sequnece of appropriate length
seq_len <- args[[1]]
out_file_slug <- args[[3]]
blank_seq <- rep("A", seq_len)


# read in the VCF (dropping unnecessary fields)
slim_vcf <- read.vcfR(args[[2]])

slim_vcf <- data.frame(pos = as.numeric(slim_vcf@fix[,2]), data.frame(slim_vcf@gt)[,-1]) %>%
  arrange(pos)

# function for dividing phased genotypes in two haplotypes
split_column <- function(x){
  
  lab <- names(slim_vcf[,-1])[x]
  x <- slim_vcf[,-1][,x]
  
  col_names <- c(paste0(lab, "_hap1"), paste0(lab, "_hap2"))
  
  data.frame(x) %>%
    separate(x, into = col_names, sep = "\\|")
  
}

# perform haplotype division
haps <- lapply(1:length(slim_vcf[,-1]), split_column) %>%
  bind_cols()

haps <- data.frame(pos = slim_vcf$pos, haps)

# function for writing haps to fasta format
# uses the blank_seq template, and applies the snp calls (i.e the 1s)
# as "T"s, then writes out to file

write_fasta_out <- function(index, haps_in, file_name){
  
  hap_df <- data.frame(pos = haps$pos, gt = haps_in[,-1][,index]) %>%
    filter(gt != 0)
  
  #apply the haplotype mask
  out_seq <- blank_seq
  out_seq[hap_df$pos] <- "T"
  out_seq <- paste0(out_seq, collapse = "")
  
  out_header <- paste0(">", names(haps_in[,-1])[index])
  
  # write to file
  cat(out_header, out_seq, file = file_name, sep = "\n",append = file.exists(file_name))
  
  
}


# total number of haplotypes
num_haps <- ncol(haps) - 1
out_file_name <- paste0("data/fasta_output/", out_file_slug, ".fasta")

# downsample haplotypes (optionally)
#num_haps <- floor(ncol(haps)/2)
num_haps <- 15
haps <- data.frame(pos = haps[,1], haps[,-1][,sample(names(haps[,-1]), num_haps)])

haps <- haps[-apply((haps[,-1] == 0), 1, function(x) all(x)),]

# process haplotypes to a fasta file
file.remove(out_file_name)
invisible(lapply(1:num_haps, write_fasta_out, haps_in = haps, file_name = out_file_name))

# output the ancestral allele priors (all anc alleles are A's)
# remove redundant snps (if downsampling)

anc_prior <- haps %>%
  select(pos) %>%
  mutate(pos = pos - 1) %>%
  mutate(A = 1.0, C = 0.0, G = 0.0, T = 0.0)

anc_out_name <- gsub("\\.fasta", "_anc_prior.txt", out_file_name)

write.table(anc_prior, file  = anc_out_name, col.names = FALSE, row.names = FALSE, quote = FALSE)
