#!/bin/bash

#SBATCH --mem=20000
#SBATCH --cpus-per-task=5
#SBATCH -J summarize
#SBATCH -o slurm-summarize.out


Rscript 04_summarize_ld.R

# #SBATCH -p noor