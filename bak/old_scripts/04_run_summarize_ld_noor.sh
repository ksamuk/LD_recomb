#!/bin/bash

#SBATCH --mem=52000
#SBATCH --cpus-per-task=14
#SBATCH -p noor
#SBATCH -J smrz_ld

Rscript 04_summarize_ld.R

