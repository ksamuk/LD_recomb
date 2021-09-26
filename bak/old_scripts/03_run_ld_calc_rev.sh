#!/bin/bash


#SBATCH --mem=52000
#SBATCH --cpus-per-task=10
#SBATCH -J ld_rev
#SBATCH -o slurm-ld_rev.out
#SBATCH -p noor

Rscript 03_compute_ld_rev.R
##SBATCH -p noor
