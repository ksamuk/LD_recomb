#!/bin/bash


#SBATCH --mem=52000
#SBATCH --cpus-per-task=10
#SBATCH -J ld_calc
#SBATCH -o slurm-ld_calc.out
#SBATCH -p noor


Rscript 03_compute_ld.R
##SBATCH -p noor
