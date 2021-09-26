#!/bin/bash

#SBATCH --mem=220000
#SBATCH --cpus-per-task=20
#SBATCH -p noor
#SBATCH -J real_nosplit
#SBATCH -o slurm-sims_realistic_no_split.out

Rscript 01_run_recomb_evolve_realistic_rates_no_split.R

