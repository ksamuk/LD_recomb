#!/bin/bash

#SBATCH -p noor
#SBATCH -J slim_launcher
#SBATCH --mail-user=ksamuk@gmail.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o tmp/slim_RE-%j.out


Rscript 01_run_recomb_evolve_realistic_rates_gene_flow.R

