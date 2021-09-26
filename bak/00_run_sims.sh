#!/bin/bash

#SBATCH --mem=220000
#SBATCH --cpus-per-task=20
#SBATCH -p noor
#SBATCH -J sims_split
#SBATCH -o slurm-sims_split.out

Rscript 01_run_recomb_evolve_simulations.R

