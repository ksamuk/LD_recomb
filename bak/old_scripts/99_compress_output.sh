#!/bin/bash

#SBATCH --mem=22000
#SBATCH --cpus-per-task=20
#SBATCH -J compress
#SBATCH -o slurm-compress.out
#SBATCH -p noor

time=$1


Rscript functions/compress_completed_slim_files.R $time
# #SBATCH -p noor
