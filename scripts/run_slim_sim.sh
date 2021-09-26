#!/bin/bash
#SBATCH -J slim_RE-%j
#SBATCH -o tmp/slim_RE/slim_RE-%j.out

# launches a single replicate of the recombination rate evolution simulation
# used by Rscript launcher to launch many replicates of each paramter combination
# individual sims can be run like:
# sh scripts/run_slim_sim 1000 1000 1.75e-06 3.075e-05 0.5 1 1e+05 1000 0.001 data/slim_output/1000_1000_1.75e-06_3.075e-05_0.5_1_1e+05_1000_0.001/

slim -d p1_si=$1 \
-d p2_si=$2 \
-d mutation_rate=$3 \
-d recomb_rate=$4 \
-d recomb_rate_cha=$5 \
-d recomb_rate_change_dir=$6 \
-d chr_length=$7 \
-d samp_interval=$8 \
-d mig_rate=$9 \
-d out_folder_ba=\'${10}\' \
slim_sim/recomb_rate_evolve.slim