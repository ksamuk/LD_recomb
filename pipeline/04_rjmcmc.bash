#!/bin/bash

sim_name=$1
window_size=$2
block_penalty=$3
num_cores=$4
out_file_slug=$5
prior=$6

time ldhelmet rjmcmc \
--num_threads $num_cores -l data/ldhelmet_output/$sim_name/${out_file_slug}/${out_file_slug}.lk \
-p data/ldhelmet_output/$sim_name/${out_file_slug}/${out_file_slug}.pade \
-m data/fasta_output/${sim_name}_mut_mat.txt \
-a data/fasta_output/${sim_name}_anc_prior.txt \
-s data/fasta_output/$sim_name.fasta \
-w $window_size \
-b $block_penalty \
--prior_rate $prior \
--burn_in 100000 \
-n 1000000 \
-o data/ldhelmet_output/$sim_name/${out_file_slug}/${out_file_slug}.post
