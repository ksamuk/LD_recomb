#!/bin/bash

sim_name=$1
num_cores=$2 
out_file_slug=$3

ldhelmet max_lk \
--num_threads $num_cores \
-l data/ldhelmet_output/$sim_name/${out_file_slug}/${out_file_slug}.lk \
-p data/ldhelmet_output/$sim_name/${out_file_slug}/${out_file_slug}.pade \
-m data/fasta_output/${sim_name}_mut_mat.txt \
-a data/fasta_output/${sim_name}_anc_prior.txt \
-s data/fasta_output/$sim_name.fasta
