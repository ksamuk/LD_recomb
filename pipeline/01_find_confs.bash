#!/bin/bash

sim_name=$1 
window_size=$2
num_cores=$3
out_file_slug=$4

ldhelmet find_confs --num_threads $num_cores \
-w $window_size \
-o data/ldhelmet_output/$sim_name/${out_file_slug}/${out_file_slug}.conf \
data/fasta_output/$sim_name.fasta
