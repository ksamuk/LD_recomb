#!/bin/bash

sim_name=$1
theta=$2
rho_lik_bins=$3
num_cores=$4
out_file_slug=$5

time ldhelmet table_gen \
--num_threads $num_cores \
-t $theta -r ${rho_lik_bins} \
-c data/ldhelmet_output/$sim_name/${out_file_slug}/${out_file_slug}.conf \
-o data/ldhelmet_output/$sim_name/${out_file_slug}/${out_file_slug}.lk 
