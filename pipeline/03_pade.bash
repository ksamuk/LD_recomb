#!/bin/bash

sim_name=$1
theta=$2
num_cores=$3 
out_file_slug=$4

ldhelmet pade \
--num_threads $num_cores \
-t $theta \
-x 12 \
--defect_threshold 40 \
-c data/ldhelmet_output/$sim_name/${out_file_slug}/${out_file_slug}.conf \
-o data/ldhelmet_output/$sim_name/${out_file_slug}/${out_file_slug}.pade
