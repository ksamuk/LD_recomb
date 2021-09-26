#!/bin/bash

sim_name=$1
out_file_slug=$2 

time ldhelmet post_to_text \
-m -p 0.025 -p 0.50 -p 0.975 \
-o data/ldhelmet_output/$sim_name/${out_file_slug}/${out_file_slug}.txt \
data/ldhelmet_output/$sim_name/${out_file_slug}/${out_file_slug}.post
