#!/bin/bash

# run like
# bash run_all.bash {sim_name} {window_size} {block_penalty} {num_cores}

sim_name=$1
window_size=$2
block_penalty=$3
num_cores=$4
out_file_slug="${sim_name}_wind_${window_size}_block_${block_penalty}"
theta=0.00072
#rho_lik_bins='0.0 0.1 10.0 1.0 100.0'
rho_lik_bins='0.0 0.0001 0.001 0.0005 0.0160 0.1 2.0 1.0 10'
prior=0.0003992727

mkdir data/ldhelmet_output/$sim_name
mkdir data/ldhelmet_output/$sim_name/${out_file_slug}

bash pipeline/01_find_confs.bash $sim_name $window_size $num_cores ${out_file_slug}  &&
bash pipeline/02_table_gen.bash $sim_name $theta "$rho_lik_bins" $num_cores ${out_file_slug}&&
bash pipeline/03_pade.bash $sim_name $theta $num_cores $out_file_slug &&
bash pipeline/04_rjmcmc.bash $sim_name $window_size $block_penalty $num_cores $out_file_slug $prior &&
bash pipeline/05_post_to_text.bash $sim_name $out_file_slug &&
bash pipeline/06_max_lk.bash $sim_name $num_cores $out_file_slug




