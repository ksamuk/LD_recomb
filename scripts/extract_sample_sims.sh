#!/bin/bash

#SBATCH -p noor
#SBATCH -J move_files
#SBATCH --mail-user=ksamuk@gmail.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o tmp/move_files-%j.out

mkdir empty_dir
rsync -a --delete empty_dir/ tmp/slim_output_sub
rm -r empty_dir

mkdir -p tmp/slim_output_sub

rsync -av --include='*_19999.vcf.gz' --include='*/' --exclude='*' data/slim_output tmp/slim_output_sub
rsync -av --include='*_20000.vcf.gz' --include='*/' --exclude='*' data/slim_output tmp/slim_output_sub
rsync -av --include='*_20100.vcf.gz' --include='*/' --exclude='*' data/slim_output tmp/slim_output_sub
rsync -av --include='*_21000.vcf.gz' --include='*/' --exclude='*' data/slim_output tmp/slim_output_sub
rsync -av --include='*_50000.vcf.gz' --include='*/' --exclude='*' data/slim_output tmp/slim_output_sub

tar czf slim_output_sub.tar.gz tmp/slim_output_sub



