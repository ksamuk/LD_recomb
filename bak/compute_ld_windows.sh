#!/bin/bash

vcf_file_slug="data/slim_output/output_vcf"

vcf-sort ${vcf_file_slug}.vcf > ${vcf_file_slug}_sorted.vcf
vcftools --hap-r2 --ld-window-bp 10000 --vcf ${vcf_file_slug}_sorted.vcf --stdout > data/vcftools_output/ld_windows.txt

