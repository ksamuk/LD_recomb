#!/bin/bash

# downsample vcfs to match lookup tables
vcf=$1
vcfout=$(echo $vcf | sed 's/\.vcf/_n25.vcf/g')
gunzip -c $vcf | bgzip > ${vcf}_tmp.gz

tabix ${vcf}_tmp.gz

bcftools view -s i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24  --output-type z ${vcf}_tmp.gz > $vcfout
rm ${vcf}_tmp.gz
rm ${vcf}_tmp.gz.tbi