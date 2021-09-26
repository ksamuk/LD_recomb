#!/bin/bash

seed=$1

cd $1

for i in `find | grep -E "\.vcf$|\.txt$"`; do gzip "$i" ; done