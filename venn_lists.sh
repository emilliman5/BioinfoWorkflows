#!/bin/bash

##This script will create feature lists based genomic interval overlaps.

files=("$@")
for f in files
do
    `bedtools -a $f -b $g -u > $out`
done

for f in files
do
    `bedtools -a $f -b $g -v > $out`
done

wc -l *.bed > $summary
