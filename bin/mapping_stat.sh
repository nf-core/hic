#!/bin/bash

bam1=$1
bam2=$2
merged=$3
tag=$4
  
tot_reads=$(samtools view -c ${merged})
map_reads=$(samtools view -c -F 4 ${merged})
gmap_reads=$(samtools view -c -F 4 ${bam1})
lmap_reads=$(samtools view -c -F 4 ${bam2})

echo -e "total_${tag}\t$tot_reads" 
echo -e "mapped_${tag}\t$map_reads"
echo -e "global_${tag}\t$gmap_reads"
echo -e "local_${tag}\t$lmap_reads" 
