#!/bin/bash

echo "\$(date)  Loading module..."
module load samtools

echo "\$(date)  Creating output directory..."
mkdir -p ${sample_id}

echo "\$(date)  Convert SAM to BAM..."
samtools view -bS ${sam_file} > ${sample_id}/${sam_file.simpleName}.bam

echo "\$(date)  COMPLETED!"