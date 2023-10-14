#!/bin/bash

echo "Loading module..."
module load samtools

echo "Creating output directory..."
mkdir -p ${sample_id}

echo "Convert SAM to BAM..."
samtools view -bS ${sam_file} > ${sample_id}/${sam_file.simpleName}.bam

echo "COMPLETED!"