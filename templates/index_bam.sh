#!/bin/bash
echo "\$(date)  Loading samtools module..."
module load samtools

echo "\$(date)  Creating directory for output..."
mkdir -p ${sample_id}

echo "\$(date)  Indexing BAM file..."
samtools index -b ${bam_file} "${sample_id}/${bam_file}.bai"

touch "ok.txt"

echo "\$(date)  COMPLETED!"