#!/bin/bash
echo "Loading GATK module..."
module load gatk/3.8

echo "Creating directory for output..."
mkdir -p ${sample_id}

java -jar \$gatk_dir/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R ${ref_path} \
    -I ${marked_dup_bam} \
    -known ${indel_db1} \
    -known ${indel_db2} \
    -targetIntervals ${realigner_intervals} \
    -o ${sample_id}/${sample_id}_realigned.bam

echo "COMPLETED!"