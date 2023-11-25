#!/bin/bash
# Ref: https://gatk.broadinstitute.org/hc/en-us/articles/360037055712-ApplyBQSR

echo "\$(date)  Loading GATK module..."
module load gatk

echo "\$(date)  Creating directory for output..."
mkdir -p ${sample_id}
cd ${sample_id}

echo "\$(date)  Apply base quality score recalibration"
gatk ApplyBQSR \
    -R "../${ref_path}" \
    -I "../${bam_file}" \
    --bqsr-recal-file "../${recal_data_table}" \
    -O "${sample_id}_recal.bam"

touch "ok.txt"

echo "\$(date)  COMPLETED!"