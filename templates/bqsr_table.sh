#!/bin/bash

echo "\$(date)  Loading GATK module..."
module load gatk

echo "\$(date)  Creating directory for output..."
mkdir -p "${sample_id}"

echo "\$(date)  Generate recalibration table..."
gatk BaseRecalibrator --java-options -Xmx4g \
    -R "${ref_path}" \
    -I "${bam_file}" \
    --known-sites "${dbSNP_vcf}" \
    -O "${sample_id}/${sample_id}_recal_data.table"

echo "\$(date)  COMPLETED!"