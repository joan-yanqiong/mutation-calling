#!/bin/bash

echo "Loading GATK module..."
module load gatk

echo "Creating directory for output..."
mkdir -p "${sample_id}"

echo "Generate recalibration table..."
gatk BaseRecalibrator \
    -R "${ref_path}" \
    -I "${bam_file}" \
    --known-sites "${dbSNP_vcf}" \
    -O "${sample_id}/${sample_id}_recal_data.table"

echo "COMPLETED!"