#!/bin/bash

echo "\$(date) Load mutect..."
module load mutect/1.1.5

echo "\$(date) Creating directory for output..."
mkdir -p ${sample_id}
cd ${sample_id}

echo "\$(date)  Run mutect..."
java -jar \$mutect_dir/muTect.jar \
    -T MuTect \
    -R "../${ref_path}" \
    --input_file:tumor "../${recal_bam}" \
    -cosmic "../${cosmic_vcf}" \
    -dbsnp "../${dbSNP_vcf}" \
    -o "${sample_id}_call_stats.txt" \
    -vcf "${sample_id}_mutect.vcf" \
    -U ALLOW_N_CIGAR_READS

echo "\$(date) COMPLETED!"