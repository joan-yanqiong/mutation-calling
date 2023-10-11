#!/bin/bash

module load mutect/1.1.5

tumor_path="*_recal.bam"

java -jar $mutect_dir/muTect.jar \
    -T MuTect \
    -R ${ref_path} \
    --input_file:tumor ${tumor_path} \
    -cosmic ${cosmic_vcf} \
    -dbsnp ${dbSNP_vcf} \
    -o ${run_id}_call_stats.txt \
    -vcf ${run_id}_mutect.vcf \
    -U ALLOW_N_CIGAR_READS