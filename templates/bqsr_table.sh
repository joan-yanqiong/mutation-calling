#!/bin/bash

module load gatk

gatk BaseRecalibrator \
    -R ${ref_path} \
    -I ${filename}_split.bam \
    --known-sites ${dbSNP_vcf} \
    -O ${run_id}_recal_data.table