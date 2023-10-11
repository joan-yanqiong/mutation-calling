#!/bin/bash

module load gatk

gatk ApplyBQSR \
    -R ${ref_path} \
    -I ${run_id}_split.bam \
    --bqsr-recal-file ${run_id}_recal_data.table \
    -O ${run_id}_recal.bam