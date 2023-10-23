#!/bin/bash

echo "\$(date)  Load oncotator..."
module load oncotator/1.9.9.0

echo "\$(date)  Creating directory for output..."
mkdir -p ${sample_id}

echo "\$(date)  Enter created folder..."
cd ${sample_id}

echo "\$(date)  Run oncotator..."
oncotator "../${mutect_vcf}" \
    ${sample_id}_onco.maf \
    hg19 \
    --db-dir "../\$DB_DIR" \
    -i VCF \
    -o TCGAMAF

echo "\$(date)  COMPLETED!"