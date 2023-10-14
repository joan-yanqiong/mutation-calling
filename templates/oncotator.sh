#!/bin/bash

echo "Load oncotator..."
module load oncotator/1.9.9.0

mkdir -p ${sample_id}
cd ${sample_id}

echo "Run oncotator..."
oncotator "../${mutect_vcf}" \
 ${sample_id}_onco.maf \
 hg19 \
 --db-dir "../\$DB_DIR" \
 -i VCF \
 -o TCGAMAF

echo "COMPLETED!"