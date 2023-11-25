#!/bin/sh

echo "\$(date)   Loading modules..."
module load gatk

echo "\$(date)  Creating output directory..."
mkdir -p ${sample_id}

echo "\$(date)  Use Funcotator for annotation..."
${params.funcotator} Funcotator \
    --variant ${input_vcf} \
    --reference ${ref_path} \
    --ref-version hg19 \
    --data-sources-path ${data_sources_path} \
    --output ${sample_id}/${sample_id}_${suffix}.funcotated.maf \
    --output-file-format MAF

touch "ok.txt"

echo "\$(date)   COMPLETED!"