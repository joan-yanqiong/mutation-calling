#!/usr/bin/env bash
echo "\$(date)  Load module..."
module load annovar

echo "\$(date)   Create output directory..."
mkdir -p ${sample_id}

echo "\$(date)   Converting VCF to Annovar input..."
perl ${params.annovar}/convert2annovar.pl \
    -format vcf4 ${mutect_vcf} \
    -outfile ${sample_id}/${sample_id}.avinput \
    -allsample \
    -withfreq

echo "\$(date)  COMPLETED!"
