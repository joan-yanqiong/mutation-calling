#!/usr/bin/env bash
echo "\$(date)  Load module..."
module load annovar

echo "\$(date)   Create output directory..."
mkdir -p ${sample_id}

echo "\$(date)   Annotating variants..."
perl ${params.annovar}/table_annovar.pl "${av_input}" \
                ${human_db} \
                -buildver hg19 \
                -out "${sample_id}/${sample_id}.annovar" \
                -remove \
                -protocol refGene \
                -operation g \
                -nastring . \
                -polish \
                --otherinfo

echo "\$(date)  COMPLETED!"
