#!/usr/bin/env bash
echo "\$(date)  Load module..."
module load annovar

echo "\$(date)   Create output directory..."
mkdir -p ${sample_id}

echo "\$(date)   Annotating variants..."
perl ${params.annovar}/table_annovar.pl \
                -vcfinput "${mutect_vcf}" \
                ${human_db} \
                -buildver hg19 \
                -remove \
                -protocol refGene \
                -operation g \
                -nastring . \
                -polish \
                --otherinfo \
                -out "${sample_id}/${sample_id}_${suffix}"
# create file as a flag to indicate successful completion
touch "ok.txt"

echo "\$(date)  COMPLETED!"
