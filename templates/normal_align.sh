#!/bin/bash

echo "Load bwa..."
module load bwa

mv $dir ${dir}_fastq

fq1="${dir}_fastq/*_1.fastq"
fq2="${dir}_fastq/*_2.fastq"

if [ -d ${sample_id} ];
then
    echo "Output directory already present. Deleting..."
    rm -r \${sample_id}
    echo "Existing output directory deleted."
fi

echo "Making directory for output..."
mkdir -p ${sample_id}

echo "Running bwa-mem on reads..."

bwa mem \
 -M \
 -R "@RG\tID:1\tLB:lib_1\tPL:ILLUMINA\tSM:${sample_id}\tPU:ILLUMINA_${sample_id}_lib1" \
 ${genome_dir} \
 \$fq1 \
 \$fq2 \
 > ${sample_id}/${sample_id}_mapped.sam

echo "COMPLETED!"