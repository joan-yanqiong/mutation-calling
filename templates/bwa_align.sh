#!/bin/bash
echo "\$(date) Load bwa..."
module load bwa

echo "\$(date) Rename folder with FASTQ files..."
mv $dir ${dir}_fastq

echo "\$(date) Setup input variables..."
fq1="${dir}_fastq/*_1.fastq"
fq2="${dir}_fastq/*_2.fastq"

if [ -d ${sample_id} ];
then
    echo "\$(date)  Output directory already present. Deleting..."
    rm -r \${sample_id}
    echo "\$(date)  Existing output directory deleted."
fi

echo "\$(date)  Create output directory..."
mkdir -p ${sample_id}

echo "\$(date)  Running bwa-mem on reads..."
bwa mem \
    -t ${task.cpus} \
    -M \
    ${index_dir}/${index_dir} \
    \$fq1 \
    \$fq2 \
    > ${sample_id}/${sample_id}_mapped.sam

echo "\$(date)  COMPLETED!"