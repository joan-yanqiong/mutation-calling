#!/bin/bash
echo "Loading STAR module..."
module load STAR

# Change name to prevent overwriting w/ nextflow
mv $fastq_dir ${fastq_dir}_fastq

echo Sample/Run ID: ${sample_id}

fq1="${fastq_dir}_fastq/*_1.fastq"
fq2="${fastq_dir}_fastq/*_2.fastq"

echo \${fq1}
echo \${fq2}

if [ -d ${sample_id} ];
then
    echo "Output directory already present. Deleting..."
    rm -r ${sample_id}
    echo "Existing output directory deleted."
fi

echo "Making directory for output..."
mkdir -p "${sample_id}/mapped"

echo "Running STAR on reads..."
# Parameters from Supp Table 12
STAR --genomeDir ${index_dir} \
--runThreadN 6 \
--readFilesIn \${fq1} \${fq2} \
--outFileNamePrefix "${sample_id}/mapped/${sample_id}_" \
--outSAMtype BAM SortedByCoordinate \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.1 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outFilterScoreMinOverLread 0.33 \
--outFilterMatchNminOverLread 0.33 \
--limitSjdbInsertNsj 1200000