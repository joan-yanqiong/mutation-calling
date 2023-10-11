#!/bin/bash
echo "Loading STAR module..."
module load STAR

fq1="${sample_fastq_dir}/*_1.fastq"
fq2="${sample_fastq_dir}/*_2.fastq"


if [ -d ${output_prefix}/$(basename ${sample_fastq_dir}) ];
then
    echo "Output directory already present. Deleting..."
    rm -r ${output_prefix}/$(basename ${sample_fastq_dir})
    echo "Existing output directory deleted."
fi

echo "Making directory for output..."
mkdir ${output_prefix}/$(basename ${sample_fastq_dir})

echo "Running STAR on reads..."
# Parameters from Supp Table 12
STAR --genomeDir ${genome_dir} \
--runThreadN 6 \
--readFilesIn ${fq1} ${fq2} \
--outFileNamePrefix ${output_prefix}/$(basename ${sample_fastq_dir}) \
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