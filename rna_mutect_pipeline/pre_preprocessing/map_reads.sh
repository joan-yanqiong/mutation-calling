#!/bin/bash

#SBATCH -p himem 
#SBATCH -t 0-10:00 
#SBATCH -c 6
#SBATCH --mem=40G 
#SBATCH --job-name mapping_reads 		
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Log%x-%j.out

genome_dir=$1 # Directory containing genome indices
output_prefix=$2 # Prefix for output files -- should include path to output directory
sample_id=$3

fq1="${sample_id}/*_1.fastq"
fq2="${sample_id}/*_2.fastq"

if [ -d ${output_prefix}/${sample_id} ];
then
    echo "Output directory already present. Deleting..."
    rm -r ${output_prefix}/${sample_id}
    echo "Existing output directory deleted."
fi

echo "Making directory for output..."
mkdir ${output_prefix}/${sample_id}

module load STAR

echo "Running STAR on reads..."
# Parameters from Supp Table 12
STAR --genomeDir ${genome_dir} \
--runThreadN 6 \
--readFilesIn ${fq1} ${fq2} \
--outFileNamePrefix ${output_prefix}/${sample_id}/ \
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