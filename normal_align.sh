#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-10:00
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH --job-name normal_align		
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

module load bwa

genome_dir=$1
sample_id=$2
output_prefix=$3

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
 
echo "Running bwa-mem on reads..."

bwa mem \
 -M \
 -R "@RG\tID:1\tLB:lib_1\tPL:ILLUMINA\tSM:${sample_id}\tPU:ILLUMINA_${sample_id}_lib1" \
 $genome_dir \
 $fq1 \
 $fq2 \
 > ${output_prefix}/${sample_id}/${sample_id}_mapped.sam