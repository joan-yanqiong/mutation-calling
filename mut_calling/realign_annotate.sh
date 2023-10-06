#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-10:00
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH --job-name realign_annotate
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

tumor_id=$1
normal_id=$2

prefix="/cluster/projects/gaitigroup/Users/Jahin"
script_path="${prefix}/rna_mutect"

cd "${prefix}/Results/Practice/STAR/${tumor_id}"

echo "writing out annotation file for upload"
echo -e "$PWD/${tumor_id}_tmp_sequence_1.fastq\n$PWD/${tumor_id}_tmp_sequence_2.fastq"   >> ${tumor_id}.rna_reads_fastq_list.list
echo "done"

cd "${prefix}/Results/bwa/${normal_id}"

echo "writing out annotation file for upload"
echo -e "$PWD/${normal_id}_vs_${tumor_id}_tmp_sequence_1.fastq\n$PWD/${normal_id}_vs_${tumor_id}_tmp_sequence_2.fastq"   >> ${normal_id}_vs_${tumor_id}.rna_reads_fastq_list.list
echo "done"