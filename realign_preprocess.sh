#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-10:00
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH --job-name realign_preprocess	
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

tumor_id=$1
normal_id=$2

prefix="/cluster/projects/gaitigroup/Users/Jahin"
script_path="${prefix}/rna_mutect"
ref_path="${prefix}/Reference/hg19-v0-Homo_sapiens_assembly19.fasta"
tumor_bam="${prefix}/Results/Practice/STAR/${tumor_id}/Aligned.sortedByCoord.out_read_groups_recal.bam"
normal_bam="${prefix}/Results/bwa/${normal_id}/${normal_id}_recal.bam"
maf="${prefix}/Results/Practice/STAR/${tumor_id}/${tumor_id}_onco.maf"

module load samtools/1.14
module load bamtools/2.4.2
module load plinkseq/0.10

cd "${prefix}/Results/Practice/STAR/${tumor_id}"

sbatch -p himem -t 7-00:00 -c 6 --mem=40G --job-name tumor_realign_preprocess \
-o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out \
${script_path}/HiSat_realign_preprocess2.sh \
 ${script_path}/jar_files \
 $maf \
 $tumor_bam \
 $ref_path \
 $tumor_id

cd "${prefix}/Results/bwa/${normal_id}"

sbatch -p himem -t 7-00:00 -c 6 --mem=40G --job-name normal_realign_preprocess \
-o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out \
${script_path}/HiSat_realign_preprocess2.sh \
 ${script_path}/jar_files \
 $maf \
 $normal_bam \
 $ref_path \
 ${normal_id}_${tumor_id}


