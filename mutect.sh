#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-10:00
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH --job-name mutect		
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

module load mutect/1.1.5

prefix="/cluster/projects/gaitigroup/Users/Jahin"

ref_path=$1
cosmic_vcf=$2
dbSNP_vcf=$3
run_id=$4
tumor_path="*_recal.bam"

java -jar $mutect_dir/muTect.jar \
    -T MuTect \
    -R ${ref_path} \
    --input_file:tumor ${tumor_path} \
    -cosmic ${cosmic_vcf} \
    -dbsnp ${dbSNP_vcf} \
    -o ${run_id}_call_stats.txt \
    -vcf ${run_id}_mutect.vcf \
    -U ALLOW_N_CIGAR_READS