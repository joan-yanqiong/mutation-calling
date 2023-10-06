#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-10:00
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH --job-name mutect		
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

module load mutect/1.1.5

prefix="/cluster/projects/gaitigroup/Users/Jahin"


ref_path="${prefix}/Reference/hg19-v0-Homo_sapiens_assembly19.fasta"
ref_gtf="${prefix}/Reference/Homo_sapiens_assembly19.gtf"
cosmic_vcf="${prefix}/Reference/b37_cosmic_v54_120711.vcf"
dbSNP_vcf="${prefix}/Reference/hg19-v0-Homo_sapiens_assembly19.dbsnp138.vcf"

tumor_id=$1
normal_id=$2
tumor_path=${tumor_id}.aligned.sorted_by_coord.hisat2.bam
normal_path="${prefix}/Results/bwa/${normal_id}/${normal_id}_${tumor_id}.aligned.sorted_by_coord.hisat2.bam"
interval_list=snp_mutations.intervals

java -jar $mutect_dir/muTect.jar \
    -T MuTect \
    -R ${ref_path} \
    --input_file:tumor ${tumor_path} \
    --input_file:normal ${normal_path} \
    -cosmic ${cosmic_vcf} \
    -dbsnp ${dbSNP_vcf} \
    -o ${tumor_id}_second_call_stats.txt \
    -vcf ${tumor_id}_second_mutect.vcf \
    -U ALLOW_N_CIGAR_READS \
    --force_output \
    -L ${interval_list}
