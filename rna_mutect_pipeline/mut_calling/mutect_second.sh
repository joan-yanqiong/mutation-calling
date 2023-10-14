#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-10:00
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH --job-name mutect
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

echo "\$(date)\tLoad mutect..."
module load mutect/1.1.5

# interval_list=snp_mutations.intervals
echo "\$(date)\tCreate output directory..."
mkdir -p ${sample_id}

echo "\$(date)\tEnter output directory..."
cd ${sample_id}

echo "\$(date)\tRun Mutect on tumor-normal pair..."
java -jar $mutect_dir/muTect.jar \
    -T MuTect \
    -R ${ref_path} \
    --input_file:tumor ${tumor_hisat_bam} \
    --input_file:normal ${normal_tumor_hisat_bam} \
    -cosmic ${cosmic_vcf} \
    -dbsnp ${dbSNP_vcf} \
    -o ${tumor_id}_second_call_stats.txt \
    -vcf ${tumor_id}_second_mutect.vcf \
    -U ALLOW_N_CIGAR_READS \
    --force_output \
    -L ${interval_list}

echo "\$(date)\tCOMPLETED!"
