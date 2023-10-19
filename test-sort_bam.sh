#!/bin/bash

#!/usr/bin/env bash
#SBATCH -J test-sort_bam
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH -c 6
#SBATCH --mem=20G
#SBATCH --time=06:00:00
##SBATCH --time=00:10:00
#SBATCH --output=slurm_out/%x_%j.out
#SBATCH --error=slurm_out/%x_%j.out

base_dir="/cluster/projects/gaitigroup/Users/Joan/"
work_dir="${base_dir}/nf_work"

project_dir="${base_dir}/h4h-mutation-calling"

read_groups_bam="/cluster/projects/gaitigroup/Users/Jahin/Results/bwa/SRR5134767/SRR5134767_read_groups.bam"
sample_id="SRR5134767"
output_dir="${project_dir}/tmp/output"

cd "${output_dir}"

echo "Loading picard module..."
module load picard

echo "Creating directory for output..."
mkdir -p ${sample_id}

echo "Sorting ${sample_id}..."
java -jar $picard_dir/picard.jar SortSam \
      I=${read_groups_bam} \
      O=${sample_id}/${sample_id}_sorted.bam \
      SORT_ORDER=queryname

echo "COMPLETED!"