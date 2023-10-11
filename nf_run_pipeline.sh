#!/usr/bin/env bash
#SBATCH -J run_nf
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=long
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=21-00:00:00
#SBATCH --output=slurm_out/%x_%j.out
#SBATCH --error=slurm_out/%x_%j.out

# module load java
# module load apptainer

# Local
base_dir="/Users/joankant/Library/CloudStorage/OneDrive-UHN/Coding/"
nf_exec="/Users/joankant/Library/CloudStorage/OneDrive-UHN/nextflow"
work_dir="${base_dir}/nf_work"
nf_profile="standard"

# Create work directory if not existing
mkdir -p $work_dir

project_dir="${base_dir}/local-mutation-calling"

# User parameters
# genome_dir=/path/to/store/genome_indices
# fasta_path=/path/to/FASTA_file
# gtf_path=/path/to/GTF_file


# RNA-ALIGNMENT
genome_dir="hg19"
fasta_path="/cluster/projects/gaitigroup/Users/Jahin/Reference/hg19-v0-Homo_sapiens_assembly19.fasta"
gtf_path="/cluster/projects/gaitigroup/Users/Jahin/Reference/Homo_sapiens_assembly19.gtf"
output_prefix="STAR"

# FASTQ files directory:
# /cluster/projects/gaitigroup/Users/Jahin/RNA_Mutect_Riaz/Riaz_FASTQfile
# Example: SRR5088817
fastq_dir="/cluster/projects/gaitigroup/Users/Jahin/RNA_Mutect_Riaz/Riaz_FASTQfile/"
sample_fastq_dir="${fastq_dir}/SRR5088817"
# file_{1,2}.fq
# input_bams="/cluster/projects/gaitigroup/Users/Jahin/Results/STAR"


# Start the pipeline
${nf_exec} run ${project_dir} \
    -resume -profile ${nf_profile} \
    -w $work_dir \
    --genome_dir ${genome_dir} \
    --fasta_path ${fasta_path} \
    --gtf_path ${gtf_path} \
    --ouput_prefix ${output_prefix} \
    --fastq_dir ${fastq_dir} \
