#!/bin/bash

#SBATCH -p himem 
#SBATCH -t 0-10:00 
#SBATCH -c 6
#SBATCH --mem=40G 
#SBATCH --job-name STAR_index		
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

genome_dir=$1 # Directory for genome indices to be stored
# Make your own directory if one doesn't exist already
fasta_path=$2 # Path to genome fasta file
# Used (anything from UCSC will do): Reference/hg19-v0-Homo_sapiens_assembly19.fasta
gtf_path=$3 # Path to genome GTF file
# Used (get from GENCODE): Reference/Homo_sapiens_assembly19.gtf

cd ${genome_dir}/.. # Not sure if this is necessary
module load STAR/2.7.9a

STAR --runThreadN 6 \
 --runMode genomeGenerate \
 --genomeDir ${genome_dir} \
 --genomeFastaFiles ${fasta_path} \
 --sjdbGTFfile ${gtf_path} \
 --sjdbOverhang 100