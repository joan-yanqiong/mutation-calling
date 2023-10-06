#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-10:00
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH --job-name hisat_index	
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

module load hisat2/2.0.4

prefix="/cluster/projects/gaitigroup/Users/Jahin"

mkdir ${prefix}/Reference/b37_hisat_index

ref_path="${prefix}/Reference/hg19-v0-Homo_sapiens_assembly19.fasta"
index_basename="${prefix}/Reference/b37_hisat_index/hg19-v0-Homo_sapiens_assembly19_index"

hisat2-build ${ref_path} ${index_basename}