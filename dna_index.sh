#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-10:00
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH --job-name dna_index	
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

module load bwa

prefix="/cluster/projects/gaitigroup/Users/Jahin/Reference"

mkdir ${prefix}/b37_bwa_index

bwa index -p ${prefix}/b37_bwa_index \
 ${prefix}/hg19-v0-Homo_sapiens_assembly19.fasta
