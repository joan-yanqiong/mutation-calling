#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-10:00
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH --job-name hisat_align_normal	
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

module load hisat2/2.0.4

prefix="/cluster/projects/gaitigroup/Users/Jahin"

index_basename="${prefix}/Reference/b37_hisat_index/hg19-v0-Homo_sapiens_assembly19_index"

filename=$1

# Read the FASTQ filenames from the FASTQ list:
n=1
while IFS= read -r "file_$n"; do
  n=$((n + 1))
done < ${filename}.rna_reads_fastq_list.list

# Parameters from Supp Table 12
hisat2 -k 5 \
 --min-intronlen 20 \
 --max-intronlen 500000 \
 -X 800 \
 -I 0 \
 -x ${index_basename} \
 -1 $file_1 \
 -2 $file_2 \
 -S $PWD/${filename}.aligned.sorted_by_coord.hisat2.sam \
 --rg-id=1 \
 --rg SM:${filename}\tPU:ILLUMINA_${filename}_lib1 \
 --rg LB:lib_1 \
 --rg PL:ILLUMINA