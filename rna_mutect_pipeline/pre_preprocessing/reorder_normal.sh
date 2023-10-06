#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-6:00
#SBATCH -c 6
#SBATCH --mem=40G 
#SBATCH --job-name reorder_normal		
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

prefix="/cluster/projects/gaitigroup/Users/Jahin"

module load picard

cd ${prefix}/Results/STAR/Normal/SRR5134767

filename="SRR5134767_recal.bam"

echo $filename

java -jar $picard_dir/picard.jar ReorderSam \
          INPUT=$filename \
          OUTPUT=reordered_normal.bam \
          REFERENCE="${prefix}/Reference/hg19-v0-Homo_sapiens_assembly19.fasta"