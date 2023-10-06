#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-10:00
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH --job-name oncotator		
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

module load oncotator/1.9.9.0

filename=$1

oncotator ${filename}_mutect.vcf \
 ${filename}_onco.maf \
 hg19 \
 --db-dir $DB_DIR \
 -i VCF \
 -o TCGAMAF