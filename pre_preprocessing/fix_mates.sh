#!/bin/bash

#SBATCH -p himem 
#SBATCH -t 0-10:00 
#SBATCH -c 6
#SBATCH --mem=40G 
#SBATCH --job-name fix_mates		
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

prefix="/cluster/projects/gaitigroup/Users/Jahin/Results"
module load picard

java -jar $picard_dir/picard.jar FixMateInformation \
       I="${prefix}/Practice/STAR/SRR5088818/Aligned.sortedByCoord.out_read_groups_recal.bam" \
       O="${prefix}/Practice/STAR/SRR5088818/SRR5088818_fixed_mate.bam"