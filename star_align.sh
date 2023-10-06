#!/bin/bash

#SBATCH -p himem 
#SBATCH -t 0-10:00 
#SBATCH -c 6
#SBATCH --mem=40G 
#SBATCH --job-name STAR_alignment 		
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Log/%x-%j.out

prefix="/cluster/projects/gaitigroup/Users/Jahin"

fq_prefix=$1 # Directory containing fastq directories for each run
output_prefix=$2 # Directory to store the directories for each individual run
genome_dir=$3
ref_path=$4
ref_gtf=$5
script_path="${prefix}/Scripts/pipeline_scripts"

sbatch ${script_path}/build_index.sh $genome_dir \
 $ref_path \
 $ref_gtf

cd ${fq_prefix}/
for sample_id in `ls`
do
    sbatch ${script_path}/map_reads.sh $genome_dir \
     $output_prefix \
     $sample_id
done

# I think this script isn't needed -- delete?