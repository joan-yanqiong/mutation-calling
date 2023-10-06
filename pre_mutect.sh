#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-10:00
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH --job-name pre-mutect-processing
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

# 1. Edit each of the files to take in relevant parameters, then pass in arguments
# to the parameters here.
#
# 2. Create dependencies so that a job doesn't begin until all the jobs it's dependent on
# are completed
#
# 3. Delete intermediate files made by each step once they're no longer needed
#
# 4. Write slurm header

prefix="/cluster/projects/gaitigroup/Users/Jahin"

ref_path="${prefix}/Reference/hg19-v0-Homo_sapiens_assembly19.fasta"
ref_gtf="${prefix}/Reference/Homo_sapiens_assembly19.gtf"
cosmic_vcf="${prefix}/Reference/b37_cosmic_v54_120711.vcf"
dbSNP_vcf="${prefix}/Reference/hg19-v0-Homo_sapiens_assembly19.dbsnp138.vcf"

tumor_fq_prefix="${prefix}/Practice_Data/FASTQ" # Directory containing fastq directories for each run
output_dir="${prefix}/Results/Practice/STAR" # Directory to store the directories for each individual run
tumor_dir="${prefix}/Results/Practice/STAR" # Directory to store the directories for each individual run


if [ -d ${output_dir}/STAR_index ];
then
    echo "Index directory already present. Deleting..."
    rm -r ${output_dir}/STAR_index
    echo "Existing index directory deleted."
fi

mkdir ${output_dir}/STAR_index

sbatch build_index.sh ${output_dir}/STAR_index \
 $ref_path \
 $ref_gtf

i=0

cd ${tumor_fq_prefix}
for run_id in SRR*
do
    cd ${run_id}

    i=`expr $i + 1`

    jid1=$(sbatch ${script_path}/map_reads.sh ${output_dir}/STAR_index ${output_dir} ${run_id} | sed 's/Submitted batch job //')
    
    jid2=$(sbatch -d=afterok:$jid1 ${script_path}/add_read_groups.sh $i ${run_id} | sed 's/Submitted batch job //')

    jid3=$(sbatch -d=afterok:$jid2 ${script_path}/mark_duplicates.sh ${run_id} | sed 's/Submitted batch job //')

    jid4=$(sbatch -d=afterok:$jid3 ${script_path}/split_cigars.sh ${ref_path} ${run_id} | sed 's/Submitted batch job //')

    jid5=$(sbatch -d=afterok:$jid4 ${script_path}/bqsr_table.sh ${ref_path} ${dbSNP_vcf} ${run_id} | sed 's/Submitted batch job //')

    jid6=$(sbatch -d=afterok:$jid5 ${script_path}/apply_bqsr.sh ${ref_path} ${run_id} | sed 's/Submitted batch job //')

    jid7=$(sbatch -d=afterok:$jid6 ${script_path}/mutect.sh ${ref_path} ${cosmic_vcf} ${dbSNP_vcf} ${run_id} | sed 's/Submitted batch job //')

    jid8=$(sbatch -d=afterok:$jid7 ${script_path}/oncotator.sh ${run_id} | sed 's/Submitted batch job //')

    cd ..
done