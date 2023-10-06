#!/bin/bash

#SBATCH -p veryhimem
#SBATCH -t 5-00:00
#SBATCH -c 6
#SBATCH --mem=180G
#SBATCH --job-name filter1
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

module load MCR/v901

tumor_id=$1

rna_mutect_dir="/cluster/projects/gaitigroup/Users/Jahin/rna_mutect"
prefix="/cluster/projects/gaitigroup/Users/Jahin"

mcr_path="/cluster/tools/software/centos7/MCR/v901"
pair_id="SRR5088817_tumor_vs_normal"
maf_path="${tumor_id}_onco_mod.maf"
call_stats_path="${tumor_id}_second_call_stats.txt"
min_alt_count=5
# pon_thr=-3
darned_mat="${rna_mutect_dir}/mat_files/Darned_mat.mat"
radar_mat="${rna_mutect_dir}/mat_files/Radar_mat.mat"
exac_mat="${rna_mutect_dir}/mat_files/Exac_mat.mat"
# pon=
# cytoband=

sh ${rna_mutect_dir}/run_FilterRNAMutationsNoPoN.sh ${mcr_path} \
 ${pair_id} \
 ${maf_path} \
 ${call_stats_paths} \
 ${min_alt_count} \
 ${darned_mat} \
 ${radar_mat} \
 ${exac_mat} 