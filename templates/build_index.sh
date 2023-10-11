#!/bin/bash


# More info:
# https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html

# TODO test if it works without entering the output folder
# cd ${genome_dir}/.. # Not sure if this is necessary

echo "Load STAR module..."
module load STAR/2.7.9a

echo "Create a directory for storing the genome indices..."
mkdir -p ${genome_dir}

echo "Build genome index..."
# --runThreadN: number of threads
# --runMode: genomeGenerate mode
# --genomeDir: /path/to/store/genome_indices
# --genomeFastaFiles: /path/to/FASTA_file
# --sjdbGTFfile: /path/to/GTF_file
# --sjdbOverhang: readlength -1

STAR --runThreadN 6 \
 --runMode genomeGenerate \
 --genomeDir ${genome_dir} \
 --genomeFastaFiles ${fasta_path} \
 --sjdbGTFfile ${gtf_path} \
 --sjdbOverhang 100