#!/bin/bash


# More info:
# https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html

echo "\$(date)  Load STAR module..."
module load STAR/2.7.9a

echo "\$(date)  Create a directory for storing the genome indices..."
mkdir -p ${index_dir}_${ref_path.simpleName}

echo "\$(date)  Build genome index..."
# --runThreadN: number of threads
# --runMode: genomeGenerate mode
# --genomeDir: /path/to/store/genome_indices
# --genomeFastaFiles: /path/to/FASTA_file
# --sjdbGTFfile: /path/to/GTF_file
# --sjdbOverhang: readlength -1

STAR --runThreadN ${task.cpus} \
    --runMode genomeGenerate \
    --genomeDir ${index_dir}_${ref_path.simpleName} \
    --genomeFastaFiles "${ref_path}" \
    --sjdbGTFfile "${gtf_path}" \
    --sjdbOverhang 100

echo "\$(date)  COMPLETED!"