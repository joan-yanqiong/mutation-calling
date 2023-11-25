#!/bin/bash

echo "\$(date)\tLoad samtools..."
module load samtools

echo "\$(date)\tCreating the FASTA index file..."
samtools faidx ${ref_genome}

touch "ok.txt"

echo "\$(date)\tCOMPLETED"
