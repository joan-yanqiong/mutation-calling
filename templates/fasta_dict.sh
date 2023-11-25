#!/bin/bash

echo "\$(date)  Load gatk..."
module load gatk

echo "\$(date)\tCreating the FASTA sequence dictionary file..."
gatk CreateSequenceDictionary -R ${ref_genome}

touch "ok.txt"

echo "\$(date)  COMPLETED"