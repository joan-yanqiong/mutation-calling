#!/bin/bash
echo "\$(date)\tLoading bwa module..."
module load bwa

echo "\$(date)\tCreating directory for output database..."
mkdir -p ${index_dir}_${ref_path.simpleName}

echo "\$(date)\tEnter created folder..."
cd ${index_dir}_${ref_path.simpleName}

echo "\$(date)\tBuild genome index..."
bwa index -p ${index_dir}_${ref_path.simpleName} "../${ref_path}"

echo "\$(date)\tCOMPLETED!"