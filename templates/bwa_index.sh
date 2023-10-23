#!/bin/bash
echo "\$(date)  Loading bwa module..."
module load bwa

echo "\$(date)  Creating directory for output database..."
mkdir -p ${index_dir}_${ref_path.simpleName}

echo "\$(date)  Enter created folder..."
cd ${index_dir}_${ref_path.simpleName}

echo "\$(date)  Build genome index..."
bwa index -p ${index_dir}_${ref_path.simpleName} "../${ref_path}"

echo "\$(date)\tCOMPLETED!"