#!/bin/bash
echo "\$(date)  Loading hisat2 module..."
module load hisat2/2.0.4

echo "\$(date)  Making directory for output..."
mkdir -p ${index_dir}_${ref_path.simpleName}

echo "\$(date)  Enter created folder..."
cd ${index_dir}_${ref_path.simpleName}

echo "\$(date)  Build Hisat2 index..."
hisat2-build -p ${task.cpus} "../${ref_path}" ${index_dir}_${ref_path.simpleName}

echo "\$(date)  COMPLETED!"