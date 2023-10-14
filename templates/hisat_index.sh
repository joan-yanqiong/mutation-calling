#!/bin/bash
echo "\$(date)\tLoading hisat2 module..."
module load hisat2/2.0.4

echo "\$(date)\tMaking directory for output..."
mkdir -p ${index_dir}_${ref_path.simpleName}

# Setup paths/directories
cd ${index_dir}_${ref_path.simpleName}
echo "\$(date)\tBuild Hisat2 index..."
hisat2-build "../${ref_path}" ${index_dir}_${ref_path.simpleName}

echo "\$(date)\tCOMPLETED!"