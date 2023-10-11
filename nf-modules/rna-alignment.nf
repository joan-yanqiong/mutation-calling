

// process INFER_CELLCHAT {
//     label 'conda_env'

//     publishDir "${projectDir}/output/${params.output_run_name}/200_cci_cellchat", mode: "copy"

//     input:
//     path input_file
//     path interactions_db

//     output:
//     path "cellchat__${input_file.baseName}.rds", emit: cellchat_obj

//     script:
//     """
//     #!/usr/bin/env bash

//     Rscript "/${params.code_dir}/scripts/200_cci_cellchat.R" \
//     --gene_expr \$PWD/${input_file} \
//     --output_dir "\$PWD" \
//     --resource \$PWD/${interactions_db} \
//     --ident_col ${params.ident_col} \

//     --n_perm ${params.n_perm}

//     """
// }

process TEST {
    input:
    path genome_dir
    path fasta_path

    output:
    path "test_output.txt"

    //  TODO: what is the output?
    script:
    template "test.sh"

}
process BUILD_INDEX {
    /* builds genome index for STAR (hg19 genome)
    inputs:
    1. genome_dir: Directory for genome indices to be stored. Make your own directory if one doesn't exist already
    2. fasta_path: Path to genome fasta file. Used (anything from UCSC will do): Reference/hg19-v0-Homo_sapiens_assembly19.fasta
    3. gtf_path: Path to genome GTF file.
    Used (get from GENCODE): Reference/Homo_sapiens_assembly19.gtf
    */

    input:
    path genome_dir
    path fasta_path
    path gtf_path

    output:
    path "${genome_dir}/*", emit: genome_indices_gen

    script:
    template build_index.sh
}

process MAP_READS {
    /*
    1. genome_dir: directory containing genome indices
    2. output_prefix: prefix for output files (output_prefix is output
       directory, will create a directory {output_prefix}/{sample_id})
    3. sample_id: directory where fastq files are stored, containing
       {sample_id}_1.fastq and {sample_id}_2.fastq files.
    */

    input:
    path genome_dir
    path output_prefix
    path sample_fastq_dir

    output:
    path "${output_prefix}/${sample_fastq_dir.baseName}", emit: aligned_bam

    script:
    template map_reads.sh
}
