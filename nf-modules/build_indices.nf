
process BWA_INDEX {
    /*
    Summary: Index database sequences in the FASTA format using bwa, for DNA
    normal samples.

    Input:
    output_db_path: Prefix of the output database [same as db filename]
    ref_path: Path to the FASTA file reference genome.

    Output:
    output_db_path: Path to the output database
    ->
    path("${output_db_path}/b37_bwa_index.amb")
    path("${output_db_path}/b37_bwa_index.ann")
    path("${output_db_path}/b37_bwa_index.bwt")
    path("${output_db_path}/b37_bwa_index.pac")
    path("${output_db_path}/b37_bwa_index.sa"

    Ref: https://bio-bwa.sourceforge.net/bwa.shtml
    */

    publishDir "${projectDir}/data/indices", mode: "copy"

    input:
    val index_dir
    tuple path(ref_path), path(ref_path_dict), path(ref_path_fai)

    output:
    path "${index_dir}_${ref_path.simpleName}", emit: indexed_dir

    script:
    template "bwa_index.sh"

}

process STAR_INDEX {
    /*
    Summary: Generating genome indexes for RNA tumor samples

    Input:
    genome_dir - path to the directory (henceforth called ”genome directory” where the genome indices are stored
    fasta_path - FASTA files with the genome reference sequences, e.g. hg19-v0- Homo_sapiens_assembly19.fasta
    gtf_path - specifies the path to the file with annotated transcripts in the standard GTF format

    Output:
    {genome_dir}/... - various files

    Ref:
    https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf
    [p. 5]
    */
    publishDir "${projectDir}/data/indices", mode: "copy"

    input:
    val index_dir
    tuple path(ref_path), path(ref_path_dict), path(ref_path_fai)
    path gtf_path

    output:
    path "${index_dir}_${ref_path.simpleName}", emit: indexed_dir

    script:
    template "star_index.sh"
}


process HISAT_INDEX {
    publishDir "${projectDir}/data/indices", mode: "copy"

    /*
    Summary: Build HISAT index

    Input:
    index_dir: name of directory to store index files
    ref_path: path to reference genome (FASTA)

    Ref: http://daehwankimlab.github.io/hisat2/howto/
    */

    input:
    val index_dir
    tuple path(ref_path), path(ref_path_dict), path(ref_path_fai)

    output:
    path "${index_dir}_${ref_path.simpleName}", emit: indexed_dir

    script:
    template "hisat_index.sh"
}
