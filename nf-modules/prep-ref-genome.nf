process CREATE_FASTA_DICT {
    label "time_30m"
    label "mem8"
    /*
    Generating the dictionary and index files

    Input:
    ref_genome: Path to the FASTA file reference genome.

    Ref:
    https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format

    https://gatk.broadinstitute.org/hc/en-us/articles/360037422891-CreateSequenceDictionary-Picard-

    */
    publishDir "${projectDir}/data/reference_genome", mode: "symlink"

    input:
    path ref_genome

    output:
    tuple path(ref_genome), path("${ref_genome.simpleName}.dict"), emit: output
    path "ok.txt"

    script:
    template "fasta_dict.sh"
}

process CREATE_FASTA_INDEX {
    label "time_30m"
    label "mem8"
    /*
    Generating the dictionary and index files

    Input:
    ref_genome: Path to the FASTA file reference genome.

    Ref:
    https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format
    */
    publishDir "${projectDir}/data/reference_genome", mode: "symlink"

    input:
    path ref_genome

    output:
    tuple path(ref_genome), path("${ref_genome.simpleName}.fasta.fai"), emit: output
    path "ok.txt"

    script:
    template "fasta_index.sh"

}