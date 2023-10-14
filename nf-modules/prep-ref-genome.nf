process CREATE_FASTA_DICT {
    label "very_short_process"
    /*
    Generating the dictionary and index files

    Input:
    ref_genome: Path to the FASTA file reference genome.

    Ref:
    https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format

    https://gatk.broadinstitute.org/hc/en-us/articles/360037422891-CreateSequenceDictionary-Picard-

    */
    publishDir "${projectDir}/data/reference_genome", mode: "copy"
    input:
    path ref_genome

    output:
    tuple path(ref_genome), path("${ref_genome.simpleName}.dict")

    script:
    """
    #!/bin/bash

    echo "\$(date)\tLoad gatk..."
    module load gatk

    echo "\$(date)\tCreating the FASTA sequence dictionary file..."
    gatk CreateSequenceDictionary -R ${ref_genome}
    echo "\$(date)\tCOMPLETED"
    """
}

process CREATE_FASTA_INDEX {
    label "very_short_process"
    /*
    Generating the dictionary and index files

    Input:
    ref_genome: Path to the FASTA file reference genome.

    Ref:
    https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format
    */
    publishDir "${projectDir}/data/reference_genome", mode: "copy"

    input:
    path ref_genome

    output:
    tuple path(ref_genome), path("${ref_genome.simpleName}.fasta.fai")

    script:
    """
    #!/bin/bash

    echo "\$(date)\tLoad samtools..."
    module load samtools

    echo "\$(date)\tCreating the FASTA index file..."
    samtools faidx ${ref_genome}

    echo "\$(date)\tCOMPLETED"

    """

}