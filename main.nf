#! /usr/bin/env nextflow

/*
Workflow for Mutation calling

*/

nextflow.enable.dsl=2

// Import processes
include { BWA_INDEX; STAR_INDEX; HISAT_INDEX } from "./nf-modules/build_indices"

// Import workflows
include {SET_REF_GENOME} from "./nf-workflows/wf_prepare_genome_ref"
include {RNA_TUMOR_PROCESSING; MUTECT_ROUND1} from "./nf-workflows/wf_rna_tumor"
include {DNA_NORMAL_PROCESSING} from "./nf-workflows/wf_dna_normal"
include {RNA_MUTECT} from "./nf-workflows/wf_rna_mutect"

workflow {
    println """\
    ---- WORKFLOW PARAMETERS ----
    PARAMS for building genome indices:
    STAR index: ${params.star_index_dir}
    HISAT index: ${params.hisat_index_dir}
    BWA index: ${params.bwa_index_dir}

    HELPER FILES:
    sample sheet: ${params.sample_sheet}
    reference genome: ${params.ref_genome}
    dbSNP vcf: ${params.dbSNP}
    COSMIC vcf: ${params.dbcosmic}
    INDEL DB1: ${params.indel_db1}
    INDEL DB2: ${params.indel_db2}
    Ref GTF: ${params.gtf_path}
    JAR FILES: ${params.jar_files}

    """.stripIndent()
    // Setup reference files/databases
    dbSNP_set = tuple(params.dbSNP, "${params.dbSNP}.idx")
    cosmic_set = tuple(params.dbcosmic, "${params.dbcosmic}.idx")

    indel_db1_set = tuple(params.indel_db1, "${params.indel_db1}.idx")
    indel_db2_set = tuple(params.indel_db2, "${params.indel_db2}.idx")

    // ---- BUILDING INDICES ---- //
    SET_REF_GENOME()

    BWA_INDEX(
        index_dir = params.bwa_index_dir,
        ref_path = SET_REF_GENOME.out.ref_genome
    )
    STAR_INDEX(
        index_dir = params.star_index_dir,
        ref_path = SET_REF_GENOME.out.ref_genome,
        params.gtf_path
    )
    HISAT_INDEX(
        index_dir = params.hisat_index_dir,
        ref_path = SET_REF_GENOME.out.ref_genome
    )

    // ---- HANDLING TUMOR SAMPLES ---- //
    RNA_TUMOR_PROCESSING(
        star_index = STAR_INDEX.out.indexed_dir,
        ref_genome = SET_REF_GENOME.out.ref_genome,
        dbSNP_set = dbSNP_set,
        cosmic_set = cosmic_set
    )

    MUTECT_ROUND1(
        rna_tumor_processed = RNA_TUMOR_PROCESSING.out,
        ref_genome = SET_REF_GENOME.out.ref_genome,
        dbSNP_set = dbSNP_set,
        cosmic_set = cosmic_set
    )

    // ---- HANDLING NORMAL SAMPLES ---- //
    DNA_NORMAL_PROCESSING(
        index_dir = BWA_INDEX.out.indexed_dir,
        indel_db1_set = indel_db1_set,
        indel_db2_set = indel_db2_set,
        ref_genome = SET_REF_GENOME.out.ref_genome,
        dbSNP_set = dbSNP_set
    )

    // ---- RNA-MUTECT ---- //
    RNA_MUTECT(
        mutect_round1 = MUTECT_ROUND1.out,
        rna_tumor_processed = RNA_TUMOR_PROCESSING.out,
        dna_normal_processed = DNA_NORMAL_PROCESSING.out,
        ref_genome = SET_REF_GENOME.out.ref_genome,
        index_dir = HISAT_INDEX.out.indexed_dir,
        dbSNP_set = dbSNP_set,
        cosmic_set = cosmic_set
    )
}

 workflow.onComplete {

    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}