#! /usr/bin/env nextflow

/*
Workflow for Mutation calling

*/

nextflow.enable.dsl=2

include {MAP_READS; BUILD_INDEX; TEST} from "./nf-modules/rna-alignment.nf"


// workflow PRE_PROCESSING {
//     main:
//     ADD_READ_GROUPS(params.picard_dir, )
// }

workflow {
    //  FASTQ FILES FOR SAMPLES
    fastq_dir=channel.fromPath(params.fastq_dir, type: 'dir', maxDepth = 0)

    // Build genome indices
    BUILD_INDEX(params.genome_dir, params.fasta_path, params.gtf_path)

    // Alignment
    MAP_READS(BUILD_INDEX.out.genome_indices_gen, params.output_prefix, fastq_dir)

    // // Pre-processing
    // ADD_READ_GROUPS(
    //     params.picard_dir,

    //     MAP_READS.out.aligned_bam.flatten())











}
