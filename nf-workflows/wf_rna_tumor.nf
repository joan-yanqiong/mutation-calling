include { STAR_ALIGN } from "../nf-modules/rna-alignment.nf"
include { ADD_READ_GROUPS; MARK_DUPLICATES; SPLIT_CIGARS; BQSR_TABLE; APPLY_BQSR} from "../nf-modules/rna-preprocessing.nf"
include { MUTECT_R1; ONCOTATOR; ANNOVAR_R1} from "../nf-modules/rna_mutect1.nf"

workflow RNA_TUMOR_PROCESSING {
    take:
    index_dir
    ref_genome
    dbSNP_set
    cosmic_set

    main:
    samples_meta = Channel.fromPath(params.sample_sheet) \
        | splitCsv(header:true) \
        | map { row -> tuple(row.ix, row.rnaseq_tumor_id, file(row.rnaseq_tumor_fastq)) }

    STAR_ALIGN(index_dir = index_dir, samples_meta)

    ADD_READ_GROUPS(STAR_ALIGN.out.output, sample_type = "tumor")

    MARK_DUPLICATES(ADD_READ_GROUPS.out.output)

    SPLIT_CIGARS(MARK_DUPLICATES.out.output, ref_genome)

    BQSR_TABLE(SPLIT_CIGARS.out.output, ref_genome, dbSNP_set)

    APPLY_BQSR(BQSR_TABLE.out.output, ref_genome)

    emit:
    APPLY_BQSR.out.output
}

workflow MUTECT_ROUND1 {
    take:
    rna_tumor_processed
    ref_genome
    dbSNP_set
    cosmic_set

    main:
    MUTECT_R1(rna_tumor_processed, ref_genome, dbSNP_set, cosmic_set, suffix = "mutect_R1")
    ANNOVAR_R1(MUTECT_R1.out.output, human_db = params.human_db, suffix = "annovar_R1")
    ONCOTATOR(MUTECT_R1.out.output)

    emit:
    ONCOTATOR.out.onco_maf

}