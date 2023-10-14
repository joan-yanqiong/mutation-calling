include { SORT_BAM; NORMAL_MARK_DUPLICATES; INDEL_REALIGN_TARGET; INDEL_REALIGNER; NORMAL_BQSR_TABLE; NORMAL_APPLY_BQSR; } from "../nf-modules/dna-preprocessing"

workflow DNA_NORMAL_PROCESSING{
      take:
         index_dir
         indel_db1_set
         indel_db2_set
         ref_genome
         dbSNP_set

      main:

//    normal_samples_meta = Channel.fromPath(params.sample_sheet) \
//    | splitCsv(header:true) \
//    | map { row -> tuple(row.normal_id, row.tumor_id, row.pair_id, file(row.normal_mapped_sam)) }

   //  TODO: testing if FASTQ files are found/downloaded again
   //  NORMAL_ALIGN(genome_dir = BUILD_INDEX.out.indexed_db, )
   // TODO: test, current files do not work.
   //  SAM_TO_BAM(normal_sample) \
   //      | SORT_BAM \
   //      | NORMAL_MARK_DUPLICATES

   // For testing purposes use read_groups.bam (Jahin) -> named now _mapped.bam
   tmp_samples_meta = Channel.fromPath(params.sample_sheet) \
   | splitCsv(header:true) \
   | map { row -> tuple(row.normal_id, row.normal_mapped_sam)}

   // NOT WORKING????
   // SORT_BAM(tmp_samples_meta)

   tmp_sort_bam_output = tuple("SRR5134767", file("/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/output/normal/SRR5134767/SRR5134767_sorted.bam"))

   NORMAL_MARK_DUPLICATES(tmp_sort_bam_output)

   INDEL_REALIGN_TARGET(NORMAL_MARK_DUPLICATES.out.output, indel_db1_set, indel_db2_set, ref_genome)

   // TODO not working errors with recognizing type of helper files
   INDEL_REALIGNER(INDEL_REALIGN_TARGET.out.output, indel_db1_set, indel_db2_set, ref_genome)

   NORMAL_BQSR_TABLE(INDEL_REALIGNER.out.output, ref_genome, dbSNP_set)

   NORMAL_APPLY_BQSR(NORMAL_BQSR_TABLE.out.output, ref_genome)

   emit:
      NORMAL_APPLY_BQSR.out.output
}