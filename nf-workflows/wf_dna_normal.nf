include { SORT_BAM; MARK_DUPLICATES_NORMAL; INDEL_REALIGN_TARGET; INDEL_REALIGNER; NORMAL_BQSR_TABLE; NORMAL_APPLY_BQSR; ADD_READ_GROUPS_NORMAL; SORT_BAM_COORD; INDEX_BAM } from "../nf-modules/dna-preprocessing"
include { BWA_ALIGN; } from "../nf-modules/dna-alignment"
include {SAM_TO_BAM} from "../nf-modules/utils.nf"
workflow DNA_NORMAL_PROCESSING{
   take:
      index_dir
      indel_db1_set
      indel_db2_set
      ref_genome
      dbSNP_set

   main:

   normal_samples_meta = Channel.fromPath(params.sample_sheet) \
   | splitCsv(header:true) \
   | map { row -> tuple(row.idx, row.normal_id, file(row.normal_fastq)) }

   BWA_ALIGN(normal_samples_meta, index_dir)

   ADD_READ_GROUPS_NORMAL(BWA_ALIGN.out.mapped_sam, sample_type = "normal", suffix = "read_groups")
   SORT_BAM(ADD_READ_GROUPS_NORMAL.out.output, sort_order = "queryname", suffix = "sorted") \
   | MARK_DUPLICATES_NORMAL

   SORT_BAM_COORD(MARK_DUPLICATES_NORMAL.out.output, sort_order = "coordinate", suffix = "marked_dup_sorted_coord") \
   | INDEX_BAM

   INDEL_REALIGN_TARGET(INDEX_BAM.out.output, indel_db1_set, indel_db2_set, ref_genome)

   INDEL_REALIGNER(INDEL_REALIGN_TARGET.out.output, indel_db1_set, indel_db2_set, ref_genome)

   NORMAL_BQSR_TABLE(INDEL_REALIGNER.out.output, ref_genome, dbSNP_set)

   NORMAL_APPLY_BQSR(NORMAL_BQSR_TABLE.out.output, ref_genome)

   emit:
      NORMAL_APPLY_BQSR.out.output
}