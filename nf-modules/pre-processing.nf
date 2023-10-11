process ADD_READ_GROUPS {
    /* ADD READ GROUPS
    1. picard_dir
    2.
    3. i
    */
    input:
    path picard_dir
    val i
    val run_id

    script:
    template add_read_groups.sh
}

process INDEX_DNA {
    input:
    path: prefix

    script:
    template dna_index.sh

}

process MARK_DUPLICATES {
    input:
    path picard_dir
    val run_id
}

process SPLIT_CIGARS {
    input:
    path ref_path
    val run_id

    script:
    template split_cigars
}


process BQSR_TABLE {
    // Make table
    input:
    path ref_path
    path file_name
    val dbSNP_vcf
    val run_id

    script:
    template bqsr_table.sh
}

process APPLY_BQSR {
    // Recalibration
    input:
    path ref_path
    val run_id

    script:
    template apply_bqsr.sh
}
