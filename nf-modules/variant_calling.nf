process MUTECT {
    /* ADD READ GROUPS
    mutect_dir
    1. ref_path
    2. comic_vcf
    3. dbSNP_vcf
    4. run_id
    */
    input:
    path ref_path
    path comic_vcf
    path dbSNP_vcf
    val run_id
    path mutect_dir

    script:
    template mutect.sh
}