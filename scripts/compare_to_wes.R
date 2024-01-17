# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
cmd_args <- commandArgs(trailingOnly = FALSE)
has_script_filepath <- startsWith(cmd_args, "--file=")
if (sum(has_script_filepath)) {
    setwd(dirname(unlist(strsplit(cmd_args[has_script_filepath], "=")))[2])
}

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
devtools::load_all("./", export_all = FALSE)

if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Get metadata",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/test_set_rerun/analysis/ad_dp")
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(args$output_dir)

# Load additional libraries
rnaseq_mutations <- readRDS("output/test_set_rerun/all_mutect_filtered_dp_ad.rds")
wes_mutations <- readRDS("/Users/joankant/Desktop/gaitigroup/Users/Joan/wes-mutation-calling/output/test_set_rerun//all_mutect_filtered_dp_ad.rds")
wes_mutations2 <- readRDS("/Users/joankant/Desktop/gaitigroup/Users/Joan/wes-mutation-calling/output/test_set_rerun/all_mutect2_filtered_dp_ad.rds")
head(rnaseq_mutations)

sample_ids <- rnaseq_mutations %>%
    pull(wes_tumor_id) %>%
    unique()

rnaseq_mutations <- rnaseq_mutations %>% mutate(exact_mut = paste(Gene.refGene, CHROM, POS, REF, ALT, sep = "_"))
wes_mutations <- wes_mutations %>% mutate(exact_mut = paste(Gene.refGene, CHROM, POS, REF, ALT, sep = "_"))
wes_mutations2 <- wes_mutations2 %>% mutate(exact_mut = paste(Gene.refGene, CHROM, POS, REF, ALT, sep = "_"))

log_info("Creating venn diagrams...")
for (sample_id in sample_ids) {
    curr_rna_mut <- rnaseq_mutations %>%
        filter(wes_tumor_id == sample_id) %>%
        pull(Gene.refGene) %>%
        unique()
    curr_wes_mut <- wes_mutations %>%
        filter(wes_tumor_id == sample_id) %>%
        pull(Gene.refGene) %>%
        unique()
    curr_wes_mut2 <- wes_mutations2 %>%
        filter(wes_tumor_id == sample_id) %>%
        pull(Gene.refGene) %>%
        unique()

    create_venndiagram(
        x = list(RNAseq = curr_rna_mut, WES_mut = curr_wes_mut, WES_mut2 = curr_wes_mut2), main = sample_id,
        category.names = c("RNA", "WES (mutect)", "WES (mutect2)"),
        filename = glue("{args$output_dir}/{sample_id}_venn_diagram.png")
    )

    curr_rna_mut <- rnaseq_mutations %>%
        filter(wes_tumor_id == sample_id) %>%
        pull(exact_mut)
    curr_wes_mut <- wes_mutations %>%
        filter(wes_tumor_id == sample_id) %>%
        pull(exact_mut)
    curr_wes_mut2 <- wes_mutations2 %>%
        filter(wes_tumor_id == sample_id) %>%
        pull(exact_mut)

    create_venndiagram(
        x = list(RNAseq = curr_rna_mut, WES_mut = curr_wes_mut, WES_mut2 = curr_wes_mut2), main = sample_id,
        category.names = c("RNA", "WES (mutect)", "WES (mutect2)"),
        filename = glue("{args$output_dir}/{sample_id}_exact_mut_venn_diagram.png")
    )
}
log_info("Remove log files...")
sapply(list.files(args$output_dir, pattern = "*.log$", full.names = TRUE), file.remove)

# Exact match
rna_tmp <- rnaseq_mutations %>%
    select(wes_tumor_id, Gene.refGene, exact_mut) %>%
    mutate(in_RNAseq = 1)

wes_tmp <- wes_mutations %>%
    select(wes_tumor_id, Gene.refGene, exact_mut) %>%
    mutate(in_WES = 1)

wes_tmp2 <- wes_mutations2 %>%
    select(wes_tumor_id, Gene.refGene, exact_mut) %>%
    mutate(in_WES2 = 1)

combined_exact <- merge(merge(rna_tmp, wes_tmp, all = TRUE), wes_tmp2, all = TRUE)
combined_exact[is.na(combined_exact)] <- 0

# Lenient match (only gene)
rna_tmp <- rnaseq_mutations %>%
    distinct(wes_tumor_id, Gene.refGene) %>%
    mutate(in_RNAseq = 1)
wes_tmp <- wes_mutations %>%
    distinct(wes_tumor_id, Gene.refGene) %>%
    mutate(in_WES = 1)
wes_tmp2 <- wes_mutations2 %>%
    distinct(wes_tumor_id, Gene.refGene) %>%
    mutate(in_WES2 = 1)

combined_lenient <- merge(merge(rna_tmp, wes_tmp, all = TRUE), wes_tmp2, all = TRUE)
combined_lenient[is.na(combined_lenient)] <- 0

# Determine overlap (lenient)
common <- combined_lenient %>%
    mutate(common_rna_wes = (in_RNAseq + in_WES == 2), common_rna_wes2 = (in_RNAseq + in_WES2 == 2)) %>%
    group_by(wes_tumor_id) %>%
    summarise(n_common_rna_wes = sum(common_rna_wes), n_common_rna_wes2 = sum(common_rna_wes2)) %>%
    ungroup()

n_genes <- combined_lenient %>%
    mutate(common_rna_wes = (in_WES > 0), common_rna_wes2 = (in_WES2 > 0)) %>%
    group_by(wes_tumor_id) %>%
    summarise(n_common_rna_wes = sum(common_rna_wes), n_common_rna_wes2 = sum(common_rna_wes2)) %>%
    ungroup() %>%
    rename(n_genes_rna_wes = n_common_rna_wes, n_genes_rna_wes2 = n_common_rna_wes2)

lenient_overlap <- merge(common, n_genes, by = "wes_tumor_id") %>%
    mutate(overlap_rna_wes = 100 * n_common_rna_wes / n_genes_rna_wes, overlap_rna_wes2 = round(100 * n_common_rna_wes2 / n_genes_rna_wes2), 2)

write.csv(lenient_overlap, glue("{args$output_dir}/overlap_perc_lenient.csv"))

# Determine overlap (exact)
common <- combined_exact %>%
    mutate(common_rna_wes = (in_RNAseq + in_WES == 2), common_rna_wes2 = (in_RNAseq + in_WES2 == 2)) %>%
    group_by(wes_tumor_id) %>%
    summarise(n_common_rna_wes = sum(common_rna_wes), n_common_rna_wes2 = sum(common_rna_wes2)) %>%
    ungroup()

n_genes <- combined_exact %>%
    mutate(common_rna_wes = (in_WES > 0), common_rna_wes2 = (in_WES2 > 0)) %>%
    group_by(wes_tumor_id) %>%
    summarise(n_common_rna_wes = sum(common_rna_wes), n_common_rna_wes2 = sum(common_rna_wes2)) %>%
    ungroup() %>%
    rename(n_genes_rna_wes = n_common_rna_wes, n_genes_rna_wes2 = n_common_rna_wes2)

exact_overlap <- merge(common, n_genes, by = "wes_tumor_id") %>%
    mutate(overlap_rna_wes = 100 * n_common_rna_wes / n_genes_rna_wes, overlap_rna_wes2 = round(100 * n_common_rna_wes2 / n_genes_rna_wes2), 2)

write.csv(exact_overlap, glue("{args$output_dir}/overlap_perc_exact.csv"))
