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
        description = "Compare mutations",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/")
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
sample_sheet <- read.csv(glue("{here::here()}/000_misc_local/sample_sheet.csv"))

rnaseq_mut <- list.files(glue("{here::here()}/output/test_set/mutations_postprocessed"))
rnaseq_samples <- str_remove(rnaseq_mut, "_mutations.txt")

wes_mut <- list.files("/Users/joankant/Desktop/gaitigroup/Users/Joan/wes-mutation-calling/output/wes_paired/mutations_gene_filtered")
wes_samples <- str_remove(wes_mut, "_mutations_final.txt")

avail_mutations <- sample_sheet %>% filter(rnaseq_condition_id %in% rnaseq_samples, wes_condition_id %in% wes_samples)

for (i in seq_len(nrow(avail_mutations))) {
    current_rna_sample <- avail_mutations$rnaseq_condition_id[i]
    current_wes_sample <- avail_mutations$wes_condition_id[i]

    current_rna_mut <- read.table(glue("{here::here()}/output/test_set/mutations_postprocessed/{current_rna_sample}_mutations.txt"), sep = "\t", header = TRUE) %>%
        mutate(mutation_in_rna = 1) %>%
        select(Chr, Start, End, Ref, Alt, mutation_in_rna)

    current_wes_mut <- read.table(glue("/Users/joankant/Desktop/gaitigroup/Users/Joan/wes-mutation-calling/output/wes_paired/mutations_gene_filtered/{current_wes_sample}_mutations_final.txt"), sep = "\t", header = TRUE) %>%
        mutate(mutation_in_wes = 1) %>%
        select(Chr, Start, End, Ref, Alt, mutation_in_wes)

    combi <- merge(current_rna_mut, current_wes_mut, by = c("Chr", "Start", "End", "Ref", "Alt"), all = TRUE)
    combi[is.na(combi)] <- 0

    rnaseq_mutations <- rownames(combi)[combi$mutation_in_rna == 1]
    wes_mutations <- rownames(combi)[combi$mutation_in_wes == 1]

    create_venndiagram(
        x = list(rnaseq_mutations, wes_mutations),
        category.names = c("RNAseq", "WES"),
        filename = glue("{args$output_dir}/venndiagram_{current_rna_sample}_{current_wes_sample}.png"),
        main = glue("RNA_ID={current_rna_sample} vs WES_ID={current_wes_sample} \n(cond={avail_mutations$condition[i]})"), main.cex = 0.15
    )
}
