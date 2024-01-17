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
        description = "Filtering + combining all samples",
    )
    parser$add_argument("--sample_sheet",
        type = "character", help = "Path to sample sheet"
    )
    parser$add_argument("--n_cores",
        type = "integer", help = "Number of cores to use", default = 8
    )
    parser$add_argument("--input_dir",
        type = "character", help = "Path to input directory"
    )
    parser$add_argument("--min_ad_alt",
        type = "integer", help = "Minimum AD alt", default = 5
    )
    parser$add_argument("--min_dp",
        type = "integer", help = "Minimum DP", default = 20
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/test_set_rerun/")
    args$sample_sheet <- glue("{here::here()}/000_misc_local/test_set.csv")
    args$n_cores <- 8
    args$input_dir <- glue("{here::here()}/output/test_set_rerun/")
    args$min_ad_alt <- 5
    args$min_dp <- 20
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
pacman::p_load(vcfR, pbapply)

log_info("Load sample sheet...")
sample_sheet <- read.csv(args$sample_sheet)

log_info("Extract sample IDs from sample sheet...")
rna_tumor_ids <- sample_sheet %>% pull(rnaseq_tumor_id)

log_info("Extract annotated VCF files...")
rna_mutect_files <- paste0(glue("{args$input_dir}/mutations_prefiltered/"), rna_tumor_ids, "/", rna_tumor_ids, "_annovar_R2.hg19_multianno.vcf")

log_info("Load mutations and format...")
cl <- parallel::makeCluster(args$n_cores)
all_rna_mutect_res <- pblapply(rna_mutect_files, format_mutations, cl = cl) %>%
    bind_rows() %>%
    rowwise() %>%
    mutate(Indiv = str_split(Indiv, "tPU", simplify = TRUE)[, 1]) %>%
    filter(!str_detect(Indiv, "_")) %>%
    left_join(sample_sheet %>% select(wes_tumor_id, rnaseq_tumor_id), by = c("Indiv" = "rnaseq_tumor_id")) %>%
    rename(rnaseq_tumor_id = Indiv) %>%
    bind_rows() %>%
    separate(gt_AD, into = c("AD_ref", "AD_alt"), sep = ",", remove = FALSE) %>%
    mutate(AD_ref = as.numeric(AD_ref), AD_alt = as.numeric(AD_alt), gt_DP = as.numeric(gt_DP), method = "RNA Mutect")

# Filtering based on AD and DP
log_info(glue("Number of mutations before filtering: {nrow(all_rna_mutect_res)}"))
all_mutect_filtered <- all_rna_mutect_res %>% filter(AD_alt >= args$min_ad_alt, gt_DP >= args$min_dp)
log_info(glue("Number of mutations after AD + DP filtering: {nrow(all_mutect_filtered)}"))

saveRDS(all_mutect_filtered, glue("{args$output_dir}/all_mutect_filtered_dp_ad.rds"))

all_mutect_filtered <- all_mutect_filtered %>% filter(Func.refGene == "exonic")
log_info(glue("Number of mutations after filtering: {nrow(all_mutect_filtered)}"))

saveRDS(all_mutect_filtered, glue("{args$output_dir}/all_mutect_filtered_dp_ad_exonic.rds"))

log_info("COMPLETED!")
