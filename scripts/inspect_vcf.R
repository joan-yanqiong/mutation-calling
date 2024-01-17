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
        description = "Inspect VCF files",
    )

    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$suffix <- "raw"
    args$output_dir <- glue("{here::here()}/output/test_set_rerun/analysis/{args$suffix}/")
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
pacman::p_load(vcfR, ggplot2, ggpubr, pbapply)

log_info("Load sample sheet...")
sample_sheet <- read.csv(glue("{here::here()}/000_misc_local/test_set.csv"))
wes_tumor_ids <- sample_sheet %>% pull(wes_tumor_id)
rna_tumor_ids <- sample_sheet %>% pull(rnaseq_tumor_id)

# Grab files
# - Raw VCF files
rna_mutect_files <- paste0(glue("/Users/joankant/Desktop/gaitigroup/Users/Joan/mutation-calling/output/test_set_rerun/tumor/"), rna_tumor_ids, "/", rna_tumor_ids, "_mutect_R2.vcf")
mutect_files <- paste0(glue("/Users/joankant/Desktop/gaitigroup/Users/Joan/wes-mutation-calling/output/test_set_rerun/wes_tumor/"), wes_tumor_ids, "/", wes_tumor_ids, "_mutect.vcf")
mutect2_files <- paste0(glue("/Users/joankant/Desktop/gaitigroup/Users/Joan/wes-mutation-calling/output/test_set_rerun/wes_tumor/"), wes_tumor_ids, "/", wes_tumor_ids, "_somatic.vcf.gz")

# - One of the filtering steps
# rna_mutect_files <- paste0(glue("/Users/joankant/Desktop/gaitigroup/Users/Joan/mutation-calling/output/test_set_rerun/tumor/"), rna_tumor_ids, "/mutations_filtered/", rna_tumor_ids, ".", args$suffix, ".recode.vcf")
# mutect_files <- paste0(glue("/Users/joankant/Desktop/gaitigroup/Users/Joan/wes-mutation-calling/output/test_set_rerun/wes_tumor/"), wes_tumor_ids, "/mutect_filtered/", wes_tumor_ids, ".", args$suffix, ".recode.vcf")
# mutect2_files <- paste0(glue("/Users/joankant/Desktop/gaitigroup/Users/Joan/wes-mutation-calling/output/test_set_rerun/wes_tumor/"), wes_tumor_ids, "/mutect2_filtered/", wes_tumor_ids, ".", args$suffix, ".recode.vcf")


# rna_mutect_files <- paste0(glue("/Users/joankant/Desktop/gaitigroup/Users/Joan/mutation-calling/output/test_set_rerun/mutations_prefiltered/"), rna_tumor_ids, "/", rna_tumor_ids, "_annovar_R2.hg19_multianno.vcf")
# mutect_files <- paste0(glue("/Users/joankant/Desktop/gaitigroup/Users/Joan/wes-mutation-calling/output/test_set_rerun/mutations_annotated/mutect/"), wes_tumor_ids, "/", wes_tumor_ids, ".hg19_multianno.vcf")
# mutect2_files <- paste0(glue("/Users/joankant/Desktop/gaitigroup/Users/Joan/wes-mutation-calling/output/test_set_rerun/mutations_annotated/mutect2/"), wes_tumor_ids, "/", wes_tumor_ids, ".hg19_multianno.vcf")

cl <- parallel::makeCluster(8)
all_rna_mutect_res <- pblapply(rna_mutect_files, format_mutations, cl = cl) %>%
    bind_rows() %>%
    mutate(method = "RNA Mutect")

# Deal with Mutect2 files
all_mutect2_res <- pblapply(mutect2_files, format_mutations) %>%
    bind_rows() %>%
    mutate(method = "Mutect2")

# Deal with Mutect files
all_mutect_res <- pblapply(mutect_files, format_mutations) %>%
    bind_rows() %>%
    mutate(method = "Mutect")

# Combine Mutect and Mutect2
all_rna_mutect_res <- all_rna_mutect_res %>%
    rowwise() %>%
    mutate(Indiv = str_split(Indiv, "tPU", simplify = TRUE)[, 1]) %>%
    filter(!str_detect(Indiv, "_")) %>%
    left_join(sample_sheet %>% select(wes_tumor_id, rnaseq_tumor_id), by = c("Indiv" = "rnaseq_tumor_id")) %>%
    select(-Indiv) %>%
    rename(Indiv = wes_tumor_id)

all_mutect_combined <- plyr::rbind.fill(all_rna_mutect_res, all_mutect2_res, all_mutect_res) %>%
    separate(gt_AD, into = c("AD_ref", "AD_alt"), sep = ",", remove = FALSE) %>%
    mutate(AD_ref = as.numeric(AD_ref), AD_alt = as.numeric(AD_alt), gt_DP = as.numeric(gt_DP))
log_info(glue("Number of mutations: {nrow(all_mutect_combined)}"))

# Distributions of AD and DPcol
# colors <- c(colorjam::rainbowJam(3))
colors <- c("blue", "red", "green")
names(colors) <- c("RNA Mutect", "Mutect2", "Mutect")

p_dp <- ggplot(all_mutect_combined, aes(x = gt_DP, fill = method)) +
    geom_histogram(binwidth = 1, alpha = 0.6, position = "identity") +
    facet_wrap(~Indiv, scales = "free", ncol = 2) +
    custom_theme() +
    scale_fill_manual(values = colors) +
    labs(x = "DP", y = "Count", title = "Distribution of DP")
ggsave(plot = p_dp, filename = glue("{args$output_dir}/dp.pdf"), width = 8, height = 30, limitsize = FALSE)

p_ad <- ggplot(all_mutect_combined, aes(x = AD_alt, fill = method)) +
    geom_histogram(binwidth = 1, alpha = 0.6, position = "identity") +
    facet_wrap(~Indiv, scales = "free", ncol = 2) +
    custom_theme() +
    scale_fill_manual(values = colors) +
    labs(x = "AD alt", y = "Count", title = "Distribution of AD alt")
ggsave(plot = p_ad, filename = glue("{args$output_dir}/ad_alt.pdf"), width = 8, height = 30, limitsize = FALSE)

# Filtering based on AD and DP
all_mutect_filtered <- all_mutect_combined %>% filter(AD_alt >= 5, gt_DP >= 20)

mutations_per_sample <- all_mutect_filtered %>%
    group_by(method, Indiv) %>%
    summarise(n = n()) %>%
    ungroup()

p_mut_sample <- ggplot(mutations_per_sample, aes(x = method, y = n, fill = method)) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    geom_boxplot() +
    geom_point() +
    custom_theme() +
    geom_line(aes(group = Indiv)) +
    labs(y = "Number of mutations", title = "Number of mutations per sample")
ggsave(glue("{args$output_dir}/mutations_per_sample.pdf"), p_mut_sample, width = 10, height = 8, limitsize = FALSE)


# Mutations per sample - prepare text file
tmp <- mutations_per_sample %>%
    pivot_wider(names_from = method, values_from = n)
sample_sheet %>%
    select(rnaseq_tumor_id, wes_tumor_id) %>%
    left_join(tmp, by = c("wes_tumor_id" = "Indiv")) %>%
    write_csv(glue("{args$output_dir}/mutations_per_sample.csv"))
