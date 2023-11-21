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
    args$output_dir <- glue("{here::here()}/data")
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
pacman::p_load(ggplot2)
# ---- Data Availability ---- #
tab_s4 <- read.csv("/Users/joankant/Library/CloudStorage/OneDrive-UHN/004_Projects/Lupus/Riaz-github/S4-genomic_data_per_case.csv")
log_info("\n\n", capture.output(tab_s4 %>% head(5)))
log_info("Total number of patients: ", tab_s4 %>% pull(Patient) %>% unique() %>% length())

head(tab_s4)

samples_oi <- tab_s4 %>% filter(Pre.treatment.RNA.Seq == 1, Pre.treatment.Exome == 1)

sample_sheet <- read.csv(glue("{here::here()}/misc/sample_sheet.csv"))
# sample_sheet %>% filter()

all_patients <- tab_s4 %>% pull(Patient)
wes_pre_patients <- tab_s4 %>%
    filter(Pre.treatment.Exome == 1) %>%
    pull(Patient)
rnaseq_pre_patients <- tab_s4 %>%
    filter(Pre.treatment.RNA.Seq == 1) %>%
    pull(Patient)
x <- list(all_patients = all_patients, rnaseq_pre = rnaseq_pre_patients, wes_pre = wes_pre_patients)
create_venndiagram(x = x, category.names = names(x), filename = glue("{here::here()}/output/figures/number_of_patients_pre_therapy.png"))

# ---- Mutations ---- #
tab_s3 <- data.frame(fread("/Users/joankant/Library/CloudStorage/OneDrive-UHN/004_Projects/Lupus/Riaz-github/S3-pre_therapy_nonsynonmous_mutations.csv"))

log_info("\n\n", capture.output(tab_s3 %>% head(5)))

mut_pp <- tab_s3 %>%
    group_by(Patient) %>%
    summarise(n_mutations = n())

mut_pp %>% ggplot() +
    geom_boxplot(aes(x = 1, n_mutations)) +
    geom_jitter()

mut_pp_distr <- mut_pp %>% ggplot() +
    geom_histogram(aes(x = n_mutations))
ggsave(plot = mut_pp_distr, file = glue("{here::here()}/output/figures/mut_pp_distribution.pdf"))

mut_patients <- tab_s3 %>%
    pull(Patient) %>%
    unique()
n_mut_patients <- mut_patients %>%
    length()
log_info("Number of patients: ", n_patients)

total_n_mutations <- mut_pp %>%
    pull(n_mutations) %>%
    sum()
total_n_mutations_avail_rna <- mut_pp %>%
    filter(Patient %in% rnaseq_pre_patients) %>%
    pull(n_mutations) %>%
    sum()
mean_n_mutations <- mut_pp %>%
    pull(n_mutations) %>%
    mean()
log_info("Total number of mutations: ", mean_n_mutations)
mean_n_mutations_avail_rna <- mut_pp %>%
    filter(Patient %in% rnaseq_pre_patients) %>%
    pull(n_mutations) %>%
    mean()
log_info("Mean number of mutations per patient with RNA-seq data: ", mean_n_mutations_avail_rna)
