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
    args$output_dir <- "/Users/joankant/Desktop/gaitigroup/Users/Joan/000_misc/Lupus"
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(args$output_dir)

# From Github - Supplementary 4
tab_s4 <- read.csv("/Users/joankant/Desktop/gaitigroup/Users/Joan/001_data/Riaz-github/S4-genomic_data_per_case.csv")

# Determine the number of samples that are available for both WES and RNAseq
pre_treatment_avail <- tab_s4 %>%
    rename(WES_available = Pre.treatment.Exome, RNAseq_available = Pre.treatment.RNA.Seq) %>%
    select(Patient, WES_available, RNAseq_available) %>%
    mutate(class = "Pre-treatment")

on_treatment_avail <- tab_s4 %>%
    rename(WES_available = On.treatment.Exome, RNAseq_available = On.treatment.RNA.Seq) %>%
    select(Patient, WES_available, RNAseq_available) %>%
    mutate(class = "On-treatment")
all_data <- rbind(pre_treatment_avail, on_treatment_avail)

avail_paired_data <- all_data %>% filter(WES_available == 1, RNAseq_available == 1)
patient_to_include <- avail_paired_data %>%
    pull(Patient) %>%
    unique()
log_info(glue("Number of samples: {nrow(avail_paired_data)}"))
log_info(glue("Number of patients with paired data: {length(patient_to_include)}"))


write.csv(avail_paired_data %>% select(Patient, class), glue("{args$output_dir}/riaz_avail_wes_and_rna_patients.csv"), row.names = FALSE)
log_info("COMPLETED!")
