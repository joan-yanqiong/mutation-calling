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
        description = "Prepare Samplesheet for pipeline", default_output = "misc"
    )
    parser$add_argument(
        "--gseo",
        type = "character",
        required = TRUE,
        help = "GEO metadata file"
    )
    parser$add_argument(
        "--wes_sra",
        type = "character",
        required = TRUE,
        help = "SRA metadata file for WES"
    )
    parser$add_argument(
        "--rnaseq_sra",
        type = "character",
        required = TRUE,
        help = "SRA metadata file for RNAseq"
    )
    parser$add_argument(
        "--fastq_dir",
        type = "character",
        required = TRUE,
        help = "Directory with fastq files"
    )
    parser$add_argument("--normal_samples_oi",
        type = "character",
        required = FALSE,
        help = "File with normal samples of interest", default = NULL
    )
    parser$add_argument("--tumor_samples_oi",
        type = "character",
        required = FALSE,
        help = "File with tumor samples of interest", default = NULL
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/misc")
    args$gseo <- glue("{here::here()}/misc/GSE91061_series_matrix.xlsx")
    args$fastq_dir <- "/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/data/Riaz"
    args$wes_sra <- glue("{here::here()}/misc/SraRunTable_WES.txt")
    args$rnaseq_sra <- glue("{here::here()}/misc/SraRunTable_RNA.txt")
    args$normal_samples_oi <- glue("/Users/joankant/Library/CloudStorage/OneDrive-UHN/004_Projects/Lupus/pilot_normal_samples_for_download.txt")
    args$tumor_samples_oi <- NULL
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
log_info("Load additional libraries...")
pacman::p_load(readxl)

log_info("Read metadata from GEO...")
rnaseq_meta <- read_excel(args$gseo, skip = 35)

# Lookup table to match RNAseq with WES
# Columns: rnasq_patient_id, GSM_id
log_info("Create lookup table...")
lookup_table <- data.frame(t(rnaseq_meta[1, ])) %>% rownames_to_column("patient_id")
colnames(lookup_table) <- c("rna_sample_id", "GSM_id")
lookup_table <- lookup_table[-1, ]
head(lookup_table)

# Downloaded from SRA selectors
log_info("Read metadata from SRA selectors...")
wes_meta_sra <- read.csv(args$wes_sra) %>%
    select(Sample.Name, Run) %>%
    rename(normal_id = Run, wes_patient_id = Sample.Name) %>%
    mutate(wes_patient_id = tolower(wes_patient_id))
rnaseq_meta_sra <- read.csv(args$rnaseq_sra) %>%
    select(Sample.Name, Run) %>%
    rename(GSM_id = Sample.Name, tumor_id = Run)


# Combine RNAseq and lookup table
log_info("Match metadata RNAseq and lookup table...")
rnaseq_meta_sra <- rnaseq_meta_sra %>%
    left_join(lookup_table, by = "GSM_id") %>%
    mutate(patient_id = str_split(rna_sample_id, simplify = TRUE, pattern = "_")[, 1]) %>%
    mutate(wes_patient_id = tolower(paste0(patient_id, "_norm")))

# Combine RNAseq and WES
log_info("Match metadata RNAseq and WES...")
sample_sheet <- rnaseq_meta_sra %>% left_join(wes_meta_sra, by = "wes_patient_id")

log_info("Remove patients without WES data...")
sample_sheet_filtered <- sample_sheet %>% filter(normal_id != "<NA>")
log_info("Number of available patients: {nrow(sample_sheet_filtered)}")

log_info("Adding additional information...")
sample_sheet_filtered <- sample_sheet_filtered %>%
    mutate(ix = row_number(), pair_id = paste0(normal_id, "_", tumor_id), tumor_fastq = paste0(args$fastq_dir, "/tumor/", tumor_id), normal_fastq = paste0(args$fastq_dir, "/normal/", normal_id)) # %>%
# select(ix, tumor_id, normal_id, pair_id, tumor_fastq, normal_fastq, )

log_info("Save sample sheet...")
write.table(sample_sheet_filtered, file = glue("{args$output_dir}/sample_sheet.csv"), sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)

if ((!is.null(args$normal_samples_oi)) || (!is.null(args$tumor_samples_oi))) {
    log_info("Subset sample sheet...")
    if (!is.null(args$normal_samples_oi)) {
        log_info("Additional filtering...")
        normal_samples_oi <- read.table(args$normal_samples_oi) %>% pull()
        sample_sheet_filtered <- sample_sheet_filtered %>% filter(normal_id %in% normal_samples_oi)
    }
    if (!is.null(args$tumor_samples_oi)) {
        log_info("Additional filtering...")
        tumor_samples_oi <- read.table(args$tumor_samples_oi) %>% pull()
        sample_sheet_filtered <- sample_sheet_filtered %>% filter(tumor_id %in% tumor_samples_oi)
    }
    log_info("Save sample sheet...")
    write.table(sample_sheet_filtered, file = glue("{args$output_dir}/sample_sheet_subset.csv"), sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

log_info("COMPLETED!")
