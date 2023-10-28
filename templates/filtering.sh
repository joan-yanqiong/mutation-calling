#!/usr/bin/env bash
<<comment

--exclude-bed Include or exclude a set of sites on the basis of a BED file. Only the first three columns (chrom, chromStart and chromEnd) are required.
The BED file is expected to have a header line.
A site will be kept or excluded if any part of any allele (REF or ALT) at a site is within the range of one of the BED entries.
--remove-filtered-all Removes all sites with a FILTER flag other than PASS.

Ref: https://vcftools.github.io/man_latest.html

comment

echo "\$(date)  Loading module..."
module load vcftools/0.1.15

echo "\$(date)  Create output directory..."
mkdir -p "${sample_id}/mutations_filtered"

echo "\$(date)  Enter output directory..."
cd "${sample_id}/mutations_filtered"

echo "\$(date) Keep only sites with a FILTER = PASS..."
vcftools --vcf "../../${mutect_vcf}" \
    --remove-filtered-all \
    --recode \
    --recode-INFO-all \
    --out ${sample_id}.pass

vcftools --vcf ${sample_id}.pass.recode.vcf \
    --recode \
    --recode-INFO-all \
    --mac ${min_alt_counts} \
    --out ${sample_id}.mac

echo "\$(date)  Remove sites from ExAC database with minor allele frequence (MAF) > 5%..."
vcftools --vcf ${sample_id}.mac.recode.vcf \
    --exclude-bed "../../${exac_bed}" \
    --recode \
    --recode-INFO-all \
    --out ${sample_id}.exac

echo "\$(date)  Remove sites from DARNED database..."
vcftools --vcf ${sample_id}.exac.recode.vcf \
    --exclude-bed "../../${darned_bed}" \
    --recode \
    --recode-INFO-all \
    --out ${sample_id}.darned

echo "\$(date)  Remove sites from RADAR database..."
vcftools --vcf ${sample_id}.darned.recode.vcf \
    --exclude-bed "../../${radar_bed}" \
    --recode \
    --recode-INFO-all \
    --out ${sample_id}.radar

echo "\$(date)  Removing sites falling into pseudogenes..."
vcftools --vcf ${sample_id}.radar.recode.vcf \
    --exclude-bed "../../${pseudo_genes_bed}" \
    --recode \
    --recode-INFO-all \
    --out ${sample_id}.a_pseudo_genes

vcftools --vcf ${sample_id}.a_pseudo_genes.recode.vcf \
    --exclude-bed "../../${pseudo_genes_bed2}" \
    --recode \
    --recode-INFO-all \
    --out ${sample_id}.b_pseudo_genes

echo "\$(date)  COMPLETED!"