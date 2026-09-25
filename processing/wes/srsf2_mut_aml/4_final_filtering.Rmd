---
title: "WES somatic variant filtering"
output: html_document
---

```{r Load the VCFS}
# Load required packages
library(vcfR)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(ggrepel)
library(qs2)

# Map sample names to the annotated Strelka2 VCFs (output of 03_strelka2_filter_annotate.sh)
vcf_list <- read.csv(text = "
Desired_name_in_R,path_to_vcf
P1,/path/to/project/output/strelka2/P1_tumor_vs_P1_germ/results/variants/merged_parsed.pass.dp50.ad5.vaf5.annotated.vcf.gz
P2,/path/to/project/output/strelka2/P2_tumor_vs_P2_germ/results/variants/merged_parsed.pass.dp50.ad5.vaf5.annotated.vcf.gz
P3,/path/to/project/output/strelka2/P3_tumor_vs_P3_germ/results/variants/merged_parsed.pass.dp50.ad5.vaf5.annotated.vcf.gz
P4,/path/to/project/output/strelka2/P4_tumor_vs_P4_germ/results/variants/merged_parsed.pass.dp50.ad5.vaf5.annotated.vcf.gz
P6,/path/to/project/output/strelka2/P6_tumor_vs_P6_germ/results/variants/merged_parsed.pass.dp50.ad5.vaf5.annotated.vcf.gz
P7,/path/to/project/output/strelka2/P7_tumor_vs_P7_germ/results/variants/merged_parsed.pass.dp50.ad5.vaf5.annotated.vcf.gz
P10,/path/to/project/output/strelka2/P10_tumor_vs_P10_germ/results/variants/merged_parsed.pass.dp50.ad5.vaf5.annotated.vcf.gz
P10_best,/path/to/project/output/strelka2/P10_tumor_vs_P10_germ/results/variants/merged_parsed.pass.dp50.ad5.vaf5.annotated.besteffect.vcf.gz
P11,/path/to/project/output/strelka2/P11_tumor_vs_P11_germ/results/variants/merged_parsed.pass.dp50.ad5.vaf5.annotated.vcf.gz
P11_best,/path/to/project/output/strelka2/P11_tumor_vs_P11_germ/results/variants/merged_parsed.pass.dp50.ad5.vaf5.annotated.besteffect.vcf.gz
P12,/path/to/project/output/strelka2/P12_tumor_vs_P12_germ/results/variants/merged_parsed.pass.dp50.ad5.vaf5.annotated.vcf.gz
P12_best,/path/to/project/output/strelka2/P12_tumor_vs_P12_germ/results/variants/merged_parsed.pass.dp50.ad5.vaf5.annotated.besteffect.vcf.gz
P13,/path/to/project/output/strelka2/P13_tumor_vs_P13_germ/results/variants/merged_parsed.pass.dp50.ad5.vaf5.annotated.vcf.gz
P13_best,/path/to/project/output/strelka2/P13_tumor_vs_P13_germ/results/variants/merged_parsed.pass.dp50.ad5.vaf5.annotated.besteffect.vcf.gz
P14,/path/to/project/output/strelka2/P14_tumor_vs_P14_germ/results/variants/merged_parsed.pass.dp50.ad5.vaf5.annotated.vcf.gz
P14_best,/path/to/project/output/strelka2/P14_tumor_vs_P14_germ/results/variants/merged_parsed.pass.dp50.ad5.vaf5.annotated.besteffect.vcf.gz
P15,/path/to/project/output/strelka2/P15_tumor_vs_P15_germ/results/variants/merged_parsed.pass.dp50.ad5.vaf5.annotated.vcf.gz
P15_best,/path/to/project/output/strelka2/P15_tumor_vs_P15_germ/results/variants/merged_parsed.pass.dp50.ad5.vaf5.annotated.besteffect.vcf.gz
P16,/path/to/project/output/strelka2/P16_tumor_vs_P16_germ/results/variants/merged_parsed.pass.dp50.ad5.vaf5.annotated.vcf.gz
P16_best,/path/to/project/output/strelka2/P16_tumor_vs_P16_germ/results/variants/merged_parsed.pass.dp50.ad5.vaf5.annotated.besteffect.vcf.gz
P17,/path/to/project/output/strelka2/P17_tumor_vs_P17_germ/results/variants/merged_parsed.pass.dp50.ad5.vaf5.annotated.vcf.gz
P17_best,/path/to/project/output/strelka2/P17_tumor_vs_P17_germ/results/variants/merged_parsed.pass.dp50.ad5.vaf5.annotated.besteffect.vcf.gz
", stringsAsFactors = FALSE)

# Holder list
all_vcf_df <- list()

# ---- Process all VCFs ----
for (i in 1:nrow(vcf_list)) {
  sample_id <- vcf_list$Desired_name_in_R[i]
  vcf_path <- vcf_list$path_to_vcf[i]

  if (file.exists(vcf_path)) {
    message("Reading: ", vcf_path)
    vcf <- read.vcfR(vcf_path, verbose = FALSE)

    fix_df <- as.data.frame(vcf@fix[, 1:7], stringsAsFactors = FALSE)
    colnames(fix_df) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")

    info_fields <- extract_info_tidy(vcf)
    fix_df$FUNCOTATION <- if ("FUNCOTATION" %in% colnames(info_fields)) info_fields$FUNCOTATION else NA

    # Assuming 2nd sample column is always tumor
    fix_df$GT     <- extract.gt(vcf, element = "GT")[, 2]
    fix_df$AD     <- extract.gt(vcf, element = "AD")[, 2]
    fix_df$DP     <- as.numeric(extract.gt(vcf, element = "DP")[, 2])
    fix_df$GQ     <- as.numeric(extract.gt(vcf, element = "GQ")[, 2])
    fix_df$MED_DP <- as.numeric(extract.gt(vcf, element = "MED_DP")[, 2])
    fix_df$MIN_DP <- as.numeric(extract.gt(vcf, element = "MIN_DP")[, 2])
    fix_df$PL     <- extract.gt(vcf, element = "PL")[, 2]
    fix_df$VAF    <- as.numeric(extract.gt(vcf, element = "VAF")[, 2])

    fix_df$SAMPLE_ID <- sample_id
    all_vcf_df[[sample_id]] <- fix_df
  } else {
    warning("VCF file not found: ", vcf_path)
  }
}

# Combine all data
combined_vcf_df <- bind_rows(all_vcf_df)

# ---- FUNCOTATION Parsing ----
selected_indices <- c(1,2,3,4,5,6,8,9,10,11,12,13,15,16,17,18,19,20,145)
simplified_names <- c(
  "gene", "build", "chrom", "start", "end", "variant_class",
  "variant_type", "ref", "allele1", "allele2", "genome_change", "transcript",
  "exon", "transcript_pos", "cdna_change", "codon_change", "protein_change",
  "gc_content", "dbsnp_id"
)

cleaned_funcotations <- combined_vcf_df %>%
  mutate(FUNCOTATION = str_remove_all(FUNCOTATION, "^\\[|\\]$")) %>%
  mutate(FUNC_SPLIT = str_split(FUNCOTATION, "\\|")) %>%
  rowwise() %>%
  mutate(FUNC_PARSED = list({
    vals <- FUNC_SPLIT
    vals <- if (length(vals) < max(selected_indices)) rep(vals, length.out = max(selected_indices)) else vals
    vals[selected_indices]
  })) %>%
  ungroup()

func_df <- as.data.frame(do.call(rbind, cleaned_funcotations$FUNC_PARSED), stringsAsFactors = FALSE)
colnames(func_df) <- simplified_names

# Final dataframe
final_df <- bind_cols(
  combined_vcf_df %>% dplyr::select(-FUNCOTATION),
  func_df
)

# Separate AD into AD_REF and AD_ALT
final_df <- final_df %>%
  separate(AD, into = c("AD_REF", "AD_ALT"), sep = ",", convert = TRUE)

message("All VCFs imported and parsed into final_df.")
rm(i,sample_id,selected_indices,simplified_names,vcf_path,vcf,vcf_list,info_fields,func_df,fix_df,combined_vcf_df,cleaned_funcotations,all_vcf_df)
```

```{r Filter out intronic variants}
final_df <- final_df %>%
  dplyr::filter(
    !variant_class %in% c("FIVE_PRIME_FLANK", "FIVE_PRIME_UTR", "IGR", "THREE_PRIME_UTR", "INTRON", "RNA")
  )

# Get all Px_best sample IDs
best_samples <- grep("_best$", unique(final_df$SAMPLE_ID), value = TRUE)
base_samples <- gsub("_best$", "", best_samples)  # e.g. P10_best -> P10

merged_df <- final_df  # start with full dataframe

for (i in seq_along(best_samples)) {
  px      <- base_samples[i]
  px_best <- best_samples[i]

  px_rows      <- final_df %>% dplyr::filter(SAMPLE_ID == px)
  px_best_rows <- final_df %>% dplyr::filter(SAMPLE_ID == px_best)

  # Find mutations in Px_best that are NOT in Px (by CHROM+POS+REF+ALT)
  px_keys      <- px_rows      %>% select(CHROM, POS, REF, ALT)
  px_best_only <- px_best_rows %>%
    anti_join(px_keys, by = c("CHROM", "POS", "REF", "ALT")) %>%
    mutate(SAMPLE_ID = px)  # relabel as Px

  # Remove the Px_best rows from merged_df, add the unique-to-best rows
  merged_df <- merged_df %>%
    dplyr::filter(SAMPLE_ID != px_best) %>%
    bind_rows(px_best_only)
}

# Verify
table(merged_df$SAMPLE_ID)
sum(final_df$protein_change == "" | is.na(final_df$protein_change))
qs_save(merged_df, "/path/to/project/objects/wes_info.qs2")
```
