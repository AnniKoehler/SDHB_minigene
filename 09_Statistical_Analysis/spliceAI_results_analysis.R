# SpliceAI / transcript summary

# Description:
#   Loads a variant table (with SpliceAI delta score columns), 
#   computes summary stats and writes a single-sheet Excel.
#
#   Key assumptions:
#     - input table has (at least) columns:
#         variant (name), hgvsc, genomic_position, ds_ag, ds_al, ds_dg, ds_dl


#   Canonical splice sites:
#     - Positive strand:  donor = exon_end + {1,2};    acceptor = exon_start - {1,2}
#     - Negative strand:  donor = exon_start - {1,2};  acceptor = exon_end + {1,2}
#


# ----------------------------
# Load necessary libraries
# ----------------------------
library(readr)
library(readxl)
library(data.table)
library(dplyr)
library(tidyr)
library(janitor)
library(stringr)
library(writexl)

# (Optional) Working directory
# setwd("~/your/path/here")

# ----------------------------
# SDHB minigene specific parameters
# ----------------------------

# 1) Input file 
INPUT_FILE <- "input_file.xlsx"

# 2) SDHB Exon coordinates

exons <- data.frame(
  exon_number = c(2, 3, 4, 5),
  exon_start  = c(17027749, 17028600, 17033060, 17044761),
  exon_end    = c(17027865, 17028736, 17033145, 17044888)
)

# 3) Gene strand and canonical window 
GENE_STRAND <- "-"           
CANONICAL_WINDOW <- 2        

# 4) Thresholds 
GAIN_THRESHOLD  <- 0.5       # for ds_ag / ds_dg "gain"
COMBO_THRESHOLD <- 0.5       # for (ds_dl > thr & ds_al > thr)

# 5) Exclusions 
EXCLUDE_HGVSC    <- c("c.300T>C", "c.225T>C", "MG_WT")  # matches hgvsc_short
EXCLUDE_VARIANTS <- character(0)                        # matches `variant`



# ----------------------------
# Data load (raw base)
# ----------------------------
df_all <- read_excel(INPUT_FILE) %>%
  mutate(
    # short HGVSc (keep MG_WT literal if present)
    hgvsc_short = ifelse(variant == "MG_WT",
                         "MG_WT",
                         stringr::str_extract(hgvsc, "(?<=:).*"))
  )

# ----------------------------
# Exclude variants (generalized)
# ----------------------------
df_filtered <- df_all %>%
  filter(
    !(hgvsc_short %in% EXCLUDE_HGVSC),
    !(variant %in% EXCLUDE_VARIANTS)
  )

# -----------------------------------------------------
# A) How many unique variants are there in total?
# -----------------------------------------------------
num_variants <- df_filtered %>%
  distinct(variant) %>%
  nrow()

# -----------------------------------------------------
# B) For each variant, compute the maximum SpliceAI score
#    across ds_ag, ds_al, ds_dg, ds_dl.
#    Then (B1) the lowest per-variant maximum,
#         (B2) the highest per-variant maximum,
#         (B3) the average per-variant maximum.
# -----------------------------------------------------
variant_max_scores <- df_filtered %>%
  group_by(variant) %>%
  summarise(
    max_score = max(ds_ag, ds_al, ds_dg, ds_dl, na.rm = TRUE),
    .groups = "drop"
  )

lowest_of_maxima  <- min(variant_max_scores$max_score, na.rm = TRUE)
highest_of_maxima <- max(variant_max_scores$max_score, na.rm = TRUE)
avg_of_maxima     <- mean(variant_max_scores$max_score, na.rm = TRUE)

# -----------------------------------------------------
# C) Identify variants on canonical splice sites (strand-aware).
#    Positive strand:
#      donor    = exon_end + {1..CANONICAL_WINDOW}
#      acceptor = exon_start - {1..CANONICAL_WINDOW}
#    Negative strand:
#      donor    = exon_start - {1..CANONICAL_WINDOW}
#      acceptor = exon_end + {1..CANONICAL_WINDOW}
# -----------------------------------------------------
# 1) Cartesian join
df_sites <- df_filtered %>%
  mutate(join_key = 1) %>%
  inner_join(exons %>% mutate(join_key = 1), by = "join_key") %>%
  select(-join_key)

# 2) Compute donor / acceptor windows once, then flag
df_sites <- df_sites %>%
  rowwise() %>%
  mutate(
    is_canonical_donor = if (GENE_STRAND == "+") {
      any(genomic_position == (exon_end + seq_len(CANONICAL_WINDOW)))
    } else {
      any(genomic_position == (exon_start - seq_len(CANONICAL_WINDOW)))
    },
    is_canonical_acceptor = if (GENE_STRAND == "+") {
      any(genomic_position == (exon_start - seq_len(CANONICAL_WINDOW)))
    } else {
      any(genomic_position == (exon_end + seq_len(CANONICAL_WINDOW)))
    }
  ) %>%
  ungroup()

# 3) Collapse to per-variant/position flags
df_sites_flagged <- df_sites %>%
  group_by(variant, genomic_position) %>%
  summarise(
    on_donor    = any(is_canonical_donor),
    on_acceptor = any(is_canonical_acceptor),
    .groups = "drop"
  )

canonical_variants <- df_sites_flagged %>%
  filter(on_donor | on_acceptor) %>%
  distinct(variant, genomic_position, on_donor, on_acceptor)

num_canonical_variants <- canonical_variants %>% distinct(variant) %>% nrow()

# -----------------------------------------------------
# D) Among canonical-site variants: split into acceptor vs donor.
# -----------------------------------------------------
num_acceptor_variants <- canonical_variants %>%
  filter(on_acceptor) %>%
  distinct(variant) %>%
  nrow()

num_donor_variants <- canonical_variants %>%
  filter(on_donor) %>%
  distinct(variant) %>%
  nrow()

# --- E) Identify all intronic variants, then subtract canonical ones ---

exons_sorted <- exons %>% arrange(exon_number)

is_intron <- function(pos) {
  if (is.na(pos) || !is.numeric(pos)) return(FALSE)
  for (i in seq_len(nrow(exons_sorted) - 1)) {
    if (pos > exons_sorted$exon_end[i] && pos < exons_sorted$exon_start[i + 1]) {
      return(TRUE)
    }
  }
  return(FALSE)
}

variant_positions <- df_filtered %>%
  distinct(variant, genomic_position)

variant_positions <- variant_positions %>%
  rowwise() %>%
  mutate(is_intron = is_intron(genomic_position)) %>%
  ungroup()

all_intronic <- variant_positions %>%
  filter(is_intron) %>%
  select(variant, genomic_position)

intronic_only <- all_intronic %>%
  anti_join(
    canonical_variants %>% select(variant, genomic_position),
    by = c("variant", "genomic_position")
  )

num_intronic_only <- intronic_only %>% distinct(variant) %>% nrow()

# -----------------------------------------------------
# F) How many variants lie inside exons?
# -----------------------------------------------------
exonic_variants <- variant_positions %>%
  rowwise() %>%
  mutate(
    is_exonic = any(
      genomic_position >= exons_sorted$exon_start &
        genomic_position <= exons_sorted$exon_end
    )
  ) %>%
  ungroup() %>%
  filter(is_exonic)

num_exonic_variants <- exonic_variants %>% distinct(variant) %>% nrow()

# -----------------------------------------------------
# G) Variants where (ds_dl > COMBO_THRESHOLD AND ds_al > COMBO_THRESHOLD) in same row
# -----------------------------------------------------
dl_al_same_row_variants <- df_filtered %>%
  filter(ds_dl > COMBO_THRESHOLD, ds_al > COMBO_THRESHOLD) %>%
  distinct(variant) %>%
  nrow()

# -----------------------------------------------------
# H) Acceptor-Gain (ds_ag > GAIN_THRESHOLD): count + min/max
# -----------------------------------------------------
df_ag_high <- df_filtered %>% filter(ds_ag > GAIN_THRESHOLD)
num_ag_variants <- df_ag_high %>% distinct(variant) %>% nrow()
min_ds_ag <- ifelse(nrow(df_ag_high) > 0, min(df_ag_high$ds_ag, na.rm = TRUE), NA)
max_ds_ag <- ifelse(nrow(df_ag_high) > 0, max(df_ag_high$ds_ag, na.rm = TRUE), NA)

# -----------------------------------------------------
# I) Donor-Gain (ds_dg > GAIN_THRESHOLD): count + min/max
# -----------------------------------------------------
df_dg_high <- df_filtered %>% filter(ds_dg > GAIN_THRESHOLD)
num_dg_variants <- df_dg_high %>% distinct(variant) %>% nrow()
min_ds_dg <- ifelse(nrow(df_dg_high) > 0, min(df_dg_high$ds_dg, na.rm = TRUE), NA)
max_ds_dg <- ifelse(nrow(df_dg_high) > 0, max(df_dg_high$ds_dg, na.rm = TRUE), NA)

# -----------------------------------------------------
# J) Assemble results into a summary tibble and write to Excel (summary only)
# -----------------------------------------------------
summary_table <- tibble::tibble(
  Metric = c(
    "Total number of unique variants",
    "Lowest per-variant maximum SpliceAI score",
    "Highest per-variant maximum SpliceAI score",
    "Average per-variant maximum SpliceAI score",
    "Number of variants on canonical splice sites",
    "  … of which are acceptor-site variants",
    "  … of which are donor-site variants",
    "Number of intronic variants (excluding canonical sites)",
    "Number of exonic variants",
    sprintf("Number of variants with (ds_dl > %.2f AND ds_al > %.2f) in same row", COMBO_THRESHOLD, COMBO_THRESHOLD),
    sprintf("Number of variants with Acceptor-Gain (ds_ag > %.2f)", GAIN_THRESHOLD),
    "  … lowest ds_ag among those rows",
    "  … highest ds_ag among those rows",
    sprintf("Number of variants with Donor-Gain (ds_dg > %.2f)", GAIN_THRESHOLD),
    "  … lowest ds_dg among those rows",
    "  … highest ds_dg among those rows"
  ),
  Value = c(
    num_variants,
    round(lowest_of_maxima, 3),
    round(highest_of_maxima, 3),
    round(avg_of_maxima, 3),
    num_canonical_variants,
    num_acceptor_variants,
    num_donor_variants,
    num_intronic_only,
    num_exonic_variants,
    dl_al_same_row_variants,
    num_ag_variants,
    ifelse(is.na(min_ds_ag), NA, round(min_ds_ag, 3)),
    ifelse(is.na(max_ds_ag), NA, round(max_ds_ag, 3)),
    num_dg_variants,
    ifelse(is.na(min_ds_dg), NA, round(min_ds_dg, 3)),
    ifelse(is.na(max_ds_dg), NA, round(max_ds_dg, 3))
  )
)

# Print to console (as before)
print(summary_table)

# Write a single-sheet Excel (summary only)
writexl::write_xlsx(
  list(summary = summary_table),
  "spliceai_summary.xlsx"
)


