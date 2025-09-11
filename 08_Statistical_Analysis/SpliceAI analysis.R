#Splice AI analysis 

# Transcript analysis 

# Load all necessary libraries
library(readr)
library(readxl)
library(data.table)
library(dplyr)
library(tidyr)
library(janitor)
library(stringr)
library(writexl)
# (ggplot2/tidyverse nicht nötig, da hier kein Plot)

# Set working directory
setwd()

# ----------------------------
# Daten laden (rohe Basis)
# ----------------------------
df_all <- read_excel("transcript_table.xlsx") %>%
  clean_names() %>%
  rename(
    hgvsg_hg38 = hgv_sg_hg38,
    hgvsc      = hgv_sc
  ) %>%
  mutate(
    # verkürzte HGVSc, außer "MG_WT"
    hgvsc_short = ifelse(variant == "MG_WT",
                         "MG_WT",
                         stringr::str_extract(hgvsc, "(?<=:).*"))
  )


# Load necessary libraries
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(writexl)

# 1. Define exon boundaries (exon number, start, end)
exons <- data.frame(
  exon_number = c(2, 3, 4, 5),
  exon_start  = c(17027749, 17028600, 17033060, 17044761),
  exon_end    = c(17027865, 17028736, 17033145, 17044888)
)



# 3. Exclude the two specified variants ("c.300T>C" and "c.225T>C")
df_filtered <- df_all %>%
  filter(!hgvsc_short %in% c("c.300T>C", "c.225T>C", "MG_WT"))

# -----------------------------------------------------
# A) How many unique variants are there in total?
# -----------------------------------------------------
num_variants <- df_filtered %>%
  distinct(variant) %>%
  nrow()

# -----------------------------------------------------
# B) For each variant, compute the maximum SpliceAI score
#    across ds_ag, ds_al, ds_dg, ds_dl.
#    Then find (B1) the lowest of those per-variant maxima,
#             (B2) the highest of those per-variant maxima.
# -----------------------------------------------------
variant_max_scores <- df_filtered %>%
  group_by(variant) %>%
  summarise(
    max_score = max(ds_ag, ds_al, ds_dg, ds_dl, na.rm = TRUE)
  )

lowest_of_maxima <- min(variant_max_scores$max_score, na.rm = TRUE)
highest_of_maxima <- max(variant_max_scores$max_score, na.rm = TRUE)
avg_of_maxima <- mean(variant_max_scores$max_score, na.rm = TRUE)

# -----------------------------------------------------
# C) Identify all variants that lie on canonical splice sites.
#    - Donor site = exon_end + 1 or exon_end + 2
#    - Acceptor site = exon_start - 1 or exon_start - 2
#    Then count how many distinct variants are on any canonical site.
# -----------------------------------------------------
# 1. Perform a cartesian join of every row with every exon
df_sites <- df_filtered %>%
  mutate(join_key = 1) %>%
  inner_join(exons %>% mutate(join_key = 1), by = "join_key") %>%
  select(-join_key) %>%
  mutate(
    is_canonical_donor = (genomic_position == exon_start - 1) | (genomic_position == exon_start - 2),
    is_canonical_acceptor = (genomic_position == exon_end + 1) | (genomic_position == exon_end + 2)
  )

# 2. For each variant & genomic_position, check if any match was donor/acceptor
df_sites_flagged <- df_sites %>%
  group_by(variant, genomic_position) %>%
  summarise(
    on_donor = any(is_canonical_donor),
    on_acceptor = any(is_canonical_acceptor)
  ) %>%
  ungroup()

# 3. Filter to those that are on any canonical site
canonical_variants <- df_sites_flagged %>%
  filter(on_donor | on_acceptor) %>%
  distinct(variant, genomic_position, on_donor, on_acceptor)

num_canonical_variants <- nrow(canonical_variants)

# -----------------------------------------------------
# D) Among those canonical‐site variants, split into
#    D1) Acceptor vs. D2) Donor.
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

# 1. Sort exons by exon_number
exons_sorted <- exons %>% arrange(exon_number)

# 2. Helper function: is the position strictly between any two consecutive exons?
is_intron <- function(pos) {
  if (is.na(pos) || !is.numeric(pos)) {
    return(FALSE)
  }
  for (i in seq_len(nrow(exons_sorted) - 1)) {
    if (pos > exons_sorted$exon_end[i] && pos < exons_sorted$exon_start[i + 1]) {
      return(TRUE)
    }
  }
  return(FALSE)
}

# 3. Build a table of one row per unique variant & position
variant_positions <- df_filtered %>%
  distinct(variant, genomic_position)

# 4. Mark every position that is intronic (without worrying about canonical yet)
variant_positions <- variant_positions %>%
  rowwise() %>%
  mutate(is_intron = is_intron(genomic_position)) %>%
  ungroup()

# 5. Keep only those marked intronic
all_intronic <- variant_positions %>%
  filter(is_intron) %>%
  select(variant, genomic_position)

# 6. Subtract out any variant/position that appeared in canonical_variants
intronic_only <- all_intronic %>%
  anti_join(
    canonical_variants %>% select(variant, genomic_position),
    by = c("variant", "genomic_position")
  )

num_intronic_only <- nrow(intronic_only)
cat("Number of intronic variants (excluding canonical splice sites): ", num_intronic_only, "\n")



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

num_exonic_variants <- nrow(exonic_variants)

# -----------------------------------------------------
# G) For how many variants do we have, in the *same row*, ds_dl > 0.5 AND ds_al > 0.5?
# -----------------------------------------------------
dl_al_same_row_variants <- df_filtered %>%
  filter(ds_dl > 0.5, ds_al > 0.5) %>%
  distinct(variant) %>%
  nrow()

# -----------------------------------------------------
# H) How many variants have an Acceptor‐Gain prediction (ds_ag > 0.5)?
#    Also report the minimum and maximum ds_ag among those rows.
# -----------------------------------------------------
df_ag_high <- df_filtered %>% filter(ds_ag > 0.5)
num_ag_variants <- df_ag_high %>% distinct(variant) %>% nrow()
min_ds_ag <- ifelse(nrow(df_ag_high) > 0, min(df_ag_high$ds_ag, na.rm = TRUE), NA)
max_ds_ag <- ifelse(nrow(df_ag_high) > 0, max(df_ag_high$ds_ag, na.rm = TRUE), NA)

# -----------------------------------------------------
# I) How many variants have a Donor‐Gain prediction (ds_dg > 0.5)?
#    Also report the minimum and maximum ds_dg among those rows.
# -----------------------------------------------------
df_dg_high <- df_filtered %>% filter(ds_dg > 0.5)
num_dg_variants <- df_dg_high %>% distinct(variant) %>% nrow()
min_ds_dg <- ifelse(nrow(df_dg_high) > 0, min(df_dg_high$ds_dg, na.rm = TRUE), NA)
max_ds_dg <- ifelse(nrow(df_dg_high) > 0, max(df_dg_high$ds_dg, na.rm = TRUE), NA)

# -----------------------------------------------------
# J) Assemble all results into a summary tibble and write to Excel
# -----------------------------------------------------
summary_table <- tibble::tibble(
  Metric = c(
    "Total number of unique variants",
    "Lowest per-variant maximum SpliceAI score",
    "Highest per-variant maximum SpliceAI score",
    "Average per-variant maximum SpliceAI score",   # <— neu
    "Number of variants on canonical splice sites",
    " … of which are acceptor-site variants",
    " … of which are donor-site variants",
    "Number of intronic variants (excluding canonical sites)",
    "Number of exonic variants",
    "Number of variants with (ds_dl > 0.5 AND ds_al > 0.5) in same row",
    "Number of variants with Acceptor-Gain (ds_ag > 0.5)",
    " … lowest ds_ag among those rows",
    " … highest ds_ag among those rows",
    "Number of variants with Donor-Gain (ds_dg > 0.5)",
    " … lowest ds_dg among those rows",
    " … highest ds_dg among those rows"
  ),
  Value = c(
    num_variants,
    round(lowest_of_maxima, 3),
    round(highest_of_maxima, 3),
    round(avg_of_maxima, 3),                        # <— neu
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

# Print the summary to console
print(summary_table)

# --- Variants mit Eintrag in "Control" erfassen ---
df_control <- df_filtered %>%
  filter(!is.na(control), stringr::str_trim(as.character(control)) != "") %>%
  mutate(control = as.character(control)) %>%
  group_by(variant, hgvsc_short) %>%
  summarise(control_note = paste(unique(control), collapse = "; "), .groups = "drop") %>%
  arrange(variant)

num_control_variants <- nrow(df_control)

# Optional: Zahl ins Summary aufnehmen
summary_table <- dplyr::bind_rows(
  summary_table,
  tibble::tibble(
    Metric = "Number of variants with a Control note",
    Value  = num_control_variants
  )
)

# Multi-Sheet-Export: Summary + Liste der Control-Varianten
writexl::write_xlsx(
  list(
    summary = summary_table,
    control_variants = df_control
  ),
  "spliceai_summary.xlsx"
)

