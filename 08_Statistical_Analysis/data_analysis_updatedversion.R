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
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Doktorarbeit/Analysis/transcript_analysis_type plot_spliceAI ")

# ----------------------------
# Daten laden (rohe Basis)
# ----------------------------
df_all <- read_excel("2025_08_18_SDHB MG transcripts shorter_ACMG_short.xlsx") %>%
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

# ----------------------------
# Kennzahlen VOR dem Filtern
# ----------------------------
# Gesamtzahl aller Transkripte (pre-filter)
total_transcripts_all <- nrow(df_all)

# No functional assessment (pre-filter)
no_func_transcripts_all <- df_all %>%
  filter(type %in% c("no functional assesment1", "no functional assesment3"))
num_no_functional_transcripts <- nrow(no_func_transcripts_all)
num_no_functional_variants    <- no_func_transcripts_all %>% distinct(variant) %>% nrow()

# Uncharacterized (pre-filter)
uncharacterized_transcripts_all <- df_all %>%
  filter(type %in% c("uncharacterized2"))
num_uncharacterized_transcripts <- nrow(uncharacterized_transcripts_all)
num_uncharacterized_variants    <- uncharacterized_transcripts_all %>% distinct(variant) %>% nrow()

# ----------------------------
# Filtern für die eigentliche Analyse
# ----------------------------
df_variants <- df_all %>%
  filter(
    variant != "MG_WT",
    !type %in% c("no functional assesment3", "no functional assesment1", "uncharacterized2")
  )

# ----------------------------
# Kennzahlen NACH dem Filtern
# ----------------------------

# Transkripte pro Variante
variant_counts <- df_variants %>%
  group_by(variant) %>%
  summarise(number_transcripts = n(), .groups = "drop")

# Varianten mit ≥2 Transkripten
variants_with_two_or_more <- variant_counts %>%
  filter(number_transcripts >= 2)
num_variants_two_or_more <- nrow(variants_with_two_or_more)

# Max und Median (nach Filter)
max_transcript_number    <- max(variant_counts$number_transcripts, na.rm = TRUE)
median_transcript_number <- median(variant_counts$number_transcripts, na.rm = TRUE)
mean_transcript_number <- mean(variant_counts$number_transcripts, na.rm = TRUE)

# Neue (nicht-kanonische) Splice-Site (nach Filter)
transcripts_with_new_splice <- df_variants %>%
  filter(!is.na(non_cannonical_splice_motif) & non_cannonical_splice_motif != "")
variants_with_new_splice    <- transcripts_with_new_splice %>% distinct(variant)
num_variants_with_new_splice <- nrow(variants_with_new_splice)

# Aberrantes Splicing (nach Filter)
aberrant_transcripts <- df_variants %>%
  filter(aberrant_splicing == 1)
num_aberrant_transcripts <- nrow(aberrant_transcripts)

# Nicht aberrant (nach Filter)
nonaberrant_transcripts <- df_variants %>%
  filter(aberrant_splicing == 0)
num_nonaberrant_transcripts <- nrow(nonaberrant_transcripts)

# Varianten mit aberrantem Splicing (distinct)
aberrant_variants <- aberrant_transcripts %>% distinct(variant)

# Varianten mit nur einem (vollen) Transkript (kein aberrantes)
variant_one_count <- df_variants %>%
  group_by(variant) %>%
  summarise(transcript_count = n(), .groups = "drop")

variants_with_one_transcript <- variant_one_count %>%
  filter(transcript_count == 1)

single_transcripts <- df_variants %>%
  semi_join(variants_with_one_transcript, by = "variant")

full_length_single_transcripts <- single_transcripts %>%
  filter(aberrant_splicing == 0)

# Anzahl Varianten, die nur ein Full-length-Transkript haben
num_full_length_single_variants <- nrow(full_length_single_transcripts)

# "Leaky" Varianten: sowohl aberrant als auch non-aberrant vorhanden
leaky_variants <- df_variants %>%
  group_by(variant) %>%
  summarise(
    has_aberrant     = any(aberrant_splicing == 1, na.rm = TRUE),
    has_non_aberrant = any(aberrant_splicing == 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(has_aberrant & has_non_aberrant)

num_leaky_variants <- nrow(leaky_variants)

# Splice-Effekt Kategorien (Transkript-Ebene, nach Filter)
intron_retention_transcripts <- df_variants %>% filter(splice_effect == "intron retention")
num_intron_retention_transcripts <- nrow(intron_retention_transcripts)

multiexon_skipping_transcripts <- df_variants %>% filter(splice_effect == "multiexon skipping")
num_multiexon_skipping_transcripts <- nrow(multiexon_skipping_transcripts)

single_exon_skipping_transcripts <- df_variants %>% filter(splice_effect == "single exon skipping")
num_single_exon_skipping_transcripts <- nrow(single_exon_skipping_transcripts)

partial_exon_skipping_transcripts <- df_variants %>% filter(splice_effect == "partial exon skipping")
num_partial_exon_skipping_transcripts <- nrow(partial_exon_skipping_transcripts)

pseudo_exon_transcripts <- df_variants %>% filter(splice_effect == "pseudo exon")
num_pseudo_exon_transcripts <- nrow(pseudo_exon_transcripts)

# Splice-Effekt Kategorien (Varianten-Ebene, nach Filter)
num_intron_retention_variants   <- intron_retention_transcripts %>% distinct(variant) %>% nrow()
num_multiexon_skipping_variants <- multiexon_skipping_transcripts %>% distinct(variant) %>% nrow()
num_single_exon_skipping_variants <- single_exon_skipping_transcripts %>% distinct(variant) %>% nrow()
num_partial_exon_skipping_variants <- partial_exon_skipping_transcripts %>% distinct(variant) %>% nrow()
num_pseudo_exon_variants        <- pseudo_exon_transcripts %>% distinct(variant) %>% nrow()

# Transcript consequences (Transkript-Ebene, nach Filter)
full_length_transcripts <- df_variants %>% filter(type == "full length")
num_full_length_transcripts <- nrow(full_length_transcripts)

PTC_transcripts <- df_variants %>% filter(type == "PTC")
num_PTC_transcripts <- nrow(PTC_transcripts)

in_frame_transcripts <- df_variants %>% filter(type == "in-frame")
num_in_frame_transcripts <- nrow(in_frame_transcripts)

# Transcript consequences (Varianten-Ebene, nach Filter)
num_full_length_variants <- full_length_transcripts %>% distinct(variant) %>% nrow()
num_PTC_variants         <- PTC_transcripts %>% distinct(variant) %>% nrow()
num_in_frame_variants    <- in_frame_transcripts %>% distinct(variant) %>% nrow()

# ----------------------------
# Summary-Tabelle
# ----------------------------
summary_table <- tibble::tibble(
  Metric = c(
    "Total transcripts (pre-filter)",
    "Variants with ≥2 transcripts",
    "Max transcripts per variant",
    "Median transcripts per variant",
    "Avergage transcripts per variant",
    "Variants with new splice site",
    "Transcripts with aberrant splicing",
    "Transcripts without aberrant splicing",
    "Variants with aberrant splicing",
    "Transcripts with no functional assessment (pre-filter)",
    "Variants with no functional assessment (pre-filter)",
    "Variants with only full-length transcript",
    "Leaky variants (WT + aberrant)",
    "Intron retention (transcripts)",
    "Multiexon skipping (transcripts)",
    "Single exon skipping (transcripts)",
    "Partial exon skipping (transcripts)",
    "Pseudo exon (transcripts)",
    "Intron retention (variants)",
    "Multiexon skipping (variants)",
    "Single exon skipping (variants)",
    "Partial exon skipping (variants)",
    "Pseudo exon (variants)",
    "Full-length transcripts",
    "PTC transcripts",
    "In-frame transcripts",
    "Uncharacterized transcripts (pre-filter)",
    "Full-length variants",
    "PTC variants",
    "In-frame variants",
    "Uncharacterized variants (pre-filter)"
  ),
  Value = c(
    total_transcripts_all,          # pre-filter
    num_variants_two_or_more,
    max_transcript_number,
    median_transcript_number,
    mean_transcript_number,
    num_variants_with_new_splice,
    num_aberrant_transcripts,
    num_nonaberrant_transcripts,
    nrow(aberrant_variants),
    num_no_functional_transcripts,  # pre-filter
    num_no_functional_variants,     # pre-filter
    num_full_length_single_variants,
    num_leaky_variants,
    num_intron_retention_transcripts,
    num_multiexon_skipping_transcripts,
    num_single_exon_skipping_transcripts,
    num_partial_exon_skipping_transcripts,
    num_pseudo_exon_transcripts,
    num_intron_retention_variants,
    num_multiexon_skipping_variants,
    num_single_exon_skipping_variants,
    num_partial_exon_skipping_variants,
    num_pseudo_exon_variants,
    num_full_length_transcripts,
    num_PTC_transcripts,
    num_in_frame_transcripts,
    num_uncharacterized_transcripts,  # pre-filter
    num_full_length_variants,
    num_PTC_variants,
    num_in_frame_variants,
    num_uncharacterized_variants      # pre-filter
  )
)

writexl::write_xlsx(summary_table, "summary_table_updated_new.xlsx")
print(summary_table)




#Splice AI analysis 

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
df_filtered <- df_variants %>%
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

# Write the summary to an Excel file named "spliceai_summary.xlsx"
write_xlsx(summary_table, "spliceai_summary.xlsx")





## plots 

# =============================
# Distribution of splice effects
# =============================

library(dplyr)
library(ggplot2)
library(scales)

# Transkript-Ebene: Anteil je Splice-Effekt
effects_transcript <- df_variants %>%
  mutate(splice_effect = tidyr::replace_na(splice_effect, "NA"),
         splice_effect = stringr::str_squish(splice_effect)) %>%
  count(splice_effect, name = "n") %>%
  mutate(percent = n / sum(n)) %>%
  arrange(desc(percent)) %>%
  mutate(
    splice_effect  = factor(splice_effect, levels = splice_effect),
    percent_label  = percent_format(accuracy = 0.1)(percent)
  )

p_effects_percent <- ggplot(effects_transcript,
                            aes(x = splice_effect, y = percent, fill = splice_effect)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = percent_label), hjust = -0.15, size = 3) +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Distribution of splice effects within assessed transcripts",
       x = "Splice effect", y = "Proportion") +
  theme_classic() +
  theme(panel.grid = element_blank()) +
  ggplot2::scale_fill_brewer(palette = "Set2") +
  expand_limits(y = max(effects_transcript$percent) * 1.12)

# Export
ggsave("splice_effects_transcript_percent.png", plot = p_effects_percent, width = 8, height = 6, dpi = 300)
ggsave("splice_effects_transcript_percent.pdf",  plot = p_effects_percent, width = 8, height = 6)
