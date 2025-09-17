# ============================================
# Transcript Results Analysis 
# ============================================
# INPUTS:  - Excel table with columns: type, variant, hgv_sc, hgv_sg_hg38, splice_effect, aberrant_splicing
#   (optional: non_cannonical_splice_motif). Set INPUT_FILE/OUTPUT_XLSX in script.
# PROCESS: 
# - Clean/rename columns; compute pre-filter counts.
# - Exclude MG_WT and non-assessable rows.
# - Per-variant transcript counts (max/median/mean).
# - Flag new splice sites; count aberrant vs non-aberrant; detect leaky variants.
# - Summarize splice-effect & consequence categories (transcript + variant level)
# OUTPUTS
# - Excel summary table (OUTPUT_XLSX) 
# - Plot of splice-effect proportions: splice_effects_transcript_percent.(png|pdf)
# ============================================
# ============================
# Load all necessary libraries
# ============================
library(readr)
library(readxl)
library(data.table)
library(dplyr)
library(tidyr)
library(janitor)
library(stringr)
library(writexl)
library(ggplot2)
library(scales)

# ============================
# User-configurable parameters
# ============================

# (Optional) Working directory
# setwd("~/your/path/here")

# Input file 
INPUT_FILE <- "input_file.xlsx"

# Column names (only change if your file uses different names)
COL_TYPE        <- "type"
COL_VARIANT     <- "variant"
COL_HGVSC       <- "hgv_sc"               # will be renamed to hgvsc
COL_HGVSG_HG38  <- "hgv_sg_hg38"          # will be renamed to hgvsg_hg38
COL_EFFECT      <- "splice_effect"
COL_ABERRANT    <- "aberrant_splicing"
COL_NONCAN_MOT  <- "non_cannonical_splice_motif"

# Category labels (adjust to your data if needed)
NO_FUNC_TYPES <- c("no functional assesment1", "no functional assesment3")
UNCHAR_TYPES  <- c("uncharacterized2")

TYPE_FULL_LENGTH <- "full length"
TYPE_PTC         <- "PTC"
TYPE_IN_FRAME    <- "in-frame"

# Aberrant coding (set to 1/0, TRUE/FALSE, "yes"/"no", etc., as in your data)
ABERRANT_TRUE  <- 1
ABERRANT_FALSE <- 0

# Output file
OUTPUT_XLSX <- "output_file.xlsx"

# ============================
# Data load (raw base)
# ============================
df_all <- read_excel(INPUT_FILE) %>%
  clean_names() %>%
  rename(
    hgvsg_hg38 = all_of(COL_HGVSG_HG38),
    hgvsc      = all_of(COL_HGVSC),
    type       = all_of(COL_TYPE),
    splice_effect = all_of(COL_EFFECT),
    aberrant_splicing = all_of(COL_ABERRANT)
  ) %>%
  mutate(
    # short HGVSc, keep literal "MG_WT" if present
    hgvsc_short = ifelse(
      !!sym(COL_VARIANT) == "MG_WT",
      "MG_WT",
      stringr::str_extract(hgvsc, "(?<=:).*")
    ),
    variant = .data[[COL_VARIANT]],
    non_cannonical_splice_motif = if (COL_NONCAN_MOT %in% names(.)) .data[[COL_NONCAN_MOT]] else NA
  )

# ============================
# Metrics BEFORE filtering
# ============================
# Total transcripts (pre-filter)
total_transcripts_all <- nrow(df_all)

# No functional assessment (pre-filter)
no_func_transcripts_all <- df_all %>%
  filter(type %in% NO_FUNC_TYPES)
num_no_functional_transcripts <- nrow(no_func_transcripts_all)
num_no_functional_variants    <- no_func_transcripts_all %>% distinct(variant) %>% nrow()

# Uncharacterized (pre-filter)
uncharacterized_transcripts_all <- df_all %>%
  filter(type %in% UNCHAR_TYPES)
num_uncharacterized_transcripts <- nrow(uncharacterized_transcripts_all)
num_uncharacterized_variants    <- uncharacterized_transcripts_all %>% distinct(variant) %>% nrow()

# ============================
# Filter for main analysis
# ============================
df_variants <- df_all %>%
  filter(
    variant != "MG_WT",
    !type %in% c(NO_FUNC_TYPES, UNCHAR_TYPES)
  )

# ============================
# Metrics AFTER filtering
# ============================

# Transcripts per variant
variant_counts <- df_variants %>%
  group_by(variant) %>%
  summarise(number_transcripts = n(), .groups = "drop")

# Variants with ≥2 transcripts
variants_with_two_or_more <- variant_counts %>%
  filter(number_transcripts >= 2)
num_variants_two_or_more <- nrow(variants_with_two_or_more)

# Max / median / mean transcripts per variant
max_transcript_number    <- max(variant_counts$number_transcripts, na.rm = TRUE)
median_transcript_number <- median(variant_counts$number_transcripts, na.rm = TRUE)
mean_transcript_number   <- mean(variant_counts$number_transcripts, na.rm = TRUE)

# New (non-canonical) splice site (post-filter)
transcripts_with_new_splice <- df_variants %>%
  filter(!is.na(non_cannonical_splice_motif) & non_cannonical_splice_motif != "")
variants_with_new_splice     <- transcripts_with_new_splice %>% distinct(variant)
num_variants_with_new_splice <- nrow(variants_with_new_splice)

# Aberrant splicing (post-filter)
aberrant_transcripts <- df_variants %>%
  filter(aberrant_splicing == ABERRANT_TRUE)
num_aberrant_transcripts <- nrow(aberrant_transcripts)

# Non-aberrant splicing (post-filter)
nonaberrant_transcripts <- df_variants %>%
  filter(aberrant_splicing == ABERRANT_FALSE)
num_nonaberrant_transcripts <- nrow(nonaberrant_transcripts)

# Variants with aberrant splicing (distinct)
aberrant_variants <- aberrant_transcripts %>% distinct(variant)

# Variants with exactly one transcript
variant_one_count <- df_variants %>%
  group_by(variant) %>%
  summarise(transcript_count = n(), .groups = "drop")

variants_with_one_transcript <- variant_one_count %>%
  filter(transcript_count == 1)

single_transcripts <- df_variants %>%
  semi_join(variants_with_one_transcript, by = "variant")

# Among those, full-length only (no aberrant)
full_length_single_transcripts <- single_transcripts %>%
  filter(aberrant_splicing == ABERRANT_FALSE)

# Number of variants that have only a full-length transcript
num_full_length_single_variants <- nrow(full_length_single_transcripts)

# “Leaky” variants: both aberrant and non-aberrant present
leaky_variants <- df_variants %>%
  group_by(variant) %>%
  summarise(
    has_aberrant     = any(aberrant_splicing == ABERRANT_TRUE,  na.rm = TRUE),
    has_non_aberrant = any(aberrant_splicing == ABERRANT_FALSE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(has_aberrant & has_non_aberrant)

num_leaky_variants <- nrow(leaky_variants)

# Splice-effect categories (transcript level, post-filter)
intron_retention_transcripts     <- df_variants %>% filter(splice_effect == "intron retention")
num_intron_retention_transcripts <- nrow(intron_retention_transcripts)

multiexon_skipping_transcripts     <- df_variants %>% filter(splice_effect == "multiexon skipping")
num_multiexon_skipping_transcripts <- nrow(multiexon_skipping_transcripts)

single_exon_skipping_transcripts     <- df_variants %>% filter(splice_effect == "single exon skipping")
num_single_exon_skipping_transcripts <- nrow(single_exon_skipping_transcripts)

partial_exon_skipping_transcripts     <- df_variants %>% filter(splice_effect == "partial exon skipping")
num_partial_exon_skipping_transcripts <- nrow(partial_exon_skipping_transcripts)

pseudo_exon_transcripts     <- df_variants %>% filter(splice_effect == "pseudo exon")
num_pseudo_exon_transcripts <- nrow(pseudo_exon_transcripts)

# Splice-effect categories (variant level, post-filter)
num_intron_retention_variants     <- intron_retention_transcripts     %>% distinct(variant) %>% nrow()
num_multiexon_skipping_variants   <- multiexon_skipping_transcripts   %>% distinct(variant) %>% nrow()
num_single_exon_skipping_variants <- single_exon_skipping_transcripts %>% distinct(variant) %>% nrow()
num_partial_exon_skipping_variants<- partial_exon_skipping_transcripts%>% distinct(variant) %>% nrow()
num_pseudo_exon_variants          <- pseudo_exon_transcripts          %>% distinct(variant) %>% nrow()

# Transcript consequences (transcript level, post-filter)
full_length_transcripts <- df_variants %>% filter(type == TYPE_FULL_LENGTH)
num_full_length_transcripts <- nrow(full_length_transcripts)

PTC_transcripts <- df_variants %>% filter(type == TYPE_PTC)
num_PTC_transcripts <- nrow(PTC_transcripts)

in_frame_transcripts <- df_variants %>% filter(type == TYPE_IN_FRAME)
num_in_frame_transcripts <- nrow(in_frame_transcripts)

# Transcript consequences (variant level, post-filter)
num_full_length_variants <- full_length_transcripts %>% distinct(variant) %>% nrow()
num_PTC_variants         <- PTC_transcripts         %>% distinct(variant) %>% nrow()
num_in_frame_variants    <- in_frame_transcripts    %>% distinct(variant) %>% nrow()

# ============================
# Summary table (Excel + print)
# ============================
summary_table <- tibble::tibble(
  Metric = c(
    "Total transcripts (pre-filter)",
    "Variants with ≥2 transcripts",
    "Max transcripts per variant",
    "Median transcripts per variant",
    "Average transcripts per variant",
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

writexl::write_xlsx(summary_table, OUTPUT_XLSX)
print(summary_table)


# ============================
# Plots (optional)
# ============================
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

ggsave("splice_effects_transcript_percent.png", plot = p_effects_percent, width = 8, height = 6, dpi = 300)
ggsave("splice_effects_transcript_percent.pdf",  plot = p_effects_percent, width = 8, height = 6)

