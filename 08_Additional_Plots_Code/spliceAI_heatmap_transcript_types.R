# =====================================================
# Combined transcript-type bar annotation + SpliceAI heatmap with region labels
# =====================================================
# Input:  Excel file "input_file.xlsx" containing (at least):
#         - hgvsg_hg38, hgvsc, variant, type, ratio, ds_al/ds_ag/ds_dl/ds_dg, genomic_position
# Process:
#   - Clean names, normalize transcript ratios per variant, exclude variants without functional assessment.
#   - Build SpliceAI heatmap (rows = DS_* channels, cols = variants).
#   - Top bar annotation shows transcript-type composition per variant.
#   - Bottom annotation labels genomic region (Exon/Intron) computed for minus strand.
# Output: "output_heatmap.pdf" (and on-screen plot)
# =====================================================

# --- Libraries ---
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(gridExtra)
library(readr)
library(readxl)
library(data.table)
library(ggplot2)
library(janitor)
library(lubridate)
library(tidyverse)

# --- (Optional) set working directory ---
# setwd("/path/to/your/folder")

# --- Read data ---
df <- read_excel("input_file.xlsx") %>%
  clean_names()

# ------------------------------------------------
# Prepare / recode
# ------------------------------------------------
df <- df %>%
  rename(hgvsg_hg38 = hgv_sg_hg38, hgvsc = hgv_sc) %>%
  mutate(
    hgvsc_short = ifelse(variant == "MG_WT", "MG_WT", str_extract(hgvsc, "(?<=:).*")),
    type = case_when(
      type == "no functional assesment3" ~ "no functional assessment³",
      type == "no functional assesment1" ~ "no functional assessment¹",
      type == "uncharacterized2"        ~ "uncharacterized²",
      TRUE ~ type
    ),
    ratio = as.numeric(ratio)
  )

# ------------------------------------------------
# Exclude variants that contain ONLY footnote types (¹ or ³)
# ------------------------------------------------
excluded_types <- c("no functional assessment¹", "no functional assessment³")

variants_to_exclude <- df %>%
  group_by(hgvsc_short) %>%
  filter(all(type %in% excluded_types)) %>%
  distinct(hgvsc_short) %>%
  pull(hgvsc_short)

df_filtered <- df %>%
  filter(!hgvsc_short %in% variants_to_exclude)

# ------------------------------------------------
# Normalize transcript ratios within variant
# ------------------------------------------------
df_filtered <- df_filtered %>%
  group_by(hgvsc_short) %>%
  mutate(ratio = ratio / sum(ratio, na.rm = TRUE)) %>%
  ungroup()

# ------------------------------------------------
# Heatmap matrix (SpliceAI scores)
# ------------------------------------------------
heatmap_matrix <- df_filtered %>%
  distinct(hgvsc_short, ds_al, ds_ag, ds_dl, ds_dg) %>%
  column_to_rownames("hgvsc_short") %>%
  t()

# ------------------------------------------------
# Column order: MG_WT first, then GENOMIC DESC; NAs to the right
# ------------------------------------------------
var_pos <- df_filtered %>%
  distinct(hgvsc_short, genomic_position)

order_tbl <- var_pos %>%
  mutate(
    is_mgwt = if_else(hgvsc_short == "MG_WT", 0L, 1L),
    gp_sort = if_else(is.na(genomic_position), Inf, genomic_position)
  ) %>%
  arrange(is_mgwt, desc(gp_sort))

col_order <- order_tbl$hgvsc_short
col_order <- intersect(col_order, colnames(heatmap_matrix))  # keep only existing columns, preserve order
heatmap_matrix <- heatmap_matrix[, col_order, drop = FALSE]

# ------------------------------------------------
# Top barplot matrix (transcript types per variant)
# ------------------------------------------------
bar_data <- df_filtered %>%
  group_by(hgvsc_short, type) %>%
  summarise(total_ratio = sum(ratio, na.rm = TRUE), .groups = "drop")

bar_wide <- tidyr::pivot_wider(
  bar_data,
  names_from = type,
  values_from = total_ratio,
  values_fill = 0
)

bar_mat <- as.matrix(bar_wide[, -1, drop = FALSE])
rownames(bar_mat) <- bar_wide$hgvsc_short

# Fill missing rows (0) so ordering matches the heatmap
missing_rows <- setdiff(col_order, rownames(bar_mat))
if (length(missing_rows) > 0) {
  zero_block <- matrix(0, nrow = length(missing_rows), ncol = ncol(bar_mat))
  rownames(zero_block) <- missing_rows
  colnames(zero_block) <- colnames(bar_mat)
  bar_mat <- rbind(bar_mat, zero_block)
}
bar_mat <- bar_mat[col_order, , drop = FALSE]

# Colors for transcript types (columns of bar_mat)
types <- colnames(bar_mat)
if (length(types) <= 12) {
  palette_colors <- RColorBrewer::brewer.pal(n = max(3, length(types)), name = "Set3")[seq_along(types)]
} else {
  palette_colors <- scales::hue_pal()(length(types))
}
names(palette_colors) <- types
bar_colors <- palette_colors[types]

# Top annotation (barplot)
top_anno <- HeatmapAnnotation(
  TranscriptTypes = anno_barplot(bar_mat, gp = gpar(fill = bar_colors)),
  annotation_name_side = "left",
  annotation_height = unit(8, "cm")
)

# Legend for transcript types (right)
type_legend <- Legend(
  title = "Transcript Type",
  labels = names(bar_colors),
  legend_gp = gpar(fill = bar_colors),
  title_gp = gpar(fontface = "bold"),
  labels_gp = gpar(fontsize = 9),
  direction = "vertical",
  ncol = 1,
  gap = unit(4, "mm")
)

# ------------------------------------------------
# Region annotation (Exon/Intron) from genomic positions (minus strand logic)
# ------------------------------------------------
exons <- data.frame(
  exon  = 1:8,
  start = c(17054032, 17044888, 17033145, 17028736, 17027865, 17024074, 17022730, 17018958),
  end   = c(17053948, 17044761, 17033060, 17028600, 17027749, 17023973, 17022608, 17018722)
)

get_region <- function(pos, exons) {
  if (is.na(pos)) return(NA_character_)
  # Exon (minus strand: start > end)
  for (i in seq_len(nrow(exons))) {
    s <- exons$start[i]; e <- exons$end[i]
    if (!is.na(s) && !is.na(e) && pos <= s && pos >= e) {
      return(paste0("Exon ", exons$exon[i]))
    }
  }
  # Intron (between two exons)
  for (i in seq_len(nrow(exons) - 1)) {
    s <- exons$start[i]
    e_next <- exons$end[i + 1]
    if (!is.na(s) && !is.na(e_next) && pos < s && pos > e_next) {
      return(paste0("Intron ", exons$exon[i]))
    }
  }
  return(NA_character_)
}

df_var_pos <- df_filtered %>%
  distinct(hgvsc_short, genomic_position) %>%
  mutate(region = sapply(genomic_position, get_region, exons = exons))

# Vector of regions in column order
region_vector <- df_var_pos$region[match(col_order, df_var_pos$hgvsc_short)]
region_vector <- factor(region_vector, levels = unique(na.omit(region_vector)))

# Colors for regions
region_colors <- RColorBrewer::brewer.pal(max(3, nlevels(region_vector)), "Pastel2")[seq_len(nlevels(region_vector))]
names(region_colors) <- levels(region_vector)

region_annotation <- HeatmapAnnotation(
  Region = region_vector,
  col = list(Region = region_colors),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  annotation_legend_param = list(title = "Genomic Region")
)

# ------------------------------------------------
# Heatmap
# ------------------------------------------------
col_fun <- colorRamp2(c(0, 0.5, 1), c("white", "#41b6c4", "#225ea8"))

ht_combined <- Heatmap(
  heatmap_matrix,
  name = "SpliceAI",
  top_annotation = top_anno,
  bottom_annotation = region_annotation,
  col = col_fun,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  column_names_rot = 90,
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  width = unit(ncol(heatmap_matrix) * 5, "mm"),
  height = unit(nrow(heatmap_matrix) * 4, "mm"),
  column_title = "Variant-wise Transcript & SpliceAI Summary",
  heatmap_legend_param = list(
    title = "Score",
    at = c(0, 0.5, 1),
    labels = c("0", "0.5", "1.0")
  )
)

# Draw to screen (with separate legend for types)
draw(ht_combined, annotation_legend_list = list(type_legend), annotation_legend_side = "right")

# ------------------------------------------------
# Export to PDF
# ------------------------------------------------
pdf("output_heatmap.pdf", width = 15, height = 8)
draw(ht_combined, annotation_legend_list = list(type_legend), annotation_legend_side = "right")
dev.off()
