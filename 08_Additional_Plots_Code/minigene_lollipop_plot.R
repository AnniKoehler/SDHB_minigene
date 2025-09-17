# =====================================================
# SDHB Lollipop-Style Schematic with Semi-Log Intron Compression
# =====================================================
# Input:  Excel file "input_file.xlsx" in the working directory, with columns:
#         - hgv_sg_hg38, hgv_sc, variant, genomic_position (and others)
# Process:
#   - Read variants, derive short HGVSc labels, and deduplicate by position.
#   - Define exons (5'→3'), compress introns (semi-log around the midpoint),
#     map genomic positions to compressed plot positions.
#   - Assign up/down y-offsets per region and draw stems + points + labels.
#   - Render exons as tiles on a baseline and flip x-axis (reverse).
# Output: "lollipop_plot.pdf" 
# =====================================================

# --- Load libraries ---
library(readr)
library(readxl)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(janitor)
library(lubridate)
library(stringr)
library(ggforce)
library(patchwork)
library(RColorBrewer)

# --- (Optional) set working directory ---
# setwd("/path/to/your/folder")

# --- Read & prepare data ---
df <- read_excel("input_file.xlsx") %>%
  clean_names() %>%
  rename(hgvsg_hg38 = hgv_sg_hg38, hgvsc = hgv_sc) %>%
  mutate(hgvsc_short = ifelse(variant == "MG_WT", "MG_WT", str_extract(hgvsc, "(?<=:).*")))

df_variants <- df %>%
  distinct(genomic_position, hgvsc_short) %>%
  arrange(genomic_position)

# --- Exon definitions (5' to 3') ---
exons <- data.frame(
  exon_start = c(17027749, 17028600, 17033060, 17044761),
  exon_end   = c(17027865, 17028736, 17033145, 17044888),
  exon_number = c(2, 3, 4, 5)
)

# --- Compute compressed positions ---
plot_pos <- 0
compressed_positions <- data.frame()
log_space <- 50

for (i in 1:nrow(exons)) {
  exon <- exons[i, ]
  exon_length <- exon$exon_end - exon$exon_start + 1
  exon_positions <- data.frame(
    genomic_position = exon$exon_start:exon$exon_end,
    plot_position = plot_pos + 1:exon_length
  )
  compressed_positions <- bind_rows(compressed_positions, exon_positions)
  plot_pos <- max(exon_positions$plot_position)
  
  if (i < nrow(exons)) {
    next_exon <- exons[i + 1, ]
    intron_start <- exon$exon_end + 1
    intron_end   <- next_exon$exon_start - 1
    intron_genomic <- intron_start:intron_end
    n <- length(intron_genomic)
    
    if (n > 1) {
      mid <- floor(n / 2)
      first_half  <- intron_genomic[1:mid]
      second_half <- intron_genomic[(mid + 1):n]
      
      df_left <- data.frame(
        genomic_position = first_half,
        scale = log1p(seq_along(first_half))
      )
      df_left$scale <- df_left$scale / max(df_left$scale)
      df_left$plot_position <- plot_pos + df_left$scale * log_space
      
      df_right <- data.frame(genomic_position = second_half)
      df_right$scale <- rev(log1p(seq_along(second_half)))
      df_right$scale <- df_right$scale / max(df_right$scale)
      df_right$plot_position <- plot_pos + log_space + (1 - df_right$scale) * log_space
      
      compressed_positions <- bind_rows(
        compressed_positions,
        df_left[, c("genomic_position", "plot_position")],
        df_right[, c("genomic_position", "plot_position")]
      )
      plot_pos <- max(df_right$plot_position)
    }
  }
}

# --- Assign region (Exon/Intron) ---
df_variants <- df_variants %>% mutate(region = NA_character_)
for (i in 1:nrow(exons)) {
  df_variants$region[df_variants$genomic_position %in% exons$exon_start[i]:exons$exon_end[i]] <- paste0("exon", exons$exon_number[i])
}
for (i in 1:(nrow(exons) - 1)) {
  intron_start <- exons$exon_end[i] + 1
  intron_end   <- exons$exon_start[i + 1] - 1
  if (intron_end > intron_start) {
    intron_mid <- floor((intron_start + intron_end) / 2)
    df_variants$region[df_variants$genomic_position %in% intron_start:intron_mid] <- paste0("intron", i, "_left")
    df_variants$region[df_variants$genomic_position %in% (intron_mid + 1):intron_end] <- paste0("intron", i, "_right")
  }
}

# --- Order & Y-values ---
region_order <- df_variants %>%
  filter(!is.na(region)) %>%
  group_by(region) %>%
  summarize(first_pos = min(genomic_position), .groups = "drop") %>%
  arrange(first_pos) %>%
  pull(region)

region_directions <- rep(c(1, -1), length.out = length(region_order))
names(region_directions) <- region_order

assign_y_values <- function(n, base = 0.3, step = 0.15) seq(from = base, by = step, length.out = n)

df_variants$y_start <- 0
df_variants$y_end <- NA
for (region_name in region_order) {
  idx <- which(df_variants$region == region_name)
  y_values <- assign_y_values(length(idx))
  df_variants$y_end[idx] <- y_values * region_directions[region_name]
}

# --- Map to compressed plot positions ---
df_variants_compressed <- df_variants %>%
  left_join(compressed_positions, by = "genomic_position") %>%
  filter(!is.na(plot_position))

# --- Exon coordinates + reverse order for labels ---
exons <- exons %>%
  rowwise() %>%
  mutate(
    xmin = min(compressed_positions$plot_position[compressed_positions$genomic_position %in% exon_start:exon_end]),
    xmax = max(compressed_positions$plot_position[compressed_positions$genomic_position %in% exon_start:exon_end])
  ) %>%
  arrange(desc(exon_start)) %>%
  ungroup()

# --- Exon labels with reversed numbering on display ---
exons$exon_label <- paste0("Exon ", rev(exons$exon_number))

# --- Build plot ---
p <- ggplot() +
  geom_segment(
    data = df_variants_compressed,
    aes(x = plot_position, xend = plot_position, y = y_start, yend = y_end),
    color = "grey"
  ) +
  geom_segment(
    data = df_variants_compressed,
    aes(x = plot_position, xend = plot_position + 5 * sign(y_end), y = y_end, yend = y_end),
    color = "grey"
  ) +
  geom_point(
    data = df_variants_compressed,
    aes(x = plot_position + 5 * sign(y_end), y = y_end),
    size = 7
  ) +
  geom_text(
    data = df_variants_compressed,
    aes(x = plot_position + 8 * sign(y_end), y = y_end, label = hgvsc_short),
    angle = 0,
    hjust = ifelse(sign(df_variants_compressed$y_end) == 1, 1, 0),
    size = 10
  ) +
  geom_rect(
    aes(
      xmin = min(compressed_positions$plot_position),
      xmax = max(compressed_positions$plot_position),
      ymin = -0.03, ymax = 0.03
    ),
    fill = "black", alpha = 0.3
  ) +
  geom_tile(
    data = exons,
    aes(x = (xmin + xmax)/2, y = 0, width = xmax - xmin, height = 0.2),
    fill = "#a8ddb5", color = "black", radius = 0.2
  ) +
  geom_text(
    data = exons,
    aes(x = (xmin + xmax)/2, y = 0, label = exon_label),
    size = 11, color = "black"
  ) +
  scale_fill_identity() +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    legend.position = "none",
    plot.margin = margin(2, 2, 2, 2, "cm")
  ) +
  coord_cartesian(clip = "off") +
  scale_x_reverse()

# --- Show and save ---
print(p)
ggsave("lollipop_plot.pdf", plot = p, device = "pdf", width = 45, height = 22, limitsize = FALSE)
