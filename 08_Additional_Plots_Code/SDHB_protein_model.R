# =====================================================
# SDHB protein linear map with exons and domains
# =====================================================
# Input:  none required (exon ranges and domains defined inline).
# Process:
#   - Draw a linear protein backbone (length = 280 aa).
#   - Overlay manually defined exon spans and domain boxes.
#   - Label exons; show a legend for domains only.
# Output: "protein_map.pdf" (vector PDF, ~15x2.5 inches).
# =====================================================

# --- Load libraries ---
library(readr)
library(readxl)
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(janitor)
library(ggplot2)
library(RColorBrewer)
library(rtracklayer)
library(R.utils)

# --- (Optional) set working directory ---
# setwd("/path/to/your/folder")

# --- Protein domains ---
domains <- data.frame(
  start  = c(40, 146, 176),
  end    = c(133, 218, 206),
  domain = c("2Fe-2S", "SDHAF1 interaction", "4Fe-4S")
)

# --- Protein length ---
protein_length <- 280

# --- Manually defined exon ranges on protein (aa) ---
exon_ranges <- data.frame(
  exon     = 1:8,
  aa_start = c(1, 26, 68, 97, 121, 159, 200, 248),
  aa_end   = c(25, 67, 96, 120, 158, 199, 247, 280)
) %>%
  mutate(
    exon_display_start = aa_start - 0.3,
    exon_display_end   = aa_end + 0.3,
    fill_color         = ifelse(exon %in% 2:5, "#a8ddb5", "#e0f3db")  # Exon 2–5 highlighted
  )

# --- Export to PDF ---
pdf("protein_map.pdf", width = 15, height = 2.5)

# --- Build plot ---
ggplot() +
  # Exons (manual fill_color vector; not mapped to legend)
  geom_rect(
    data = exon_ranges,
    aes(
      xmin = exon_display_start,
      xmax = exon_display_end,
      ymin = 0.42,
      ymax = 0.58
    ),
    fill = exon_ranges$fill_color
  ) +
  # Protein backbone
  geom_rect(aes(xmin = 0.01, xmax = protein_length, ymin = 0.48, ymax = 0.52), fill = "gray") +
  # Domains (categorical fill -> legend)
  geom_rect(
    data = domains,
    aes(xmin = start, xmax = end, ymin = 0.47, ymax = 0.53, fill = domain)
  ) +
  # Exon labels
  geom_text(
    data = exon_ranges,
    aes(x = (aa_start + aa_end) / 2, y = 0.6, label = paste0("Exon ", exon)),
    size = 3
  ) +
  # Domain colors
  scale_fill_manual(values = c(
    "2Fe-2S"             = "#b2182b",
    "SDHAF1 interaction" = "#fdae61",
    "4Fe-4S"             = "#3288bd"
  )) +
  # Layout
  theme_minimal() +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid   = element_blank(),
    legend.position = "right",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 10),
    axis.line.x     = element_line(color = "black")
  ) +
  labs(
    title = "SDHB protein",
    x = "Amino acid position",
    y = NULL
  ) +
  scale_x_continuous(breaks = seq(0, protein_length, by = 20))

# --- Close PDF device ---
dev.off()
