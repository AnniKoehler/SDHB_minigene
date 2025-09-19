# =====================================================
# Analysis of ACMG classification per transcript
# =====================================================
# Input:  Excel file "input_file.xlsx" with columns including:
#         - variant, type, splice_effect, transcript_acmg_code
# Process:
#   - Filter out MG_WT and unassessed/uncharacterized rows
#   - Order ACMG codes by frequency
#   - Summarize counts and percentages
#   - Plot stacked proportions of splice effects per ACMG code
#     (with special plotmath label for "PVS1_N/A5")
# Output:
#   - "transcript_summary.xlsx" (counts + % per ACMG code)
#   - "splice_plot.png" and "splice_plot.pdf"
# =====================================================

# Packages
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(ggplot2)
library(scales)
library(janitor)
library(writexl)
library(patchwork)

# ---- Load & filter data ----
acmg_transcripts <- read_excel("input_file.xlsx") %>%
  clean_names() %>%
  filter(
    variant != "MG_WT",
    !type %in% c("no functional assesment3", "no functional assesment1", "uncharacterized2")
  )

# ---- Plot labels: only PVS1_N/A5 uses plotmath (N/A^5) ----
lab_plotmath <- function(x) {
  ifelse(x == "PVS1_N/A5",
         'PVS1~"N/A"^5',
         sprintf('"%s"', x))
}

# ---- Order categories by frequency (ascending, as in original) ----
order_levels <- acmg_transcripts %>%
  count(transcript_acmg_code, name = "n") %>%
  arrange(n) %>%
  pull(transcript_acmg_code)

# ---- Prepare data ----
acmg_transcripts <- acmg_transcripts %>%
  mutate(
    splice_effect = replace_na(splice_effect, "NA"),
    transcript_acmg_code = factor(transcript_acmg_code, levels = order_levels),
    transcript_acmg_code_lab = lab_plotmath(as.character(transcript_acmg_code))
  )

# ---- Overview table (counts + %) ----
overview <- acmg_transcripts %>%
  count(transcript_acmg_code, name = "n") %>%
  mutate(percent = round(100 * n / sum(n), 1)) %>%
  arrange(desc(n))

# Export overview
write_xlsx(overview, "transcript_summary.xlsx")

# ---- Plot (stacked proportions per ACMG code) ----
gg <- acmg_transcripts %>%
  count(transcript_acmg_code, transcript_acmg_code_lab, splice_effect) %>%
  group_by(transcript_acmg_code, transcript_acmg_code_lab) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = transcript_acmg_code_lab, y = prop, fill = splice_effect)) +
  geom_col(width = 0.6) +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Transcript ACMG Code",
    y = "Proportion",
    fill = "Transcript splice effect",
    title = "Distribution of transcript splice effects across ACMG evidence codes"
  ) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        legend.position = "right") +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  ggplot2::scale_fill_brewer(palette = "Set2")
# Colorblind-safe alternative:
# + ggplot2::scale_fill_viridis_d(option = "D")

# ---- Footnote under the plot (superscript 5 via plotmath) ----
gg_with_caption <- gg + plot_annotation(
  caption = expression({}^5*"The transcript corresponds to the alternative transcript ENST00000485515.6")
)

# ---- Export ----
ggsave("splice_plot.png", plot = gg_with_caption, width = 10, height = 6, dpi = 300)
ggsave("splice_plot.pdf", plot = gg_with_caption, width = 10, height = 6)
