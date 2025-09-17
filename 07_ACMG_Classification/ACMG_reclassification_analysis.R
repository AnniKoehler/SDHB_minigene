# =====================================================
# ACMG points + class changes (minimal, generic I/O)
# =====================================================
# Input:  Excel file "input_file.xlsx" in the working directory.
#         Must contain columns "without_mg" and "with_mg" holding ACMG/AMP codes.
# Process:
#   - Parse ACMG codes → sum points (pathogenic positive, benign negative).
#   - Map point totals to classes; BA1 → special "Benign (BA1)".
#   - Compute signed/absolute point deltas and class-change flags.
#   - Build summary/aggregation tables and write everything to Excel.
# Output: Excel file "output_file.xlsx" with sheets:
#         Summary, Class Switches, Same Classes, Same Points,
#         Variants Detail, Extremes max–min, Classes with MG,
#         (optional) Infinite Changes, and RNA strength.
# =====================================================

# Load all necessary libraries
library(readr)
library(readxl)
library(data.table)
library(writexl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(janitor)
library(lubridate)
library(stringr)
library(purrr)
library(tidyverse)

# Set working directory (optional)
# setwd("/path/to/your/folder")

# ---------------------------
# Load data
# ---------------------------
acmg_variants <- read_excel("input_file.xlsx") %>%
  clean_names() %>%
  mutate(
    without_mg = str_replace_all(without_mg, "[\r\n]{2,}", "\n"),
    with_mg    = str_replace_all(with_mg,    "[\r\n]{2,}", "\n")
  )

# ---------------------------
# Scoring function
# ---------------------------
score_acmg <- function(x) {
  if (is.na(x) || str_trim(x) == "") return(NA_real_)
  
  s <- x |>
    str_replace_all("[,;\\n\\r\\t]+", " ") |>
    str_squish()
  
  # BA1 standalone → treat as -Inf (special handling downstream)
  if (str_detect(s, "\\bBA1\\b")) return(-Inf)
  
  re <- "(BA1)|(PVS1|PS\\d*|PM\\d*|PP\\d*|BS\\d*|BP\\d*|BM)(?:[_ ]?(Very\\s+Strong|Strong|Moderate|Supporting))?"
  m  <- str_match_all(s, re)[[1]]
  if (nrow(m) == 0) return(NA_real_)
  
  df <- tibble(
    code     = m[,3],
    strength = m[,4]
  ) |>
    mutate(
      family = case_when(
        str_detect(code, "^PVS1$")            ~ "path",
        str_detect(code, "^PS\\d*$")          ~ "path",
        str_detect(code, "^PM\\d*$")          ~ "path",
        str_detect(code, "^PP\\d*$")          ~ "path",
        str_detect(code, "^BS\\d*$")          ~ "benign",
        str_detect(code, "^BP\\d*$")          ~ "benign",
        code == "BM"                          ~ "benign",
        TRUE                                  ~ NA_character_
      ),
      strength_norm = case_when(
        !is.na(strength)                      ~ strength,
        code == "PVS1"                        ~ "Very Strong",
        str_detect(code, "^PS\\d*$")          ~ "Strong",
        str_detect(code, "^PM\\d*$")          ~ "Moderate",
        str_detect(code, "^PP\\d*$")          ~ "Supporting",
        str_detect(code, "^BS\\d*$")          ~ "Strong",
        str_detect(code, "^BP\\d*$")          ~ "Supporting",
        code == "BM"                          ~ "Moderate",
        TRUE                                  ~ NA_character_
      )
    )
  
  pts <- df |>
    mutate(
      strength_norm = str_replace_all(strength_norm, "\\s+", " "),
      raw_points = case_when(
        family == "path"   & strength_norm == "Very Strong" ~  8,
        family == "path"   & strength_norm == "Strong"      ~  4,
        family == "path"   & strength_norm == "Moderate"    ~  2,
        family == "path"   & strength_norm == "Supporting"  ~  1,
        family == "benign" & strength_norm == "Strong"      ~ -4,
        family == "benign" & strength_norm == "Moderate"    ~ -2,
        family == "benign" & strength_norm == "Supporting"  ~ -1,
        TRUE ~ NA_real_
      )
    ) |>
    pull(raw_points)
  
  if (all(is.na(pts))) return(NA_real_)
  sum(pts, na.rm = TRUE)
}

# ---------------------------
# Class from points
# ---------------------------
class_from_points <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (is.infinite(p) && p < 0) return("Benign (BA1)")
  
  if (p >= 10)  return("Pathogenic")
  if (p >= 6)   return("Likely Pathogenic")
  if (p >= 1)   return("VUS")
  if (p <= -6)  return("Benign")
  if (p <= -1)  return("Likely Benign")
  return("VUS")
}

# ---------------------------
# Apply scoring & classification
# ---------------------------
acmg_variants <- acmg_variants |>
  mutate(
    points_without_mg = map_dbl(without_mg, score_acmg),
    points_with_mg    = map_dbl(with_mg,    score_acmg),
    class_without_mg  = vapply(points_without_mg, class_from_points, character(1)),
    class_with_mg     = vapply(points_with_mg,    class_from_points, character(1))
  )

# ---------------------------
# Changes (absolute and signed)
# ---------------------------
acmg_variants <- acmg_variants %>%
  mutate(
    point_change_signed = points_with_mg - points_without_mg,               # signed delta
    point_change_abs    = case_when(                                        # absolute delta
      is.na(point_change_signed) ~ NA_real_,
      is.infinite(points_with_mg) | is.infinite(points_without_mg) ~ Inf,
      TRUE ~ abs(point_change_signed)
    ),
    class_change        = if_else(class_with_mg != class_without_mg, TRUE, FALSE)
  )

# ---------------------------
# Aggregates & tables
# ---------------------------
n_class_switches <- sum(acmg_variants$class_change, na.rm = TRUE)
n_point_switches <- sum(!is.na(acmg_variants$point_change_abs) & acmg_variants$point_change_abs > 0)

class_switch_table <- acmg_variants %>%
  filter(class_change) %>%
  count(class_without_mg, class_with_mg, name = "n") %>%
  arrange(desc(n))

mean_abs_point_change <- acmg_variants %>%
  filter(is.finite(point_change_abs)) %>%
  summarise(mean = mean(point_change_abs, na.rm = TRUE)) %>%
  pull(mean)

same_class_table <- acmg_variants %>%
  filter(!is.na(class_without_mg), class_without_mg == class_with_mg) %>%
  count(class_without_mg, name = "n") %>%
  arrange(desc(n))

same_points_table <- acmg_variants %>%
  filter(!is.na(point_change_abs), point_change_abs == 0) %>%
  count(class_without_mg, name = "n") %>%
  arrange(desc(n))

max_abs_point_change <- suppressWarnings(max(acmg_variants$point_change_abs, na.rm = TRUE))
min_abs_point_change <- suppressWarnings(min(acmg_variants$point_change_abs, na.rm = TRUE))

min_abs_point_change_nonzero <- acmg_variants %>%
  filter(is.finite(point_change_abs), point_change_abs > 0) %>%
  summarise(min_abs_nz = ifelse(n() > 0, min(point_change_abs, na.rm = TRUE), NA_real_)) %>%
  pull(min_abs_nz)

extreme_max_df <- acmg_variants %>%
  filter(
    (!is.infinite(max_abs_point_change) & is.finite(point_change_abs) & point_change_abs == max_abs_point_change) |
      (is.infinite(max_abs_point_change)  & is.infinite(point_change_abs))
  )

extreme_min_df <- acmg_variants %>%
  filter(is.finite(point_change_abs), point_change_abs == min_abs_point_change)

inf_changes_df <- acmg_variants %>%
  filter(is.infinite(point_change_abs))

has_infinite_changes <- any(is.infinite(acmg_variants$point_change_abs), na.rm = TRUE)
n_infinite_changes   <- sum(is.infinite(acmg_variants$point_change_abs))

summary_df <- tibble::tibble(
  metric = c("n_class_switches",
             "n_point_switches_abs_gt_0",
             "mean_abs_point_change_finite",
             "max_abs_point_change_incl_inf",
             "min_abs_point_change_incl_0",
             "min_abs_point_change_nonzero_finite",
             "has_infinite_changes",
             "n_infinite_changes"),
  value  = c(n_class_switches,
             n_point_switches,
             mean_abs_point_change,
             max_abs_point_change,
             min_abs_point_change,
             min_abs_point_change_nonzero,
             has_infinite_changes,
             n_infinite_changes)
)

# ---------------------------
# New classes (with MG)
# ---------------------------
new_classes_list <- acmg_variants %>%
  distinct(class_with_mg) %>%
  arrange(class_with_mg)

new_classes_count <- acmg_variants %>%
  count(class_with_mg, name = "n") %>%
  arrange(desc(n))

# ---------------------------
# Excel export (safe sheet names)
# ---------------------------
wb <- createWorkbook()

addWorksheet(wb, "Summary")
addWorksheet(wb, "Class Switches")
addWorksheet(wb, "Same Classes")
addWorksheet(wb, "Same Points")
addWorksheet(wb, "Variants Detail")
addWorksheet(wb, "Extremes max–min")    # avoid "/" in sheet names
addWorksheet(wb, "Classes with MG")
if (n_infinite_changes > 0) addWorksheet(wb, "Infinite Changes")

writeData(wb, "Summary", summary_df)
writeData(wb, "Class Switches", class_switch_table)
writeData(wb, "Same Classes", same_class_table)
writeData(wb, "Same Points", same_points_table)
writeData(wb, "Variants Detail", acmg_variants)

# Extremes sheet
writeData(wb, "Extremes max–min", tibble::tibble(note = "MAX absolute change"))
writeData(wb, "Extremes max–min", extreme_max_df, startRow = 3)
start_row_min <- nrow(extreme_max_df) + 6
writeData(wb, "Extremes max–min", tibble::tibble(note = "MIN absolute change (incl. 0)"),
          startRow = start_row_min)
writeData(wb, "Extremes max–min", extreme_min_df, startRow = start_row_min + 2)

# New classes sheet
writeData(wb, "Classes with MG", tibble::tibble(note = "List of classes (with MG)"))
writeData(wb, "Classes with MG", new_classes_list, startRow = 3)
start_row_counts <- nrow(new_classes_list) + 6
writeData(wb, "Classes with MG", tibble::tibble(note = "Frequencies of classes (with MG)"),
          startRow = start_row_counts)
writeData(wb, "Classes with MG", new_classes_count, startRow = start_row_counts + 2)

if (n_infinite_changes > 0) {
  writeData(wb, "Infinite Changes", inf_changes_df)
}

saveWorkbook(wb, "output_file.xlsx", overwrite = TRUE)

# --- Count overall_rna_strength types (simple) ---
rna_counts <- acmg_variants %>%
  transmute(type_raw = str_squish(as.character(overall_rna_strength))) %>%
  mutate(
    type = case_when(
      is.na(type_raw) | type_raw == "" ~ "Missing",
      str_to_upper(type_raw) %in% c("NA", "N/A") ~ "Missing",
      TRUE ~ type_raw
    )
  ) %>%
  count(type, name = "n") %>%
  arrange(desc(n))



# --- Write to Excel (append / replace sheet) ---
out_xlsx <- "output_file.xlsx"

wb <- if (file.exists(out_xlsx)) openxlsx::loadWorkbook(out_xlsx) else openxlsx::createWorkbook()
sheet_name <- "RNA strength"

if (sheet_name %in% openxlsx::sheets(wb)) openxlsx::removeWorksheet(wb, sheet_name)
openxlsx::addWorksheet(wb, sheet_name)
openxlsx::writeData(wb, sheet_name, rna_counts, startRow = 1, startCol = 1)

openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
