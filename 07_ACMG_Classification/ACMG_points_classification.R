# =====================================================
# ACMG Scoring Script
# =====================================================
# Input:  Excel file with two text columns (without_mg / with_mg)
#         containing ACMG/AMP criteria codes (e.g. "PVS1, PS3, BP7").
# Process: 
#   - Parse criteria codes and map them to point values 
#     (pathogenic evidence positive, benign evidence negative).
#   - Sum up points for each variant.
#   - Assign a classification based on thresholds:
#       Pathogenic, Likely Pathogenic, VUS, Likely Benign, Benign, or "Benign (BA1)".
# Output: Excel file with added columns:
#         points_without_mg, points_with_mg,
#         class_without_mg, class_with_mg
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
library(tidyverse)  # optional – includes many of the above packages as well

# (Optional) Working directory
# setwd("~/your/path/here")

# Read input table and normalize multi-line cells
acmg_variants <- read_excel("Input_file.xlsx") %>%
  clean_names() %>%
  mutate(
    without_mg = str_replace_all(without_mg, "[\r\n]{2,}", "\n"),
    with_mg    = str_replace_all(with_mg,    "[\r\n]{2,}", "\n")
  )

# ---------------------------
# Scoring function
# ---------------------------
score_acmg <- function(x) {
  if (is.na(x) || str_trim(x) == "") return(NA_real_)   # keep empty as NA
  
  s <- x |>
    str_replace_all("[,;\\n\\r\\t]+", " ") |>
    str_squish()
  
  if (str_detect(s, "\\bBA1\\b")) return(-Inf)  # BA1 special case
  
  re <- "(BA1)|(PVS1|PS\\d*|PM\\d*|PP\\d*|BS\\d*|BP\\d*|BM)(?:[_ ]?(Very\\s+Strong|Strong|Moderate|Supporting))?"
  m <- str_match_all(s, re)[[1]]
  if (nrow(m) == 0) return(NA_real_)
  
  df <- tibble(
    code = m[,3],
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
# Classification function
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
# Apply functions to dataframe
# ---------------------------
acmg_variants <- acmg_variants |>
  mutate(
    points_without_mg = map_dbl(without_mg, score_acmg),
    points_with_mg    = map_dbl(with_mg,    score_acmg),
    class_without_mg  = vapply(points_without_mg, class_from_points, character(1)),
    class_with_mg     = vapply(points_with_mg,    class_from_points, character(1))
  )

# Write output
write_xlsx(acmg_variants, "Output_file.xlsx")
