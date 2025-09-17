# ===============================================================
# Transcript pooling, Variant classification & ACMG Table Builder
# ===============================================================
# Input:    - main variant sheet
#           - sheet that provides HGVSp
# Process:  - Loads a curated variant sheet with RNA/protein HGVS and per-variant ratios
#           - Builds a compact per-variant info string (e.g., "0.42 (r.c.XXX, p.YYY)")
#           - Optionally tags special cases as FL_WT or N/A
#           - Appends Unicode superscripts from transcript_acmg_code (e.g., N/A¹²)
#           - Normalizes ACMG codes (collapsing any N/A* to plain "N/A")
#           - Pools per-variant evidence into a wide ACMG table (PVS1, PVS1_S, PVS1_M, BP7_S, N/A)
#           - Computes a combined strength string (e.g., "PVS1 (0.32) + BP7_S (1.00)")
#           - Assigns overall ACMG strength (for clear cases)
#           - Joins a second sheet to add HGVSp and writes the final Excel table
# Output:   - "ACMG_table_extended.xlsx"  (final, wide table)
# ===============================================================

# -------------------------------
# Load libraries
# -------------------------------
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

# -------------------------------
# Working directory (edit to your path)
# -------------------------------
# setwd("~/")

# -------------------------------
# Helpers: superscripts & parsing
# -------------------------------

# Map digits to Unicode superscripts
to_superscript <- function(x) {
  nums <- c("0"="⁰", "1"="¹", "2"="²", "3"="³", "4"="⁴", "5"="⁵",
            "6"="⁶", "7"="⁷", "8"="⁸", "9"="⁹", ","=",")
  chars <- unlist(strsplit(x, ""))
  paste0(map_chr(chars, ~ ifelse(is.na(nums[.x]), .x, nums[.x])), collapse = "")
}

# Extract e.g. "N/A12,34" → "¹²³⁴" from transcript_acmg_code
extract_na_superscripts <- function(code_column) {
  matches <- stringr::str_extract_all(code_column, "N/A[0-9,]+")[[1]]
  if (length(matches) == 0) return("")
  numbers <- unlist(stringr::str_extract_all(matches, "\\d"))
  if (length(numbers) == 0) return("")
  to_superscript(paste0(numbers, collapse = ""))
}

# Sum all numeric ratio tokens in a text like: "0.32 (…); 0.58 (…)"
sum_ratios_from_text <- function(text) {
  nums <- stringr::str_extract_all(text, "\\d+[,.]\\d+")[[1]]
  if (length(nums) == 0) return(0)
  sum(as.numeric(stringr::str_replace(nums, ",", ".")))
}

# Build combined strength: e.g., "PVS1 (0.32) + BP7_S (1.00)"
build_combined_strength <- function(row, codes) {
  parts <- c()
  for (code in codes) {
    val <- row[[code]]
    if (!is.na(val) && val != "") {
      sum_val <- sum_ratios_from_text(val)
      if (sum_val > 0) {
        parts <- c(parts, paste0(code, " (", format(round(sum_val, 2), nsmall = 2), ")"))
      }
    }
  }
  paste(parts, collapse = " + ")
}

# Classify PVS1 overall strength from combined_strength
get_pvs1_strength_class <- function(strength_string) {
  if (is.na(strength_string) || strength_string == "") return(NA)
  
  extract_val <- function(code) {
    pattern <- paste0(code, " \\((\\d+\\.\\d+)\\)")  # dot-decimal
    val <- stringr::str_match(strength_string, pattern)[, 2]
    if (!is.na(val)) as.numeric(val) else 0
  }
  
  PVS1   <- extract_val("PVS1")
  PVS1_S <- extract_val("PVS1_S")
  PVS1_M <- extract_val("PVS1_M")
  sum_supporting <- PVS1_S + PVS1_M
  
  if (PVS1 > 0.9) {
    "PVS1_Strong (RNA)"
  } else if (PVS1_S > 0.9 || (PVS1 > 0.8 && PVS1 <= 0.9)) {
    "PVS1_Moderate (RNA)"
  } else if (PVS1_M > 0.9 || (sum_supporting > 0.8 && sum_supporting <= 0.9)) {
    "PVS1_Supporting (RNA)"
  } else if (sum_supporting >= 0.7 && sum_supporting <= 0.8) {
    "PVS1_N/A (RNA)"
  } else {
    NA
  }
}

# Classify BP7 from combined_strength
check_bp7_strength <- function(strength_string) {
  if (is.na(strength_string) || strength_string == "") return(NA)
  val <- stringr::str_match(strength_string, "BP7_S \\((\\d+\\.\\d+)\\)")[,2]
  if (!is.na(val) && as.numeric(val) > 0.9) "BP7 (RNA)" else NA
}

# -------------------------------
# Load & prepare main variant table
# -------------------------------
df_variants <- readxl::read_excel("2025_08_07_SDHB MG transcripts shorter_ACMG.xlsx") %>%
  janitor::clean_names() %>%
  dplyr::rename(
    hgvsg_hg38   = hgv_sg_hg38,
    hgvsc        = hgv_sc,
    rna_hgvs     = rna_hgvs_nc_000001_11_nm_003000_3_3,
    protein_hgvs = protein_hgvs_np_002991_22
  ) %>%
  dplyr::mutate(
    hgvsc_short        = ifelse(variant == "MG_WT", "MG_WT", stringr::str_extract(hgvsc, "(?<=:).*")),
    rna_hgvs_short     = stringr::str_extract(rna_hgvs, "(?<=:).*"),
    protein_hgvs_short = stringr::str_extract(protein_hgvs, "(?<=:).*"),
    ratio_clean        = as.numeric(stringr::str_replace(ratio, ",", ".")),
    acmg_info          = paste0(sprintf("%.2f", ratio_clean), " (", rna_hgvs_short, ", ", protein_hgvs_short, ")")
  ) %>%
  # keep analysis on variants only
  dplyr::filter(variant != "MG_WT")

# Optional: force FL_WT / N/A labels when both RNA & protein are NA
df_variants <- df_variants %>%
  dplyr::mutate(
    acmg_info = dplyr::case_when(
      is.na(rna_hgvs_short) & is.na(protein_hgvs_short) & transcript_acmg_code == "BP7_S" ~
        paste0(sprintf("%.2f", ratio_clean), " (FL_WT)"),
      is.na(rna_hgvs_short) & is.na(protein_hgvs_short) ~
        paste0(sprintf("%.2f", ratio_clean), " (N/A)"),
      TRUE ~ acmg_info
    )
  )

# Superscripts: derive from transcript_acmg_code and append to acmg_info tail
df_variants <- df_variants %>%
  dplyr::mutate(na_superscripts = purrr::map_chr(transcript_acmg_code, extract_na_superscripts)) %>%
  dplyr::mutate(
    na_superscripts = dplyr::coalesce(na_superscripts, ""),
    acmg_info       = dplyr::if_else(na_superscripts != "", paste0(acmg_info, na_superscripts), acmg_info)
  )

# Clean ACMG classes (collapse any "N/A…" to plain "N/A")
df_variants <- df_variants %>%
  dplyr::mutate(
    acmg_code_clean = dplyr::case_when(
      stringr::str_detect(transcript_acmg_code, "N/A") ~ "N/A",
      TRUE ~ transcript_acmg_code
    )
  )

# -------------------------------
# Aggregate to wide ACMG table
# -------------------------------
acmg_codes <- c("PVS1", "PVS1_S", "PVS1_M", "BP7_S", "N/A")

df_filtered <- df_variants %>%
  dplyr::filter(acmg_code_clean %in% acmg_codes) %>%
  dplyr::mutate(acmg_code_clean = factor(acmg_code_clean, levels = acmg_codes))

df_acmg_table <- df_filtered %>%
  dplyr::group_by(variant, acmg_code_clean) %>%
  dplyr::summarise(
    combined_info = paste(na.omit(acmg_info), collapse = "; "),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from  = acmg_code_clean,
    values_from = combined_info,
    values_fill = ""
  )

# Keep one hgvsc_short per variant
df_hgvsc <- df_filtered %>%
  dplyr::group_by(variant) %>%
  dplyr::summarise(hgvsc_short = dplyr::first(hgvsc_short), .groups = "drop")

df_acmg_table_extended <- df_acmg_table %>%
  dplyr::left_join(df_hgvsc, by = "variant") %>%
  dplyr::relocate(hgvsc_short, .after = variant)

# -------------------------------
# Join HGVSp from second sheet
# -------------------------------
all_variants <- readxl::read_excel("Variants 0,7 filtered.xlsx") %>%
  dplyr::mutate(HGVSp = stringr::str_extract(HGVSp, "(?<=:).*")) %>%
  dplyr::mutate(variant = paste0("V", V))

df_acmg_table_extended <- df_acmg_table_extended %>%
  dplyr::left_join(dplyr::select(all_variants, variant, HGVSp), by = "variant") %>%
  dplyr::select(variant, hgvsc_short, HGVSp, dplyr::all_of(acmg_codes))

# -------------------------------
# Combined strength & classifications
# -------------------------------
df_acmg_table_extended <- df_acmg_table_extended %>%
  dplyr::rowwise() %>%
  dplyr::mutate(combined_strength = build_combined_strength(dplyr::pick(dplyr::all_of(acmg_codes)), acmg_codes)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    pvs1_strength_classification = sapply(combined_strength, get_pvs1_strength_class),
    bp7_strength_classification  = sapply(combined_strength, check_bp7_strength)
  )

# Append consolidated N/A superscripts to combined_strength (if present)
na_superscripts_all <- df_variants %>%
  dplyr::filter(stringr::str_detect(transcript_acmg_code, "N/A")) %>%
  dplyr::select(variant, na_superscripts) %>%
  dplyr::filter(!is.na(na_superscripts), na_superscripts != "") %>%
  dplyr::group_by(variant) %>%
  dplyr::summarise(
    na_superscripts = paste0(
      sort(unique(unlist(strsplit(paste0(na_superscripts, collapse = ""), "")))),
      collapse = ""
    ),
    .groups = "drop"
  )

df_acmg_table_extended <- df_acmg_table_extended %>%
  dplyr::left_join(na_superscripts_all, by = "variant") %>%
  dplyr::mutate(na_superscripts = tidyr::replace_na(na_superscripts, "")) %>%
  dplyr::mutate(
    combined_strength = ifelse(
      na_superscripts != "" & stringr::str_detect(combined_strength, "N/A \\(\\d+\\.\\d+\\)"),
      stringr::str_replace(combined_strength, "(N/A \\(\\d+\\.\\d+\\))", paste0("\\1", na_superscripts)),
      combined_strength
    )
  )

# -------------------------------
# Write Excel output
# -------------------------------
openxlsx::write.xlsx(df_acmg_table_extended,
                     file = "ACMG_table_extended.xlsx",
                     asTable = TRUE)

# End of script
