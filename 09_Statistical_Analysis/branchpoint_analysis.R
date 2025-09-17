# -----------------------------------------------------------
# Intron features analysis — PPT & branch point finder
# -----------------------------------------------------------
# Input:   - seq_list — a named list of DNAString objects (adjust sequences below).
# Process:
#   Analyzes intron sequences (as Biostrings::DNAString) to locate:
#     - the strongest polypyrimidine tract (PPT) within A-40..A-5
#     - a candidate branch point (BP) within A-60..A-10 using consensus motifs
#   Assumes each input intron ends with the canonical acceptor "AG";
#   the position of the 'A' in that terminal "AG" is used as reference (c.(Exon)-2).
# Output:  - Excel file "intron_analysis_results.xlsx" with one row per intron.
# -----------------------------------------------------------

# -------------------------------
# 1) Load packages
# -------------------------------
library(Biostrings)
library(IRanges)  
library(writexl)

# -------------------------------
# 2) Define three intron sequences
#    (Adjust the strings to your actual sequences.)
# -------------------------------
seq_list <- list(
  "Intron 2" = DNAString("ATGGATAAGCTAATACATCCAGGTGTCTCCGATTATATTATGATAAAGTGTAGGGAGGTTGAACTTTACATAAATACCACTGGATATTTTTCTTTGTTAG"),
  "Intron 3" = DNAString("TGCCCCTGATGGCACAGCAAGGAGGATCCAGAAGAAAGTATTTGGGGCAGGACTGATTCCGGATATGGGTGAGGATGTGTTAAATGTGTGTCTCTTTCAG"),
  "Intron 4" = DNAString("CGAGTAGTCAGTGTCCAAGAAATGGGGTAAATAAAGCTGAGGTGATGATGGAATCTGATCCTTTTCTTCTTCTTCTTCTTCTTCTTCCTCTTAACCACAG")
)

# -------------------------------
# 3) Function to analyze a single intron
# -------------------------------
analyze_intron <- function(intron_seq, seq_name = "") {
  # L = intron length (including terminal "AG")
  L     <- length(intron_seq)
  # pos_A = 1-based position of the 'A' in the terminal "AG"
  pos_A <- L - 1
  
  # --- Polypyrimidine tract (PPT) ---
  # Window: pos_A-40 .. pos_A-5 (1-based, inclusive)
  ppt_start_idx <- max(1, pos_A - 40)
  ppt_end_idx   <- max(1, pos_A - 5)
  ppt_window    <- subseq(intron_seq, start = ppt_start_idx, end = ppt_end_idx)
  
  # Convert to vector of single bases
  ppt_chars <- strsplit(as.character(ppt_window), "")[[1]]
  
  # Run-length encode pyrimidine runs (C/T)
  mask       <- ppt_chars %in% c("C", "T")
  rle_mask   <- Rle(mask)
  run_vals   <- runValue(rle_mask)    # TRUE/FALSE
  run_lens   <- runLength(rle_mask)   # run lengths
  pyr_runs   <- which(run_vals == TRUE)
  
  if (length(pyr_runs) == 0) {
    ppt_start_global <- NA_integer_
    ppt_end_global   <- NA_integer_
    ppt_cDNA         <- NA_character_
  } else {
    best_run_idx      <- pyr_runs[which.max(run_lens[pyr_runs])]
    run_ends_cumsum   <- cumsum(run_lens)
    run_starts_cumsum <- run_ends_cumsum - run_lens + 1
    ppt_rel_start     <- run_starts_cumsum[best_run_idx]
    ppt_rel_end       <- run_ends_cumsum[best_run_idx]
    ppt_start_global  <- ppt_start_idx + ppt_rel_start - 1
    ppt_end_global    <- ppt_start_idx + ppt_rel_end   - 1
    # cDNA notation relative to the exon acceptor: A of AG = c.(Exon)-2
    n_old_start  <- pos_A - ppt_start_global
    n_old_end    <- pos_A - ppt_end_global
    n_cdna_start <- n_old_start + 2
    n_cdna_end   <- n_old_end   + 2
    ppt_cDNA     <- paste0("c.(Exon)-", n_cdna_start, "_-", n_cdna_end)
  }
  
  # --- Branch point (BP) ---
  # Extended window: pos_A-60 .. pos_A-10
  bp_start_idx <- max(1, pos_A - 60)
  bp_end_idx   <- max(1, pos_A - 10)
  bp_window    <- subseq(intron_seq, start = bp_start_idx, end = bp_end_idx)
  bp_string    <- as.character(bp_window)
  
  # Branch-point consensus motif regexes (liberal set)
  motifs <- list(
    YNCURAY = "[CT][ACGT]C[CT][AG]A[CT]",
    CTRAY   = "C[CT][AG]A[CT]",
    NAGR    = "[ACGT]AG[AG]",
    TCTAAC  = "TCTAAC",
    YTNCTR  = "[CT]T[ACGT]CT[AG]"
  )
  
  # Scan for motifs; pick the one closest to the acceptor (largest global position)
  bp_candidates <- list()
  for (motif_name in names(motifs)) {
    pat <- motifs[[motif_name]]
    matches <- gregexpr(pattern = pat, text = bp_string, perl = TRUE)[[1]]
    if (all(matches == -1)) next
    ends <- matches + attr(matches, "match.length") - 1
    best_idx <- which.max(ends)  # closest to acceptor within this motif group
    rel_start <- matches[best_idx]
    rel_end   <- ends[best_idx]
    global_pos <- bp_start_idx + rel_start - 1
    bp_candidates[[length(bp_candidates) + 1]] <- data.frame(
      motif      = motif_name,
      rel_start  = rel_start,
      rel_end    = rel_end,
      global_pos = global_pos,
      length     = rel_end - rel_start + 1,
      stringsAsFactors = FALSE
    )
  }
  
  if (length(bp_candidates) > 0) {
    bp_df <- do.call(rbind, bp_candidates)
    best_bp_row <- bp_df[which.max(bp_df$global_pos), ]  # overall closest to acceptor
    bp_global <- best_bp_row$global_pos
    n_bp_old  <- pos_A - bp_global
    n_bp_cdna <- n_bp_old + 2
    bp_cDNA   <- paste0("c.(Exon)-", n_bp_cdna)
    bp_motif  <- best_bp_row$motif
  } else {
    # Fallback: choose the rightmost 'A' in the window
    chars <- strsplit(bp_string, "")[[1]]
    all_A_positions <- which(chars == "A")
    if (length(all_A_positions) > 0) {
      relA      <- max(all_A_positions)
      bp_global <- bp_start_idx + relA - 1
      n_bp_old  <- pos_A - bp_global
      n_bp_cdna <- n_bp_old + 2
      bp_cDNA   <- paste0("c.(Exon)-", n_bp_cdna)
      bp_motif  <- "Fallback_A"
    } else {
      bp_global <- NA_integer_
      bp_cDNA   <- NA_character_
      bp_motif  <- NA_character_
    }
  }
  
  # Result row
  list(
    Sequence      = seq_name,
    Intron_Length = L,
    pos_A         = pos_A,
    PPT_Start     = ppt_start_global,
    PPT_End       = ppt_end_global,
    PPT_cDNA      = ppt_cDNA,
    BP_motif      = bp_motif,
    BP_Global     = bp_global,
    BP_cDNA       = bp_cDNA
  )
}

# -------------------------------
# 4) Analyze all introns
# -------------------------------
all_results <- lapply(names(seq_list), function(nm) {
  analyze_intron(seq_list[[nm]], seq_name = nm)
})

# Convert to data.frame
df_results <- do.call(rbind, lapply(all_results, as.data.frame, stringsAsFactors = FALSE))

# -------------------------------
# 5) Write Excel
# -------------------------------
write_xlsx(df_results, path = "intron_analysis_results.xlsx")
cat("Results saved to 'intron_analysis_results.xlsx' in:\n", getwd(), "\n")
