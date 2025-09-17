# ===============================================================
# Sashimi plot batch runner (ggsashimi, custom unique+cutoff)
# ===============================================================
# What this does
# - Loops over subfolders in INPUT_DIR
# - For each sample, finds one BAM and the corresponding STAR junctions (*_SJ.out.tab)
# - Computes a per-sample dynamic junction cutoff = 1% of the max junction count (min=1)
# - Calls adjusted ggsashimi script (my_ggsashimi.py) with that cutoff
# - Produces one PDF per sample in OUT_DIR
#
# Notes
# - This script assumes STAR-style *_SJ.out.tab where column 7 = junction read count.
# - Unique-read handling and cutoff application are implemented in my_ggsashimi.py.
# - Defaults below are your original example values — adjust as needed.
# ===============================================================

# --- User settings (defaults = your example) ---
INPUT_DIR="./data"
OUT_DIR="./sashimi_junction_plots"
REGION="MG_SDHB_3.1:688-3593"
GTF="sdhb_minigene_plasmid_spliced_fixed.gtf"

# Python + patched ggsashimi (adjust these paths for your machine)
PYTHON="/Users/.venv/bin/python"
GGSASHIMI_PY="/Users/my_ggsashimi.py"

# Visual parameters (passed through to ggsashimi)
COLORSET="-C 2"
ALPHA="--alpha 0.6"

# ---------------------------------------------------------------

mkdir -p "${OUT_DIR}"

shopt -s nullglob
for dir in "${INPUT_DIR}"/*/; do
  # pick the first BAM and STAR junction table in the subfolder
  bam_files=("$dir"/*.bam)
  sj_tabs=("$dir"/*_SJ.out.tab)

  if [[ ${#bam_files[@]} -eq 0 || ${#sj_tabs[@]} -eq 0 ]]; then
    echo "[SKIP] No BAM or *_SJ.out.tab in: $dir"
    continue
  fi

  bam="${bam_files[0]}"
  sj_tab="${sj_tabs[0]}"
  base=$(basename "$bam" .bam)
  out_pref="${OUT_DIR}/sashimi_${base}"

  # 1) Compute max junction support (col 7) and a 1% cutoff, rounded; min=1
  #    (awk version avoids external bc/sort dependencies)
  read -r max_count cutoff < <(
    awk 'max<$7{max=$7} END{
      if(max=="") max=0;
      c = int(max*0.01 + 0.5);  # round to nearest
      if(c < 1) c = 1;
      printf("%d %d\n", max, c)
    }' "$sj_tab"
  )
  echo "[${base}] max_junction=${max_count}, 1% cutoff=${cutoff}"

  # 2) Run ggsashimi (coverage + junctions) with shrink and dynamic cutoff
  cmd=( "$PYTHON" "$GGSASHIMI_PY"
        -b "$bam"
        -j "$sj_tab"
        -c "$REGION"
        -g "$GTF"
        -M "$cutoff"
        -o "$out_pref"
        $COLORSET
        $ALPHA )

  # Add grouping file if it exists
  if [[ -n "$GROUPS_FILE" && -f "$GROUPS_FILE" ]]; then
    cmd+=( -P "$GROUPS_FILE" )
  fi

  echo "[RUN] ${cmd[*]}"
  "${cmd[@]}"

  echo "-> ${out_pref}.pdf generated"
done
