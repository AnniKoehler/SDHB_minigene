# Custom ggsashimi (my_ggsashimi.py)

## Overview
This repository contains a customized version of [ggsashimi](https://github.com/guigolab/ggsashimi) for creating Sashimi plots.  
The base is the original ggsashimi script (version 1.1.5).  

The customized version is provided as:
- `my_ggsashimi.py`

For convenient usage, a shell script is included:
- `sashimi_plot_general.sh`

---

## Modifications compared to the original
- **Arcs (splice junctions)** in the plot are no longer derived from BAM files.  
  Instead, they are read directly from STAR `SJ.out.tab` files.  
- Only **unique read counts** (column 7 of `SJ.out.tab`) are considered.  
- All other functionality, parameters, and visualizations remain unchanged from the original script.

---

## Requirements
- **Python 3** with the package `pysam`
- **R** with the following packages:
  - `ggplot2`
  - `grid`
  - `gridExtra`
  - `data.table`
  - `gtable`

---

## Usage

### 1. Prepare your files
- Place your BAM files and the corresponding STAR `*_SJ.out.tab` files into subfolders under `./data/`.  
- Provide a GTF file with exon annotations (e.g. `sdhb_minigene_plasmid_spliced_fixed.gtf`).

### 2. Run the shell script
The included script `sashimi_plot_general.sh` will process all samples automatically:

```bash
bash run_sashimi.sh
