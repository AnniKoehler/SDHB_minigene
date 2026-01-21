# Minigene-based characterization and classification of splice-associated variants in succinate dehydrogenase B


## Overview
This repository contains scripts and reference files used for the functional assessment of **SDHB splice-associated variants** using a **minigene-based splicing assay** in HEK293T cells.

The workflow integrates:
- variant prioritization  
- minigene-specific reference construction  
- targeted NGS based RNA sequencing and quantitative transcript analysis  
- HGVS annotation of aberrant transcripts  
- integration of RNA analysis into ACMG/AMP criteria  
- assignment of ACMG/AMP points and final variant classification 

This repository accompanies a manuscript submitted to **npj Precision Oncology** and is intended to ensure analytical transparency and reproducibility.

## Experimental Context

| Item | Description |
|----|----|
| **Gene** | SDHB |
| **Reference transcript** | MANE Select NM_003000.3 |
| **Cell line** | HEK293T |
| **Assay type** | Minigene splicing assay |
| **Primary readout** | Targeted NGS based RNA sequencing |


## Repository Structure

SDHB_minigene/\
  \t01_Synthetic_Variant_Table/\
  \t02_SpliceAI_Run/\
  03_Variant_Prioritization/\
  04_Minigene_Reference/\
  05_NGS_Workflow/\
  06_Automated_HGVS_Nomenclature/\
  07_ACMG_Classification/\
  08_Additional_Plots_Code/\
  09_Statistical_Analysis/\
  README.md\

## Folder Descriptions

### 01_Synthetic_Variant_Table
Generation of a synthetic variant table creating all possible single nucleotide substitutions in SDHB. Variants were annotated using the Ensembl Variant Effect Predictor (VEP) (accessed 30.06.2023) and subsequently used for variant prioritization.

### 02_SpliceAI_Run
In silico annotation of synthetic variants located within the SDHB minigene region (exons 2–5, including 325 bp flanking intronic sequence) using SpliceAI, reporting both ALT and REF scores.

### 03_Variant_Prioritization
Prioritization of variants for functional testing based on predicted splice-altering potential using a SpliceAI Δ-score cutoff ≥ 0.25 to maximize sensitivity.

### 04_Minigene_Reference
Construction of custom minigene-specific FASTA and GTF references, including exon–intron boundaries and flanking plasmid sequence corresponding to the SDHB minigene construct. These references were used for read alignment.

### 05_NGS_Workflow
Targeted RNA sequencing preprocessing and alignment pipeline using STAR (Dobin et al., 2013), adapted to a custom minigene reference and optimized parameter settings. BAM files and splice-junction files will be deposited at the German Human Genome Archive (GHGA).

### 06_Automated_HGVS_Nomenclature
Automated generation of HGVS-compliant RNA- and protein-level variant nomenclature for insertions (intron retention) and deletions (partial or complete exon skipping). RNA-level variants (r.( )) were reconstructed from splice-junction coordinates, and corresponding protein consequences (p.( )) were inferred by translation of the altered mRNA sequence (Baumann, 2021). Complex splicing events—including pseudo-exons, multi-exon skipping, indels, and double intron retention—were annotated manually.

### 07_ACMG_Classification
Scripts used to integrate functional splicing results into ACMG/AMP-compliant variant classification, including: pooling of transcript-level code strengths, overall RNA evidence assignment, calculation of ACMG/AMP points per variant analysis and reclassification impact

### 08_Additional_Plots_Code
Supplementary plotting scripts used to generate:
Figure 1: SDHB protein model, lollipop plot, branchpoint analysis
Figure 2: SpliceAI heatmap and transcript analysis
Figure 4: Sashimi plots

### 09_Statistical_Analysis
Statistical analysis of transcript-level readouts and SpliceAI scores for variants selected for functional testing.

## Citation
If you use this repository or adapt its workflow, please cite: Köhler A. et al.
Minigene-based characterization and classification of splice-associated variants in succinate dehydrogenase B.
npj Precision Oncology (in submission).
A DOI will be added upon publication.
