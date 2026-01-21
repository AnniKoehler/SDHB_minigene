Minigene-based characterization and classification of splice-associated variants in succinate dehydrogenase B 

Overview
This repository contains scripts and reference files used for the functional assessment of SDHB splice-associated variants using a minigene-based splicing assay in HEK293T cells. The workflow integrates variant prioritization, minigene-specific reference construction, targeted RNA mapping, quantitative transcript analysis, HGVS annotation of aberrant transcripts, integration of RNA analysis into ACMG/AMP criteria and assigment of ACMG/AMP points and classes per variant for overall clinical classification. 
This repository accompanies a manuscript submitted to npj Precision Oncology and is intended to ensure analytical transparency and reproducibility.
 
Experimental Context
Gene: SDHB
Reference transcript: MANE Select NM_003000.3
Cell line: HEK293T
Assay type: Minigene splicing assay
Primary readout: targeted NGS of RNA  
Repository Structure
SDHB_minigene/
├── 01_Synthetic_Variant_Table/
├── 02_SpliceAI_Run/
├── 03_Variant_Prioritization/
├── 04_Minigene_Reference/
├── 05_NGS_Workflow/
├── 06_Automated_HGVS_Nomenclature/
├── 07_ACMG_Classification/
├── 08_Additional_Plots_Code/
├── 09_Statistical_Analysis/
└── README.md
 
Folder Descriptions
01_Synthetic_Variant_Table
Generation of a synthetic variant table creating all possible single nucleotide changes in SDHB. The variants were annotated with using the Ensembl Variant Effect Predictor (VEP; accessed 30.06.2023). The table was further used for variant prioritization. 
 
02_SpliceAI_Run
Annotation of synthetic variants in genomic region of SDHB minigene (exon 2 to exon 5 and flanking 325 bp intronic region) using in-silico SpliceAI predicition with ALT and REF scores. 

 
03_Variant_Prioritization
Variants were prioritized for potential splice-altering effect with a SpliceAI Δ-score cut-off ≥ 0.25 to maximize sensitivity. 
 
04_Minigene_Reference
Construction of custom minigene-specific FASTA and GTF references, including exon–intron boundaries and flanking plasmid sequence corresponding to the SDHB minigene construct. These references were used for read alignment. 
 
05_NGS_Workflow
Targeted RNA sequencing preprocessing and alignment pipeline using the Spliced Transcripts Alignment to a Reference (STAR) algorithm (Dobin et al., 2013), adapted to a custom minigene reference and optimized parameter settings. BAM- and splice-junction files will be deposited at the German Human Genome Archive (GHGA).
 
06_Automated_HGVS_Nomenclature
Automated generation of HGVS-compliant RNA- and protein-level variant nomenclature for insertions (intron retention) and deletions (partial or complete exon skipping), ensuring consistent reporting across figures, tables, and ACMG interpretation. The algorithm reconstructed r.( ) RNA variants from the flanking splice-junction coordinates and inferred the corresponding p.( ) protein consequences by translating the altered mRNA sequence (Baumann, 2021). Complex splicing effects - including pseudo-exons, multi-exon skipping, indels and double intron retentions - were annotated manually. 
 
07_ACMG_Classification
Scripts used to integrate functional splicing results into ACMG/AMP-compliant variant classification, including pooling of code strengths assigned per transcript for overall RNA code assignment, calculation of ACMG/AMP points per variants derived from assigned criteria and analysis of reclassification impact. 
 
08_Additional_Plots_Code
Supplementary plotting scripts used to generate Fig. 1 (SDHB protein model, lollipop plot, branchpoint analysis), Fig. 2 (Splice AI heatmap and transcript analysis) as well as Fig. 4 (Sashimi plots).
 
09_Statistical_Analysis
Analysis of transcript readout and SpliceAI scores of variants selected for testing. 
 
Citation
If you use this repository or adapt its workflow, please cite:
Köhler A. et al. Minigene-based characterization and classification of splice-associated variants in succinate dehydrogenase B. npj Precision Oncology (in submission).
A DOI will be added upon publication.

