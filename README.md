# RNA Sequence Analysis

## RNA-seq Analysis of Reprogrammed and non-Reprogrammed Somatic Cell

This directory contains bioinformatic analysis of RNA-seq data from a mouse embryonic stem cell (ESC) reprogramming experiment. The analysis focuses on the role of the TRIM28 gene in reprogramming, based on the article "TRIM28 Is an Epigenetic Barrier to Induced Pluripotent Stem Cell Reprogramming" by Myles et al. ([http://dx.doi.org/10.1002/stem.2453](http://dx.doi.org/10.1002/stem.2453)).

This project was submitted as part of the Bioinformatics Department's OMICS task 1 at the Moscow Institute of Physics and Technology (MIPT).

### Data Description

- **Date downloaded**: March 16th, 2024
- **Experiment**: RNA-sequencing (RNA-seq)
- **Organism**: Mouse
- **Cell lines**:
  - Reprogrammed (iPSCs):
    - SRR3414629
    - SRR3414630
    - SRR3414631
  - Non-reprogrammed (control) - mouse ESCs:
    - SRR3414635
    - SRR3414636
    - SRR3414637
- **Library type**: Stranded RNA-seq
- **Sequencing platform**: Illumina HiSeq 2500
- **Source**: Transcriptomic
- **Selection**: Poly(A) selected RNA
- **Layout**: Single-end
- **Construction protocol**:
  - Cells were harvested using standard procedures and RNA was extracted using Trizol reagent.
  - Illumina TruSeq RNA Sample Prep Kit (Cat#FC-122-1001) was used with 1 ug of total RNA for the construction of sequencing libraries.
  - RNA libraries were prepared for sequencing using standard Illumina protocols.

### Purpose

This analysis aims to investigate the differential gene expression between the reprogrammed and non-reprogrammed dataset.

### Target Audience

This project is intended for researchers with a basic understanding of RNA-seq analysis.

### Software

The analysis pipeline utilizes tools like HISAT2 for reads alignment, HTcount for counts aligned reads, DESeq2 for Differential Expression Analysis, and fgsea for Gene Set Expression Analysis.

### Scripts

Analysis scripts can be found in the 'scripts' subdirectory.

