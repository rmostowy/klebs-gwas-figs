# GWAS and capsular specificity in *Klebsiella pnaeumoniae*

## Overview

This repository contains the code used to generate all figures in the manuscript titled **“Capsular specificity in temperate phages of *Klebsiella pneumoniae* is driven by diverse receptor-binding enzymes”** by **Aleksandra Otwinowska, Janusz Koszucki, Vyshakh R. Panicker, Jade Leconte, Sebastian Olejniczak, Kathryn E. Holt, Edward J. Feil, Eduardo P. C. Rocha, Bogna Smug, Barbara Maciejewska, Zuzanna Drulis-Kawa and Rafal J. Mostowy** (submitted to *PLOS Biology*). This includes both main-text figures and supplementary figures. The analyses integrate genome-wide association studies (GWAS) for prophage-encoded proteins in *Klebsiella pneumoniae* with experimental validation of selected candidate enzymes. The latest version of the preprint can be accessed [here](https://www.biorxiv.org/content/10.1101/2025.06.25.661490v2).

## Repository Structure
```
├── scripts/				# Scripts generating Figures 1–7 and all Supplementary Figures
   ├── data/				# Scripts to process some data and save locally
   ├── figure1/				# All files pertaining to figure 1
   ├── figure2/				# All files pertaining to figure 2
   ├── ...					# etc.
   ├── other/				# folders containing other supplementary figures
   ├── list-of-figures.txt	# Translation between names of files and actual figures and supplementary figures in the paper
├── config/					# Config file with all paths and paramteres used
└── README.md
```

Each script is named according to the figure it generates (e.g. `fig2_*.R`, `supp_fig12_*.R`).

## Reproducibility

All scripts were used to generate the figures included in the submitted manuscript version.  The data required to reproduce the figures are available via [Figshare](https://doi.org/10.6084/m9.figshare.29181188) as described in the manuscript Data Availability statement.

This repository has been archived in Zenodo to ensure a permanent, citable version corresponding to the submitted manuscript. Please refer to the manuscript or the Figshare repository to obtain an up-to-date DOI.

## Contact

For questions regarding the code or figure generation, please contact:

Rafal Mostowy  
rafal.mostowy@uj.edu.pl
