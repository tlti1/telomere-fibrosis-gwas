# telomere-fibrosis-gwas
Genome wide analysis of genetic overlap between telomere length and fibrosis across 12 organ systems using GWAS summary stats, LDSC, colocolisation and bioinformatic follow up. 

# Genetic Overlap Between Telomere Length and Fibrosis Across 12 Organ Systems

## Overview
This repository contains all code used for the Individual Research Project (IRP) 
component of the MSc Bioinformatics programme at the University of Leicester (2025-2026).

The project investigates the genetic overlap between telomere length and fibrosis 
across 12 organ systems using GWAS summary statistics, with analyses including 
LD Score Regression, colocalisation and bioinformatic follow-up.

## Data Sources
- Telomere length GWAS: Burren et al. 2024, Nature Genetics (GCST90435144)
- Fibrosis GWAS: Joof et al. 2026, medRxiv

Data is not stored in this repository.

## Analysis Pipeline
1. Quality control and variant filtering
2. Variant lookup across 12 fibrosis organ systems
3. Allele alignment and direction of effect analysis
4. LD Score Regression (LDSC) genetic correlation
5. Bonferroni correction and significant variant identification
6. Colocalisation analysis (coloc.abf)
7. Bioinformatic follow-up (GTEx, Open Targets)
8. Visualisation (mirror plots, forest plots)

## Requirements
### R packages
- coloc (v6.0.1)
- ggplot2
- writexl

### Command line tools
- PLINK (v1.9)
- LDSC (v1.0.1)
- Python 2.7 (for LDSC)

## Author
Teresa Ironman
MSc Bioinformatics, University of Leicester
Supervisors: Dr Richard Allen, Dr Beatriz Guillen Guio
