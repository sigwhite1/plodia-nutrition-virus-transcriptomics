# plodia-nutrition-virus-transcriptomics
Analysis code for transcriptomic responses to baculovirus infection in Plodia interpunctella under contrasting nutritional environments, examining how resource availability reorganizes host antiviral transcriptional programs.

## Study overview
Host-pathogen interactions are strongly influenced by environmental conditions. This project analyzes a microarray dataset examining how nutritional environment shapes the transcriptional response of *Plodia interpunctella* larvae following infection with *Plodia interpunctella granulosis virus* (PiGV).

The analysis identifies coordinated transcriptional modules and evaluates how resource availability reorganizes host physiological programs during infection.

## Analysis workflow
The analysis pipeline is organized into numbered scripts that should be run in the following order:
1. `R_QC_microarray.R`  
   Import raw Agilent two-color microarray data, perform quality control, background correction, and normalization. Excludes low-quality arrays based on A-value correlation.
2. `02_normalization_limma.R`  
   Normalize expression data and fit differential expression models using limma.
3. `03_define_deg_sets.R`  
   Identify differentially expressed genes across infection and nutritional treatments.
4. `04_build_trajectory_matrix.R`  
   Construct temporal expression trajectories across the infection time course.
5. `05_k_sweep_clustering.R`  
   Evaluate clustering solutions across a range of k values.
6. `06_cluster_stability.R`  
   Assess robustness of clusters using bootstrap resampling.
7. `07_cluster_summaries.R`  
   Summarize transcriptional modules and visualize temporal patterns.
8. `08_functional_enrichment.R`  
   Perform functional enrichment and immune gene set analyses.

## Computational environment
All analyses were conducted in R.
Key packages include:

- limma
- tidyverse
- cluster
- pheatmap
- fgsea
- gprofiler2

## Data availability
The microarray dataset analyzed in this project is available through the
European Bioinformatics Institute ArrayExpress/BioStudies repository.

Accession: E-MTAB-5868
https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-5868

