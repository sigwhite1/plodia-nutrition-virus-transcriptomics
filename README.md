# plodia-nutrition-virus-transcriptomics
Analysis code for transcriptomic responses to baculovirus infection in Plodia interpunctella under contrasting nutritional environments, examining how resource availability reorganizes host antiviral transcriptional programs.

## Study overview
Host-pathogen interactions are strongly influenced by environmental conditions. This project analyzes a microarray dataset examining how nutritional environment shapes the transcriptional response of *Plodia interpunctella* larvae following infection with *Plodia interpunctella granulosis virus* (PiGV).

The analysis identifies coordinated transcriptional modules and evaluates how resource availability reorganizes host physiological programs during infection.

## Analysis workflow
The analysis pipeline is organized into numbered scripts that should be run in the following order:
1. `R_QC_microarray.R`  
   Import raw Agilent two-color microarray data, perform quality control, background correction, and normalization. Excludes low-quality arrays based on A-value correlation.
2. `R_00_config.R`  
   Define shared configuration, file paths, analysis constants, and plotting defaults used by all downstream scripts.
3. `R_01_load_inputs.R`  
   Load and validate QC-normalized data and sample metadata.
4. `R_02_design_matrix.R`  
   Construct the loop design matrix for two-color microarray analysis using limma's modelMatrix approach.
5. `R_03_contrasts.R`  
   Define virus-vs-control contrasts for each diet × time combination.
6. `R_04_limma_fit.R`  
   Fit limma linear models, apply empirical Bayes moderation, and extract differential expression results across all contrasts.
7. `R_05_reporting_barplot.R`  
   Summarize DEG counts per contrast, visualize diet overlap, and compute Jaccard similarity between dietary conditions.
8. `R_06_select_DEGs_for_clustering.R`  
   Define the DEG universe for clustering as the union of significant probes across all contrasts.
9. `R_07_cluster_k_sweep.R`
    Evaluate hierarchical clustering solutions across k = 2–5 using silhouette width and within-cluster SSE, and assess bootstrap stability of the chosen k = 4 solution.
10. `R_08_clustering_all_contrasts.R`
    Assign DEGs to four expression clusters using hierarchical clustering with stabilized cluster labels anchored to a biological reference contrast.
11. `R_08b_annotate_clusters.R`
    Map cluster assignments to Drosophila melanogaster orthologs via DIAMOND BLASTP results and build the master probe annotation table.
12. `R_09_enrichment_GO_clusters.R`
    Perform GO enrichment analysis for each expression cluster using g:Profiler against a custom microarray background.
13. `R_10_cluster_summaries.R`
    Compute cluster-level trajectory summaries and generate mean ± 95% CI expression plots across diet × time.
14. `R_11_immune_gene_evaluation.R`
    Identify immune-annotated probes using GO:BP immune categories and generate heatmaps of immune gene expression across diet × time contrasts.
15. `R_12_immune_pathway_dynamics.R`
    Summarize immune gene dynamics using gene-blocked linear models and test for diet-dependent effects across GO immune categories.
16. `R_13_deg_cluster_enrichment.R`
    Test whether DEGs at each diet × time contrast are non-randomly distributed across expression clusters using Fisher's exact tests.

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

