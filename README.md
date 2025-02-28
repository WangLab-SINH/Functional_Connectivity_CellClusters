# Cell Types Associated with Human Brain Functional Connectomes and Their Implications in Psychiatric Diseases
The human whole-brain single-cell data used in this study were obtained from https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443 (Siletti K, et al.). The data were further processed using the script preprocess_cortical_seuratdata.R, resulting in the retention of 151 cell clusters from the expression profiles of 1,033,759 cells for analysis.

# System Requirements
Before running the analysis, ensure that you have the following prerequisites:
- Matlab R2021a
- R 4.0+ (for statistical modules)
- Seurat v4.3.0.1 package
- NeuronChat v1.0.0 package

# Code
Included are the codes necessary to replicate the analyses. All the codes takes approximately two days to run completely.
clusterPls.m and genePls.m are MATLAB functions used to calculate the cell clusters and genes associated with the connectivity strength of resting-state functional network nodes.
preprocess_cortical_seuratdata.R is an R script used for preprocessing the downloaded single-cell data from different brain regions and calculating the proportion of different cell clusters in each brain region.
corticalregion_cluster_pseudobulk.R is an R script used to calculate the pseudobulk gene expression levels of different cell clusters in different brain regions and perform normalization.
neuronchat_ana.R is an R script used to calculate the cell communication strength between different cell clusters across different brain regions and within regions.

# Data
Included are the proportions of different cell clusters in various brain regions (cluster_density_df.csv), the strength of resting-state functional networks at different cut-offs (degree_centrality.csv), and normalized pseudobulk data of different cell clusters in cortical regions (cluster_pseudobulk).
