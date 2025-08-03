# AEnet
Alternative splicing (AS), as a crutical promoter of proteomic diversity, constitutes a core source of cellular heterogeneity together with gene expression levels. AS is closely linked to various physiological and pathological states including aging, embryonic development, and tumor progression. Recently, rapid developed single-cell transcriptomics technology has greatly advanced the process of deeply analyzing cellular heterogeneity from gene expression level. However, due to technical limitations such as sequencing depth, dropout events, and batch effects, we have not yet found an ideal approach to comprehensively analyze cellular heterogeneity and its regulatory mechanisms at alternative splicing level. Therefore, this study aims to construct a novel method named alternative splicing network(AEN) combining gene expression levels with AS patterns to analyze cellular heterogeneity, identify novel cell types, and infer the regulatory mechanisms of alternative splicing. By applying the AEN method to gastrula data, pan-cancer T-cell data, and lung cancer tumor cell data, we have successfully revealed the hidden alternative splicing heterogeneity in these cells and identified key alternative splicing events and related splicing factors during cell transformation. The application of AEN method provides new insights into the understanding of cellular heterogeneity and its related physiological and pathological processes.

## The principle of AEnet
AEnet is a computational method designed to uncover cellular alternative splicing (AS) heterogeneity and regulatory mechanisms by integrating splicing patterns with gene expression in single-cell RNA-seq (scRNA-seq) data. It operates in three major phases:

Quantification of Alternative Splicing Patterns (ASPs):
AEnet computes the percent spliced-in (PSI) values of ASPs using junction reads. PSI is calculated per cell, with missing values (NaNs) arising from undetected junctions due to data sparsity. To address the high dropout rates in scRNA-seq, AEnet focuses only on cells with valid PSI and expression values when analyzing ASP–gene expression (ASP-EXP) relationships.

Construction of ASP–Expression Links:
Spearman correlations between ASPs and gene expression levels are computed across valid cells. Only statistically significant and reproducible ASP-EXP links are retained—especially those consistent across multiple samples, thereby minimizing batch effects. These links reflect potential regulatory relationships, such as splicing factor-driven AS regulation or splicing-associated expression dynamics.

Clustering and Functional Interpretation:
AEnet ranks ASPs based on the number of significant ASP-EXP links and selects top-ranking ASPs as "anchors." Using Jaccard similarity between their associated gene sets, AEnet clusters ASPs and associated genes into splicing modules and co-expression programs. These clusters enable the identification of:

Cell subpopulations with distinct splicing profiles,

Key splicing factors influencing ASP usage,

Biological pathways modulated by specific ASPs.

In summary, AEnet systematically links AS variation with gene expression at single-cell resolution to detect splicing-driven cell states and regulatory mechanisms, providing insights into transcriptomic complexity beyond conventional expression-based clustering.

![Scheme of AEnet methods.](https://github.com/liushang17/AEnet/blob/version1.1/image.png)

## The steps for AEnet
AEnet comprises six analytical steps: (1) Alternative Splicing Pattern (ASP) Determination, (2) ASP–Gene Expression Network Construction, (3) Multi-Sample Integration, (4) Key ASP/Gene Identification, (5) ASP/Gene Cluster Inference, and (6) Cell Type Annotation Assistance and Regulatory Mechanism Prediction.


## The output of AEnet


## The demo datasets
The demonstration datasets can be downloaded from the following links:

iPSC dataset: https://drive.google.com/file/d/1Qkg4De3DER4Qs5V_vP-M7GgwprwjwkLA/view?usp=drive_link

T cell dataset: https://drive.google.com/file/d/1zuut5OlYgFeYXytU5kUZsYKbx07CKb56/view?usp=drive_link

For browser-based visualization of the analysis results:

Single-sample dataset (iPSC): https://liushang17.github.io/ipsc.html

Multi-sample dataset (T cells): https://liushang17.github.io/tcells.html

Two example R scripts, demo.iPSC.R and demo.Tcells.R, are provided to demonstrate the usage of AEnet with the iPSC and T cell datasets, respectively. 

## install
devtools::install_github("https://github.com/liushang17/AEnet")

