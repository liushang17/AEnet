AEnet
Alternative splicing (AS), as a crutical promoter of proteomic diversity, constitutes a core source of cellular heterogeneity together with gene expression levels. AS is closely linked to various physiological and pathological states including aging, embryonic development, and tumor progression. Recently, rapid developed single-cell transcriptomics technology has greatly advanced the process of deeply analyzing cellular heterogeneity from gene expression level. However, due to technical limitations such as sequencing depth, dropout events, and batch effects, we have not yet found an ideal approach to comprehensively analyze cellular heterogeneity and its regulatory mechanisms at alternative splicing level. Therefore, this study aims to construct a novel method named alternative splicing network(AEN) combining gene expression levels with AS patterns to analyze cellular heterogeneity, identify novel cell types, and infer the regulatory mechanisms of alternative splicing. By applying the AEN method to gastrula data, pan-cancer T-cell data, and lung cancer tumor cell data, we have successfully revealed the hidden alternative splicing heterogeneity in these cells and identified key alternative splicing events and related splicing factors during cell transformation. The application of AEN method provides new insights into the understanding of cellular heterogeneity and its related physiological and pathological processes.

# The steps for AEnet
AENet consists of six steps: ASP Determination, ASP-Expression Network Construction, Multiple-Sample Integration, Key ASP Identification, ASP Cluster Determination, Cell Type Identification Assistance, and Regulatory Mechanism Prediction. Detailed information for each step is provided below.

## ASP Determination


# install
devtools::install_github("https://github.com/liushang17/AEN")

# Demo
We showed the application of a demo in the demo.r. The dataset used in the demo.r is shown in the link (https://pan.baidu.com/s/14ZMS2e6fMk3P1CLdOCK6tw . Extraction code ：fyS9).
