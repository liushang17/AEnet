# AEnet
Alternative splicing (AS), as a crutical promoter of proteomic diversity, constitutes a core source of cellular heterogeneity together with gene expression levels. AS is closely linked to various physiological and pathological states including aging, embryonic development, and tumor progression. Recently, rapid developed single-cell transcriptomics technology has greatly advanced the process of deeply analyzing cellular heterogeneity from gene expression level. However, due to technical limitations such as sequencing depth, dropout events, and batch effects, we have not yet found an ideal approach to comprehensively analyze cellular heterogeneity and its regulatory mechanisms at alternative splicing level. Therefore, this study aims to construct a novel method named alternative splicing network(AEN) combining gene expression levels with AS patterns to analyze cellular heterogeneity, identify novel cell types, and infer the regulatory mechanisms of alternative splicing. By applying the AEN method to gastrula data, pan-cancer T-cell data, and lung cancer tumor cell data, we have successfully revealed the hidden alternative splicing heterogeneity in these cells and identified key alternative splicing events and related splicing factors during cell transformation. The application of AEN method provides new insights into the understanding of cellular heterogeneity and its related physiological and pathological processes.

## The steps for AEnet
AENet consists of six steps: ASP Determination, ASP-Expression Network Construction, Multiple-Sample Integration, Key ASP Identification, ASP Cluster Determination, Cell Type Identification Assistance, and Regulatory Mechanism Prediction. Detailed information for each step is provided below.

### ASP Determination
The goal of the ASP Determination step is to identify alternative splicing patterns within the dataset.

Inputs: The input required for this step: A series of junctions in the dataset from DESJ-detection (https://github.com/liushang17/DESJ-detection). 

Step1: ASP Determination. Exon-exon junctions with the same starting point or endpoint were defined as alternative splicing patterns (ASP).

Outputs: The output of ASP determination is one metadata, including ASPs and the junction composition of ASPs. 

Below is an example command: 

asp <- asp(annj)  # annj is a dataframe with the first column as junctions.

### ASP-Expression Network
At the ASP-Expression Network step, AEN would detect the relationship between ASP usage preference and gene expression dynamic and construct ASP-Expression Network for each sample. 

Inputs: The Gene expression matrix with cell as column, gene as row, and expression as value, as well as the cell metadata with the cell resource information. The Cell-junction CPM matrix and the ASPs metadata from the above step are also needed. 

Step1: detect the relationship between ASPs usage preference and gene expression dynamic. Spearman correlation was calculated between the PSI values of ASPs and gene expression when the number of cells expressing non-missing values of PSI exceeded 10.

Step2: integrate the relationship across all genes and ASPs. All ASPs and genes were represented as nodes in the AEN graph (Alternative Splicing Pattern-Expression Network). ASPs were connected to their significantly associated genes (P value < 0.01), with the line width indicating the strength and direction (positive or negative) of the correlation.

Output: The output of ASP-Expression Network is the AEN network in a file with the four columns: ASP, Gene symbol, correlation coefficient, and p value.

Below is an example command:

cor_list <- step1_cor(exp,mat,asp,ann) # Construct the AEN network for one sample 

cor_list <- mutli(exp,mat,asp,ann) # Construct the AEN network for each sample

### Multiple-Samples Integration
At the multiple samples integration, we would integrate the AEN network from all samples to form the final pan-samples AEN network. 

Inputs: Multiple inputs are required for this step: the AEN network from all samples.

Step1: filter the ASP-Gene pairs in the given batch. AEN filters out ASP-gene pairs that are present in only one sample and excludes those with opposing correlations across different samples. For the remaining pairs, the strength is determined by the number of samples in which the ASP-gene pair is present. The direction of A-E links aligns with the direction observed in single samples.

Output: The output of the steps is the high quality ASP-Gene pairs, with three columns: ASP, Gene Symbol, number of samples, and direction.

Below is an example command:
corm <- merge_cor(cor_list) # Merge the AEN network across multiple samples

### The determination of anchor ASP
At the anchor ASP determination, AEN would identify the anchor ASP events based on the degree of ASPs in the network in the high quality AEN network.

Inputs: The high quality AEN network from the above steps.

Step1: The identification of the key ASP. The degree of ASP nodes in the AEN network is determined by the number of links associated with the ASP. The top 1,500 ASPs with the highest degree were identified as anchor ASPs.

Step2: the final AEN network. The anchor ASPs and their associated links were retained to construct a high-quality AEN network.

Output: The output of the steps is the final ASP-Gene pairs, with three columns: the anchor ASP, Gene Symbol, number of samples, and direction.

Below is an example command:

corm2 <- asp_selection(corm) # To construct the high quality AEN

### ASP cluster Determination
The next step aims to identify ASP clusters. The ASPs within each cluster are associated with similar phenotypes, corresponding to a specific gene set. 

Inputs: The AEN network of the anchor ASPs from the above steps.

Step1: The transformation of AEN. To begin, the graph is transformed into a matrix where gene features are represented as rows, ASPs as columns, and the correlation direction (with 1 indicating a positive correlation and -1 indicating a negative correlation) as the matrix elements.

Step2: The similarity across ASPs. The Jaccard similarity coefficient between ASPs was calculated based on their associated genes. ASPs with a higher number of shared associated genes exhibit a higher Jaccard similarity coefficient.

Step3: Clustering on the ASPs node using hierarchical clustering. Hierarchical clustering, was then performed on the jaccard similarity coefficient matrix from the previous step to group ASPs into multiple sub-clusters. ASPs within the same sub-cluster exhibit similar correlations with the phenotype. The number of sub-clusters was determined as one-fourth of the total ASP count (25 by default). The number of sub-clusters is adjustable.

Step4: ASP sub-clusters filtering. Sub-clusters are retained when more than 10% of ASP pairs within the sub-cluster have a similarity greater than 0.1.

Step5: ASP sub-clusters merging. Clusters with fewer than 10 ASPs are merged with the most similar clusters, provided that the similarity between clusters is above 0.1. The similarity between two ASP clusters is defined as the proportion of ASP pairs with a similarity score greater than 0.1, where one ASP in the pair belongs to one cluster and the other ASP belongs to the other cluster.

Output: The ASP clusters as well as the Jaccard similarity coefficient matrix.

Below is an example command:

ann_junc_list <- junction_clustering(corm2) # ASP clusters

ann_junc <- ann_junc_list$asp_clusters

asp_simm <- ann_junc_list$asp_simm

### Assistance of cell type identification
The aim is to uncover cellular splicing heterogeneity based on the ASP clusters. 

Inputs: The ASP clusters and the Cell-junction CPM (counts per million) matrix.

Step1: The enrichment score of the ASP clusters for each cell. AEN constructs an enrichment score matrix, using the enrichment score as an indicator of the extent to which a cell is influenced by the ASP clusters. The enrichment score is calculated by determining the mean PSI of ASPs within the ASP clusters across the cells.

Step2: the cell heterogeneity. AEN will apply clustering on the Cell-ASP classes enrichment matrix. This clustering is performed using the FindNeighbors and FindClusters functions in Seurat, with the ASP clusters serving as the principal component (PC). The desired level of resolution can be set as an optional parameter during the clustering process.

Step3: Visualization. The enrichment score matrix is used as input for UMAP (Uniform Manifold Approximation and Projection) using the UMAP package.

Output: The output is the cell annotation determined by alternative splicing.

Below is an example command:

asp_mat <- asp_score(mat,ann_junc) # To construct the Cell-ASP classes Enrichment Score Matrix

cell_info <- cell_clus(asp_mat,resolution = 1,min.dist = 0.1) # To clustering based on the matrix

### Regulatory mechanism prediction
The aim is to elucidate potential regulatory mechanisms between alternative splicing and gene expression. A total of 404 splicing factors was downloaded from the SpliceAidF database (http://srv00.recas.ba.infn.it/SpliceAidF/). 

Input: The ASP clusters, as well as the AEN network.

Step1: The extraction of SFs sub-AEN. AEN extracted the splicing factors as well as the associated ASPs in the AEN to form the SFs sub-AEN.

Step2: The key splicing factors. For a given ASP cluster, we further extracted the SFs-ASPs sub-AEN from the SFs sub-AEN. The degree of the SF nodes within the SFs-ASPs sub-AEN network is determined by the number of links associated with each splicing factor. The top SFs with the highest degrees were identified as the most critical factors. Additionally, the regulatory direction of splicing factors on the specified set of ASPs was assessed based on the proportion of positive and negative correlations. Specifically, if 75% or more of the splicing factor–ASP pairs were positive (or negative), the splicing factor was classified as positively (or negatively) regulating the ASP sets.

Step3: The identification of the different pathways affected by ASPs. For a single ASP event, we identified gene sets positively and negatively correlated with the event. Functional enrichment analysis was performed on the ASP gene and the positively and negatively correlated gene sets. The pathways to which the ASP gene belongs is affected by the ASP events. Metascape was used for the pathway analysis.

Output: The ranking and associated splicing factors with the specific ASP cluster.

Below is an example command:
key_sf_C1 <- key_sf(corm2,sfname,ann_junc,"C_1") # To detect the key splicing factors for each ASP classes



## install
devtools::install_github("https://github.com/liushang17/AEN")

## Demo
We showed the application of a demo in the demo.r. The dataset used in the demo.r is shown in the link (https://pan.baidu.com/s/14ZMS2e6fMk3P1CLdOCK6tw . Extraction code ：fyS9).
