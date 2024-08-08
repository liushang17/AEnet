library(AEN)

setwd("/Users/shangliu/01.terms/02.AEN/02.result/11.method/AEN/AEN/tests/")
############# input ###############
mat <- readRDS("/Users/shangliu/01.terms/02.AEN/02.result/11.method/AEN/AEN/tests/Normalize.step3.rds")
exp <- readRDS("/Users/shangliu/01.terms/02.AEN/02.result/11.method/AEN/AEN/tests/Neu.exp.newname.rds")
ann <- read.table("/Users/shangliu/01.terms/02.AEN/02.result/11.method/AEN/AEN/tests/Neu.cluster.xls",sep = "\t")
sf <- read.table("/Users/shangliu/01.terms/02.AEN/00.basic/splicing.factor.v2.txt",sep = "\t")

############ Transformation ########
colnames(ann) <- c("cell","Type")
ann$Patient <- "Cell_line"
annj <- data.frame(V1 = rownames(mat),type = "junction")
sfname <- sf$V1
############ Running ###############
asp <- asp(annj) # First steps to detect the ASP patterns within the dataset
cor_list <- step1_cor(exp,mat,asp,ann) # Construct the AEN network for each sample
corm <- merge_cor(cor_list) # Merge the AEN network across multiple samples
corm2 <- asp_selection(corm) # To construct the high quality AEN
ann_junc_list <- junction_clustering(corm2) # to find ASP clusters
ann_junc <- ann_junc_list$asp_clusters
asp_simm <- ann_junc_list$asp_simm
asp_mat <- asp_score(mat,ann_junc) # To construct the Cell-ASP Clusters Enrichment Score Matrix
cell_info <- cell_clus(asp_mat,resolution = 1,min.dist = 0.1) # To clustering based on the matrix
key_sf_C1 <- key_sf(corm2,sfname,ann_junc,"C_1") # To detect the key splicing factors for each ASP Clusters

########## Save ####################
saveRDS(asp,file = "./Result/ASP.list.rds")
saveRDS(corm,file = "./Result/AEN.network.rds")
saveRDS(corm2,file = "./Result/AEN.network.high.rds")
saveRDS(asp_mat,file = "./Result/ASP.Enrichment.rds")
saveRDS(cell_info,file = "./Result/Cell.info.rds")
saveRDS(ann_junc_list,file = "./Result/ASP.cluster.rds")

####### Plot ######################
### UMAP
cell_info$cell <- rownames(cell_info)
cell_info1 <- merge(cell_info,ann,by = "cell")
library(ggplot2)
ggplot(cell_info1,aes(x=umap_1,y=umap_2,colour = cluster))+geom_point()+theme_classic()
ggplot(cell_info1,aes(x=umap_1,y=umap_2,colour = Type))+geom_point()+theme_classic()

### ASP Clusters
ann_junc$pattern <- ann_junc$PC
pos <- which(ann_junc$pattern %in% "Other")
if(length(pos) > 0){ann_junc <- ann_junc[-pos,]}

ann_junc <- ann_junc[order(ann_junc$pattern),]
ann_junc1 <- data.frame(clus = ann_junc$pattern)
rownames(ann_junc1) <- ann_junc$symbol
asp_simm1 <- asp_simm[rownames(ann_junc1),rownames(ann_junc1)]
pheatmap::pheatmap(asp_simm1,annotation_row = ann_junc1,annotation_col = ann_junc1,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F)

### ASP Clusters-Cell Enrichment Matirx
cell_info2 <- cell_info1[order(cell_info1$cluster,cell_info1$Type),]
cell_info3 <- data.frame(asp_cluster = cell_info2$cluster,gene_cluster = cell_info2$Type)
rownames(cell_info3) <- cell_info2$cell
asp_mat1 <- asp_mat[,2:ncol(asp_mat)]
rownames(asp_mat1) <- asp_mat$Cell
asp_mat2 <- asp_mat1[rownames(cell_info3),]

paletteLength <- 100
myColor <- colorRampPalette(c("lightblue", "white","red"))(paletteLength)
myBreaks <- c(seq(min(asp_mat2), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(asp_mat2)/paletteLength, max(asp_mat2), length.out=floor(paletteLength/2)))

pheatmap::pheatmap((asp_mat2), color=myColor, breaks=myBreaks,annotation_row = cell_info3,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = T)

### ASP Clusters
corm3 <- ann_junc
corm3$junction <- corm3$symbol
corm3$junction1 <- 0
corm3$junction2 <- 0
corm3$junction1.length <- 0
corm3$junction2.length <- 0
corm3$gene <- 0
for(i in 1:nrow(corm3)){
  tmp <- strsplit(corm3$junction[i],split = "\\.")[[1]]
  tmp1 <- strsplit(tmp[1],split = "_")[[1]]
  tmp2 <- as.numeric(as.character(tmp1[3])) - as.numeric(as.character(tmp1[2]))
  tmp3 <- strsplit(tmp[2],split = "_")[[1]]
  tmp4 <- as.numeric(as.character(tmp3[3])) - as.numeric(as.character(tmp3[2]))

  corm3$junction1[i] <- tmp[1]
  corm3$junction2[i] <- tmp[2]
  corm3$junction1.length[i] <- tmp2
  corm3$junction2.length[i] <- tmp4
}
pos <- which(corm3$PC %in% c("C_11","C_12","C_5","C_8","C_3"))
corm4 <- corm3[pos,]
t.test((corm4$junction1.length),(corm4$junction2.length),alternative = "greater")
library(ggplot2)
ann <- data.frame(type = rep(c("iPSC","Other"),each = nrow(corm4)),value = c(corm4$junction1.length,corm4$junction2.length))
ggplot(ann,aes(x=type,y=log2(value)))+geom_boxplot()+theme_classic()

###############
pos <- which(ann_junc$pattern %in% "C_5")
ann_junc1 <- ann_junc[pos,]

i <- 32
juncs <- strsplit(as.character(ann_junc1$symbol[i]),split = "\\.")[[1]]
pos <- which(rownames(mat) %in% juncs)
mat1 <- mat[pos,]
tmp <- colSums(mat1)
tmp1 <- as.numeric(as.character(as.matrix(mat1[1,])))
psi <- data.frame(cell = colnames(mat),value = tmp1 / tmp)
pos <- which(psi$value %in% NaN)
if(length(pos) > 0){psi <- psi[-pos,]}
psi2 <- merge(psi,cell_info1,by = "cell")
library(ggplot2)
ggplot(psi2,aes(x=Type,y = value))+geom_boxplot()

