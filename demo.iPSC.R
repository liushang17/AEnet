library(AEN)
library(ggsci)
library(ggrepel)
library(ggpubr)
col <- c(pal_npg("nrc")(10)[1:6],"#FAFD7CFF","#FF6F00FF",pal_lancet("lanonc")(9)[c(1,3,7)],"#660099FF","#B5CF6BFF","#B24745FF","#CCFF00FF",
         "#FFCD00FF","#800000FF","#20854EFF","#616530FF","#FF410DFF","#EE4C97FF","#FF1463FF","#00FF00FF","#990080FF","#00FFFFFF",
         "#666666FF","#CC33FFFF","#00D68FFF","#4775FFFF","#C5B0D5FF","#FDAE6BFF","#79CC3DFF","#996600FF","#FFCCCCFF","#0000CCFF",
         "#7A65A5FF","#1A5354FF","#24325FFF")

############# input ###############
indir <- "/Users/shangliu/01.terms/02.AEN/04.result/iPSC/input/"

mat <- readRDS(paste0(indir,"/ASP.rds"))
exp <- readRDS(paste0(indir,"/exp.rds"))
ann <- read.table(paste0(indir,"/meta_data.xls"),sep = "\t")
aspann <- read.table(paste0(indir,"/junction_anno.xls"),sep = "\t")

outdir <- "/Users/shangliu/01.terms/02.AEN/04.result/iPSC/Result/"

############ Transformation ########
colnames(ann) <- c("cell","Type")
ann$Patient <- "Cell_line"
sfname <- sf$V1

########### Pre-processing ##########
annj <- asp(mat, outdir,n_cores=2) ##### ASP detection
psi_matrix <- calculate_psi(mat,annj,outdir) ##### PSI calculation

pos <- which(colnames(exp) %in% colnames(psi_matrix))
exp <- exp[,pos]

multi(as.matrix(psi_matrix),as.matrix(exp),ann,outdir,parallel = T,inner_cores = 20) ###### network construction for each sample
corm <- process_cor_results(outdir) ####### final AEN network

pos <- which(corm$p.value < 0.00001) ####### The high quality AEN network
corm <- corm[pos,]
saveRDS(corm,file = paste0(outdir,"step3.network.rds"))

corm2 <- asp_selection(corm,sample_cutoff = 0,link_cutoff = 15) ######### key ASP filter
junc_list <- junction_clustering(corm2) ####### ASP cluster
saveRDS(junc_list,file = paste0(outdir,"step4.asp.cluster.rds"))

junc_list <- readRDS(paste0(outdir,"step4.asp.cluster.rds"))

ann_junc <- junc_list$asp_clusters
corm2 <- gene_selection(corm,aspann,sample_cutoff = 0) ######### key Gene filter
exp_list <- gene_clustering(corm2) ############ Gene Cluster
saveRDS(exp_list,file = paste0(outdir,"step4.gene.cluster.rds"))

exp_list <- readRDS(paste0(outdir,"step4.gene.cluster.rds"))

exp_junc <- exp_list$asp_clusters

########### Cell Clustering ##########

asp_matrix <- asp_score(psi_matrix,ann_junc) ######## enrichment score for ASP
exp_matrix <- gene_score(as.matrix(exp),exp_junc) ########### enrichment score for Gene

res <- cell_clus(asp_matrix,exp_matrix,scale = FALSE) ########### Cell clustering
colnames(ann)[1] <- "Cell"
ann1 <- merge(ann,res,by = "Cell")
ggplot(ann1,aes(x=as.numeric(umap_1),y=as.numeric(umap_2),color = Type))+geom_point()+
  scale_color_manual(values = col) + labs(x="UMAP_1",y="UMAP_2") + theme_classic()
ggplot(ann1,aes(x=as.numeric(umap_1),y=as.numeric(umap_2),color = cluster))+geom_point()+
  scale_color_manual(values = col) + labs(x="UMAP_1",y="UMAP_2") + theme_classic()

siz <- NULL
pin <- c(0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55)
for(i in pin){
  res <- cell_clus(asp_matrix,exp_matrix,resolution = i,scale = FALSE)
  colnames(ann)[1] <- "Cell"
  ann1 <- merge(ann,res,by = "Cell")

  if (require("mclust")) {
    tmp <- adjustedRandIndex(ann1$Type,ann1$cluster)
  }
  siz <- c(siz,tmp)
}
siz

siz <- NULL
pin <- c(0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55)
for(i in pin){
  res <- cell_clus(NULL,exp_matrix,resolution = i,scale = FALSE)
  colnames(ann)[1] <- "Cell"
  ann1 <- merge(ann,res,by = "Cell")

  if (require("aricode")) {
    tmp <- NMI(ann1$Type,ann1$cluster)
  }
  siz <- c(siz,tmp)
}
siz

########### Key ASP ##########
res1 <- res[order(res$cluster),]
res2 <- data.frame(clus = res1$cluster)
rownames(res2) <- res1$Cell

rownames(asp_matrix) <- asp_matrix$Cell
asp_mat <- asp_matrix[, setdiff(colnames(asp_matrix), "Cell"), drop = FALSE]
asp_mat <- asp_mat[rownames(res2),]
pheatmap::pheatmap(asp_mat,cluster_rows = F,cluster_cols = F,annotation_row = res2,show_rownames = F,scale = "column")

###### C_11 is the key feature for Cell Cluster 0
sf <- read.table(paste0(indir,"/splicing.factor.v2.txt"),sep = "\t")
sf <- key_sf(corm,sfname,ann_junc,cluster = "C_11")

sf1 <- sf[1:10,]
sf1$SF <- factor(sf1$SF,levels = as.character(sf1$SF))
ggplot(sf1,aes(x=SF,y=importance))+geom_bar(stat = "identity")+theme_classic()

