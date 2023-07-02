args <- commandArgs(T)

library(pheatmap)
library("reshape2")
library("tidyr")
library(data.table)
library(Matrix)
library(umap)
library(ggsci)
col <- c(pal_npg("nrc")(10)[1:6],"#FAFD7CFF","#FF6F00FF",pal_lancet("lanonc")(9)[c(3,7)],"#660099FF","#B5CF6BFF","#B24745FF","#CCFF00FF",
         "#FFCD00FF","#800000FF","#20854EFF","#616530FF","#FF410DFF","#EE4C97FF","#FF1463FF","#00FF00FF","#990080FF","#00FFFFFF",
         "#666666FF","#CC33FFFF","#00D68FFF","#4775FFFF","#C5B0D5FF","#FDAE6BFF","#79CC3DFF","#996600FF","#FFCCCCFF","#0000CCFF",
         "#7A65A5FF","#1A5354FF","#24325FFF")

########### loading parameters
infile <- args[1]
phenotype_cluster_cutoff <- as.numeric(as.character(args[2])) # default 30
asp_cluster_cutoff <- as.numeric(as.character(args[3])) # default 30
phenotype_num_cutoff <- as.numeric(as.character(args[4])) # default 30
asp_num_cutoff <- as.numeric(as.character(args[5])) # default 30
outdir <- agrs[6]

############
ase <- read.table(infile,sep = "\t")
ase$V7 <- paste0(ase$V6,"_",ase$V4)
pos <- which(ase$V2 > 0)
ase$V7[pos] <- paste0(ase$V6,"_",ase$V5)[pos]
ase$V8 <- paste0(ase$V6,"_",ase$V4,"_",ase$V5)
colnames(ase) <- c("symbol","FC","P.value","cluster1","cluster2","Gene","asp","aspa")

pos <- which(abs(ase$FC) >= 1 & ase$P.value < 0.01)
ase1 <- ase[pos,]
ase1$FC <- abs(ase1$FC)

########### filter gene with as and de
pos <- which(ase1$Gene == ase1$symbol)
ase2 <- ase1[pos,]
pos <- which(ase1$aspa %in% ase2$aspa)
ase3 <- ase1[-pos,]

########### filter symbol
sui <- data.frame(table(ase3$symbol))
sui$num <- sui$Freq
pos <- which(sui$num >= 100)
sui$num[pos] <- 100
#plot(density(sui$num)) + abline(v= 20)
pos <- which(sui$num >= phenotype_num_cutoff)
sui1 <- sui[pos,]
pos <- which(ase3$symbol %in% sui1$Var1)
ase3 <- ase3[pos,]

######## filter asp
sui <- data.frame(table(ase3$asp))
sui$num <- sui$Freq
pos <- which(sui$num >= 100)
sui$num[pos] <- 100
#plot(density(sui$num)) + abline(v= 50)
pos <- which(sui$num >= asp_num_cutoff)
sui1 <- sui[pos,]
pos <- which(ase3$asp %in% sui1$Var1)
ase3 <- ase3[pos,]

###### trans
ase3$FC <- abs(ase3$FC)
ase3$FC[ase3$FC >= 3] <- 3
ase3$cellID<-ase3$asp
ase3$geneID <- ase3$symbol
pos <- which(ase3$symbol %in% NA)
if(length(pos) > 0){ase3 <- ase3[-pos,]}
gene=unique(ase3$geneID)
cell=unique(ase3$cellID)
gene_idx=c(1:length(gene))
cell_idx=c(1:length(cell))
names(gene_idx)=gene
names(cell_idx)=cell
ase3$num <- 1
mat=sparseMatrix(i=gene_idx[ase3$geneID],j=cell_idx[ase3$cellID],x=ase3$FC)
mit=sparseMatrix(i=gene_idx[ase3$geneID],j=cell_idx[ase3$cellID],x=ase3$num)
mit[mit == 0] <- 1

mat1 <- mat / mit
mat1 <- as.matrix(mat1)
mat1[mat1 == NA] <- 0

rownames(mat1)=gene
colnames(mat1)=cell

######## asp first clustering
p1 <- pheatmap::pheatmap(mat1,clustering_distance_rows = "correlation",clustering_distance_cols = "correlation")

num <- length(integer(nrow(mat1) / phenotype_cluster_cutoff))
rowinfo <- data.frame(cutree(p1$tree_row,k = num))
rowinfo$name <- rownames(rowinfo)
colnames(rowinfo) <- c("cluster","symbol")
sui <- data.frame(table(rowinfo$cluster))
sui$asp <- ""
sui$final <- sui$Var1

for(i in 1:nrow(sui)){
  pos <- which(rowinfo$cluster %in% sui$Var1[i])
  rowinfo1 <- rowinfo[pos,]
  pos <- which(rownames(mat1) %in% rowinfo1$symbol)
  mat3 <- mat1[pos,]
  if(length(pos) > 1){
    mat3[mat3 > 0] <- 1
    tmp <- colSums(mat3) / nrow(mat3)
    pos <- which(tmp >= 0.5)
    tmp1 <- colnames(mat3)[pos]
    sui$asp[i] <- paste(tmp1,collapse = ",")
  }else{
    pos <- which(mat3 > 0)
    tmp1 <- colnames(mat1)[pos]
    sui$asp[i] <- paste(tmp1,collapse = ",")
  }
}


pos <- which(sui$Freq >= cutoff_num)
sui1 <- sui[pos,]
sui2 <- sui[-pos,]

for(i in 1:nrow(sui2)){
  ari <- NULL
  tmp <- as.character(strsplit(sui2$asp[i],split = ",")[[1]])
  for(j in 1:nrow(sui1)){
    tmp1 <- as.character(strsplit(sui1$asp[j],split = ",")[[1]])
    pos <- which(tmp %in% tmp1)
    tmp2 <- length(pos) / (length(unique(c(as.character(tmp1),as.character(tmp)))))
    ari <- c(ari,tmp2)
  }
  
  pos <- which(ari == max(ari))
  tmp <- as.character(sui1$Var1[pos])[1]
  pos <- which(sui$Var1 %in% sui2$Var1[i])
  sui$final[pos] <- tmp
}

rowinfo$final_cluster <- rowinfo$cluster
for(i in 1:nrow(sui)){
  pos <- which(rowinfo$cluster %in% sui$Var1[i])
  rowinfo$final_cluster[pos] <- sui$final[i]
}

######## first final 
rowinfo <- rowinfo[order(as.numeric(rowinfo$final_cluster)),]
rowinfo1 <- data.frame(clus = as.character(rowinfo$final_cluster))
rownames(rowinfo1) <- rowinfo$symbol
mat2 <- mat1[rownames(rowinfo1),]
p2 <- pheatmap::pheatmap(mat2,annotation_row = rowinfo1,cluster_rows = F,cluster_cols = T,clustering_distance_cols = "correlation")

######## asp second clustering
num <- length(integer(ncol(mat2) / asp_cluster_cutoff))
colinfo <- data.frame(cutree(p2$tree_col,k = num))
colinfo$name <- rownames(colinfo)
colnames(colinfo) <- c("cluster","symbol")
colinfo1 <- colinfo[order(colinfo$cluster),]
colinfo2 <- data.frame(cluster = as.character(colinfo1$cluster))
rownames(colinfo2) <- colinfo1$symbol
mat2 <- mat2[,rownames(colinfo2)]

clus <- unique(as.character(rowinfo$final_cluster))
matinfo <- matrix(nrow = length(clus),ncol = ncol(mat2))
siz <- NULL
for(i in 1:length(clus)){
  pos <- which(rowinfo$final_cluster %in% clus[i])
  rowinfo1 <- rowinfo[pos,]
  pos <- which(rownames(mat2) %in% rowinfo1$symbol)
  mat3 <- mat2[pos,]
  mat3[mat3 > 0] <- 1
  tmp <- colSums(mat3) / nrow(mat3)
  matinfo[i,] <- tmp
  siz <- c(siz,tmp)
}
colnames(matinfo) <- colnames(mat2)
rownames(matinfo) <- clus

matinfo2 <- matinfo[,rownames(colinfo2)]

sui <- data.frame(table(colinfo2$cluster))
sui$Var1 <- as.character(sui$Var1)
pos <- which(sui$Freq <= 3)
sui1 <- sui[pos,]
sui2 <- sui[-pos,]

colinfo2$newcluster <- colinfo2$cluster
for(i in 1:nrow(sui1)){
  pos <- which(colinfo2$cluster %in% sui1$Var1[i])
  celltmp <- rownames(colinfo2)[pos]
  if(length(pos) > 1){
    pos <- which(colnames(matinfo2) %in% celltmp)
    matinfo3 <- matinfo2[,pos]
    tmp <- rowMeans(matinfo3)
    
  }else{
    pos <- which(colnames(matinfo2) %in% celltmp)
    tmp <- matinfo2[,pos]
    
  }
  tmpall <- NULL
  for(j in 1:nrow(sui2)){
    pos <- which(colinfo2$cluster %in% sui2$Var1[j])
    celltmp <- rownames(colinfo2)[pos]
    
      pos <- which(colnames(matinfo2) %in% celltmp)
      matinfo3 <- matinfo2[,pos]
      tmp1 <- rowMeans(matinfo3)
      tmp2 <- cor(tmp,tmp1)
      tmpall <- c(tmpall,tmp2)
    
  }
  pos1 <- which(tmpall == max(tmpall))
  pos2 <- which(colinfo2$cluster %in% sui1$Var1[i])
  colinfo2$newcluster[pos2] <- (sui2$Var1)[pos1]
}

clus_col <- unique(colinfo2$newcluster)
matinfo_clus <- matrix(nrow = nrow(matinfo2),ncol = length(clus_col))
for(i in 1:length(clus_col)){
  pos <- which(colinfo2$newcluster %in% clus_col[i])
  celltmp <- rownames(colinfo2)[pos]
  
    pos <- which(colnames(matinfo2) %in% celltmp)
    matinfo3 <- matinfo2[,pos]
    tmp <- rowMeans(matinfo3)
    matinfo_clus[,i] <- tmp
  
}
rownames(matinfo_clus) <- rownames(matinfo2)
colnames(matinfo_clus) <- clus_col

matinfo_clus1 <- cor(matinfo_clus)
pheatmap::pheatmap(matinfo_clus1,display_numbers = T)

matinfo_clus1[matinfo_clus1 >= 0.5] <- 1
matinfo_clus1[matinfo_clus1 < 0.5] <- 0

colinfo2$final_cluster <- colinfo2$newcluster
siztmp <- NULL
for(i in 1:nrow(matinfo_clus1)){
  tmp <- as.numeric(as.character(matinfo_clus1[i,]))
  pos <- which(tmp > 0)
  tmp1 <- paste(colnames(matinfo_clus1)[pos], collapse= ",")
  siztmp <- c(siztmp,tmp1)
}
tmpinfo <- data.frame(cluster = rownames(matinfo_clus1),type = siztmp)
tmpinfo1 <- data.frame(table(tmpinfo$type))
pos <- which(tmpinfo1$Freq >= 2)
tmpinfo2 <- tmpinfo1[pos,]
for(i in 1:nrow(tmpinfo2)){
  pos <- which(tmpinfo$type %in% tmpinfo2$Var1[i])
  tmpinfo3 <- tmpinfo[pos,]
  pos <- which(colinfo2$newcluster %in% tmpinfo3$cluster)
  colinfo2$final_cluster[pos] <- tmpinfo3$cluster[1]
}

colinfo2 <- colinfo2[order(colinfo2$final_cluster),]
colinfo3 <- data.frame(cluster = as.character(colinfo2$final_cluster))
rownames(colinfo3) <- rownames(colinfo2)
mat2 <- mat2[,rownames(colinfo3)]
rowinfo1 <- rowinfo[order(rowinfo$final_cluster),]
rowinfo2 <- data.frame(clus = as.character(rowinfo1$final_cluster))
rownames(rowinfo2) <- rowinfo1$symbol
mat2 <- mat2[rownames(rowinfo2),]

pdf(paste0(outdir,"/Phenotype.ASP.pdf"))
pheatmap::pheatmap(mat2,annotation_row = rowinfo2,annotation_col = colinfo3,cluster_rows = F,cluster_cols = F)
dev.off()

rowinfo2$features <- rownames(rowinfo2)
colinfo2$features <- rownames(colinfo2)

colnames(rowinfo2) <- c("merge_clusters","features")
colnames(colinfo2) <- c("merge_clusters","features")

write.table(rowinfo2,file = paste0(outdir,"/Phenotype.v1.xls"),sep = "\t",quote = F)
write.table(colinfo2,file = paste0(outdir,"/ASP.v1.xls",sep = "\t",quote=F))

tmpmit <- data.frame(asp = rownames(colinfo2),gene = gsub("_.*","",rownames(colinfo2)))
genename <- unique(as.character(tmpmit$gene))
write.table(genename,file = paste0(outdir,"/Key.gene.xls"),sep = "\t",quote=F,row.names = F,col.names = F)


