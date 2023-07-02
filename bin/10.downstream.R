args = commandArgs(T)

##################
library(Seurat)
library(ggplot2)
library(umap)
library(ggplot2)
library(ggrepel)
library(ggsci)
col <- c(pal_npg("nrc")(10)[1:6],"#FAFD7CFF","#FF6F00FF",pal_lancet("lanonc")(9)[c(3,7)],"#660099FF","#B5CF6BFF","#B24745FF","#CCFF00FF",
         "#FFCD00FF","#800000FF","#20854EFF","#616530FF","#FF410DFF","#EE4C97FF","#FF1463FF","#00FF00FF","#990080FF","#00FFFFFF",
         "#666666FF","#CC33FFFF","#00D68FFF","#4775FFFF","#C5B0D5FF","#FDAE6BFF","#79CC3DFF","#996600FF","#FFCCCCFF","#0000CCFF",
         "#7A65A5FF","#1A5354FF","#24325FFF")

######### loading parameters
in_file <- args[1]
sf_file <- args[2]
phenotype_file <- args[3]
asp_file <- args[4]
outdir <- args[5]
########## loading file
ase <- read.table(in_file,sep = "\t")
sf <- read.table(sf_file,sep = ",",head = T)
rowinfo <- read.table(phenotype_file,sep = "\t",header = T,row.names = 1)
colinfo <- read.table(asp_file,sep = "\t")

####### pre-processing
ase$V7 <- paste0(ase$V6,"_",ase$V4)
pos <- which(ase$V2 > 0)
ase$V7[pos] <- paste0(ase$V6,"_",ase$V5)[pos]
ase$V8 <- paste0(ase$V6,"_",ase$V4,"_",ase$V5)
colnames(ase) <- c("symbol","FC","P.value","cluster1","cluster2","Gene","asp","aspa")

pos <- which(abs(ase$FC) >= 1 & ase$P.value < 0.01)
ase1 <- ase[pos,]
ase1$FC <- abs(ase1$FC)

pos <- which(ase1$Gene == ase1$symbol)
ase2 <- ase1[pos,]
pos <- which(ase1$aspa %in% ase2$aspa)
ase3 <- ase1[-pos,]

##########
pos <- which(ase3$symbol %in% sf$Gene)
ase3 <- ase3[pos,]
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

##########
pos <- which(rownames(colinfo) %in% ase3$asp)
colinfo <- colinfo[pos,]
colinfo1 <- colinfo[order(colinfo$merge_cluster),]
colinfo2 <- data.frame(asp_cluster = as.character(colinfo1$merge_cluster))
rownames(colinfo2) <- rownames(colinfo1)
mat2 <- mat1[,rownames(colinfo2)]

pdf(paste0(outdir,"/SF.heatmap.pdf"))
pheatmap::pheatmap(mat2,cluster_rows = F,cluster_cols = F,show_rownames = T,show_colnames = F,annotation_col = colinfo2)
dev.off()


############

aspmat <- readRDS("/Users/shangliu/01.terms/02.AEN/02.result/02.Neu/Neu.cell.chose.isoform.new.rds")

clusters <- unique(as.character(colinfo$merge_cluster))
aspmat1 <- matrix(nrow = nrow(aspmat),ncol = length(clusters))
for(i in 1:length(clusters)){
  pos <- which(colinfo$merge_cluster %in% clusters[i])
  colinfo3 <- colinfo[pos,]
  pos <- which(colnames(aspmat) %in% rownames(colinfo3))
  aspmat_tmp <- aspmat[,pos]
  aspmat_tmp[aspmat_tmp < 2] <- 0
  aspmat_tmp[aspmat_tmp == 2] <- 1
  patt_tmp <- data.frame(asp = colnames(aspmat_tmp),gene = gsub("_.*","",colnames(aspmat_tmp)))
  genetmp <- unique(as.character(patt_tmp$gene))
  tmp <- rowSums(aspmat_tmp)
  tmp <- tmp / length(genetmp)
  aspmat1[,i] <- tmp
}
rownames(aspmat1) <- rownames(aspmat)
colnames(aspmat1) <- clusters


#######
aspmat2 <- aspmat1
colnames(aspmat2) <- paste0("PC_",colnames(aspmat2))
write.table(aspmat2,file = paste0(outdir,"/Enrichment.xls"),sep = "\t",quote=F)

###########
tms <- FindNeighbors(aspmat2)
tms1 <- FindClusters(tms$snn,resolution = 2)
colnames(tms1) <- "cluster"

tsne_out <- umap((aspmat2))
umapinfo <- data.frame(tsne_out$layout)
umapinfo <- umapinfo[rownames(tms1),]
tms1$umap_1 <- umapinfo$X1
tms1$umap_2 <- umapinfo$X2
write.table(tms1,file = paste0(outdir,"/UMAP.meta.xls"),sep = "\t",quote=F)


pdf(paste0(outdir,"/UMAP.pdf"))
ggplot(tms1,aes(x=as.numeric(umap_1),y=as.numeric(umap_2),color = cluster))+geom_point()+
  scale_color_manual(values = col) + labs(x="UMAP_1",y="UMAP_2") + theme_classic()
dev.off()


