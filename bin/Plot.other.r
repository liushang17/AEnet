args = commandArgs(T)

mat <- readRDS(args[1])
exp <- read.table(args[2],header = T,sep = "\t")
clus <- as.character(strsplit(args[3],split = ",")[[1]])
gene <- as.character(args[4])
ann <- read.table(args[5],sep = "\t")
outdir <- args[6]

ann$V7 <- paste0(ann$V6,"_",ann$V3)
pos <- which(ann$V6 %in% gene)
ann1 <- ann[pos,]
junc <- NULL
for(i in 1:nrow(ann1)){
  tmp <- as.character(strsplit(ann1$V1[i],split = ",")[[1]])
  junc <- c(junc,tmp)
}
pos <- which(rownames(mat) %in% junc)
mat1 <- mat[pos,]

pa <- rep("Pattern",nrow(exp1))
for(i in 1:nrow(ann1)){
  tmp <- exp[,ann1$V7[i]]
  pos <- which(tmp == 2)
  pa[pos] <- ann1$V7[i]
}
cluinfo <- data.frame(cell = rownames(exp1),clus = pa)
pos <- which(cluinfo$clus %in% "Pattern")
cluinfo1 <- cluinfo[-pos,]
cluinfo2 <- cluinfo1[order(cluinfo1$clus),]
cluinfo3 <- data.frame(clus = cluinfo2$clus)
rownames(cluinfo3) <- cluinfo2$cell

mat2 <- mat1[,cluinfo2$cell]

pdf(paste0(outdir,"/",gene,".pattern.pdf"))
pheatmap::pheatmap(log2(mat2+1),cluster_row = T,cluster_col = F, show_rownames = T,show_colnames=F,annotation_col = cluinfo3)
dev.off()

