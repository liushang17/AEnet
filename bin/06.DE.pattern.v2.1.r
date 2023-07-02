args = commandArgs(T)
library(edgeR)
library(limma)

####################### loading data #########################
mat <- readRDS(args[1])
colann <- readRDS(args[2])
outdir <- args[3]
st <- as.numeric(args[4])
en <- as.numeric(args[5])

######################## filter junctions #####################
ann <- data.frame(genename = gsub("_.*","",colnames(colann)),asp = colnames(colann))

tes <- data.frame(table(as.character(ann$genename)))
pos <- which(tes$Freq > 1)
tes1 <- tes[pos,]
geneall <- as.character(tes1$Var1)


########################## DE analysis ########################
for(i in st:en){
	genename <- geneall[i]
	pos <- which(ann$genename %in% genename)
	ann1 <- ann[pos,]
	num <- nrow(ann1)
	clus <- ann1$asp
	clustype <-  gsub(".*_","",clus)

  	for(n in 1:(num-1)){
    		for(m in (n+1):num){
			filename <- paste0(outdir,"/",genename)
			pos1 <- which(colnames(colann) %in% clus[n])
			pos2 <- which(as.numeric(as.character(as.matrix(colann[,pos1]))) == 2)
			ex <- as.character(rownames(colann))[pos2]

			pos1 <- which(colnames(colann) %in% clus[m])
                        pos2 <- which(as.numeric(as.character(as.matrix(colann[,pos1]))) == 2)
                        nex <- as.character(rownames(colann))[pos2]

			pos <- which(colnames(mat) %in% ex)
			matex=mat[,pos]
			pos <- which(colnames(mat) %in% nex)
			matnex=mat[,pos]
       			exl=length(colnames(matex))
       			nexl=length(colnames(matnex))
			if(exl >= 20 & nexl >= 20){
       				matexnex=cbind(matex,matnex)
       				matexnex=matexnex[rowSums(matexnex)!=0,]
       				keep <- rowSums(matexnex>0) >=0.1*(exl+nexl)
       				matexnex<-matexnex[keep,]
       				alll=length(colnames(matexnex))
       				dge <- DGEList(counts=matexnex)
       				group1=1
       				for (j in 1:exl){group1[j]="EX"}
       				for (j in (exl+1):alll){group1[j]="NEX"}
       				group<-factor(group1)
       				design <- model.matrix(~group)
       				v <- voom(dge,design,normalize="quantile")
       				fit <- lmFit(v,design)
       				fit2 <- eBayes(fit)
       				tempOutput = topTable(fit2, coef=ncol(design), n=Inf)
    
   				filename <- paste0(filename,".",clustype[n],"_",clustype[m],".xls")
       				write.table(tempOutput,file=filename,sep="\t",quote=F)
				}
			}
		}
}
