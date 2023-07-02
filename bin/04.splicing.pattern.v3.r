args = commandArgs(T)

############# loading library and parameters
library(pheatmap)
library(amap)
library(edgeR)
library(limma)
source('../filter.cell.r')

mincell <- 10
maxsd <- 0.2
maxmean <- 0.2
fc_cut <- 0.1
sim_cut <- 0.1

## loading data
mat <- readRDS(args[1])
rowann <- read.table(args[2],sep = "\t")


############# gene names
genese <- data.frame(table(as.character(rowann$V3)))
pos <- which(genese$Freq >= 2)
genese_f <- genese[pos,]
geneall <- genese_f$Var1

############# pattern 
for(n in 1:length(geneall)){
  genename <- geneall[n]
  print(as.character(genename))
  pos <- which(rowann$V3 %in% genename)
  rowann1 <- rowann[pos,]
  pos <- which(rownames(mat) %in% rowann1$V1)
  mattmp <- mat[pos,]
  
  ########## Step1: QC
  tes <- filter(mattmp,mincell,maxsd,maxmean)
  sui <- tes
  numfilter <- sui$cl
  cluallinfo <- sui$cluinfo
  cluallinfo$clu = 'filter'
  pos <- which(cluallinfo$suitmp.cluster != numfilter)
  cluallinfo$clu[pos] <- cluallinfo$suitmp.cluster[pos]
  mit <- cluallinfo
  mit1 <- mit[order(mit$suitmp.cluster),]
  mit2 <- data.frame(clu = mit1$clu)
  rownames(mit2) <- rownames(mit1)
  mattmp1 <- mattmp[,rownames(mit2)]
  
  mit2$cluster <- "A"
  pos <- which(mit2$clu == "filter")
  mit2 <- mit2[-pos,]
  pos <- which(colnames(mattmp1) %in% rownames(mit2))
  mattmp2 <- mattmp1[,pos]
  
  ########## Step2: PSI value
  if(length(pos) > 3){
  tmp <- colSums(mattmp2)
  pos <- which(tmp == 0)
  if(length(pos) > 0){
    mattmp2 <- mattmp2[,-pos]
    tmp <- tmp[-pos]
  }
  
  mattmp3 <- t(t(mattmp2) / tmp)
  
  pos <- which(rownames(mit2) %in% colnames(mattmp3))
  mit2 <- mit2[pos,]
  
  mit2 <- mit2[order(mit2$cluster),]
  mattmp3 <- mattmp3[,rownames(mit2)]
  
  ######### Step3: Clustering based on hierarchical clustering
  cell_dist <- Dist(t(mattmp3),method = "correlation")
  cell_clus <- hclust(cell_dist)
  
  cor_cutoff <- floor(ncol(mattmp3) / 30) + 1
  if(cor_cutoff > 20){cor_cutoff <- 20}
  
  cell_clus_info <- data.frame(cutree(cell_clus,k=cor_cutoff))
  cell_clus_info <- data.frame(cell = rownames(cell_clus_info),cluster = cell_clus_info[,1])
  
  mit <- cell_clus_info
  mit1 <- mit[order(mit$cluster),]
  mit2 <- data.frame(clu = as.character(mit1$cluster),cluster = "A")
  rownames(mit2) <- mit1$cell
  
  sui <- data.frame(table(mit2$clu))
  pos <- which(sui$Freq > 3)
  sui <- sui[pos,]
  pos <- which(mit2$clu %in% sui$Var1)
  mit2 <- mit2[pos,]
  
  mattmp4 <- mattmp3[,rownames(mit2)]
  
  ######### Step4: similarity analysis
  mit2$new_clus <- mit2$clu
  sui <- data.frame(table(mit2$clu))
  sui$new_clus <- sui$Var1
  pos <- which(sui$Freq < 10)
  if(length(pos) > 0){
    sui1 <- sui[pos,]
    
    matmerge <- matrix(ncol = nrow(sui),nrow = nrow(mattmp4))
    for(i in 1:nrow(sui)){
      pos1 <- which(mit2$clu %in% sui$Var1[i])
      pos2 <- which(colnames(mattmp4) %in% rownames(mit2)[pos1])
      mattmp5 <- mattmp4[,pos2]
      if(length(pos2) > 1){
        matmerge[,i] <- rowMeans(mattmp5)
      }else{
        matmerge[,i] <- (mattmp5)
      }
    }
    colnames(matmerge) <- sui$Var1
    matcor <- cor((matmerge),method = "pearson")
    matcor[matcor %in% NA] <- 1
    
    rownames(matcor) <- sui$Var1
    colnames(matcor) <- sui$Var1
    
    matcor[matcor == 1] <- 0
    for(i in 1:nrow(sui1)){
      pos <- which(rownames(matcor) %in% sui1$Var1[i])
      tmp <- matcor[pos,]
      pos1 <- which(tmp == max(tmp))
      
      if(length(pos1)  > 0 & max(tmp) > 0.1){
        pos <- which(mit2$clu %in% sui1$Var1[i])
        pos2 <- which(mit2$clu %in% rownames(matcor)[pos1])
        mit2$new_clus[pos] <- mit2$new_clus[pos2[1]]
      }
    }
  }
  
  ######### Step5: merge clusters by difference
  mitsui <- data.frame(table(mit2$new_clus))
  mit2$final_clus <- mit2$new_clus
  pos <- which(mitsui$Freq >= 10)
  clus <- as.character(mitsui$Var1)[pos]
  if(length(pos) >= 2){
    pos <- which(mit2$new_clus %in% mitsui$Var1[pos])
    mit2 <- mit2[pos,]
    mitz <- matrix(nrow = length(clus),ncol = length(clus))
    mitz[,] <- 0
  
    for(i in 1:length(clus)){
      for(j in 1:length(clus)){
        if(i != j){
        pos1 <- which(mit2$new_clus %in% clus[i])
        pos2 <- which(colnames(mattmp4) %in% rownames(mit2)[pos1])
        mattmp5 <- mattmp4[,pos2]
        
        pos1 <- which(mit2$new_clus %in% clus[j])
        pos2 <- which(colnames(mattmp4) %in% rownames(mit2)[pos1])
        mattmp6 <- mattmp4[,pos2]
        
        matex <- mattmp5
        matnex <- mattmp6
        exl <- ncol(mattmp5)
        nexl <- ncol(mattmp6)
        
        tempOutput <- matrix(nrow = nrow(mattmp5),ncol = 3)
      	for(t in 1:nrow(mattmp5)){
          tmps1 <- mattmp5[t,]
          tmps2 <- mattmp6[t,]
	  if(sd(c(tmps1,tmps2))== 0){
	    tmps4 = 1
          }else{
	    if(sd(tmps1) == 0 & sd(tmps2) == 0){
		tmps4 = 0
	    }else{
		tmps3 <- t.test(tmps1,tmps2)
		tmps4 = tmps3$p.value
	    }
	  }
          tempOutput[t,1] <- rownames(mattmp5)[t]
          tempOutput[t,2] <- mean(tmps2) - mean(tmps1)
          tempOutput[t,3] <- tmps4
        }
        tempOutput <- data.frame(tempOutput)
        colnames(tempOutput) <- c("chr","logFC","adj.P.Val")
        tempOutput$logFC <- as.numeric(as.character(tempOutput$logFC))
        tempOutput$adj.P.Val <- as.numeric(as.character(tempOutput$adj.P.Val)) 
        
        pos <- which(abs(tempOutput$logFC) >= fc_cut & tempOutput$adj.P.Val < 0.001 )
        if(length(pos) > 0){
          mitz[i,j] <- 0
        }else{
          mitz[i,j] <- 1
        }
      }else{
        mitz[i,j] <- 1
      }
    }
  }
  
  rownames(mitz) <- clus
  colnames(mitz) <- clus
  }
  
  if(length(clus) >= 2){
    for(i in 2:length(clus)){
      for(j in 1:(i-1)){
        if(mitz[i,j] == 1){
          tmp <- mitz[i,] - mitz[j,]
          tmp1 <- sum(abs(tmp)) / length(clus)
          
          if(tmp1 < sim_cut){
            pos <- which(mit2$new_clus %in% clus[j])
            pos1 <- which(mit2$new_clus %in% clus[i])
            mit2$final_clus[pos1] <- mit2$final_clus[pos[1]]
          }
        }
      }
    }
  }
  
  mit2 <- mit2[,c("final_clus","new_clus","clu")]
  pos <- which(colnames(mattmp) %in% rownames(mit2))
  cellname <- colnames(mattmp)[-pos]
  mit3 <- data.frame(final_clus = rep(0,length(cellname)),new_clus = rep(0,length(cellname)), clu = rep(0,length(cellname)))
  rownames(mit3) <- cellname
  mit4 <- rbind(mit2,mit3)
  mit4 <- mit4[colnames(mattmp),]

  write.table(mit4,file = paste0(args[3],"/",genename,".pattern.xls"),sep = "\t",quote=F)

  
  annmit <- data.frame(table(mit4$final_clus))
  pos1 <- which(annmit$Freq >= 10 & annmit$Var1 != 0)
  if(length(pos1) > 0){
    mat1 <- matrix(nrow=nrow(mit4),ncol = length(pos1)+1)
    mat1[,1] <- rownames(mit4)
    for(s in 1:length(pos1)){
        mat1[,s+1] <- 0
	pos <- which(rownames(mit4) %in% rownames(mit3))
	mat1[-pos,s+1] <- 1
	pos <- which(mit4$final_clus %in% annmit$Var1[pos1[s]])
	mat1[pos,s+1] <- 2
    }
    colnames(mat1) <- c("cell",paste0(genename,"_",annmit$Var1[pos1]))
    write.table(mat1,file = paste0(args[3],"/",genename,".pattern.v1.xls"),sep = "\t",quote=F,row.names=F,col.names=T)
  }
  }
}

  








