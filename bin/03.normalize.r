args = commandArgs(T)

exp_mat <- read.table(args[1],header = T)
cellinfo <- read.table(args[2],sep = "\t")
outdir <- args[3]

colnames(cellinfo) <- c("cell","Frag")

exp_mat1 <- merge(exp_mat,cellinfo,by = "cell")

for(i in 2:(ncol(exp_mat1)-1)){
	exp_mat1[,i] <- exp_mat1[,i] / exp_mat1[,ncol(exp_mat1)] * 1000000
}

saveRDS(exp_mat1,file = paste0(outdir,"/Normalize.step1.rds"))

exp_mat2 <- exp_mat1[,2:(ncol(exp_mat1)-1)]
rownames(exp_mat2) <- exp_mat1$cell
saveRDS(exp_mat2,file = paste0(outdir,"/Normalize.step2.rds"))

exp_mat3 <- t(exp_mat2)
saveRDS(exp_mat3,file = paste0(outdir,"/Normalize.step3.rds"))
