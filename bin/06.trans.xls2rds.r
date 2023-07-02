args = commandArgs(T)

infile <- args[1]

mat <- read.table(inflile,sep = "\t",header = T,row.names= 1)
outfile <- gsub(".csv|.xls",".rds",infile)
saveRDS(mat,file = outfile)