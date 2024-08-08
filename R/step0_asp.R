#' @title ASP
#' @description The function to detect the ASP (Alternative Splicing Pattern)
#'
#' @param desj a data frame containing the junction used in the samples with junction as the first column
#'
#' @return a data frame containing the information of ASP as well as the composition of ASP
#' @export
#'
#' @examples asp <- asp(desj)

asp <- function(desj){
  annj <- desj
  colnames(annj)[1] <- "V1"
  annj <- annj[order(annj$V1),]
  annj$chr <- 0
  annj$st <- 0
  annj$en <- 0
  annj$site <- 0

  annj <- unique(annj)
  for(i in 1:nrow(annj)){
    tmp <- as.character(strsplit(annj$V1[i],split = "_")[[1]])
    annj$chr[i] <- tmp[1]
    annj$st[i] <- tmp[2]
    annj$en[i] <- tmp[3]
    annj$site[i] <- tmp[4]
  }

  annj$chr_st <- paste0(annj$chr,"_",annj$st,"_",annj$site)
  annj$chr_en <- paste0(annj$chr,"_",annj$en,"_",annj$site)

  sui <- data.frame(table(annj$chr_st))
  pos <- which(sui$Freq > 1)
  sui <- sui[pos,]

  mit <- NULL
  for(i in 1:nrow(sui)){
    pos <- which(annj$chr_st %in% sui$Var1[i])
    annj1 <- annj[pos,]
    for(m in 1:(nrow(annj1)-1)){
      for(n in (m+1):nrow(annj1)){
        mitt <- data.frame(all = paste0(annj1$V1[m],".",annj1$V1[n]),junction1 = annj1$V1[m],junction2 = annj1$V1[n])
        mit <- rbind(mit,mitt)
      }
    }
  }

  sui <- data.frame(table(annj$chr_en))
  pos <- which(sui$Freq > 1)
  sui <- sui[pos,]

  for(i in 1:nrow(sui)){
    pos <- which(annj$chr_en %in% sui$Var1[i])
    annj1 <- annj[pos,]
    for(m in 1:(nrow(annj1)-1)){
      for(n in (m+1):nrow(annj1)){
        mitt <- data.frame(all = paste0(annj1$V1[m],".",annj1$V1[n]),junction1 = annj1$V1[m],junction2 = annj1$V1[n])
        mit <- rbind(mit,mitt)
      }
    }
  }
  return(mit)
}

