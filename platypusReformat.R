library(vcfR)
library(stringr)
library(dplyr)


args = commandArgs(trailingOnly=TRUE)

vcf_platypus <- read.vcfR(args[1])

dat_fam = as.data.frame(vcf_platypus@gt)
fix = data.frame(vcf_platypus@fix)

platypus_dg=data.frame(vcf_platypus@fix, stringsAsFactors = FALSE)
platypus_body=data.frame(vcf_platypus@fix[,1:6], stringsAsFactors=FALSE)
platypus_dg[is.na(platypus_dg)] <- '.'
fname <- gsub('.vcf',"",args[1])

##taking TC and TR and also calculation AD, DP, RD and AF

for(i in 1:dim(platypus_dg)[1]){
  L5=strsplit(as.character(platypus_dg[i,8]), ";")

  TC = as.numeric(gsub('TC=',"",L5[[1]][15]))
  TR = gsub(',[0-9]*',"",L5[[1]][18])
  TR = as.numeric(gsub('TR=',"",TR[1]))
  DP = TC
  AD = TR
  RD = DP - AD

  AF=round((AD/DP)*100, digits = 2)

  platypus_dg[i,9:12] = c(DP,AD,RD,AF)
}

colnames(platypus_dg)[9:12] <- c("DP","AD","RD","AF")
platypus_dg[,12] = paste0(platypus_dg[,12], "%")

##we have to take the gt from datfram

for(i in 1:dim(dat_fam)[1]){
  L1=strsplit(as.character(dat_fam[i,2]), ":")
  L2=strsplit(as.character(dat_fam[i,3]), ":")
  L3=strsplit(as.character(dat_fam[i,4]), ":")
  L4=strsplit(as.character(dat_fam[i,5]), ":")

  GT = unique(c(L1[[1]][1], L2[[1]][1], L3[[1]][1], L4[[1]][1]))
  GT = paste(GT,collapse=",")

  platypus_body[i,7] = (GT[length(GT)])
}

  platypus_dg = platypus_dg %>%
    select(1:5,9:12 )

  df2 <- cbind(platypus_dg, GT = platypus_body$V7)
  df2$qual = sample(".", replace = TRUE)
  df2[,12] <- paste("GT=",df2[,10], ":DP=",df2[,6], ":RD=",df2[,7], ":AD=",df2[,8], ":AF=",df2[,9])
  df2 = df2 %>%
    select(1:5,10,11,12)
  names(df2) = NULL
  write.csv(df2, paste(fname,'rf.vcf',sep = "_"), row.names = F, quote = F, col.names = F, sep = "\t")
