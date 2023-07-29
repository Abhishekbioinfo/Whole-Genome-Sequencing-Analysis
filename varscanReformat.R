library(vcfR)
library(stringr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

vcf_varscan <- read.vcfR(args[1])
dat_fam = as.data.frame(vcf_varscan@gt)

varscan_dg=data.frame()
varscan_dg=data.frame(vcf_varscan@fix[,1:5])
varscan_dg[is.na(varscan_dg)] <- "."
fname <- gsub('.vcf',"",args[1])

for(i in 1:dim(dat_fam)[1]){
  L1=strsplit(as.character(dat_fam[i,3]), ":")
  
  GT = L1[[1]][1]
  
  RD = as.numeric(L1[[1]][4])
  
  AD = as.numeric(L1[[1]][5])
  
  DP = as.numeric(L1[[1]][3])
  
  VAF=(AD/DP)*100
  
  format <- paste("GT=",GT,";DP=",DP,";RD=",RD,";AD=",AD,";AF=",VAF, sep = "")
  
  varscan_dg[i,6:8] = c(GT,".",format)
}

write.table(varscan_dg, paste(fname,'rf.vcf',sep = "_") , row.names = F,
            quote = F, col.names = F, sep = "\t")

