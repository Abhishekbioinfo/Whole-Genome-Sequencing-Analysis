library(vcfR)
library(stringr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
vcf_lofreq <- read.vcfR(args[1])
dat_fam = as.data.frame(vcf_lofreq@fix)

lofreq_dg=data.frame()
lofreq_dg=data.frame(vcf_lofreq@fix[,1:5])
lofreq_dg[is.na(lofreq_dg)] <- "."
fname <- gsub('.vcf',"",args[1])

for(i in 1:dim(dat_fam)[1]){
  L1=strsplit(as.character(dat_fam[i,8]), ";")
  
  GT = "."
  
  DP4 = strsplit(gsub('DP4=',"",L1[[1]][4]), ",")

  RD = sum(as.numeric(DP4[[1]][1]),as.numeric(DP4[[1]][2]))
  
  AD = sum(as.numeric(DP4[[1]][3]),as.numeric(DP4[[1]][4]))
  
  DP = sum(RD,AD)
  
  VAF=(AD/DP)*100
  
  format <- paste("GT=",GT,";DP=",DP,";RD=",RD,";AD=",AD,";AF=",VAF, sep = "")
  
  lofreq_dg[i,6:8] = c(GT,".",format)
}

write.table(lofreq_dg, paste(fname,'rf.vcf',sep = "_"), row.names = F, quote = F, col.names = F, sep = "\t")

