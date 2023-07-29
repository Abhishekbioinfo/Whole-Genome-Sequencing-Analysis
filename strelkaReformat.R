library(vcfR)
library(stringr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

vcf_strelka <- read.vcfR("test.vcf")
dat_fam = as.data.frame(vcf_strelka@gt)

strelka_dg=data.frame()
strelka_dg=data.frame(vcf_strelka@fix[,1:5])
strelka_dg[is.na(strelka_dg)] <- "."
fname <- gsub('.vcf.gz',"","test.vcf")



for(i in 1:dim(dat_fam)[1]){
  
  if(vcf_strelka@fix[i,7] == "LowGQX;NoPassedVariantGTs") {
  
  L1=strsplit(as.character(dat_fam[i,2]), ":")
  
  GT = L1[[1]][1]
  
  RD = as.numeric(gsub(',[0-9]*',"",L1[[1]][6]))
  
  AD = as.numeric(gsub('[0-9]*,',"",L1[[1]][6]))
  
  DP = RD+AD
  
  VAF=(AD/DP)*100
  
  format <- paste("GT=",GT,";DP=",DP,";RD=",RD,";AD=",AD,";AF=",VAF, sep = "")
  
  strelka_dg[i,6:8] = c(GT,".",format)

  }
  
  else {
    
    L1=strsplit(as.character(dat_fam[i,2]), ":")
    
    GT = L1[[1]][1]
    
    RD = as.numeric(gsub(',[0-9]*',"",L1[[1]][5]))
    
    AD = as.numeric(gsub('[0-9]*,',"",L1[[1]][5]))
    
    DP = AD+RD
    
    VAF=(AD/DP)*100
    
    format <- paste("GT=",GT,";DP=",DP,";RD=",RD,";AD=",AD,";AF=",VAF, sep = "")
    
    strelka_dg[i,6:8] = c(GT,".",format)
    
  }
  
}
  
write.table(strelka_dg, paste(fname,'rf.vcf',sep = "_"), row.names = F, quote = F, col.names = F, sep = "\t")

