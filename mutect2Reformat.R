library(vcfR)
library(stringr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

vcf_mutect <- read.vcfR(args[1])
dat_fam = as.data.frame(vcf_mutect@gt)

mutect_dg=data.frame()
mutect_dg=data.frame(vcf_mutect@fix[,1:5])
mutect_dg[is.na(mutect_dg)] <- "."
fname <- gsub('.vcf',"",args[1])

for(i in 1:dim(dat_fam)[1]){
    L1=strsplit(as.character(dat_fam[i,2]), ":")
    L2=strsplit(as.character(dat_fam[i,3]), ":")
    L3=strsplit(as.character(dat_fam[i,4]), ":")
    L4=strsplit(as.character(dat_fam[i,5]), ":")

    GT = unique(c(L1[[1]][1], L2[[1]][1], L3[[1]][1], L4[[1]][1]))

    RD = sum(as.numeric(gsub(',[0-9]*',"",L1[[1]][2])),
             as.numeric(gsub(',[0-9]*',"",L2[[1]][2])),
             as.numeric(gsub(',[0-9]*',"",L3[[1]][2])),
             as.numeric(gsub(',[0-9]*',"",L4[[1]][2])))

    AD = sum(as.numeric(gsub('[0-9]*,',"",L1[[1]][2])),
             as.numeric(gsub('[0-9]*,',"",L2[[1]][2])),
             as.numeric(gsub('[0-9]*,',"",L3[[1]][2])),
             as.numeric(gsub('[0-9]*,',"",L4[[1]][2])))

    DP = sum(as.numeric(L1[[1]][4]),
             as.numeric(L2[[1]][4]),
             as.numeric(L3[[1]][4]),
             as.numeric(L4[[1]][4]))

    VAF=(AD/DP)*100

    format <- paste("GT=",GT,";DP=",DP,";RD=",RD,";AD=",AD,";AF=",VAF, sep = "")

    mutect_dg[i,6:8] = c(GT,".",format)
}

write.table(mutect_dg, paste(fname,'rf.vcf',sep = "_"), row.names = F, quote = F, col.names = F, sep = "\t")

