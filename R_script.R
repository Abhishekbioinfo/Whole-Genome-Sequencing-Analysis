library(tidyverse)
args = commandArgs(trailingOnly=TRUE)
#setwd("/Users/debdutta.chatterjee_/Desktop/CNV_pytor_SS_LUNG_ELP/")
cnv_raw_data <- read_delim(args[1], col_names = F)

#print(cnv_raw_data)
#/mnt/data/CNV_pytor/KH-22RS70-FT_S4
#cnv_raw_data <- read_delim("/mnt/data/CNV_pytor/KH-22RS70-FT_S4/KH-22RS70-FT_S4_hard_filtered_calls_10k.tsv", col_names = F)

colnames(cnv_raw_data) = c("Sample","Method","CNV_type", "CHR","Start","End", "CNV_size", "CNV_level", 
                                "e-val1", "e-val2", "e-val3", "e-val4", "q0", "pN",
                                "dG", "Gene"
)
# cnv_raw_data = cnv_raw_data %>%
#   mutate(CN = log(10 ^ CNV_level))

cnv_raw_data = cnv_raw_data %>%
  mutate(Copy_Number = (CNV_level*2))
# for(i in 1:nrow(cnv_raw_data)) {
#   if ((cnv_raw_data$CN[i] - floor(cnv_raw_data$CN[i]))*10 >= 7.5) {
#     cnv_raw_data$CN[i] = ceiling(cnv_raw_data$CN[i])
#   } else{
#     cnv_raw_data$CN[i] = floor(cnv_raw_data$CN[i])
#   }
# }

map = read_delim(args[2], col_names = F)
#map = read_delim("/mnt/data/CNV_pytor/KH-22RS70-FT_S4/cyto.txt", col_names = F)
map = map[,1:4]
colnames(map) = c("CHROM", "START", "END", "ARM")
for(i in 1:nrow(cnv_raw_data)){
  for(j in 1:nrow(map)){
    if(cnv_raw_data[i,4] == map[j,1]){
      if(cnv_raw_data[i,5] >= map[j,2] & cnv_raw_data[i,6] <= map[j,3]){
        cnv_raw_data[i,18] = map[j,4]
      }else if(cnv_raw_data[i,5] > map[j,2] & map[j,3] - cnv_raw_data[i,5] > cnv_raw_data[i,6] - map[j+1,2]){
        cnv_raw_data[i,18] = map[j,4]
      }
      else if(cnv_raw_data[i,5] > map[j,2] & cnv_raw_data[i,6] - map[j+1,2] > map[j,3] - cnv_raw_data[i,5]){
        cnv_raw_data[i,18] = map[j+1,4]
      }
    } 
  } 
}

#cnv_raw_data_2 = cnv_raw_data %>%
#  mutate(CHR = str_c(CHR, ARM, sep = "_")) %>%
#  select(CNV_TYPE = CNV_type, CHR, START = Start, END = End, GENE = Gene, CN)


cnv_raw_data_2 = cnv_raw_data %>%
  mutate(ARM_pq = str_extract(ARM, "[a-z]*")) %>%
  mutate(CHR_pq = str_c(CHR,ARM_pq, sep = "")) %>%
  select(-ARM_pq) %>%
  mutate(CHR = str_c(CHR, ARM, sep = "")) %>%
  select(1,2,3,4,5,6,CHR_pq, Copy_Number,Gene, everything()) %>%
  mutate(Gene = str_replace_all(Gene, "\\((.*?)\\)", ""))






fname <- gsub('.tsv', "",args[1])
write.table(cnv_raw_data_2, paste(fname,'CNV_details.csv', sep = "_") , row.names = F, quote = F, sep = "\t")
#write.table(cnv_raw_data_2, paste('CNV_details.csv') , row.names = F, quote = F, sep = "\t")

