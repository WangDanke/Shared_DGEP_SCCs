####
#step1: read expression data
#######GSE87410，ESCC(6个批次的数据集)，HNSC(六个批次的数据集)，CESC(两个批次的数据集)，CSCC(三个批次的数据集)，LUSC(两个批次的数据集)
setwd("D:\\SCC_DGEP_shared\\SCC_ENA_seq\\LUSC\\GSE87410")

# 
ID_list <- pheno_GSE87410$ID_SRR
GSE87410 <- data.frame(matrix(nrow = 62754, ncol = 22))
colnames(GSE87410) = ID_list
# 
for (file_id in ID_list) {
  # 
  file_path <- paste("D:\\SCC_DGEP_shared\\SCC_ENA_seq\\LUSC\\GSE87410\\Counts\\", file_id,"\\",file_id,".txt", sep = "")
  # 
  data<- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  GSE87410[,file_id] = data[,7]
  rownames(GSE87410) = data$Geneid
  GSE87410$gene_length = data$Length
}

save(GSE87410,pheno_GSE87410,file = "GSE87410.RData")

##step2: anno
library(org.Hs.eg.db)
columns(org.Hs.eg.db)
library('clusterProfiler')
library(biomaRt)

# 
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# 
gene_id <- rownames(GSE87410)
# 
gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = gene_id,
                   mart = ensembl)

GSE87410$ensembl_gene_id = rownames(GSE87410)
GSE87410 = merge(GSE87410,gene_info,by = "ensembl_gene_id",all = FALSE)

save(GSE87410,pheno_GSE87410,file = "GSE87410_anno.RData")


######step3: TPM transform
#
setwd("D:\\鳞癌分析\\泛鳞癌seq数据_标准参考基因组\\TPM_SCCs")
load("D:/鳞癌分析/泛鳞癌seq数据_标准参考基因组/Counts_SCCs/GSE144293_anno.RData")
########
library(dplyr)
##去重
GSE144293 <- GSE144293 %>% distinct(external_gene_name, .keep_all = TRUE)

kb = GSE144293$gene_length/1000
counts <- GSE144293 %>%
  select(-ensembl_gene_id, -gene_length, -external_gene_name)
counts = as.matrix(counts)

rpk = counts/kb
tpm = t(t(rpk)/colSums(rpk)*1000000)

GSE144293_EXPR = as.data.frame(log2(tpm+1))
rownames(GSE144293_EXPR) = GSE144293$external_gene_name

save(GSE144293_EXPR,pheno_GSE144293,file = "GSE144293_TPM_log.RData")






























