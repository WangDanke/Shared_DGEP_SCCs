# the SCC gene expression data were downloaded from TCGA data base (including TCGA-ESCA, TCGA-CESC, TCGA-HNSC and TCGA-LUSC projects)
main_folder <- "D:\\SCC_DGEP_shared\\SCC_TCGA\\SCCs_gene_expression_primary_tumor"
sub_folders <- list.dirs(main_folder, recursive = FALSE)
tsv_files <- list.files(path = sub_folders, pattern = "\\.tsv$", full.names = TRUE)

basename(sub_folders)
basename(tsv_files)

merge_data <- read.table("D:\\SCC_DGEP_shared\\SCC_TCGA\\SCCs_gene_expression_primary_tumor\\000b6b94-572d-4d06-a8f4-2e43829f83d4\\e74f321c-217f-4bdc-ad17-f132501b5157.rna_seq.augmented_star_gene_counts.tsv", header=T, sep="\t",fill = T)
merge_data <- merge_data[-c(1:4),]
#choose the TPM data
merge_data <- merge_data[,c(1:3,7)]
#
colnames(merge_data)[4] <- basename(tsv_files[1])
n <- length(tsv_files)
for (i in 2:n){
  new_data = read.table(file = tsv_files[i], header=T, sep="\t",fill = T)
  new_data <- new_data[-c(1:4),]
  colnames(new_data)[7] <- basename(tsv_files[i])
  merge_data=cbind(merge_data,new_data[,7])
}

##
colnames(merge_data)[4:1172] = basename(sub_folders[1:1169])
save(merge_data,file = "TCGA_SCC.RData")

###prepare the metadata
setwd("D:\\SCC_DGEP_shared\\SCC_TCGA\\")
pheno_SCC_TCGA = read.csv("clinical.tsv", header=T, sep="\t",fill = T)
pheno_SCC_TCGA <- pheno_SCC_TCGA[!duplicated(pheno_SCC_TCGA$case_submitter_id, fromLast=T), ] 
pheno_SCC_TCGA <- pheno_SCC_TCGA[!duplicated(pheno_SCC_TCGA$case_id, fromLast=T), ]
pheno_SCC_TCGA <- pheno_SCC_TCGA[,c('case_submitter_id',
                                           'case_id',
                                           'project_id',
                                           'age_at_index',
                                           'days_to_birth',
                                           'days_to_death',
                                           'gender',
                                           'vital_status',
                                           'age_at_diagnosis',
                                           'ajcc_clinical_t',
                                           'days_to_last_follow_up'
)]

save(pheno_SCC_TCGA,file = "pheno_SCC_TCGA.RData")

###match metadata and expression data
sample_sheet= read.table("sample_sheet.tsv", header=T, sep="\t",fill = T)
matchSN = match(colnames(merge_data)[4:1172], sample_sheet$File.ID)
rownames(sample_sheet) = sample_sheet$File_ID
sample_sheet = sample_sheet[matchSN,]
colnames(merge_data)[4:1172] = sample_sheet$Case.ID

A = intersect(colnames(merge_data)[4:1172], pheno_SCC_TCGA$case_submitter_id)
rownames(pheno_SCC_TCGA)=pheno_SCC_TCGA$case_submitter_id
pheno_SCC_TCGA = pheno_SCC_TCGA[A,]
TCGA_SCC_EXPR = merge_data[,A]
TCGA_SCC_EXPR = cbind(merge_data[,c(1,2,3)],TCGA_SCC_EXPR)

##
save(pheno_SCC_TCGA,TCGA_SCC_EXPR,file = 'TCGA_4_SCCs.RData')













