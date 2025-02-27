#Drug sensitive analysis
install.packages("data.table")
if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("MultiAssayExperiment")
install.packages("RCurl")
install.packages("SummarizedExperiment")
install.packages('GenomicRanges')
install.packages("GenomeInfoDb")
library(GenomeInfoDb)
library(GenomicRanges)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(data.table)
library(oncoPredict)
library(dplyr)
library(ggsci)
library(ggpubr)
library(gtools)
library(oncoPredict)
library(reshape2)
library(tidyverse)


load("D:\\SCC_DGEP_shared\\SCCs_array\\Subtype_analysis\\Prognostic_model\\TCGA_Group_Risk.RData")
load("D:\\SCC_DGEP_shared\\SCC_TCGA\\TCGA_4_SCCs_8816x1169.RData")
pheno_SCC_TCGA = pheno_SCC_TCGA[!pheno_SCC_TCGA$project_id=="TCGA-LUSC",]
TCGA_SCC_EXPR = TCGA_SCC_EXPR[,rownames(pheno_SCC_TCGA)]


setwd("D:\\SCC_DGEP_shared\\SCCs_array\\14_Drug\\")
GDSC2_Expr <- readRDS("GDSC2_Expr (RMA Normalized and Log Transformed).rds")
GDSC2_Drug <- readRDS("GDSC2_Res.rds")
GDSC2_Drug = exp(GDSC2_Drug)
identical(rownames(GDSC2_Drug),colnames(GDSC2_Expr))


calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Drug,
              testExprData = as.matrix(TCGA_SCC_EXPR), 
              batchCorrect = 'eb',  
              powerTransformPhenotype = F,
              minNumSamples = 20,
              printOutput = T,
              removeLowVaryingGenes = 0.2,
              removeLowVaringGenesFrom = "homogenizeData"
)

Drug_prediction = read.csv("./calcPhenotype_Output/DrugPredictions.csv",header = T)
colnames(Drug_prediction)[1] = "case_submitter_id"

TCGA_data$case_submitter_id = rownames(TCGA_data)
phe = TCGA_data[,c("case_submitter_id","risk_group")]
phe = phe[!phe$risk_group=="Medium",]
phe$risk_group = factor(phe$risk_group,levels = c("High","Low"))


##Figure 7E
Sub_Drug_prediction = cbind(Drug_prediction[,1],
                            Drug_prediction[,c(
                              "Cisplatin_1005", "Afatinib_1032","Vinblastine_1004",
                              "Gemcitabine_1190","Irinotecan_1088","Vinorelbine_2048"
                                               )])
colnames(Sub_Drug_prediction)[1] = "case_submitter_id"
rownames(Sub_Drug_prediction) = Sub_Drug_prediction$case_submitter_id
Sub_Drug_prediction %>% 
  select(1:7) %>% 
  inner_join(phe, "case_submitter_id") %>% 
  pivot_longer(2:7, names_to = "drugs", values_to = "ic50") %>%
  mutate(drugs = factor(drugs, levels = c(
    "Cisplatin_1005", "Afatinib_1032", "Gemcitabine_1190",
    "Irinotecan_1088","Vinblastine_1004","Vinorelbine_2048"
  ))) %>%
  ggplot(., aes(risk_group, ic50)) +
  geom_boxplot(aes(fill = risk_group),
               notch = TRUE, notchwidth = 0.9) +
  scale_fill_jama() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none") +
  facet_wrap(vars(drugs), scales = "free_y", nrow = 1) +
  stat_compare_means()


####
###
gene_names <- c("COL1A1", "MMP1", "SERPINE1", "KRT6A", "IGF2BP3", "SPP1")
gene_expr = as.data.frame(t(TCGA_SCC_EXPR[gene_names,]))
gene_expr = gene_expr[rownames(Sub_Drug_prediction),]
gene_expr = gene_expr[rownames(TCGA_data),]
gene_expr = log2(gene_expr+1)
gene_expr$risk_score = as.numeric(TCGA_data$risk_score)

View(gene_expr)

##Correlation analysis between gene expression and drug sensitivity
# 
Sub_Drug_prediction = Sub_Drug_prediction[,-1]
Sub_Drug_prediction = Sub_Drug_prediction[rownames(gene_expr),]
# 
cor_results <- data.frame(Column1 = character(), Column2 = character(),
                          Correlation = numeric(), PValue = numeric(), stringsAsFactors = FALSE)

for (col1 in colnames(gene_expr)) {
  for (col2 in colnames(Sub_Drug_prediction)) {
    
    # 
    data1 <- gene_expr[[col1]]
    data2 <- Sub_Drug_prediction[[col2]]
    
    # 
    test_result <- cor.test(data1, data2, method="spearman")
    
    # 
    cor_results <- rbind(cor_results, data.frame(Column1 = col1, Column2 = col2,
                                                 Correlation = test_result$estimate,
                                                 PValue = test_result$p.value))
  }
}

# 
print(cor_results)

library(ggplot2)

cor_results$significance <- ifelse(cor_results$PValue < 0.001, "***",
                                   ifelse(cor_results$PValue < 0.01, "**",
                                          ifelse(cor_results$PValue < 0.05, "*", "")))


cor_results$Column2 = factor(cor_results$Column2,levels = c("Cisplatin_1005", "Afatinib_1032", "Gemcitabine_1190",
                                                            "Irinotecan_1088","Vinblastine_1004","Vinorelbine_2048"))
cor_results$Column1 = factor(cor_results$Column1,levels = c("COL1A1", "MMP1", "SERPINE1", "KRT6A", "IGF2BP3", "SPP1","risk_score"))
# Figure 7F
ggplot(cor_results, aes(Column2,Column1,  fill = Correlation)) +
  geom_tile(color = "black") +  
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text.y = element_text(angle = 0, hjust = 1)) +  
  labs(title = "Correlation Heatmap", x = "IC50", y = "Expression level and risk score")+
  geom_text(aes(label = significance), color = "black", size = 3)




