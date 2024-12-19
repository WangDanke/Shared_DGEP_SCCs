setwd("D:\\SCC_DGEP_shared\\SCCs_array\\Subtype_analysis")
load("D:\\SCC_DGEP_shared\\SCCs_array\\WGCNA\\SCC_combined_combat_for_WGCNA.RData")
load("D:D:\\SCC_DGEP_shared\\SCC_TCGA\\TCGA_4_SCCs.RData")
###Use the TCGA datasets to identify Subtypes, and GEO microarray dataset to validate the Subtypes
TCGA_SCC_EXPR <- TCGA_SCC_EXPR[!duplicated(TCGA_SCC_EXPR$gene_name, fromLast=T), ]
rownames(TCGA_SCC_EXPR) = TCGA_SCC_EXPR$gene_name
TCGA_SCC_EXPR = TCGA_SCC_EXPR[rownames(SCC_EXPR),]
save(TCGA_SCC_EXPR, pheno_SCC_TCGA,file = "TCGA_4_SCCs_8816x1169.RData")

#####use the hub genes found in WGCNA,
hub_gene = read.csv("9_module_hub_gene_10%.csv",header = T)
##hub_gene = read.csv("9_module_hub_gene_5%.csv",header = T)使用top5%的样本
library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)
library(dplyr)
library(Hmisc)

####identified subtypes of cross-organ (ESCC, HNSC, CESC) SCCs.

library("ConsensusClusterPlus")
######
##聚类鉴别分型_TCGA
TCGA_SCC_EXPR = TCGA_SCC_EXPR[,-c(1,2,3)]
SCC_TCGA_hub = TCGA_SCC_EXPR[hub_gene$gene,]
SCC_TCGA_hub = log2(SCC_TCGA_hub+1)
pheno_SCC_TCGA = subset(pheno_SCC_TCGA,!pheno_SCC_TCGA$project_id=="TCGA-LUSC")##Removal of LUSC samples
SCC_TCGA_hub = SCC_TCGA_hub[,rownames(pheno_SCC_TCGA)]

##Normalization
SCC_TCGA_hub = sweep(SCC_TCGA_hub,1, apply(SCC_TCGA_hub,1,median,na.rm=T))
SCC_TCGA_hub = as.matrix(SCC_TCGA_hub)
SCC_TCGA_hub = SCC_TCGA_hub[hub_gene$gene,]
##identified the Subtypes
title="maxK_10_TCGA_3SSC_hub_10"
Cluster_TCGA = ConsensusClusterPlus(SCC_TCGA_hub,maxK=10,reps=50,pItem=0.8,pFeature=1,
                                    title = title,clusterAlg="hc",distance="pearson",seed=1262118,plot="png")


##Prognosis of the Subtypes. Fig. 5A
library(survival)
library(survminer)
library(preprocessCore)

BB = as.data.frame(Cluster_TCGA[[4]][["consensusClass"]])
pheno_SCC_TCGA$Subtype_3SCC_10= BB$`Cluster_TCGA[[4]][["consensusClass"]]`

pheno_SCC_TCGA$time=as.numeric(pheno_SCC_TCGA$days_to_death)+as.numeric(pheno_SCC_TCGA$days_to_last_follow_up)
pheno_SCC_TCGA$time = pheno_SCC_TCGA$time/365
pheno_SCC_TCGA$event=ifelse(pheno_SCC_TCGA$vital_status=='Alive',0,1)

sfit <- survfit(Surv(time, event)~Subtype_3SCC_10, data=pheno_SCC_TCGA)
print(sfit)
colors <- c("#C6B3D3","#ED9F9B","#80BA8A","#9CD1C8")  # 对应的颜色
labels <- c("1","2","3","4")  
ggsurvplot(sfit, 
           conf.int = FALSE, 
           pval = TRUE, 
           palette = colors,    
           legend.labs = labels 
)

p


###########
###########Table S4.
########description o f the Subtypes
annotation_col = pheno_SCC_TCGA[,c("gender",'age','HPV_status','project_id',"ajcc_clinical_t","Subtype_3SCC_10")]
library(compareGroups)
annotation_col$gender = as.factor(annotation_col$gender)
annotation_col$project_id = as.factor(annotation_col$project_id)
annotation_col$HPV_status = as.factor(annotation_col$HPV_status)
annotation_col$ajcc_clinical_t = as.factor(annotation_col$ajcc_clinical_t)
annotation_col$Subtype_3SCC_10 = as.factor(annotation_col$Subtype_3SCC_10)

pheno_SCC_TCGA = pheno_SCC_TCGA[rownames(annotation_col),]
annotation_col$age = pheno_SCC_TCGA$age_at_index
annotation_col$age = as.numeric(annotation_col$age)
save(annotation_col,file = "annotation_col.RData")

annotation_col = annotation_col[,c(5,1,2,3,4)]
res1 <- compareGroups(Subtype_3SCC_10 ~ ., data = annotation_col, ref = 1)
restab = createTable(res1, show.ratio = TRUE)
restab
export2word(restab, file='table1.docx')



####
##Validation of the Subtypes in GEO microarray dataset
pheno_SCC_T = subset(pheno_SCC,pheno_SCC$sample_type=="T")
pheno_SCC_4 = subset(pheno_SCC_T,!pheno_SCC_T$organ=="lung")
pheno_SCC_3 = subset(pheno_SCC_4,!pheno_SCC_4$organ=="skin")
SCC_EXPR_T_3 = SCC_EXPR[,rownames(pheno_SCC_3)]

SCC_GEO_hub = SCC_EXPR_T_3[hub_gene$gene,]

##normailzation
SCC_GEO_hub = sweep(SCC_GEO_hub,1, apply(SCC_GEO_hub,1,median,na.rm=T))
SCC_GEO_hub = as.matrix(SCC_GEO_hub)
title="maxK_10_GEO_3SCC_hub_10"
Cluster_GEO = ConsensusClusterPlus(SCC_GEO_hub,maxK=10,reps=50,pItem=0.8,pFeature=1,
                               title = title,clusterAlg="hc",distance="pearson",seed=1262118,plot="png")
###match metadata and Subtypes
AA = as.data.frame(Cluster_GEO[[4]][["consensusClass"]])
pheno_SCC_3$Subtype = AA$`Cluster_GEO[[4]][["consensusClass"]]`

pheno_SCC_3<-pheno_SCC_3 %>%
  mutate(Subtype = case_when(
    Subtype == 1 ~ "subtype1",
    Subtype == 2 ~ "subtype2",
    Subtype == 3 ~ "subtype3",
    Subtype == 4 ~ "subtype4"
  ))
pheno_SCC_3 <- pheno_SCC_3[order(pheno_SCC_3$Subtype), ]
annotation_col = pheno_SCC_3[,c(4,5,7)]

annotation_col <- annotation_col[order(annotation_col$Subtype), ]
annotation_row =  read.csv("9_module_hub_gene_10%.csv",header = T,row.names = 1)
SCC_GEO_hub = SCC_GEO_hub[,rownames(annotation_col)]

annotation_colors = list(
  module = c("CD1"= "turquoise","CD2"="blue","CD3"="brown","CD4"="yellow","CD5"="green","CD6"="red",
             "CD7" = "black","CD8"="pink","CD9"="magenta"),
  gender = c("female"="#FFC1C1","male"="#87CEEB"),
  project_id = c("TCGA-CESC"="#778899","TCGA-ESCA" = "#CD5555","TCGA-HNSC"="#CDAA7D"),
  ajcc_clinical_t = c("NA" = "#E8E8E8","T1"="#CFCFCF","T2" = "#B5B5B5","T2" = "#9C9C9C","T3" = "#828282",
                      "T4" = "#696969","T4a" = "#4F4F4F","T4b"="#363636"),
  Subtype_3SCC_10 = c("subtype1"="#C6B3D3","subtype2"="#ED9F9B","subtype3"="#80BA8A","subtype4"="#9CD1C8")
  
)
#Fig. S9D
library(pheatmap)
pheatmap(SCC_GEO_hub ,cluster_rows = F,cluster_cols = F,
         color=colorRampPalette(c("navy","white","firebrick3"))(100),
         show_colnames = F,border_color = NA,show_rownames =T,
         annotation_col = annotation_col,annotation_row = annotation_row)


save(pheno_SCC_3,SCC_EXPR_T_3,file = "GEO_4subtypes_without_CSCC_LUSC.RData")
save(pheno_SCC_3,SCC_GEO_hub,file = "GEO_4subtypes_3_SCC.RData")










                     