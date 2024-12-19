#Combined the datasets from the same type of SCC (Microarray)
##use the Combat function 
library(affy)
library(limma)
library(biomaRt)
library(sva)
library(stringi)
library(stringr)
library(GEOquery)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(WGCNA); library(lumi); library(ggplot2); library(nlme);library(affy)
library(ggplot2); library(Cairo); library(GEOquery);library(nlme);
library(biomaRt); library(sva);library(limma);library(illuminaio)
#####
#####
####CESC, three datasets
setwd("D:\\SCC_DGEP_shared\\SCCs_array\\After_preprocess\\CESC")

pheno_GSE166466$batch = "GSE166466"
pheno_GSE7803$batch = "GSE7803"
pheno_GSE9750$batch = "GSE9750"
sel_phe = c("Sample_title","ID","sample_type","batch","organ")
##select the genes that appears in all datasets
AA = intersect(rownames(GSE166466_anno),rownames(GSE7803_anno))
BB = intersect(rownames(GSE9750_anno),AA)

##combined three datasets
CESC_EXPR = cbind(GSE166466_anno[BB,],GSE7803_anno[BB,],GSE9750_anno[BB,])
pheno_CESC = rbind(pheno_GSE166466[,sel_phe],pheno_GSE7803[,sel_phe],pheno_GSE9750[,sel_phe])

##match
matchSN = match(colnames(CESC_EXPR), rownames(pheno_CESC))
CESC_EXPR = CESC_EXPR[,matchSN]
save(file = "CESC_EXPR_3_datasets.RData",CESC_EXPR,pheno_CESC)
##batch correct
##boxplot before batch correct

pheno_CESC$sample_type = factor(pheno_CESC$sample_type,levels = c("N","T"))
pheno_CESC$Group = as.numeric(ifelse(pheno_CESC$sample_type=="T","1","0"))

boxplot(CESC_EXPR, range = 0, col= as.numeric(pheno_CESC$sample_type), main ="Boxplot", ylab = "Intensity")
ddb.pca <- PCA(t(CESC_EXPR), graph = F )
fviz_pca_ind(ddb.pca,
             geom.ind = "point",     # show points only (but not "text")
             col.ind = pheno_CESC$batch, # color by groups
             addEllipses = TRUE,     # Concentration ellipses
             legend.title = "Groups",
             range=0
)
#correct
mod = model.matrix(~sample_type, data=pheno_CESC)
batch = factor(pheno_CESC$batch)
CESC_EXPR = ComBat(CESC_EXPR, batch=batch, mod=mod)
CESC_EXPR = as.data.frame(CESC_EXPR)
##PCA plot, boxplot and histogram after batch correct
boxplot(CESC_EXPR, range = 0, col= as.numeric(pheno_CESC$sample_type), main ="Boxplot", ylab = "Intensity")
ddb.pca <- PCA(t(CESC_EXPR), graph = F )
fviz_pca_ind(ddb.pca,
             geom.ind = "point",     # show points only (but not "text")
             col.ind = pheno_CESC$batch, # color by groups
             addEllipses = TRUE,     # Concentration ellipses
             legend.title = "Groups",
             range=0
)
# Histogram
i = 1; plot(density((CESC_EXPR[,i]), na.rm=T), col = i, main="Hist of Log2 Exp", xlab = "log2 exp")
for(i in 2:dim(CESC_EXPR)[2])
  lines(density((CESC_EXPR[,i]), na.rm=T), col =as.numeric(pheno_CESC$sample_type))
legend("topright", levels(pheno_CESC$sample_type), cex=0.7, text.col = 1:2)
save(CESC_EXPR,pheno_CESC,file = "CESC_3_combat.RData")

#####
#####
###CSCC,two datasets
setwd("D:\\SCC_DGEP_shared\\SCCs_array\\After_preprocess\\CSCC")

pheno_GSE108008$batch = "GSE108008"
pheno_GSE98780_com$batch = "GSE98780_com"

sel_phe = c("Sample_title","ID","sample_type", "organ","batch","Group")

pheno_GSE108008 = pheno_GSE108008[,sel_phe]
pheno_GSE98780_com = pheno_GSE98780_com[,sel_phe]

CC = intersect(rownames(GSE108008_anno),rownames(GSE98780_com_anno))

GSE108008_anno = GSE108008_anno[CC,]
GSE98780_com_anno = GSE98780_com_anno[CC,]

#combined
CSCC_EXPR = cbind(GSE108008_anno,GSE98780_com_anno)
pheno_CSCC = rbind(pheno_GSE108008,pheno_GSE98780_com)

#match
matchSN = match(colnames(CSCC_EXPR), rownames(pheno_CSCC))
CSCC_EXPR = CSCC_EXPR[,matchSN]
save(file = "CSCC_EXPR_2_datasets.RData",CSCC_EXPR,pheno_CSCC)
##before batch correct
boxplot(CSCC_EXPR, range = 0, col= as.numeric(pheno_CSCC$sample_type), main ="Boxplot", ylab = "Intensity")
ddb.pca <- PCA(t(CSCC_EXPR), graph = F )
fviz_pca_ind(ddb.pca,
             geom.ind = "point",     # show points only (but not "text")
             col.ind = pheno_CSCC$batch, # color by groups
             addEllipses = TRUE,     # Concentration ellipses
             legend.title = "Groups",
             range=0
)
#####batch correction
mod = model.matrix(~sample_type, data=pheno_CSCC)
batch = factor(pheno_CSCC$batch)
CSCC_EXPR = ComBat(CSCC_EXPR, batch=batch, mod=mod)
CSCC_EXPR = as.data.frame(CSCC_EXPR)
#after batch correct
boxplot(CSCC_EXPR, range = 0, col= as.numeric(pheno_CSCC$sample_type), main ="Boxplot", ylab = "Intensity")
ddb.pca <- PCA(t(CSCC_EXPR), graph = F )
fviz_pca_ind(ddb.pca,
             geom.ind = "point",     # show points only (but not "text")
             col.ind = pheno_CSCC$sample_type, # color by groups,batch or sample_type
             addEllipses = TRUE,     # Concentration ellipses
             legend.title = "Groups",
             range=0
)
# Histogram
i = 1; plot(density((CSCC_EXPR[,i]), na.rm=T), col = i, main="Hist of Log2 Exp", xlab = "log2 exp")
for(i in 2:dim(CSCC_EXPR)[2])
  lines(density((CSCC_EXPR[,i]), na.rm=T), col =as.numeric(pheno_CSCC$sample_type))
legend("topright", levels(pheno_CSCC$sample_type), cex=0.7, text.col = 1:2)
save(file = "CSCC_combat_2.RData",CSCC_EXPR,pheno_CSCC)

#####
#####
#####HNSC, five datasets
setwd("D:\\SCC_DGEP_shared\\SCCs_array\\After_preprocess\\HNSC")

pheno_GSE138206$batch = "GSE138206"
pheno_GSE172120$batch = "GSE172120"
pheno_GSE23558$batch = "GSE23558"
pheno_GSE31056$batch = "GSE31056"
pheno_GSE201777$batch = "GSE201777"

sel_phe = c("Sample_title","ID","sample_type","organ","batch")

pheno_GSE138206 = pheno_GSE138206[,sel_phe]
pheno_GSE172120 = pheno_GSE172120[,sel_phe]
pheno_GSE23558 = pheno_GSE23558[,sel_phe]

pheno_GSE201777 = pheno_GSE201777[,sel_phe]
pheno_GSE31056 = pheno_GSE31056[,sel_phe]

AA = intersect(rownames(GSE138206_anno),rownames(GSE172120_anno))
BB = intersect(rownames(GSE31056_anno),rownames(GSE23558_anno))
EE = intersect(AA,BB)
GG = intersect(EE,rownames(GSE201777_anno))

GSE138206_anno = GSE138206_anno[GG,]
GSE172120_anno = GSE172120_anno[GG,]
GSE23558_anno = GSE23558_anno[GG,]
GSE201777_anno = GSE201777_anno[GG,]
GSE31056_anno = GSE31056_anno[GG,]


##combined
HNSC_EXPR = cbind(GSE138206_anno,GSE172120_anno,GSE23558_anno,
                  GSE201777_anno,GSE31056_anno)
pheno_HNSC = rbind(pheno_GSE138206,pheno_GSE172120,pheno_GSE23558,
                   pheno_GSE201777,pheno_GSE31056)

##match
matchSN = match(colnames(HNSC_EXPR), rownames(pheno_HNSC))
HNSC_EXPR = HNSC_EXPR[,matchSN]
save(file = "HNSC_EXPR_5_datasets.RData",HNSC_EXPR,pheno_HNSC)


#before batch correct
###clean the metadata
pheno_HNSC$sample_type = factor(pheno_HNSC$sample_type,levels = c("N","T"))
pheno_HNSC$Group = as.numeric(ifelse(pheno_HNSC$sample_type=="T","1","0"))

boxplot(HNSC_EXPR, range = 0, col= as.numeric(pheno_HNSC$Group), main ="Boxplot", ylab = "Intensity")
library("FactoMineR")
library("factoextra")
ddb.pca <- PCA(t(HNSC_EXPR), graph = F )
fviz_pca_ind(ddb.pca,
             geom.ind = "point",     # show points only (but not "text")
             col.ind = pheno_HNSC$batch, # color by groups
             addEllipses = TRUE,     # Concentration ellipses
             legend.title = "Groups",
             range=0
)

#batch correction
mod = model.matrix(~sample_type, data=pheno_HNSC)
batch = factor(pheno_HNSC$batch)
HNSC_EXPR = ComBat(HNSC_EXPR, batch=batch, mod=mod)
HNSC_EXPR = as.data.frame(HNSC_EXPR)

##after batch correct
boxplot(HNSC_EXPR, range = 0, col= as.numeric(pheno_HNSC$sample_type), main ="Boxplot", ylab = "Intensity")
legend("topright", legend=c("N", "T"), col = 1:2, pch=19)
ddb.pca <- PCA(t(HNSC_EXPR), graph = F )
fviz_pca_ind(ddb.pca,
             geom.ind = "point",     # show points only (but not "text")
             col.ind = pheno_HNSC$batch, # color by groups
             addEllipses = TRUE,     # Concentration ellipses
             legend.title = "Groups",
             range=0
)
# Histogram
i = 1; plot(density((HNSC_EXPR[,i]), na.rm=T), col = i, main="Hist of Log2 Exp", xlab = "log2 exp")
for(i in 2:dim(HNSC_EXPR)[2])
  lines(density((HNSC_EXPR[,i]), na.rm=T), col =as.numeric(pheno_HNSC$sample_type))
legend("topright", levels(pheno_HNSC$sample_type), cex=0.7, text.col = 1:2)

save(HNSC_EXPR,pheno_HNSC,file = "HNSC_5_combat.RData")

######
######
###### ESCC, seven datasets
setwd("D:\\SCC_DGEP_shared\\SCCs_array\\After_preprocess\\ESCC")

pheno_GSE17351$batch = "GSE17351"
pheno_GSE20347$batch = "GSE20347"
pheno_GSE38129$batch = "GSE38129"
pheno_GSE23400_96$batch = "GSE23400"
pheno_GSE77861$batch = "GSE77861"
pheno_GSE161533$batch = "GSE161533"
pheno_GSE100942$batch = "GSE100942"

sel_phe = c("Sample_title","ID","sample_type","organ","batch")

pheno_GSE17351 = pheno_GSE17351[,sel_phe]
pheno_GSE20347 = pheno_GSE20347[,sel_phe]
pheno_GSE38129 = pheno_GSE38129[,sel_phe]
pheno_GSE23400 = pheno_GSE23400_96[,sel_phe]
pheno_GSE100942 = pheno_GSE100942[,sel_phe]
pheno_GSE161533 = pheno_GSE161533[,sel_phe]
pheno_GSE77861 = pheno_GSE77861[,sel_phe]

AA = intersect(rownames(GSE17351_anno),rownames(GSE20347_anno))
BB = intersect(rownames(GSE38129_anno),rownames(GSE23400_96_anno))
CC = intersect(rownames(GSE100942_anno),rownames(GSE161533_anno))
EE = intersect(AA,BB)
FF = intersect(CC,rownames(GSE77861_anno))
GG = intersect(EE,FF)

GSE17351_anno = GSE17351_anno[GG,]
GSE20347_anno = GSE20347_anno[GG,]
GSE38129_anno = GSE38129_anno[GG,]
GSE23400_anno = GSE23400_96_anno[GG,]
GSE161533_anno = GSE161533_anno[GG,]
GSE100942_anno = GSE100942_anno[GG,]
GSE77861_anno = GSE77861_anno[GG,]

##combined
ESCC_EXPR = cbind(GSE17351_anno,GSE20347_anno,GSE38129_anno,GSE23400_anno,
                  GSE161533_anno, GSE100942_anno,GSE77861_anno)
pheno_ESCC = rbind(pheno_GSE17351,pheno_GSE20347,pheno_GSE38129,pheno_GSE23400,
                   pheno_GSE100942,pheno_GSE161533,pheno_GSE77861)

##match
matchSN = match(colnames(ESCC_EXPR), rownames(pheno_ESCC))
ESCC_EXPR = ESCC_EXPR[,matchSN]
save(file = "ESCC_EXPR_7_datasets.RData",ESCC_EXPR,pheno_ESCC)


#before batch correct
###clean the metadata
pheno_ESCC$sample_type = factor(pheno_ESCC$sample_type,levels = c("N","T"))
pheno_ESCC$Group = as.numeric(ifelse(pheno_ESCC$sample_type=="T","1","0"))

boxplot(ESCC_EXPR, range = 0, col= as.numeric(pheno_ESCC$Group), main ="Boxplot", ylab = "Intensity")
ddb.pca <- PCA(t(ESCC_EXPR), graph = F )
fviz_pca_ind(ddb.pca,
             geom.ind = "point",     # show points only (but not "text")
             col.ind = pheno_ESCC$batch, # color by groups
             addEllipses = TRUE,     # Concentration ellipses
             legend.title = "Groups",
             range=0
)


#batch correction
mod = model.matrix(~sample_type, data=pheno_ESCC)
batch = factor(pheno_ESCC$batch)
ESCC_EXPR = ComBat(ESCC_EXPR, batch=batch, mod=mod)
ESCC_EXPR = as.data.frame(ESCC_EXPR)
##
boxplot(ESCC_EXPR, range = 0, col= as.numeric(pheno_ESCC$sample_type), main ="Boxplot", ylab = "Intensity")
ddb.pca <- PCA(t(ESCC_EXPR), graph = F )
fviz_pca_ind(ddb.pca,
             geom.ind = "point",     # show points only (but not "text")
             col.ind = pheno_ESCC$batch, # color by groups
             addEllipses = TRUE,     # Concentration ellipses
             legend.title = "Groups",
             range=0
)

#####
#####
#####LUSC four datasets
setwd("D:\\SCC_DGEP_shared\\SCCs_array\\After_preprocess\\LUSC")

pheno_GSE18842$batch = "GSE18842"
pheno_GSE126533$batch = "GSE126533"
pheno_GSE74706$batch = "GSE74706"
pheno_GSE67061_com$batch = "GSE67061_com"

sel_phe = c("Sample_title","ID","sample_type","batch")

pheno_GSE18842 = pheno_GSE18842[,sel_phe]
pheno_GSE126533 = pheno_GSE126533[,sel_phe]
pheno_GSE74706 = pheno_GSE74706[,sel_phe]
pheno_GSE67061_com = pheno_GSE67061_com[,sel_phe]

AA = intersect(rownames(GSE18842_anno),rownames(GSE126533_anno))
BB = intersect(rownames(GSE74706_anno),rownames(GSE67061_com_anno))
CC = intersect(AA,BB)

GSE18842_anno = GSE18842_anno[CC,]
GSE126533_anno = GSE126533_anno[CC,]
GSE74706_anno = GSE74706_anno[CC,]
GSE67061_com_anno = GSE67061_com_anno[CC,]

#combined
LUSC_EXPR = cbind(GSE18842_anno,GSE126533_anno,GSE74706_anno,GSE67061_com_anno)
pheno_LUSC = rbind(pheno_GSE18842,pheno_GSE126533,pheno_GSE74706,pheno_GSE67061_com)

#match
matchSN = match(colnames(LUSC_EXPR), rownames(pheno_LUSC))
LUSC_EXPR = LUSC_EXPR[,matchSN]
save(file = "LUSC_EXPR_4_datasets.RData",LUSC_EXPR,pheno_LUSC)

#before batch correct
boxplot(LUSC_EXPR, range = 0, col= as.numeric(pheno_LUSC$sample_type), main ="Boxplot", ylab = "Intensity")
ddb.pca <- PCA(t(LUSC_EXPR), graph = F )
fviz_pca_ind(ddb.pca,
             geom.ind = "point",     # show points only (but not "text")
             col.ind = pheno_LUSC$batch, # color by groups
             addEllipses = TRUE,     # Concentration ellipses
             legend.title = "Groups",
             range=0
)
#####batch correction
mod = model.matrix(~sample_type, data=pheno_LUSC)
batch = factor(pheno_LUSC$batch)
LUSC_EXPR = ComBat(LUSC_EXPR, batch=batch, mod=mod)
LUSC_EXPR = as.data.frame(LUSC_EXPR)
#after batch correct
boxplot(LUSC_EXPR, range = 0, col= as.numeric(pheno_LUSC$sample_type), main ="Boxplot", ylab = "Intensity")
ddb.pca <- PCA(t(LUSC_EXPR), graph = F )
fviz_pca_ind(ddb.pca,
             geom.ind = "point",     # show points only (but not "text")
             col.ind = pheno_LUSC$batch, # color by groups,batch or sample_type
             addEllipses = TRUE,     # Concentration ellipses
             legend.title = "Groups",
             range=0
)
# Histogram
i = 1; plot(density((LUSC_EXPR[,i]), na.rm=T), col = i, main="Hist of Log2 Exp", xlab = "log2 exp")
for(i in 2:dim(LUSC_EXPR)[2])
  lines(density((LUSC_EXPR[,i]), na.rm=T), col =as.numeric(pheno_LUSC$sample_type))
legend("topright", levels(pheno_LUSC$sample_type), cex=0.7, text.col = 1:2)
save(file = "LUSC_combat_4.RData",LUSC_EXPR,pheno_LUSC)

#####
#####
#####EAC, two datasets
setwd("D:\\SCC_DGEP_shared\\SCCs_array\\After_preprocess\\EAC")

pheno_GSE74553$batch = "GSE74553"
pheno_GSE26886$batch = "GSE26886"

sel_phe = c("Sample_title","ID","sample_type","batch","Group")

pheno_GSE74553 = pheno_GSE74553[,sel_phe]
pheno_GSE26886 = pheno_GSE26886[,sel_phe]

CC = intersect(rownames(GSE74553_anno),rownames(GSE26886_anno))

GSE74553_anno = GSE74553_anno[CC,]
GSE26886_anno = GSE26886_anno[CC,]

#combined
EAC_EXPR = cbind(GSE74553_anno,GSE26886_anno)
pheno_EAC = rbind(pheno_GSE74553,pheno_GSE26886)

#match
matchSN = match(colnames(EAC_EXPR), rownames(pheno_EAC))
EAC_EXPR = EAC_EXPR[,matchSN]
save(file = "EAC_EXPR_2_datasets.RData",EAC_EXPR,pheno_EAC)
##before batch correct
boxplot(EAC_EXPR, range = 0, col= as.numeric(pheno_EAC$sample_type), main ="Boxplot", ylab = "Intensity")
ddb.pca <- PCA(t(EAC_EXPR), graph = F )
fviz_pca_ind(ddb.pca,
             geom.ind = "point",     # show points only (but not "text")
             col.ind = pheno_EAC$batch, # color by groups
             addEllipses = TRUE,     # Concentration ellipses
             legend.title = "Groups",
             range=0
)
#####batch correction
mod = model.matrix(~sample_type, data=pheno_EAC)
batch = factor(pheno_EAC$batch)
EAC_EXPR = ComBat(EAC_EXPR, batch=batch, mod=mod)
EAC_EXPR = as.data.frame(EAC_EXPR)
#after batch correct
boxplot(EAC_EXPR, range = 0, col= as.numeric(pheno_EAC$sample_type), main ="Boxplot", ylab = "Intensity")
ddb.pca <- PCA(t(EAC_EXPR), graph = F )
fviz_pca_ind(ddb.pca,
             geom.ind = "point",     # show points only (but not "text")
             col.ind = pheno_EAC$sample_type, # color by groups,batch or sample_type
             addEllipses = TRUE,     # Concentration ellipses
             legend.title = "Groups",
             range=0
)
# Histogram
i = 1; plot(density((EAC_EXPR[,i]), na.rm=T), col = i, main="Hist of Log2 Exp", xlab = "log2 exp")
for(i in 2:dim(EAC_EXPR)[2])
  lines(density((EAC_EXPR[,i]), na.rm=T), col =as.numeric(pheno_EAC$sample_type))
legend("topright", levels(pheno_EAC$sample_type), cex=0.7, text.col = 1:2)
save(file = "EAC_combat_2.RData",EAC_EXPR,pheno_EAC)

#####
####
####LUAD, four datasets
setwd("D:\\SCC_DGEP_shared\\SCCs_array\\After_preprocess\\LUAD")

pheno_GSE118370$batch = "GSE118370"
pheno_GSE130779$batch = "GSE130779"
pheno_GSE136043$batch = "GSE136043"
pheno_GSE140797$batch = "GSE140797"

sel_phe = c("Sample_title","ID","sample_type","batch","Group")

pheno_GSE118370 = pheno_GSE118370[,sel_phe]
pheno_GSE130779 = pheno_GSE130779[,sel_phe]
pheno_GSE136043 = pheno_GSE136043[,sel_phe]
pheno_GSE140797 = pheno_GSE140797[,sel_phe]

AA = intersect(rownames(GSE118370_anno),rownames(GSE130779_anno))
BB = intersect(rownames(GSE136043_anno),rownames(GSE140797_anno))
CC = intersect(AA,BB)

GSE118370_anno = GSE118370_anno[CC,]
GSE130779_anno = GSE130779_anno[CC,]
GSE136043_anno = GSE136043_anno[CC,]
GSE140797_anno = GSE140797_anno[CC,]

#combined
LUAD_EXPR = cbind(GSE118370_anno,GSE130779_anno,GSE136043_anno,GSE140797_anno)
pheno_LUAD = rbind(pheno_GSE118370,pheno_GSE130779,pheno_GSE136043,pheno_GSE140797)
save(file = "LUAD_EXPR_4_datasets.RData",LUAD_EXPR,pheno_LUAD)
#match
matchSN = match(colnames(LUAD_EXPR), rownames(pheno_LUAD))
LUAD_EXPR = LUAD_EXPR[,matchSN]
#before batch correct
boxplot(LUAD_EXPR, range = 0, col= as.numeric(pheno_LUAD$sample_type), main ="Boxplot", ylab = "Intensity")
ddb.pca <- PCA(t(LUAD_EXPR), graph = F )
fviz_pca_ind(ddb.pca,
             geom.ind = "point",     # show points only (but not "text")
             col.ind = pheno_LUAD$batch, # color by groups
             addEllipses = TRUE,     # Concentration ellipses
             legend.title = "Groups",
             range=0
)
#####batch correct 
mod = model.matrix(~sample_type, data=pheno_LUAD)
batch = factor(pheno_LUAD$batch)
LUAD_EXPR = ComBat(LUAD_EXPR, batch=batch, mod=mod)
LUAD_EXPR = as.data.frame(LUAD_EXPR)
##after batch correct
boxplot(LUAD_EXPR, range = 0, col= as.numeric(pheno_LUAD$sample_type), main ="Boxplot", ylab = "Intensity")
ddb.pca <- PCA(t(LUAD_EXPR), graph = F )
fviz_pca_ind(ddb.pca,
             geom.ind = "point",     # show points only (but not "text")
             col.ind = pheno_LUAD$batch, # color by groups,batch or sample_type
             addEllipses = TRUE,     # Concentration ellipses
             legend.title = "Groups",
             range=0
)

# Histogram
i = 1; plot(density((LUAD_EXPR[,i]), na.rm=T), col = i, main="Hist of Log2 Exp", xlab = "log2 exp")
for(i in 2:dim(LUAD_EXPR)[2])
  lines(density((LUAD_EXPR[,i]), na.rm=T), col =as.numeric(pheno_LUAD$sample_type))
legend("topright", levels(pheno_LUAD$sample_type), cex=0.7, text.col = 1:2)
save(file = "LUAD_combat_4.RData",LUAD_EXPR,pheno_LUAD)


















