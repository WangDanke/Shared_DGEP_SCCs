#ESCC，GSE23400_96，clinical message and raw data process.
###loading packages
BiocManager::install("biomaRt")
BiocManager::install("FactoMineR")
library(WGCNA)
library(affy)
library(limma)
library(biomaRt)
library(sva)
library(PCA)
library(ggplot2)
library(FactoMineR)

setwd("D:\\SCC_DGEP_shared\\SCC_array\\ESCC\\GSE23400_RAW")
getwd()
# list the ".CEL" files
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)
basename(cel.files)

library(affy)
GSE23400_96_raw<- ReadAffy(filenames = cel.files)#
sampleNames(GSE23400_96_raw)

##rename samples
library(stringi)
sampleNames(GSE23400_96_raw)<-stri_sub(sampleNames(GSE23400_96_raw),1,9)
sampleNames(GSE23400_96_raw)

####function "RMA" to preprocess
GSE23400_96_raw_rma <- rma(GSE23400_96_raw)
GSE23400_96 <- exprs(GSE23400_96_raw_rma)

###anno the gene symbol
library(GEOquery)
GSE23400_96_anno <- as.data.frame(GSE23400_96)
GSE23400_96_anno$ID<-rownames(GSE23400_96_anno)
gpl571<-getGEO("GPL96",destdir=".")
iddata<-as.data.frame(Table(gpl571)[c('ID','Gene Symbol')])
GSE23400_96_anno<-merge(x=GSE23400_96_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE23400_96_anno <- GSE23400_96_anno[,-1]
###The expression levels of the same gene names were averaged
GSE23400_96_anno<- aggregate(GSE23400_96_anno,by = list(GSE23400_96_anno$`Gene Symbol`),FUN = mean)
head(GSE23400_96_anno)

GSE23400_96_anno <- GSE23400_96_anno[-1,]##blank gene name was dropped
rownames(GSE23400_96_anno) <- GSE23400_96_anno$Group.1
GSE23400_96_anno <- GSE23400_96_anno[,-c(1,108)]

save(GSE23400_96_anno,file = "GSE23400_96_NEW.RData")

###整理临床信息
library(tidyr)
library(dplyr)
library(stringr)
###
###GSE23400_96
phe_GSE23400_96 = read.csv("GSE23400-GPL96_series_matrix.csv", header = F)
phe_GSE23400_96 = phe_GSE23400_96[c(35:42),]
phe_GSE23400_96 = phe_GSE23400_96 %>% t()
colnames(phe_GSE23400_96) = phe_GSE23400_96[1,]
phe_GSE23400_96 = phe_GSE23400_96[-1,]
rownames(phe_GSE23400_96) = phe_GSE23400_96[,2]
phe_GSE23400_96 = as.data.frame(phe_GSE23400_96)

colnames(phe_GSE23400_96)
pheno_GSE23400_96 = phe_GSE23400_96[,c(1,2,8)]
View(pheno_GSE23400_96)
pheno_GSE23400_96$sample_type = ifelse(pheno_GSE23400_96$`!Sample_source_name_ch1` == "adjacent normal tissue","N","T")
pheno_GSE23400_96 = pheno_GSE23400_96[,c(1,2,4)]
colnames(pheno_GSE23400_96) = c("Sample_title","ID","sample_type")
pheno_GSE23400_96$organ = "esophagus"

####clean the datMet
pheno_GSE23400_96$sample_type = factor(pheno_GSE23400_96$sample_type,levels = c("N","T"))
pheno_GSE23400_96$Group = as.numeric(ifelse(pheno_GSE23400_96$sample_type=="T","1","0"))

save(file = "GSE23400_96_new.RData",GSE23400_96_anno,pheno_GSE23400_96)

##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE23400_96_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE23400_96$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomright", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE23400_96_anno)[outliers]); print(table(outliers))
GSE23400_96_anno = GSE23400_96_anno[,!outliers]
pheno_GSE23400_96 = pheno_GSE23400_96[!outliers,]

#clean the datMet
pheno_GSE23400_96$sample_type = factor(pheno_GSE23400_96$sample_type,levels = c("N","T"))
pheno_GSE23400_96$Group = as.numeric(ifelse(pheno_GSE23400_96$sample_type=="T","1","0"))

####Adjusting for unknown covariates
mod = model.matrix(~as.factor(Group),data = pheno_GSE23400_96)
mod0 = model.matrix(~1,data = pheno_GSE23400_96)
n.sv = num.sv(GSE23400_96_anno, mod, method="be")
GSE23400_96_anno = as.matrix(GSE23400_96_anno)

svobj = sva(GSE23400_96_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE23400_96$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 102)

X = svobj$sv
Y = GSE23400_96_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE23400_96_anno = GSE23400_96_anno-t(to_regress)
boxplot(GSE23400_96_anno)

save(file = "GSE23400_96.RData",GSE23400_96_anno,pheno_GSE23400_96)


