#ESCC,GSE161533,clinical message and raw data process.
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

setwd("D:\\SCC_DGEP_shared\\SCC_array\\ESCC\\GSE161533_RAW")
getwd()
# list the ".CEL" files
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)
basename(cel.files)

library(affy)
GSE161533_raw<- ReadAffy(filenames = cel.files)
sampleNames(GSE161533_raw)

##rename samples
library(stringi)
sampleNames(GSE161533_raw)<-stri_sub(sampleNames(GSE161533_raw),1,10)
sampleNames(GSE161533_raw)

####function "RMA" to preprocess
GSE161533_raw_rma <- rma(GSE161533_raw)
GSE161533 <- exprs(GSE161533_raw_rma)



###anno the gene symbol
library(GEOquery)
GSE161533_anno <- as.data.frame(GSE161533)
GSE161533_anno$ID<-rownames(GSE161533_anno)
gpl571=getGEO (filename = 'GPL570.soft')
iddata<-as.data.frame(Table(gpl571)[c('ID','Gene Symbol')])
GSE161533_anno<-merge(x=GSE161533_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE161533_anno <- GSE161533_anno[,-1]
###The expression levels of the same gene names were averaged
GSE161533_anno<- aggregate(GSE161533_anno,by = list(GSE161533_anno$`Gene Symbol`),FUN = mean)
head(GSE161533_anno)

GSE161533_anno <- GSE161533_anno[-1,]##blank gene name was dropped
rownames(GSE161533_anno) <- GSE161533_anno$Group.1
GSE161533_anno <- GSE161533_anno[,-c(1,58)]


save(GSE161533_anno,file = "GSE161533_NEW.RData")

#prepare the metadata
library(tidyr)
library(dplyr)
library(stringr)
###
###GSE161533
phe_GSE161533 = read.csv("GSE161533_series_matrix.csv", header = F)
phe_GSE161533 = phe_GSE161533[c(28:42),]
phe_GSE161533 = phe_GSE161533 %>% t()
colnames(phe_GSE161533) = phe_GSE161533[1,]
phe_GSE161533 = phe_GSE161533[-1,]
rownames(phe_GSE161533) = phe_GSE161533[,2]
phe_GSE161533 = as.data.frame(phe_GSE161533)

colnames(phe_GSE161533)
pheno_GSE161533 = phe_GSE161533[,c(1,2,10,11,12,13,14,15)]
View(pheno_GSE161533)
pheno_GSE161533$sample_type = ifelse(pheno_GSE161533$`!Sample_characteristics_ch1`=="tissue: normal tissue","N","T")
pheno_GSE161533$gender = ifelse(pheno_GSE161533$`!Sample_characteristics_ch1.3`=="gender: Male","M","F")
AA = as.data.frame(str_split_fixed(pheno_GSE161533$`!Sample_characteristics_ch1.2`,":",2))
BB = as.data.frame(str_split_fixed(pheno_GSE161533$`!Sample_characteristics_ch1.4`,":",2))
CC = as.data.frame(str_split_fixed(pheno_GSE161533$`!Sample_characteristics_ch1.5`,":",2))

pheno_GSE161533= cbind(pheno_GSE161533[,c(1,2,9,10)],AA[,2],BB[,2],CC[,2])

colnames(pheno_GSE161533) = c("Sample_title","ID","sample_type","gender","age","smoke","drink")
pheno_GSE161533$organ = "esophagus"

####clean the datMet
pheno_GSE161533$sample_type = factor(pheno_GSE161533$sample_type,levels = c("N","T"))
pheno_GSE161533$Group = as.numeric(ifelse(pheno_GSE161533$sample_type=="T","1","0"))

save(file = "GSE161533_new.RData",GSE161533_anno,pheno_GSE161533)

##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE161533_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE161533$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomright", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE161533_anno)[outliers]); print(table(outliers))
GSE161533_anno = GSE161533_anno[,!outliers]
pheno_GSE161533 = pheno_GSE161533[!outliers,]

#clean the datMet
pheno_GSE161533$sample_type = factor(pheno_GSE161533$sample_type,levels = c("N","T"))
pheno_GSE161533$Group = as.numeric(ifelse(pheno_GSE161533$sample_type=="T","1","0"))

####Adjusting for unknown covariates
mod = model.matrix(~as.factor(Group),data = pheno_GSE161533)
mod0 = model.matrix(~1,data = pheno_GSE161533)
n.sv = num.sv(GSE161533_anno, mod, method="be")
GSE161533_anno = as.matrix(GSE161533_anno)

svobj = sva(GSE161533_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE161533$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 54)

X = svobj$sv
Y = GSE161533_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE161533_anno = GSE161533_anno-t(to_regress)

save(file = "GSE161533.RData",GSE161533_anno,pheno_GSE161533)


