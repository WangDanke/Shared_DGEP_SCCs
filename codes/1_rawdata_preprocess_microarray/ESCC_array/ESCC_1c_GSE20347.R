#ESCC，GSE20347，clinical message and raw data process.
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

setwd("D:\\SCC_DGEP_shared\\SCC_array\\ESCC\\GSE20347_RAW")
getwd()
# list the ".CEL" files
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)
basename(cel.files)

library(affy)
GSE20347_raw<- ReadAffy(filenames = cel.files)#
sampleNames(GSE20347_raw)

##rename samples
library(stringi)
sampleNames(GSE20347_raw)<-stri_sub(sampleNames(GSE20347_raw),1,9)##
sampleNames(GSE20347_raw)

####function "RMA" to preprocess
GSE20347_raw_rma <- rma(GSE20347_raw)
GSE20347 <- exprs(GSE20347_raw_rma)



###anno the gene symbol
library(GEOquery)
GSE20347_anno <- as.data.frame(GSE20347)
GSE20347_anno$ID<-rownames(GSE20347_anno)
gpl571<-getGEO("GPL571",destdir=".")
iddata<-as.data.frame(Table(gpl571)[c('ID','Gene Symbol')])
GSE20347_anno<-merge(x=GSE20347_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE20347_anno <- GSE20347_anno[,-1]
###The expression levels of the same gene names were averaged
GSE20347_anno<- aggregate(GSE20347_anno,by = list(GSE20347_anno$`Gene Symbol`),FUN = mean)
head(GSE20347_anno)

GSE20347_anno <- GSE20347_anno[-1,]##blank gene name was dropped
rownames(GSE20347_anno) <- GSE20347_anno$Group.1
GSE20347_anno <- GSE20347_anno[,-c(1,36)]


save(GSE20347_anno,file = "GSE20347_NEW.RData")

###prepare the metadata

library(tidyr)
library(dplyr)
library(stringr)
###
###GSE20347
phe_GSE20347 = read.csv("GSE20347_series_matrix.csv", header = F)
phe_GSE20347 = phe_GSE20347[c(35:44),]
phe_GSE20347 = phe_GSE20347 %>% t()
colnames(phe_GSE20347) = phe_GSE20347[1,]
phe_GSE20347 = phe_GSE20347[-1,]
rownames(phe_GSE20347) = phe_GSE20347[,2]
phe_GSE20347 = as.data.frame(phe_GSE20347)

colnames(phe_GSE20347)
pheno_GSE20347 = phe_GSE20347[,c(1,2,10)]
View(pheno_GSE20347)
pheno_GSE20347$sample_type = ifelse(pheno_GSE20347$`!Sample_characteristics_ch1` == "tissue: normal adjacent esophageal tissue","N","T")
pheno_GSE20347 = pheno_GSE20347[,c(1,2,4)]
colnames(pheno_GSE20347) = c("Sample_title","ID","sample_type")
pheno_GSE20347$organ = "esophagus"

####clean the datMet
pheno_GSE20347$sample_type = factor(pheno_GSE20347$sample_type,levels = c("N","T"))
pheno_GSE20347$Group = as.numeric(ifelse(pheno_GSE20347$sample_type=="T","1","0"))

save(file = "GSE20347_new.RData",GSE20347_anno,pheno_GSE20347)

##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE20347_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE20347$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomright", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE20347_anno)[outliers]); print(table(outliers))
GSE20347_anno = GSE20347_anno[,!outliers]
pheno_GSE20347 = pheno_GSE20347[!outliers,]

#clean the datMet
pheno_GSE20347$sample_type = factor(pheno_GSE20347$sample_type,levels = c("N","T"))
pheno_GSE20347$Group = as.numeric(ifelse(pheno_GSE20347$sample_type=="T","1","0"))

####Adjusting for unknown covariates 
mod = model.matrix(~as.factor(Group),data = pheno_GSE20347)
mod0 = model.matrix(~1,data = pheno_GSE20347)
n.sv = num.sv(GSE20347_anno, mod, method="be")
GSE20347_anno = as.matrix(GSE20347_anno)

svobj = sva(GSE20347_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE20347$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 32)

X = svobj$sv
Y = GSE20347_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE20347_anno = GSE20347_anno-t(to_regress)
boxplot(GSE20347_anno)

save(GSE20347_anno,pheno_GSE20347,file = "GSE20347.RData")


























