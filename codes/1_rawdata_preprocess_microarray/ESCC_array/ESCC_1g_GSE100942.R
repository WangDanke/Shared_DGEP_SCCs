#ESCC,GSE100942,clinical message and raw data process.
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

setwd("D:\\SCC_DGEP_shared\\SCC_array\\ESCC\\GSE100942_RAW")
getwd()
# list the ".CEL" files
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)
basename(cel.files)

library(affy)
GSE100942_raw<- ReadAffy(filenames = cel.files)#
sampleNames(GSE100942_raw)

##rename samples
library(stringi)
sampleNames(GSE100942_raw)<-stri_sub(sampleNames(GSE100942_raw),1,10)##
sampleNames(GSE100942_raw)

####function "RMA" to preprocess
GSE100942_raw_rma <- rma(GSE100942_raw)
GSE100942 <- exprs(GSE100942_raw_rma)



###anno the gene symbol
library(GEOquery)
GSE100942_anno <- as.data.frame(GSE100942)
GSE100942_anno$ID<-rownames(GSE100942_anno)
gpl571=getGEO (filename = 'GPL570.soft')
iddata<-as.data.frame(Table(gpl571)[c('ID','Gene Symbol')])
GSE100942_anno<-merge(x=GSE100942_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE100942_anno <- GSE100942_anno[,-1]
###The expression levels of the same gene names were averaged
GSE100942_anno<- aggregate(GSE100942_anno,by = list(GSE100942_anno$`Gene Symbol`),FUN = mean)
head(GSE100942_anno)

GSE100942_anno <- GSE100942_anno[-1,]##blank gene name was dropped
rownames(GSE100942_anno) <- GSE100942_anno$Group.1
GSE100942_anno <- GSE100942_anno[,-c(1,10)]


save(GSE100942_anno,file = "GSE100942_NEW.RData")

###prepare the metadata
setwd("D:\\ESCC_Subtype_research\\GSE100942_RAW")
getwd()
library(tidyr)
library(dplyr)
library(stringr)
###
###GSE100942
phe_GSE100942 = read.csv("GSE100942_series_matrix.csv", header = F)
phe_GSE100942 = phe_GSE100942[c(27:36),]
phe_GSE100942 = phe_GSE100942 %>% t()
colnames(phe_GSE100942) = phe_GSE100942[1,]
phe_GSE100942 = phe_GSE100942[-1,]
rownames(phe_GSE100942) = phe_GSE100942[,2]
phe_GSE100942 = as.data.frame(phe_GSE100942)

colnames(phe_GSE100942)
pheno_GSE100942 = phe_GSE100942[,c(1,2,10)]
View(pheno_GSE100942)
pheno_GSE100942$sample_type = ifelse(pheno_GSE100942$`!Sample_characteristics_ch1`=="tissue: Esophagus tumor","T","N")
pheno_GSE100942 = pheno_GSE100942[,c(1,2,4)]

colnames(pheno_GSE100942) = c("Sample_title","ID","sample_type")
pheno_GSE100942$organ = "esophagus"

pheno_GSE100942 = pheno_GSE100942[colnames(GSE100942_anno),]

####clean the datMet
pheno_GSE100942$sample_type = factor(pheno_GSE100942$sample_type,levels = c("N","T"))
pheno_GSE100942$Group = as.numeric(ifelse(pheno_GSE100942$sample_type=="T","1","0"))

save(file = "GSE100942_new.RData",GSE100942_anno,pheno_GSE100942)

##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE100942_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE100942$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomright", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE100942_anno)[outliers]); print(table(outliers))
GSE100942_anno = GSE100942_anno[,!outliers]
pheno_GSE100942 = pheno_GSE100942[!outliers,]

#clean the datMet
pheno_GSE100942$sample_type = factor(pheno_GSE100942$sample_type,levels = c("N","T"))
pheno_GSE100942$Group = as.numeric(ifelse(pheno_GSE100942$sample_type=="T","1","0"))

####Adjusting for unknown covariates
mod = model.matrix(~as.factor(Group),data = pheno_GSE100942)
mod0 = model.matrix(~1,data = pheno_GSE100942)
n.sv = num.sv(GSE100942_anno, mod, method="be")
GSE100942_anno = as.matrix(GSE100942_anno)

svobj = sva(GSE100942_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE100942$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 8)

X = svobj$sv
Y = GSE100942_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE100942_anno = GSE100942_anno-t(to_regress)
boxplot(GSE100942_anno)
save(pheno_GSE100942,GSE100942_anno,file = "GSE100942.RData")

