#ESCC，GSE38129，clinical message and raw data process.
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

setwd("D:\\SCC_DGEP_shared\\SCC_array\\ESCC\\GSE38129_RAW")
getwd()
# list the ".CEL" files
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)
basename(cel.files)
library(affy)
GSE38129_raw<- ReadAffy(filenames = cel.files)#读入文件
sampleNames(GSE38129_raw)

##rename samples
library(stringi)
sampleNames(GSE38129_raw)<-stri_sub(sampleNames(GSE38129_raw),1,9)
sampleNames(GSE38129_raw)

####function "RMA" to preprocess
GSE38129_raw_rma <- rma(GSE38129_raw)
GSE38129 <- exprs(GSE38129_raw_rma)

###anno the gene symbol
library(GEOquery)
GSE38129_anno <- as.data.frame(GSE38129)
GSE38129_anno$ID<-rownames(GSE38129_anno)
gpl571<-getGEO("GPL571",destdir=".")
iddata<-as.data.frame(Table(gpl571)[c('ID','Gene Symbol')])
GSE38129_anno<-merge(x=GSE38129_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE38129_anno <- GSE38129_anno[,-1]
###The expression levels of the same gene names were averaged
GSE38129_anno<- aggregate(GSE38129_anno,by = list(GSE38129_anno$`Gene Symbol`),FUN = mean)
head(GSE38129_anno)

GSE38129_anno <- GSE38129_anno[-1,]##blank gene name was dropped
rownames(GSE38129_anno) <- GSE38129_anno$Group.1
GSE38129_anno <- GSE38129_anno[,-c(1,62)]


save(GSE38129_anno,file = "GSE38129_NEW.RData")

###prepare the meta data
library(tidyr)
library(dplyr)
library(stringr)
###
###GSE38129
phe_GSE38129 = read.csv("GSE38129_series_matrix.csv", header = F)
phe_GSE38129 = phe_GSE38129[c(30:40),]
phe_GSE38129 = phe_GSE38129 %>% t()
colnames(phe_GSE38129) = phe_GSE38129[1,]
phe_GSE38129 = phe_GSE38129[-1,]
rownames(phe_GSE38129) = phe_GSE38129[,2]
phe_GSE38129 = as.data.frame(phe_GSE38129)

colnames(phe_GSE38129)
pheno_GSE38129 = phe_GSE38129[,c(1,2,8)]
View(pheno_GSE38129)
pheno_GSE38129$sample_type = ifelse(pheno_GSE38129$`!Sample_source_name_ch1` == "adjacent normal tissue","N","T")
pheno_GSE38129 = pheno_GSE38129[,c(1,2,4)]
colnames(pheno_GSE38129) = c("Sample_title","ID","sample_type")
pheno_GSE38129$organ = "esophagus"

####clean the datMet
pheno_GSE38129$sample_type = factor(pheno_GSE38129$sample_type,levels = c("N","T"))
pheno_GSE38129$Group = as.numeric(ifelse(pheno_GSE38129$sample_type=="T","1","0"))

save(file = "GSE38129_new.RData",GSE38129_anno,pheno_GSE38129)

##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE38129_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE38129$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomright", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE38129_anno)[outliers]); print(table(outliers))
GSE38129_anno = GSE38129_anno[,!outliers]
pheno_GSE38129 = pheno_GSE38129[!outliers,]

#clean the datMet
pheno_GSE38129$sample_type = factor(pheno_GSE38129$sample_type,levels = c("N","T"))
pheno_GSE38129$Group = as.numeric(ifelse(pheno_GSE38129$sample_type=="T","1","0"))

####Adjusting for unknown covariates 
mod = model.matrix(~as.factor(Group),data = pheno_GSE38129)
mod0 = model.matrix(~1,data = pheno_GSE38129)
n.sv = num.sv(GSE38129_anno, mod, method="be")
GSE38129_anno = as.matrix(GSE38129_anno)

svobj = sva(GSE38129_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE38129$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 55)

X = svobj$sv
Y = GSE38129_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE38129_anno = GSE38129_anno-t(to_regress)
boxplot(GSE38129_anno)
save(file = "GSE38129.RData",GSE38129_anno,pheno_GSE38129)




























