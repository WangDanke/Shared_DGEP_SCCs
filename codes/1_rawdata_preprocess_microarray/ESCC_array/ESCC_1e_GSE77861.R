#ESCC，GSE77861，clinical message and raw data process.
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

setwd("D:\\SCC_DGEP_shared\\SCC_array\\ESCC\\GSE77861_RAW")
getwd()
# list the ".CEL" files
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)
basename(cel.files)

library(affy)
GSE77861_raw<- ReadAffy(filenames = cel.files)#
sampleNames(GSE77861_raw)

##rename samples
library(stringi)
sampleNames(GSE77861_raw)<-stri_sub(sampleNames(GSE77861_raw),1,10)##
sampleNames(GSE77861_raw)

####function "RMA" to preprocess
GSE77861_raw_rma <- rma(GSE77861_raw)
GSE77861 <- exprs(GSE77861_raw_rma)



###anno the gene symbol
library(GEOquery)
GSE77861_anno <- as.data.frame(GSE77861)
GSE77861_anno$ID<-rownames(GSE77861_anno)
gpl571=getGEO (filename = 'GPL570.soft')
iddata<-as.data.frame(Table(gpl571)[c('ID','Gene Symbol')])
GSE77861_anno<-merge(x=GSE77861_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE77861_anno <- GSE77861_anno[,-1]
###The expression levels of the same gene names were averaged
GSE77861_anno<- aggregate(GSE77861_anno,by = list(GSE77861_anno$`Gene Symbol`),FUN = mean)
head(GSE77861_anno)

GSE77861_anno <- GSE77861_anno[-1,]##blank gene name was dropped
rownames(GSE77861_anno) <- GSE77861_anno$Group.1
GSE77861_anno <- GSE77861_anno[,-c(1,16)]


save(GSE77861_anno,file = "GSE77861_NEW.RData")

###prepare the metadata
library(tidyr)
library(dplyr)
library(stringr)
###
###GSE77861
phe_GSE77861 = read.csv("GSE77861_series_matrix.csv", header = F)
phe_GSE77861 = phe_GSE77861[c(39:59),]
phe_GSE77861 = phe_GSE77861 %>% t()
colnames(phe_GSE77861) = phe_GSE77861[1,]
phe_GSE77861 = phe_GSE77861[-1,]
rownames(phe_GSE77861) = phe_GSE77861[,2]
phe_GSE77861 = as.data.frame(phe_GSE77861)

colnames(phe_GSE77861)
pheno_GSE77861 = phe_GSE77861[,c(1,2,21)]
View(pheno_GSE77861)
AA = as.data.frame(str_split_fixed(pheno_GSE77861$`!Sample_description`,"-",2))
pheno_GSE77861= cbind(pheno_GSE77861[,c(1,2)],AA[,1])
pheno_GSE77861$sample_type= ifelse(pheno_GSE77861$`AA[, 1]`=="Healthy Tissue ","N","T")
pheno_GSE77861=pheno_GSE77861[,c(1,2,4)]

colnames(pheno_GSE77861) = c("Sample_title","ID","sample_type")
pheno_GSE77861$organ = "esophagus"

####clean the datMet
pheno_GSE77861$sample_type = factor(pheno_GSE77861$sample_type,levels = c("N","T"))
pheno_GSE77861$Group = as.numeric(ifelse(pheno_GSE77861$sample_type=="T","1","0"))

save(file = "GSE77861_new.RData",GSE77861_anno,pheno_GSE77861)

##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE77861_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE77861$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomright", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE77861_anno)[outliers]); print(table(outliers))
GSE77861_anno = GSE77861_anno[,!outliers]
pheno_GSE77861 = pheno_GSE77861[!outliers,]

#clean the datMet
pheno_GSE77861$sample_type = factor(pheno_GSE77861$sample_type,levels = c("N","T"))
pheno_GSE77861$Group = as.numeric(ifelse(pheno_GSE77861$sample_type=="T","1","0"))

####Adjusting for unknown covariates
mod = model.matrix(~as.factor(Group),data = pheno_GSE77861)
mod0 = model.matrix(~1,data = pheno_GSE77861)
n.sv = num.sv(GSE77861_anno, mod, method="be")
GSE77861_anno = as.matrix(GSE77861_anno)

svobj = sva(GSE77861_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE77861$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 13)

X = svobj$sv
Y = GSE77861_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE77861_anno = GSE77861_anno-t(to_regress)
boxplot(GSE77861_anno)
save(file = "GSE77861.RData",GSE77861_anno,pheno_GSE77861)


