#LUSC_1a_GSE18842 91samples，###Affymetrix
library(affy)
library(limma)
library(biomaRt)
library(sva)
library(stringi)
library(stringr)
library(GEOquery)
library(ggplot2)
library(FactoMineR)
library(WGCNA); library(lumi); library(ggplot2); library(nlme);library(affy)
library(ggplot2); library(Cairo); library(GEOquery);library(nlme);
library(biomaRt); library(sva);library(limma);library(illuminaio)

##选择工作路径
setwd("D:\\SCC_DGEP_shared\\SCC_array\\LUSC\\GSE18842_RAW")
dir = getwd()

##list end of .cel.gz files
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)

basename(cel.files) 
GSE18842_raw<- ReadAffy(filenames = cel.files) 
sampleNames(GSE18842_raw) 
#the length of name maybe 8,9 or 10
sampleNames(GSE18842_raw)<-stri_sub(sampleNames(GSE18842_raw),1,9)
##rma(Robust Multi-Array Average expression measure)：Background correcting，Normalizing，Calculating Expression
GSE18842_raw_rma <- rma(GSE18842_raw)

#
GSE18842 <- exprs(GSE18842_raw_rma) 

###anno the gene symbol
GSE18842_anno <- as.data.frame(GSE18842)
GSE18842_anno$ID<-rownames(GSE18842_anno)
gpl<-getGEO("GPL570",destdir=".")
iddata<-as.data.frame(Table(gpl)[c('ID','Gene Symbol')])
GSE18842_anno<-merge(x=GSE18842_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE18842_anno <- GSE18842_anno[,-1]
###The expression levels of the same gene names were averaged
GSE18842_anno<- aggregate(GSE18842_anno,by = list(GSE18842_anno$`Gene Symbol`),FUN = mean)
head(GSE18842_anno)

GSE18842_anno <- GSE18842_anno[-c(1),]##blank gene name was dropped
rownames(GSE18842_anno) <- GSE18842_anno$Group.1
GSE18842_anno <- GSE18842_anno[,-c(1,93)]
GSE18842_anno[1:2,1:2]
boxplot(GSE18842_anno)

#####prepare the metadata
library(tidyr)
library(dplyr)
library(stringr)
###
###GSE18842
phe_GSE18842 = read.csv("GSE18842_series_matrix.csv", header = F)
phe_GSE18842 = phe_GSE18842[c(38:48),]
phe_GSE18842 = phe_GSE18842 %>% t()
colnames(phe_GSE18842) = phe_GSE18842[1,]
phe_GSE18842 = phe_GSE18842[-1,]
rownames(phe_GSE18842) = phe_GSE18842[,2]
phe_GSE18842 = as.data.frame(phe_GSE18842)

colnames(phe_GSE18842)
pheno_GSE18842 = phe_GSE18842[,c(1,2,11)]
View(pheno_GSE18842)
pheno_GSE18842$sample_type = ifelse(pheno_GSE18842$`!Sample_characteristics_ch1` == "sample type: tumor","T","N")
pheno_GSE18842$organ = "lung"
pheno_GSE18842$batch = "GSE18842"
pheno_GSE18842 = pheno_GSE18842[,c(1,2,4,5,6)]
colnames(pheno_GSE18842) = c("Sample_title","ID","sample_type","organ","batch")
save(pheno_GSE18842,file = "pheno_GSE18842.RData")

#####match metadata and expression data
matchSN = match(colnames(GSE18842_anno), rownames(pheno_GSE18842))
GSE18842_anno = GSE18842_anno[,matchSN]

#clean the datMet
pheno_GSE18842$sample_type = factor(pheno_GSE18842$sample_type,levels = c("N","T"))
pheno_GSE18842$Group = as.numeric(ifelse(pheno_GSE18842$sample_type=="T","1","0"))

save(file = "GSE18842_match.RData",GSE18842_anno,pheno_GSE18842)

##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE18842_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE18842$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomright", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE18842_anno)[outliers]); print(table(outliers))
GSE18842_anno = GSE18842_anno[,!outliers]
pheno_GSE18842 = pheno_GSE18842[!outliers,]

##Regressing Unknown Covariates
mod = model.matrix(~as.factor(Group),data = pheno_GSE18842)
mod0 = model.matrix(~1,data = pheno_GSE18842)
n.sv = num.sv(GSE18842_anno, mod, method="be")
GSE18842_anno = as.matrix(GSE18842_anno)

svobj = sva(GSE18842_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE18842$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 87)

X = svobj$sv
Y = GSE18842_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE18842_anno = GSE18842_anno-t(to_regress)
boxplot(GSE18842_anno)
save(file = "GSE18842.RData",GSE18842_anno,pheno_GSE18842)



















