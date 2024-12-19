#ESCC，GSE17351，clinical message and raw data process.
###loading packages
BiocManager::install("biomaRt")
BiocManager::install("FactoMineR")
library(WGCNA)
library(affy)
library(limma)
library(biomaRt)
library(sva)
library(stringi)
library(GEOquery)
library(PCA)
library(ggplot2)
library(FactoMineR)

setwd("D:\\SCC_DGEP_shared\\SCC_array\\ESCC\\GSE17351_RAW")
getwd()
# list the ".CEL" files
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)
basename(cel.files)
GSE17351_raw<- ReadAffy(filenames = cel.files)#读入文件
sampleNames(GSE17351_raw)

##rename samples
sampleNames(GSE17351_raw)<-stri_sub(sampleNames(GSE17351_raw),1,9)##
sampleNames(GSE17351_raw)

####function "RMA" to preprocess
GSE17351_raw_rma <- rma(GSE17351_raw)
GSE17351 <- exprs(GSE17351_raw_rma)

###gene annotation
GSE17351_anno <- as.data.frame(GSE17351)
GSE17351_anno$ID<-rownames(GSE17351_anno)
gpl570<-getGEO("GPL570",destdir=".")
iddata<-as.data.frame(Table(gpl570)[c('ID','Gene Symbol')])
GSE17351_anno<-merge(x=GSE17351_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE17351_anno <- GSE17351_anno[,-1]
###The expression levels of genes with the same gene name were averaged
GSE17351_anno<- aggregate(GSE17351_anno,by = list(GSE17351_anno$`Gene Symbol`),FUN = mean)
head(GSE17351_anno)
#Group.1 GSM433786 GSM433787 GSM433788 GSM433789 GSM433790 GSM433791 GSM433792 GSM433793 GSM433794
#1           5.057084  5.015055  4.923569  5.103173  4.908568  4.823771  4.857703  4.858489  4.884393
#2     A1BG  5.969950  5.974494  5.924251  7.201177  6.272978  5.803666  5.890439  6.503881  6.157236
#3 A1BG-AS1  5.304607  5.007162  4.770723  5.618097  5.301334  4.663768  5.073829  5.102902  4.967380
#4     A1CF  4.760489  4.571400  4.829012  4.745067  4.801011  4.558535  4.837814  4.739964  5.018903
#5      A2M  8.270933  7.301832  8.776535  8.291396  8.742367  8.263867  8.628180  8.722209  8.686100
#6  A2M-AS1  5.266746  4.252222  7.447606  6.074372  5.916305  5.997484  5.759067  5.351163  5.499208
GSE17351_anno <- GSE17351_anno[-1,]
rownames(GSE17351_anno) <- GSE17351_anno$Group.1
GSE17351_anno <- GSE17351_anno[,-c(1,12)]

save(GSE17351_anno,file = "GSE17351_NEW.RData")

###prepare the metadata
library(tidyr)
library(dplyr)
library(stringr)
###
###GSE17351
phe_GSE17351 = read.csv("GSE17351_series_matrix.csv", header = F)
phe_GSE17351 = phe_GSE17351[c(30:41),]
phe_GSE17351 = phe_GSE17351 %>% t()
colnames(phe_GSE17351) = phe_GSE17351[1,]
phe_GSE17351 = phe_GSE17351[-1,]
rownames(phe_GSE17351) = phe_GSE17351[,2]
phe_GSE17351 = as.data.frame(phe_GSE17351)

colnames(phe_GSE17351)
pheno_GSE17351 = phe_GSE17351[,c(1,2,10,11,12)]
View(pheno_GSE17351)
pheno_GSE17351$gender = ifelse(pheno_GSE17351$`!Sample_characteristics_ch1.1` == "gender: male, Japanese","M","F")
pheno_GSE17351$sample_type = ifelse(pheno_GSE17351$`!Sample_characteristics_ch1` == "tissue: esophagus, normal","N","T")

BB = as.data.frame(str_split_fixed(pheno_GSE17351$`!Sample_characteristics_ch1.2`,":",2))
pheno_GSE17351 = cbind(pheno_GSE17351[,c(1,2,6,7)],BB[,2])
colnames(pheno_GSE17351) = c("Sample_title","ID","gender","sample_type","age")
pheno_GSE17351$organ = "esophagus"

####clean the datMet
pheno_GSE17351$sample_type = factor(pheno_GSE17351$sample_type,levels = c("N","T"))
pheno_GSE17351$Group = as.numeric(ifelse(pheno_GSE17351$sample_type=="T","1","0"))

save(file = "GSE17351_new.RData",GSE17351_anno,pheno_GSE17351)

##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE17351_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE17351$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomright", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE17351_anno)[outliers]); print(table(outliers))
GSE17351_anno = GSE17351_anno[,!outliers]
pheno_GSE17351 = pheno_GSE17351[!outliers,]

###Adjusting for unknown covariates 
mod = model.matrix(~as.factor(Group),data = pheno_GSE17351)
mod0 = model.matrix(~1,data = pheno_GSE17351)
n.sv = num.sv(GSE17351_anno, mod, method="be")
GSE17351_anno = as.matrix(GSE17351_anno)

svobj = sva(GSE17351_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE17351$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 10)

X = svobj$sv
Y = GSE17351_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE17351_anno = GSE17351_anno-t(to_regress)
boxplot(GSE17351_anno)

save(file = "GSE17351.RData",GSE17351_anno,pheno_GSE17351)





























