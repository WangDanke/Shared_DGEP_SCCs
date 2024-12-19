#LUAD_1d_GSE118370
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

setwd("D:\\鳞癌分析\\泛鳞癌array数据\\腺癌\\LUAD\\GSE118370_RAW")
getwd()
# list the ".CEL" files
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)
basename(cel.files)
library(affy)
GSE118370_raw<- ReadAffy(filenames = cel.files)#读入文件
sampleNames(GSE118370_raw)

##rename samples
library(stringi)
sampleNames(GSE118370_raw)<-stri_sub(sampleNames(GSE118370_raw),1,10)##8或9或10
sampleNames(GSE118370_raw)

####function "RMA" to preprocess
GSE118370_raw_rma <- rma(GSE118370_raw)
GSE118370 <- exprs(GSE118370_raw_rma)

###anno the gene symbol
library(GEOquery)
GSE118370_anno <- as.data.frame(GSE118370)
GSE118370_anno$ID<-rownames(GSE118370_anno)
gpl571<-getGEO("GPL570",destdir=".")
iddata<-as.data.frame(Table(gpl571)[c('ID','Gene Symbol')])
GSE118370_anno<-merge(x=GSE118370_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE118370_anno <- GSE118370_anno[,-1]
###The expression levels of the same gene names were averaged
GSE118370_anno<- aggregate(GSE118370_anno,by = list(GSE118370_anno$`Gene Symbol`),FUN = mean)
head(GSE118370_anno)

GSE118370_anno <- GSE118370_anno[-1,]##blank gene name was dropped
rownames(GSE118370_anno) <- GSE118370_anno$Group.1
GSE118370_anno <- GSE118370_anno[,-c(1,14)]

save(GSE118370_anno,file = "GSE118370_anno.RData")

###prepare the mertdata

library(tidyr)
library(dplyr)
library(stringr)
###
###GSE118370
phe_GSE118370 = read.csv("GSE118370_series_matrix.csv", header = F)
phe_GSE118370 = phe_GSE118370[c(27:36),]
phe_GSE118370 = phe_GSE118370 %>% t()
colnames(phe_GSE118370) = phe_GSE118370[1,]
phe_GSE118370 = phe_GSE118370[-1,]
rownames(phe_GSE118370) = phe_GSE118370[,2]
phe_GSE118370 = as.data.frame(phe_GSE118370)

colnames(phe_GSE118370)
pheno_GSE118370 = phe_GSE118370[,c(1,2,10)]
View(pheno_GSE118370)
pheno_GSE118370$sample_type = ifelse(pheno_GSE118370$`!Sample_characteristics_ch1` == "tissue: normal tissue","N","T")
pheno_GSE118370 = pheno_GSE118370[,c(1,2,4)]
colnames(pheno_GSE118370) = c("Sample_title","ID","sample_type")


####clean the datMet
pheno_GSE118370$sample_type = factor(pheno_GSE118370$sample_type,levels = c("N","T"))
pheno_GSE118370$Group = as.numeric(ifelse(pheno_GSE118370$sample_type=="T","1","0"))
#####match metadata and expression data
matchSN = match(colnames(GSE118370_anno), rownames(pheno_GSE118370))
GSE118370_anno = GSE118370_anno[,matchSN]
save(file = "GSE118370_match.RData",GSE118370_anno,pheno_GSE118370)

##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE118370_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE118370$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomright", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE118370_anno)[outliers]); print(table(outliers))
GSE118370_anno = GSE118370_anno[,!outliers]
pheno_GSE118370 = pheno_GSE118370[!outliers,]

#clean the datMet
pheno_GSE118370$sample_type = factor(pheno_GSE118370$sample_type,levels = c("N","T"))
pheno_GSE118370$Group = as.numeric(ifelse(pheno_GSE118370$sample_type=="T","1","0"))

####Regressing Unknown Covariates 
mod = model.matrix(~as.factor(Group),data = pheno_GSE118370)
mod0 = model.matrix(~1,data = pheno_GSE118370)
n.sv = num.sv(GSE118370_anno, mod, method="be")
GSE118370_anno = as.matrix(GSE118370_anno)

svobj = sva(GSE118370_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE118370$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 11)

X = svobj$sv
Y = GSE118370_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE118370_anno = GSE118370_anno-t(to_regress)
boxplot(GSE118370_anno)
save(file = "GSE118370.RData",GSE118370_anno,pheno_GSE118370)