#LUAD_1c_GSE130779
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

setwd("D:\\SCC_DGEP_shared\\Adenocarcinoma\\EAC\\GSE130779_RAW")
dir = getwd()
# 
txt_files <- list.files(dir, pattern = "\\.txt.gz", full.names = TRUE)

GSE130779_raw <- read.maimages(txt_files,
                               source="agilent", 
                               green.only=TRUE,
                               other.columns = "gIsWellAboveBG")
dim(GSE130779_raw)

##background correction
GSE130779_1 <- backgroundCorrect(GSE130779_raw,method = "normexp") 
##normalization
GSE130779_1 <- normalizeBetweenArrays(GSE130779_1, method="quantile")
class(GSE130779_1)

##Find the control probe
Control <- GSE130779_1$genes$ControlType==1L
table(Control)

##Found that the probe name is blank
NoSymbol <- is.na(GSE130779_1$genes$ProbeName)
table(NoSymbol)

##Find the probes with low expression level, and the expression level of more than half of the samples is greater than the background noise.
IsExpr <- rowSums(GSE130779_1$other$gIsWellAboveBG > 0) >= 8
table(IsExpr)

#finally
GSE130779_1_filt <- GSE130779_1[!Control & !NoSymbol & IsExpr, ]
dim(GSE130779_1_filt)

GSE130779 = GSE130779_1_filt@.Data[[1]]
boxplot(GSE130779) 


#Name the column Sample Number
colnames(GSE130779) = str_extract(colnames(GSE130779),"GSM\\d*")
GSE130779[1:2,1:2]
#Give the row names the probe numbers
rownames(GSE130779)= GSE130779_1_filt$genes$ProbeName

###anno the gene symbol
GSE130779_anno <- as.data.frame(GSE130779)
GSE130779_anno$ID<-rownames(GSE130779_anno)
gpl<-getGEO("GPL20115",destdir=".")
iddata<-as.data.frame(Table(gpl)[c('ID','GeneSymbol')])
GSE130779_anno<-merge(x=GSE130779_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE130779_anno <- GSE130779_anno[,-1]
###The expression levels of the same gene names were averaged
GSE130779_anno<- aggregate(GSE130779_anno,by = list(GSE130779_anno$`GeneSymbol`),FUN = mean)
head(GSE130779_anno)

GSE130779_anno <- GSE130779_anno[-c(1),]##blank gene name was dropped
rownames(GSE130779_anno) <- GSE130779_anno$Group.1
GSE130779_anno <- GSE130779_anno[,-c(1,18)]
GSE130779_anno[1:2,1:2]
boxplot(GSE130779_anno)



#####prepare the metadata
library(tidyr)
library(dplyr)
library(stringr)
###
###GSE130779
phe_GSE130779 = read.csv("GSE130779_series_matrix.csv", header = F)
phe_GSE130779 = phe_GSE130779[c(25:34),]
phe_GSE130779 = phe_GSE130779 %>% t()
colnames(phe_GSE130779) = phe_GSE130779[1,]
phe_GSE130779 = phe_GSE130779[-1,]
rownames(phe_GSE130779) = phe_GSE130779[,2]
phe_GSE130779 = as.data.frame(phe_GSE130779)

colnames(phe_GSE130779)
pheno_GSE130779 = phe_GSE130779[,c(1,2,10)]
View(pheno_GSE130779)

pheno_GSE130779$sample_type = ifelse(pheno_GSE130779$`!Sample_characteristics_ch1`=="tissue: lung adenocarcinoma","T","N")
pheno_GSE130779 = pheno_GSE130779[,c(1,2,4)]
colnames(pheno_GSE130779)=c("Sample_title","ID","sample_type")
pheno_GSE130779$batch = "GSE130779"
save(pheno_GSE130779,file = "pheno_GSE130779.RData")

#####match metadata and expression data
matchSN = match(colnames(GSE130779_anno), rownames(pheno_GSE130779))
GSE130779_anno = GSE130779_anno[,matchSN]


#clean the datMet
pheno_GSE130779$sample_type = factor(pheno_GSE130779$sample_type,levels = c("N","T"))
pheno_GSE130779$Group = as.numeric(ifelse(pheno_GSE130779$sample_type=="T","1","0"))
save(file = "GSE130779_match.RData",GSE130779_anno,pheno_GSE130779)


##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE130779_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE130779$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE130779_anno)[outliers]); print(table(outliers))
GSE130779_anno = GSE130779_anno[,!outliers]
pheno_GSE130779 = pheno_GSE130779[!outliers,]


#Regressing Unknown Covariates
pheno_GSE130779$sample_type = factor(pheno_GSE130779$sample_type,levels = c("N","T"))
pheno_GSE130779$Group = as.numeric(ifelse(pheno_GSE130779$sample_type=="T","1","0"))

mod = model.matrix(~as.factor(Group),data = pheno_GSE130779)
mod0 = model.matrix(~1,data = pheno_GSE130779)
n.sv = num.sv(GSE130779_anno, mod, method="be")
GSE130779_anno = as.matrix(GSE130779_anno)

svobj = sva(GSE130779_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE130779$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 16)

X = svobj$sv
Y = GSE130779_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE130779_anno = GSE130779_anno-t(to_regress)
boxplot(GSE130779_anno)
save(file = "GSE130779.RData",GSE130779_anno,pheno_GSE130779)
