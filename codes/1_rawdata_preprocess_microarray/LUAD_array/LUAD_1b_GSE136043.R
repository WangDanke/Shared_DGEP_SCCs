#LUAD_1b_GSE136043
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
setwd("D:\\SCC_DGEP_shared\\Adenocarcinoma\\EAC\\GSE136043_RAW")
dir = getwd()
#
txt_files <- list.files(dir, pattern = "\\.txt.gz", full.names = TRUE)

GSE136043_raw <- read.maimages(txt_files,
                               source="agilent", 
                               green.only=TRUE,
                               other.columns = "gIsWellAboveBG")
dim(GSE136043_raw)
##background correction
GSE136043_1 <- backgroundCorrect(GSE136043_raw,method = "normexp") 
##normalization
GSE136043_1 <- normalizeBetweenArrays(GSE136043_1, method="quantile")
class(GSE136043_1)

###Find the control probe
Control <- GSE136043_1$genes$ControlType==1L
table(Control)

##Found that the probe name is blank
NoSymbol <- is.na(GSE136043_1$genes$ProbeName)
table(NoSymbol)

##Find the probes with low expression level, and the expression level of more than half of the samples is greater than the background noise.
IsExpr <- rowSums(GSE136043_1$other$gIsWellAboveBG > 0) >= 5
table(IsExpr)

#finally
GSE136043_1_filt <- GSE136043_1[!Control & !NoSymbol & IsExpr, ]
dim(GSE136043_1_filt)
GSE136043 = GSE136043_1_filt@.Data[[1]]
boxplot(GSE136043) #


#Name the column Sample Number
colnames(GSE136043) = str_extract(colnames(GSE136043),"GSM\\d*")
GSE136043[1:2,1:2]
#Give the row names the probe numbers
rownames(GSE136043)= GSE136043_1_filt$genes$ProbeName


###anno the gene symbol
GSE136043_anno <- as.data.frame(GSE136043)
GSE136043_anno$ID<-rownames(GSE136043_anno)
gpl<-getGEO("GPL13497",destdir=".")
iddata<-as.data.frame(Table(gpl)[c('ID','GENE_SYMBOL')])
GSE136043_anno<-merge(x=GSE136043_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE136043_anno <- GSE136043_anno[,-1]
###The expression levels of the same gene names were averaged
GSE136043_anno<- aggregate(GSE136043_anno,by = list(GSE136043_anno$`GENE_SYMBOL`),FUN = mean)
head(GSE136043_anno)

GSE136043_anno <- GSE136043_anno[-c(1),]##blank gene name was dropped
rownames(GSE136043_anno) <- GSE136043_anno$Group.1
GSE136043_anno <- GSE136043_anno[,-c(1,12)]
GSE136043_anno[1:2,1:2]
boxplot(GSE136043_anno)



#####prepare the metadata
library(tidyr)
library(dplyr)
library(stringr)

###GSE136043
phe_GSE136043 = read.csv("GSE136043_series_matrix .csv", header = F)
phe_GSE136043 = phe_GSE136043[c(27:41),]
phe_GSE136043 = phe_GSE136043 %>% t()
colnames(phe_GSE136043) = phe_GSE136043[1,]
phe_GSE136043 = phe_GSE136043[-1,]
rownames(phe_GSE136043) = phe_GSE136043[,2]
phe_GSE136043 = as.data.frame(phe_GSE136043)

colnames(phe_GSE136043)
pheno_GSE136043 = phe_GSE136043[,c(1,2,10,13,14,15)]
View(pheno_GSE136043)
pheno_GSE136043$gender = ifelse(pheno_GSE136043$`!Sample_characteristics_ch1.2`=="gender: Female","F","M")
pheno_GSE136043$sample_type = ifelse(pheno_GSE136043$`!Sample_characteristics_ch1`=="tissue: non-tumor lung tissue","N","T")
pheno_GSE136043 = pheno_GSE136043[,c(1,2,4,6,7,8)]
colnames(pheno_GSE136043)=c("Sample_title","ID","age","t_stage","gender","sample_type")
pheno_GSE136043$batch = "GSE136043"
save(pheno_GSE136043,file = "pheno_GSE136043.RData")

#####match metadata and expression data
matchSN = match(colnames(GSE136043_anno), rownames(pheno_GSE136043))
GSE136043_anno = GSE136043_anno[,matchSN]


#clean the datMet
pheno_GSE136043$sample_type = factor(pheno_GSE136043$sample_type,levels = c("N","T"))
pheno_GSE136043$Group = as.numeric(ifelse(pheno_GSE136043$sample_type=="T","1","0"))
save(file = "GSE136043_match.RData",GSE136043_anno,pheno_GSE136043)


##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE136043_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE136043$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE136043_anno)[outliers]); print(table(outliers))
GSE136043_anno = GSE136043_anno[,!outliers]
pheno_GSE136043 = pheno_GSE136043[!outliers,]


##Regressing Unknown Covariates
pheno_GSE136043$sample_type = factor(pheno_GSE136043$sample_type,levels = c("N","T"))
pheno_GSE136043$Group = as.numeric(ifelse(pheno_GSE136043$sample_type=="T","1","0"))

mod = model.matrix(~as.factor(Group),data = pheno_GSE136043)
mod0 = model.matrix(~1,data = pheno_GSE136043)
n.sv = num.sv(GSE136043_anno, mod, method="be")
GSE136043_anno = as.matrix(GSE136043_anno)

svobj = sva(GSE136043_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE136043$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 9)

X = svobj$sv
Y = GSE136043_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE136043_anno = GSE136043_anno-t(to_regress)
boxplot(GSE136043_anno)
save(file = "GSE136043.RData",GSE136043_anno,pheno_GSE136043)
