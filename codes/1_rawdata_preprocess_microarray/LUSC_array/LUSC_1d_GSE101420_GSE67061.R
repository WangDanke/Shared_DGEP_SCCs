###GSE102420 and GSE67061 are source from the same Study, so combined them
#LUSC_GSE101420 60samples，###Agilent
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
setwd("D:\\SCC_DGEP_shared\\SCC_array\\LUSC\\GSE101420_RAW")
dir = getwd()

# 
txt_files <- list.files(dir, pattern = "\\.txt.gz", full.names = TRUE)

GSE101420_raw <- read.maimages(txt_files,
                               source="agilent", 
                               green.only=TRUE,
                               other.columns = "gIsWellAboveBG")
dim(GSE101420_raw)


GSE101420_1 <- backgroundCorrect(GSE101420_raw,method = "normexp") 

GSE101420_1 <- normalizeBetweenArrays(GSE101420_1, method="quantile")
class(GSE101420_1)

Control <- GSE101420_1$genes$ControlType==1L
table(Control)

NoSymbol <- is.na(GSE101420_1$genes$ProbeName)
table(NoSymbol)

IsExpr <- rowSums(GSE101420_1$other$gIsWellAboveBG > 0) >= 30
table(IsExpr)

GSE101420_1_filt <- GSE101420_1[!Control & !NoSymbol & IsExpr, ]
dim(GSE101420_1_filt)

GSE101420 = GSE101420_1_filt@.Data[[1]]
boxplot(GSE101420) 

colnames(GSE101420) = str_extract(colnames(GSE101420),"GSM\\d*")
GSE101420[1:2,1:2]

rownames(GSE101420)= GSE101420_1_filt$genes$ProbeName

###anno the gene symbol
GSE101420_anno <- as.data.frame(GSE101420)
GSE101420_anno$ID<-rownames(GSE101420_anno)
gpl<-getGEO("GPL6480",destdir=".")
iddata<-as.data.frame(Table(gpl)[c('ID','GENE_SYMBOL')])
GSE101420_anno<-merge(x=GSE101420_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE101420_anno <- GSE101420_anno[,-1]
###The expression levels of the same gene names were averaged
GSE101420_anno<- aggregate(GSE101420_anno,by = list(GSE101420_anno$`GENE_SYMBOL`),FUN = mean)
head(GSE101420_anno)

GSE101420_anno <- GSE101420_anno[-c(1),]##blank gene name was dropped
rownames(GSE101420_anno) <- GSE101420_anno$Group.1
GSE101420_anno <- GSE101420_anno[,-c(1,62)]
GSE101420_anno[1:2,1:2]
boxplot(GSE101420_anno)



#####prepare the metadata
library(tidyr)
library(dplyr)
library(stringr)
###
###GSE101420
phe_GSE101420 = read.csv("GSE101420_series_matrix.csv", header = F)
phe_GSE101420 = phe_GSE101420[c(28:40),]
phe_GSE101420 = phe_GSE101420 %>% t()
colnames(phe_GSE101420) = phe_GSE101420[1,]
phe_GSE101420 = phe_GSE101420[-1,]
rownames(phe_GSE101420) = phe_GSE101420[,2]
phe_GSE101420 = as.data.frame(phe_GSE101420)

colnames(phe_GSE101420)
pheno_GSE101420 = phe_GSE101420[,c(1,2,11)]
View(pheno_GSE101420)
AA = as.data.frame(str_split_fixed(pheno_GSE101420$`!Sample_characteristics_ch1`,":",2))
BB = as.data.frame(str_split_fixed(pheno_GSE101420$`!Sample_characteristics_ch1.1`,":",2))
CC = as.data.frame(str_split_fixed(pheno_GSE101420$`!Sample_characteristics_ch1.2`,":",2))
pheno_GSE101420 = cbind(pheno_GSE101420[,c(1,2)],AA[,2],BB[,2],CC[,2])
colnames(pheno_GSE101420)=c("Sample_title","ID","gender","age","t_stage")

pheno_GSE101420$sample_type ="N"
pheno_GSE101420$organ = "lung"
pheno_GSE101420$batch = "GSE101420"
save(pheno_GSE101420,file = "pheno_GSE101420.RData")

#####match metadata and expression data
matchSN = match(colnames(GSE101420_anno), rownames(pheno_GSE101420))
GSE101420_anno = GSE101420_anno[,matchSN]

#clean the datMet
pheno_GSE101420$sample_type = factor(pheno_GSE101420$sample_type,levels = c("N","T"))
pheno_GSE101420$Group = as.numeric(ifelse(pheno_GSE101420$sample_type=="T","1","0"))

save(file = "GSE101420_match.RData",GSE101420_anno,pheno_GSE101420)

####GSE67061整理
setwd("D:\\SCC_DGEP_shared\\SCC_array\\LUSC\\GSE67061_RAW")
dir = getwd()

#####prepare the metadata
library(tidyr)
library(dplyr)
library(stringr)
###
###GSE67061
phe_GSE67061 = read.csv("GSE67061-GPL6480_series_matrix.csv", header = F)
phe_GSE67061 = phe_GSE67061[c(39:54),]
phe_GSE67061 = phe_GSE67061 %>% t()
colnames(phe_GSE67061) = phe_GSE67061[1,]
phe_GSE67061 = phe_GSE67061[-1,]
rownames(phe_GSE67061) = phe_GSE67061[,2]
phe_GSE67061 = as.data.frame(phe_GSE67061)

colnames(phe_GSE67061)
pheno_GSE67061 = phe_GSE67061[-c(1:8),]##Only the squamous cell carcinoma samples are selected
pheno_GSE67061 = pheno_GSE67061[,c(1,2,11,12,13,16)]

View(pheno_GSE67061)
pheno_GSE67061$gender = ifelse(pheno_GSE67061$`!Sample_characteristics_ch1.1`=="gender: Male","M","F")
pheno_GSE67061$sample_type = ifelse(pheno_GSE67061$`!Sample_characteristics_ch1`=="tissue: immortalized human bronchial epithelial cell lines","cell_line","T")
AA = as.data.frame(str_split_fixed(pheno_GSE67061$`!Sample_characteristics_ch1.2`,":",2))
BB = as.data.frame(str_split_fixed(pheno_GSE67061$`!Sample_characteristics_ch1.3`,":",2))
pheno_GSE67061 = cbind(pheno_GSE67061[,c(1,2,7,8)],AA[,2],BB[,2])
colnames(pheno_GSE67061)=c("Sample_title","ID","gender","sample_type","age","t_stage")

pheno_GSE67061$organ = "lung"
pheno_GSE67061$batch = "GSE67061"
save(pheno_GSE67061,file = "pheno_GSE67061.RData")

#list txt files
txt_files <- list.files(dir, pattern = "\\.txt.gz", full.names = TRUE)

GSE67061_raw <- read.maimages(txt_files,
                               source="agilent", 
                               green.only=TRUE,
                               other.columns = "gIsWellAboveBG")
dim(GSE67061_raw)

GSE67061_1 <- backgroundCorrect(GSE67061_raw,method = "normexp") 

GSE67061_1 <- normalizeBetweenArrays(GSE67061_1, method="quantile")
class(GSE67061_1)

Control <- GSE67061_1$genes$ControlType==1L
table(Control)

NoSymbol <- is.na(GSE67061_1$genes$ProbeName)
table(NoSymbol)

IsExpr <- rowSums(GSE67061_1$other$gIsWellAboveBG > 0) >= 34
table(IsExpr)

GSE67061_1_filt <- GSE67061_1[!Control & !NoSymbol & IsExpr, ]
dim(GSE67061_1_filt)

GSE67061 = GSE67061_1_filt@.Data[[1]]
boxplot(GSE67061) 

colnames(GSE67061) = str_extract(colnames(GSE67061),"GSM\\d*")
GSE67061[1:2,1:2]

rownames(GSE67061)= GSE67061_1_filt$genes$ProbeName

###anno the gene symbol
GSE67061_anno <- as.data.frame(GSE67061)
GSE67061_anno$ID<-rownames(GSE67061_anno)
gpl<-getGEO("GPL6480",destdir=".")
iddata<-as.data.frame(Table(gpl)[c('ID','GENE_SYMBOL')])
GSE67061_anno<-merge(x=GSE67061_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE67061_anno <- GSE67061_anno[,-1]
###The expression levels of the same gene names were averaged
GSE67061_anno<- aggregate(GSE67061_anno,by = list(GSE67061_anno$`GENE_SYMBOL`),FUN = mean)
head(GSE67061_anno)

GSE67061_anno <- GSE67061_anno[-c(1),]##blank gene name was dropped
rownames(GSE67061_anno) <- GSE67061_anno$Group.1
GSE67061_anno <- GSE67061_anno[,-c(1,79)]
GSE67061_anno[1:2,1:2]
boxplot(GSE67061_anno)



#####match metadata and expression data
matchSN = match(colnames(GSE67061_anno), rownames(pheno_GSE67061))
GSE67061_anno = GSE67061_anno[,matchSN]

save(file = "GSE67061_match.RData",GSE67061_anno,pheno_GSE67061)


###combined GSE101420 and GSE67610
gene = intersect(rownames(GSE101420_anno),rownames(GSE67061_anno))
sel = c("Sample_title","ID","gender","age","t_stage","sample_type","organ","batch")
pheno_GSE67061_com = rbind(pheno_GSE101420[,sel],pheno_GSE67061[,sel])
GSE67061_com_anno = cbind(GSE67061_anno[gene,],GSE101420_anno[gene,])
boxplot(GSE67061_com_anno)


##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE67061_com_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE67061_com$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE67061_com_anno)[outliers]); print(table(outliers))
GSE67061_com_anno = GSE67061_com_anno[,!outliers]
pheno_GSE67061_com = pheno_GSE67061_com[!outliers,]


##Regression  unknown Covariates
pheno_GSE67061_com$sample_type = factor(pheno_GSE67061_com$sample_type,levels = c("N","T"))
pheno_GSE67061_com$Group = as.numeric(ifelse(pheno_GSE67061_com$sample_type=="T","1","0"))

mod = model.matrix(~as.factor(Group),data = pheno_GSE67061_com)
mod0 = model.matrix(~1,data = pheno_GSE67061_com)
n.sv = num.sv(GSE67061_com_anno, mod, method="be")
GSE67061_com_anno = as.matrix(GSE67061_com_anno)

svobj = sva(GSE67061_com_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE67061_com$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 123)

X = svobj$sv
Y = GSE67061_com_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE67061_com_anno = GSE67061_com_anno-t(to_regress)
boxplot(GSE67061_com_anno)
save(file = "GSE67061_com.RData",GSE67061_com_anno,pheno_GSE67061_com)



















