##LUAD_1a_GSE140797
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
setwd("D:\\SCC_DGEP_shared\\Adenocarcinoma\\EAC\\GSE140797_RAW")
dir = getwd()
# 
txt_files <- list.files(dir, pattern = "\\.txt.gz", full.names = TRUE)

GSE140797_raw <- read.maimages(txt_files,
                               source="agilent", 
                               green.only=TRUE,
                               other.columns = "gIsWellAboveBG")
dim(GSE140797_raw)

##background correction
GSE140797_1 <- backgroundCorrect(GSE140797_raw,method = "normexp") 
##normalization
GSE140797_1 <- normalizeBetweenArrays(GSE140797_1, method="quantile")
class(GSE140797_1)

##Find the control probe
Control <- GSE140797_1$genes$ControlType==1L
table(Control)

##Found that the probe name is blank
NoSymbol <- is.na(GSE140797_1$genes$ProbeName)
table(NoSymbol)

##Find the probes with low expression level, and the expression level of more than half of the samples is greater than the background noise.
IsExpr <- rowSums(GSE140797_1$other$gIsWellAboveBG > 0) >= 7
table(IsExpr)

#finally
GSE140797_1_filt <- GSE140797_1[!Control & !NoSymbol & IsExpr, ]
dim(GSE140797_1_filt)

#
GSE140797 = GSE140797_1_filt@.Data[[1]]
boxplot(GSE140797) 


#Name the column Sample Number
colnames(GSE140797) = str_extract(colnames(GSE140797),"GSM\\d*")
GSE140797[1:2,1:2]
#Give the row names the probe numbers
rownames(GSE140797)= GSE140797_1_filt$genes$ProbeName

###anno the gene symbol
GSE140797_anno <- as.data.frame(GSE140797)
GSE140797_anno$ID<-rownames(GSE140797_anno)
gpl<-getGEO("GPL13497",destdir=".")
iddata<-as.data.frame(Table(gpl)[c('ID','GENE_SYMBOL')])
GSE140797_anno<-merge(x=GSE140797_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE140797_anno <- GSE140797_anno[,-1]
###The expression levels of the same gene names were averaged
GSE140797_anno<- aggregate(GSE140797_anno,by = list(GSE140797_anno$`GENE_SYMBOL`),FUN = mean)
head(GSE140797_anno)

GSE140797_anno <- GSE140797_anno[-c(1),]##blank gene name was dropped
rownames(GSE140797_anno) <- GSE140797_anno$Group.1
GSE140797_anno <- GSE140797_anno[,-c(1,16)]
GSE140797_anno[1:2,1:2]
boxplot(GSE140797_anno)

#prepare the metadata
###GSE140797
phe_GSE140797 = read.csv("GSE140797_series_matrix.csv", header = F)
phe_GSE140797 = phe_GSE140797[c(29:40),]
phe_GSE140797 = phe_GSE140797 %>% t()
colnames(phe_GSE140797) = phe_GSE140797[1,]
phe_GSE140797 = phe_GSE140797[-1,]
rownames(phe_GSE140797) = phe_GSE140797[,2]
phe_GSE140797 = as.data.frame(phe_GSE140797)

colnames(phe_GSE140797)
pheno_GSE140797 = phe_GSE140797[,c(1,2,10,11,12)]
View(pheno_GSE140797)
pheno_GSE140797$gender = ifelse(pheno_GSE140797$`!Sample_characteristics_ch1.1`=="gender: female","F","M")
pheno_GSE140797$sample_type = ifelse(pheno_GSE140797$`!Sample_characteristics_ch1`=="tissue: lung adenocarcinoma","T","N")
AA = as.data.frame(str_split_fixed(pheno_GSE140797$`!Sample_characteristics_ch1.2`,":",2))

pheno_GSE140797 = cbind(pheno_GSE140797[,c(1,2,6,7)],AA[,2])
colnames(pheno_GSE140797)=c("Sample_title","ID","gender","sample_type","age")
pheno_GSE140797$batch = "GSE140797"
save(pheno_GSE140797,file = "pheno_GSE140797.RData")

#####match metadata and expression data
matchSN = match(colnames(GSE140797_anno), rownames(pheno_GSE140797))
GSE140797_anno = GSE140797_anno[,matchSN]


#clean the datMet
pheno_GSE140797$sample_type = factor(pheno_GSE140797$sample_type,levels = c("N","T"))
pheno_GSE140797$Group = as.numeric(ifelse(pheno_GSE140797$sample_type=="T","1","0"))
save(file = "GSE140797_match.RData",GSE140797_anno,pheno_GSE140797)


##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE140797_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE140797$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE140797_anno)[outliers]); print(table(outliers))
GSE140797_anno = GSE140797_anno[,!outliers]
pheno_GSE140797 = pheno_GSE140797[!outliers,]


##Regressing Unknown Covariates
pheno_GSE140797$sample_type = factor(pheno_GSE140797$sample_type,levels = c("N","T"))
pheno_GSE140797$Group = as.numeric(ifelse(pheno_GSE140797$sample_type=="T","1","0"))

mod = model.matrix(~as.factor(Group),data = pheno_GSE140797)
mod0 = model.matrix(~1,data = pheno_GSE140797)
n.sv = num.sv(GSE140797_anno, mod, method="be")
GSE140797_anno = as.matrix(GSE140797_anno)

svobj = sva(GSE140797_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE140797$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 14)

X = svobj$sv
Y = GSE140797_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE140797_anno = GSE140797_anno-t(to_regress)
boxplot(GSE140797_anno)
save(file = "GSE140797.RData",GSE140797_anno,pheno_GSE140797)

















