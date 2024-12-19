#HNSC_1a_GSE23558 32samples，###Agilent
library(WGCNA); library(lumi); library(ggplot2); library(nlme)
library(ggplot2); library(Cairo); library(GEOquery);library(nlme);
library(biomaRt); library(sva);library(limma);library(illuminaio)
####Agilent 
raw_dir = "D:\\SCC_DGEP_shared\\SCC_array\\HNSC\\GSE23558"
raw_datas = paste0(raw_dir,"/",dir(raw_dir))
setwd("D:\\SCC_DGEP_shared\\SCC_array\\HNSC\\GSE23558")
GSE23558_raw <- read.maimages(raw_datas,
                              source="agilent", 
                              green.only=TRUE,
                              other.columns = "gIsWellAboveBG")
dim(GSE23558_raw)

##background correct
GSE23558_1 <- backgroundCorrect(GSE23558_raw,method = "normexp") 
##quantile normalization
GSE23558_1 <- normalizeBetweenArrays(GSE23558_1, method="quantile")
class(GSE23558_1)

##Find the control probe
Control <- GSE23558_1$genes$ControlType==1L
table(Control)

##Found that the probe name is empty白
NoSymbol <- is.na(GSE23558_1$genes$ProbeName)
table(NoSymbol)

##Find probes with low expression levels, and keep those with more than half of the samples whose expression levels are greater than the background noise
IsExpr <- rowSums(GSE23558_1$other$gIsWellAboveBG > 0) >= 8
table(IsExpr)

#finally
GSE23558_1_filt <- GSE23558_1[!Control & !NoSymbol & IsExpr, ]
dim(GSE23558_1_filt)

#
GSE23558 = GSE23558_1_filt@.Data[[1]]
boxplot(GSE23558) 

#Name the column Sample Number
colnames(GSE23558) = str_extract(colnames(GSE23558),"GSM\\d*")
GSE23558[1:2,1:2]
#Give the row names the probe numbers
rownames(GSE23558)= GSE23558_1_filt$genes$ProbeName

save(GSE23558,file = "GSE23558.RData")



###anno the gene symbol
GSE23558_anno <- as.data.frame(GSE23558)
GSE23558_anno$ID<-rownames(GSE23558_anno)
gpl<-getGEO("GPL10526",destdir=".")
iddata<-as.data.frame(Table(gpl)[c('ID','SYMBOL')])
GSE23558_anno<-merge(x=GSE23558_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE23558_anno <- GSE23558_anno[,-1]
###The expression levels of the same gene names were averaged
GSE23558_anno<- aggregate(GSE23558_anno,by = list(GSE23558_anno$`SYMBOL`),FUN = mean)
head(GSE23558_anno)

GSE23558_anno <- GSE23558_anno[-c(1),]##blank gene name was dropped
rownames(GSE23558_anno) <- GSE23558_anno$Group.1
GSE23558_anno <- GSE23558_anno[,-c(1,98)]
GSE23558_anno[1:2,1:2]
boxplot(GSE23558_anno)
save(GSE23558_anno,file = "GSE23558_anno.RData")

##match the metadata with the expression data
pheno_GSE23558$sample_type = c("T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","N","N","N","N","N")
matchSN = match(colnames(GSE23558_anno), rownames(pheno_GSE23558))
GSE23558_anno = GSE23558_anno[,matchSN]
save(file = "GSE23558.RData",GSE23558_anno,pheno_GSE23558)

##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE23558_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE23558$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE23558_anno)[outliers]); print(table(outliers))
GSE23558_anno = GSE23558_anno[,!outliers]
pheno_GSE23558 = pheno_GSE23558[!outliers,]

#clean the datMet
pheno_GSE23558$sample_type = factor(pheno_GSE23558$sample_type,levels = c("N","T"))
pheno_GSE23558$Group = as.numeric(ifelse(pheno_GSE23558$sample_type=="T","1","0"))
pheno_GSE23558$gender = as.factor(pheno_GSE23558$gender)
pheno_GSE23558$age = as.numeric(pheno_GSE23558$age)
pheno_GSE23558$t_stage = as.numeric(pheno_GSE23558$stage)

#clean the datMet
pheno_GSE23558$sample_type = factor(pheno_GSE23558$sample_type,levels = c("N","T"))
pheno_GSE23558$Group = as.numeric(ifelse(pheno_GSE23558$sample_type=="T","1","0"))

####Adjustment for unknown covariates
mod = model.matrix(~as.factor(Group),data = pheno_GSE23558)
mod0 = model.matrix(~1,data = pheno_GSE23558)
n.sv = num.sv(GSE23558_anno, mod, method="be")
GSE23558_anno = as.matrix(GSE23558_anno)

svobj = sva(GSE23558_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE23558$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 30)

X = svobj$sv
Y = GSE23558_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE23558_anno = GSE23558_anno-t(to_regress)
boxplot(GSE23558_anno)

save(file = "GSE23558.RData",GSE23558_anno,pheno_GSE23558)
















