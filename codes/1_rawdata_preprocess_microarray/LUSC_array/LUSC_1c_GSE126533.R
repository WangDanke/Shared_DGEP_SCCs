#LUSC_1c_GSE126533, Agilent
library(stringi)
library(stringr)
library(FactoMineR)
library(WGCNA); library(lumi); library(ggplot2); library(nlme);library(affy)
library(ggplot2); library(Cairo); library(GEOquery);library(nlme);
library(biomaRt); library(sva);library(limma);library(illuminaio);


raw_dir = "D:\\SCC_DGEP_shared\\SCC_array\\LUSC\\GSE126533"
raw_datas = paste0(raw_dir,"/",dir(raw_dir))
setwd("D:\\SCC_DGEP_shared\\SCC_array\\LUSC\\GSE126533")



library(affy)
# 
txt_files <- list.files("D:\\鳞癌分析\\泛鳞癌array数据\\LUSC\\GSE126533", pattern = "\\.txt.gz", full.names = TRUE)

GSE126533_raw <- read.maimages(txt_files,
                               source="agilent", 
                               green.only=TRUE,
                               other.columns = "gIsWellAboveBG")
dim(GSE126533_raw)


GSE126533_1 <- backgroundCorrect(GSE126533_raw,method = "normexp") 

GSE126533_1 <- normalizeBetweenArrays(GSE126533_1, method="quantile")
class(GSE126533_1)

Control <- GSE126533_1$genes$ControlType==1L
table(Control)

NoSymbol <- is.na(GSE126533_1$genes$ProbeName)
table(NoSymbol)

IsExpr <- rowSums(GSE126533_1$other$gIsWellAboveBG > 0) >= 5
table(IsExpr)

GSE126533_1_filt <- GSE126533_1[!Control & !NoSymbol & IsExpr, ]
dim(GSE126533_1_filt)

GSE126533 = GSE126533_1_filt@.Data[[1]]
boxplot(GSE126533) #

colnames(GSE126533) = str_extract(colnames(GSE126533),"GSM\\d*")
GSE126533[1:2,1:2]

rownames(GSE126533)= GSE126533_1_filt$genes$ProbeName

###anno the gene symbol
GSE126533_anno <- as.data.frame(GSE126533)
GSE126533_anno$ID<-rownames(GSE126533_anno)

#use the annotation files download from the website
gpl<-read.csv("GPL22120-25936.csv")
gpl = gpl[-c(1:32),]
colnames(gpl) = gpl[1,]
gpl = gpl[-1,]
iddata = gpl[,c("ID","GENE_SYMBOL")]

GSE126533_anno<-merge(x=GSE126533_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE126533_anno <- GSE126533_anno[,-1]
###The expression levels of the same gene names were averaged
GSE126533_anno<- aggregate(GSE126533_anno,by = list(GSE126533_anno$`GENE_SYMBOL`),FUN = mean)
head(GSE126533_anno)

GSE126533_anno <- GSE126533_anno[-c(1),]##blank gene name was dropped
rownames(GSE126533_anno) <- GSE126533_anno$Group.1
GSE126533_anno <- GSE126533_anno[,-c(1,12)]
GSE126533_anno[1:2,1:2]
boxplot(GSE126533_anno)

save(GSE126533_anno,file = "GSE126533_anno.RData")

##match the metadata with the expression data
#clean the datMet
pheno_GSE126533$sample_type = factor(pheno_GSE126533$sample_type,levels = c("N","T"))
pheno_GSE126533$Group = as.numeric(ifelse(pheno_GSE126533$sample_type=="T","1","0"))

matchSN = match(colnames(GSE126533_anno), rownames(pheno_GSE126533))
GSE126533_anno = GSE126533_anno[,matchSN]
save(file = "GSE126533.RData",GSE126533_anno,pheno_GSE126533)


##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE126533_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE126533$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE126533_anno)[outliers]); print(table(outliers))
GSE126533_anno = GSE126533_anno[,!outliers]
pheno_GSE126533 = pheno_GSE126533[!outliers,]

###Regressing Unknown Covariates
mod = model.matrix(~as.factor(Group),data = pheno_GSE126533)
mod0 = model.matrix(~1,data = pheno_GSE126533)
n.sv = num.sv(GSE126533_anno, mod, method="be")
GSE126533_anno = as.matrix(GSE126533_anno)

svobj = sva(GSE126533_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE126533$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 9)

X = svobj$sv
Y = GSE126533_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE126533_anno = GSE126533_anno-t(to_regress)
boxplot(GSE126533_anno) 

save(file = "GSE126533.RData",GSE126533_anno,pheno_GSE126533)











