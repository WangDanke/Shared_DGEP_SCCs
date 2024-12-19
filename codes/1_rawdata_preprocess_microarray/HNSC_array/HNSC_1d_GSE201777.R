###HNSC GSE201777 31samples
##选择工作路径
library(affy)
library(limma)
library(biomaRt)
library(sva)
library(stringi)
library(GEOquery)
library(ggplot2)
library(FactoMineR)
setwd("D:\\SCC_DGEP_shared\\SCC_array\\HNSC\\GSE201777")
dir = getwd()
##
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)

basename(cel.files) 
GSE201777_raw<- ReadAffy(filenames = cel.files) 
sampleNames(GSE201777_raw)
#the length of name maybe 8,9 or 10 
sampleNames(GSE201777_raw)<-stri_sub(sampleNames(GSE201777_raw),1,10)
##rma(Robust Multi-Array Average expression measure)：Background correcting，Normalizing，Calculating Expression
GSE201777_raw_rma <- rma(GSE201777_raw)
#
GSE201777 <- exprs(GSE201777_raw_rma) 


##annotation
###anno the gene symbol
GSE201777_anno <- as.data.frame(GSE201777)
GSE201777_anno$ID<-rownames(GSE201777_anno)
gpl571<-getGEO("GPL15207",destdir=".")
iddata<-as.data.frame(Table(gpl571)[c('ID','Gene Symbol')])
GSE201777_anno<-merge(x=GSE201777_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE201777_anno <- GSE201777_anno[,-1]
###The expression levels of the same gene names were averaged
GSE201777_anno<- aggregate(GSE201777_anno,by = list(GSE201777_anno$`Gene Symbol`),FUN = mean)
head(GSE201777_anno)

GSE201777_anno <- GSE201777_anno[-c(1,2),]##blank gene name was dropped
rownames(GSE201777_anno) <- GSE201777_anno$Group.1
GSE201777_anno <- GSE201777_anno[,-c(1,49)]
GSE201777_anno[1:2,1:2]

##match the metadata with the expression data
matchSN = match(colnames(GSE201777_anno), rownames(pheno_GSE201777))
GSE201777_anno = GSE201777_anno[,matchSN]
save(file = "GSE201777.RData",GSE201777_anno,pheno_GSE201777)

##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE201777_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE201777$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE201777_anno)[outliers]); print(table(outliers))
GSE201777_anno = GSE201777_anno[,!outliers]
pheno_GSE201777 = pheno_GSE201777[!outliers,]

#clean the datMet
pheno_GSE201777$sample_type = factor(pheno_GSE201777$sample_type,levels = c("N","T"))
pheno_GSE201777$Group = as.numeric(ifelse(pheno_GSE201777$sample_type=="T","1","0"))

#### Adjustment for covariates
mod = model.matrix(~as.factor(Group),data = pheno_GSE201777)
mod0 = model.matrix(~1,data = pheno_GSE201777)
n.sv = num.sv(GSE201777_anno, mod, method="be")
GSE201777_anno = as.matrix(GSE201777_anno)

svobj = sva(GSE201777_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE201777$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 31)

X = svobj$sv
Y = GSE201777_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE201777_anno = GSE201777_anno-t(to_regress)

boxplot(GSE201777_anno)

save(file = "GSE201777.RData",GSE201777_anno,pheno_GSE201777)
