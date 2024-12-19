##OSCC_1e_GSE31056 
###Affymetrix
library(WGCNA); library(lumi); library(ggplot2); library(nlme)
library(ggplot2); library(Cairo); library(GEOquery);library(nlme);library(stringi)
library(biomaRt); library(sva);library(limma);library(illuminaio);library(affy)
setwd("D:\\SCC_DGEP_shared\\SCC_array\\HNSC\\GSE31056")
dir = getwd()

##
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)

basename(cel.files) 
GSE31056_raw<- ReadAffy(filenames = cel.files)
sampleNames(GSE31056_raw) 
#the length of name maybe 8,9 or 10 
sampleNames(GSE31056_raw)<-stri_sub(sampleNames(GSE31056_raw),1,9)
##rma(Robust Multi-Array Average expression measure)：Background correcting，Normalizing，Calculating Expression
GSE31056_raw_rma <- rma(GSE31056_raw)
#
GSE31056 <- exprs(GSE31056_raw_rma) 
save(GSE31056,file = "GSE31056.RData")


###anno the gene symbol
GSE31056_anno <- as.data.frame(GSE31056)
GSE31056_anno$ID<-rownames(GSE31056_anno)
library(hgu133plus2.db)
ids = toTable(hgu133plus2SYMBOL)##
colnames(ids) = c("ID","Gene Symbol")
GSE31056_anno<-merge(x=GSE31056_anno,y=ids,by='ID',all.x=T,all.y=F)
GSE31056_anno <- GSE31056_anno[,-1]
###The expression levels of the same gene names were averaged
GSE31056_anno<- aggregate(GSE31056_anno,by = list(GSE31056_anno$`Gene Symbol`),FUN = mean)
head(GSE31056_anno)
##blank gene name was dropped
rownames(GSE31056_anno) <- GSE31056_anno$Group.1
GSE31056_anno <- GSE31056_anno[,-c(1,98)]
GSE31056_anno[1:2,1:2]
boxplot(GSE31056_anno)
save(GSE31056_anno,file = "GSE31056_anno.RData")


##match the metadata with the expression data
matchSN = match(colnames(GSE31056_anno), rownames(pheno_GSE31056))
GSE31056_anno = GSE31056_anno[,matchSN]
save(file = "GSE31056.RData",GSE31056_anno,pheno_GSE31056)

##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE31056_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE31056$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE31056_anno)[outliers]); print(table(outliers))
GSE31056_anno = GSE31056_anno[,!outliers]
pheno_GSE31056 = pheno_GSE31056[!outliers,]

#clean the datMet
pheno_GSE31056$sample_type = factor(pheno_GSE31056$sample_type,levels = c("N","T"))
pheno_GSE31056$Group = as.numeric(ifelse(pheno_GSE31056$sample_type=="T","1","0"))

#### Adjustment for unknown covariates
mod = model.matrix(~as.factor(Group),data = pheno_GSE31056)
mod0 = model.matrix(~1,data = pheno_GSE31056)
n.sv = num.sv(GSE31056_anno, mod, method="be")
GSE31056_anno = as.matrix(GSE31056_anno)

svobj = sva(GSE31056_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE31056$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 43)

X = svobj$sv
Y = GSE31056_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE31056_anno = GSE31056_anno-t(to_regress)
boxplot(GSE31056_anno)

save(file = "GSE31056.RData",GSE31056_anno,pheno_GSE31056)

























