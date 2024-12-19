###HNSC_1e_GSE172120, 6samples
library(affy)
library(limma)
library(biomaRt)
library(sva)
library(stringi)
library(GEOquery)
library(ggplot2)
library(FactoMineR)
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\HNSC\\GSE172120")
dir = getwd()

##
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)

basename(cel.files) 
GSE172120_raw<- ReadAffy(filenames = cel.files) 
sampleNames(GSE172120_raw) 
#the length of name maybe 8,9 or 10 
sampleNames(GSE172120_raw)<-stri_sub(sampleNames(GSE172120_raw),1,10)
##rma(Robust Multi-Array Average expression measure)：Background correcting，Normalizing，Calculating Expression
GSE172120_raw_rma <- rma(GSE172120_raw)
#
GSE172120 <- exprs(GSE172120_raw_rma) 


###anno the gene symbol
GSE172120_anno <- as.data.frame(GSE172120)
GSE172120_anno$ID<-rownames(GSE172120_anno)
gpl571<-getGEO("GPL15207",destdir=".")
iddata<-as.data.frame(Table(gpl571)[c('ID','Gene Symbol')])
GSE172120_anno<-merge(x=GSE172120_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE172120_anno <- GSE172120_anno[,-1]
###The expression levels of the same gene names were averaged
GSE172120_anno<- aggregate(GSE172120_anno,by = list(GSE172120_anno$`Gene Symbol`),FUN = mean)
head(GSE172120_anno)

GSE172120_anno <- GSE172120_anno[-c(1,2),]##blank gene name was dropped
rownames(GSE172120_anno) <- GSE172120_anno$Group.1
GSE172120_anno <- GSE172120_anno[,-c(1,49)]
GSE172120_anno[1:2,1:2]

##match the metadata with the expression data
matchSN = match(colnames(GSE172120_anno), rownames(pheno_GSE172120))
GSE172120_anno = GSE172120_anno[,matchSN]
save(file = "GSE172120.RData",GSE172120_anno,pheno_GSE172120)

##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE172120_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE172120$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE172120_anno)[outliers]); print(table(outliers))
GSE172120_anno = GSE172120_anno[,!outliers]
pheno_GSE172120 = pheno_GSE172120[!outliers,]

#clean the datMet
pheno_GSE172120$sample_type = factor(pheno_GSE172120$sample_type,levels = c("N","T"))
pheno_GSE172120$Group = as.numeric(ifelse(pheno_GSE172120$sample_type=="T","1","0"))

#### Adjustment for covariates
mod = model.matrix(~as.factor(Group),data = pheno_GSE172120)
mod0 = model.matrix(~1,data = pheno_GSE172120)
n.sv = num.sv(GSE172120_anno, mod, method="be")
GSE172120_anno = as.matrix(GSE172120_anno)

svobj = sva(GSE172120_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE172120$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 6)

X = svobj$sv
Y = GSE172120_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE172120_anno = GSE172120_anno-t(to_regress)

boxplot(GSE172120_anno)

save(file = "GSE172120.RData",GSE172120_anno,pheno_GSE172120)
