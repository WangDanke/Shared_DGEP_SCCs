#OSCC_1c_GSE138206, Affymetrix
library(WGCNA); library(lumi); library(ggplot2); library(nlme)
library(ggplot2); library(Cairo); library(GEOquery);library(nlme);
library(biomaRt); library(sva);library(limma);library(illuminaio)
setwd("D:\\SCC_DGEP_shared\\SCC_array\\HNSC\\GSE138206")
dir = getwd()
#
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)
basename(cel.files) #
GSE138206_raw<- ReadAffy(filenames = cel.files) #
sampleNames(GSE138206_raw) 
#the length of name maybe 8,9 or 10 
sampleNames(GSE138206_raw)<-stri_sub(sampleNames(GSE138206_raw),1,10)
##rma(Robust Multi-Array Average expression measure)：Background correcting，Normalizing，Calculating Expression
GSE138206_raw_rma <- rma(GSE138206_raw)
#
GSE138206 <- exprs(GSE138206_raw_rma) 
save(GSE138206,file = "GSE138206.RData")

###anno the gene symbol
GSE138206_anno <- as.data.frame(GSE138206)
GSE138206_anno$ID<-rownames(GSE138206_anno)
gpl<-getGEO("GPL570",destdir=".")
iddata<-as.data.frame(Table(gpl)[c('ID','Gene Symbol')])
GSE138206_anno<-merge(x=GSE138206_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE138206_anno <- GSE138206_anno[,-1]
###The expression levels of the same gene names were averaged
GSE138206_anno<- aggregate(GSE138206_anno,by = list(GSE138206_anno$`Gene Symbol`),FUN = mean)
head(GSE138206_anno)

GSE138206_anno <- GSE138206_anno[-c(1),]##blank gene name was dropped
rownames(GSE138206_anno) <- GSE138206_anno$Group.1
GSE138206_anno <- GSE138206_anno[,-c(1,20)]
GSE138206_anno[1:2,1:2]
boxplot(GSE138206_anno)
save(GSE138206_anno,file = "GSE138206_anno.RData")


##match the metadata with the expression data
matchSN = match(colnames(GSE138206_anno), rownames(pheno_GSE138206))
GSE138206_anno = GSE138206_anno[,matchSN]
save(file = "GSE138206.RData",GSE138206_anno,pheno_GSE138206)

##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE138206_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE138206$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE138206_anno)[outliers]); print(table(outliers))
GSE138206_anno = GSE138206_anno[,!outliers]
pheno_GSE138206 = pheno_GSE138206[!outliers,]


#clean the datMet
pheno_GSE138206$sample_type = factor(pheno_GSE138206$sample_type,levels = c("N","T"))
pheno_GSE138206$Group = as.numeric(ifelse(pheno_GSE138206$sample_type=="T","1","0"))

####Adjustment for covariates 
mod = model.matrix(~as.factor(Group),data = pheno_GSE138206)
mod0 = model.matrix(~1,data = pheno_GSE138206)
n.sv = num.sv(GSE138206_anno, mod, method="be")
GSE138206_anno = as.matrix(GSE138206_anno)

svobj = sva(GSE138206_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE138206$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 11)

X = svobj$sv
Y = GSE138206_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE138206_anno = GSE138206_anno-t(to_regress)
boxplot(GSE138206_anno)

save(file = "GSE138206.RData",GSE138206_anno,pheno_GSE138206)



















