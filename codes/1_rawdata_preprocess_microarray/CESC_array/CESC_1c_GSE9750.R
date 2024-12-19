#cervical_SCC_1c_GSE9750, Affymetrix
library(stringi)
library(stringr)
library(FactoMineR)
library(WGCNA); library(lumi); library(ggplot2); library(nlme);library(affy)
library(ggplot2); library(Cairo); library(GEOquery);library(nlme);
library(biomaRt); library(sva);library(limma);library(illuminaio);

setwd("D:\\SCC_DGEP_shared\\SCCs_array\\CESC\\GSE9750_RAW")

dir = getwd()

##List files ending with .cel.gz
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)

basename(cel.files) #
GSE9750_raw<- ReadAffy(filenames = cel.files) #
sampleNames(GSE9750_raw) #
#the length of name maybe 8,9 or 10 
sampleNames(GSE9750_raw)<-stri_sub(sampleNames(GSE9750_raw),1,9)

##rma(Robust Multi-Array Average expression measure)：Background correcting，Normalizing，Calculating Expression
GSE9750_raw_rma <- rma(GSE9750_raw)
#
GSE9750 <- exprs(GSE9750_raw_rma)


###anno the gene symbol
GSE9750_anno <- as.data.frame(GSE9750)
GSE9750_anno$ID<-rownames(GSE9750_anno)
gpl<-getGEO("GPL96",destdir=".")
iddata<-as.data.frame(Table(gpl)[c('ID','Gene Symbol')])
GSE9750_anno<-merge(x=GSE9750_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE9750_anno <- GSE9750_anno[,-1]
###The expression levels of the same gene names were averaged
GSE9750_anno<- aggregate(GSE9750_anno,by = list(GSE9750_anno$`Gene Symbol`),FUN = mean)
head(GSE9750_anno)

GSE9750_anno <- GSE9750_anno[-c(1),]##blank gene name was dropped
rownames(GSE9750_anno) <- GSE9750_anno$Group.1
GSE9750_anno <- GSE9750_anno[,-c(1,68)]
GSE9750_anno[1:2,1:2]
boxplot(GSE9750_anno)

#prepare the metadata
phe_GSE9750 = read.csv("GSE9750_series_matrix.csv", header = F)
phe_GSE9750 = phe_GSE9750[c(35:44),]
phe_GSE9750 = phe_GSE9750 %>% t()
colnames(phe_GSE9750) = phe_GSE9750[1,]
phe_GSE9750 = phe_GSE9750[-1,]
rownames(phe_GSE9750) = phe_GSE9750[,2]
pheno_GSE9750 = as.data.frame(phe_GSE9750)
pheno_GSE9750 = pheno_GSE9750[-c(1:9),]   ##remove the cell samples
pheno_GSE9750 = pheno_GSE9750[,c(1,2,8,10)]
write.csv(pheno_GSE9750, file = "pheno_GSE9750.csv")
pheno_GSE9750 = read.csv("pheno_GSE9750.csv", header = F)
pheno_GSE9750 = pheno_GSE9750[,c(2,3,6,7)]
pheno_GSE9750 = pheno_GSE9750[-1,]
rownames(pheno_GSE9750) = pheno_GSE9750$V3
colnames(pheno_GSE9750) = c("Sample_title","ID","sample_type","age")
pheno_GSE9750$organ = "ovary"
save(pheno_GSE9750,file = "pheno_GSE9750.RData")

##match the metadata with the expression data
GSE9750_anno = GSE9750_anno[,rownames(pheno_GSE9750)]
matchSN = match(colnames(GSE9750_anno), rownames(pheno_GSE9750))
GSE9750_anno = GSE9750_anno[,matchSN]
save(file = "GSE9750.RData",GSE9750_anno,pheno_GSE9750)

##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE9750_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE9750$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE9750_anno)[outliers]); print(table(outliers))
GSE9750_anno = GSE9750_anno[,!outliers]
pheno_GSE9750 = pheno_GSE9750[!outliers,]
save(file = "GSE9750_remove.RData",pheno_GSE9750,GSE9750_anno)
#clean the datMet
pheno_GSE9750$sample_type = factor(pheno_GSE9750$sample_type,levels = c("N","T"))
pheno_GSE9750$Group = as.numeric(ifelse(pheno_GSE9750$sample_type=="T","1","0"))

###Regressing Unknown Covariates
mod = model.matrix(~as.factor(Group),data = pheno_GSE9750)
mod0 = model.matrix(~1,data = pheno_GSE9750)
n.sv = num.sv(GSE9750_anno, mod, method="be")
GSE9750_anno = as.matrix(GSE9750_anno)

svobj = sva(GSE9750_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE9750$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 53)

X = svobj$sv
Y = GSE9750_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE9750_anno = GSE9750_anno-t(to_regress)
boxplot(GSE9750_anno)
save(file = "GSE9750.RData",pheno_GSE9750,GSE9750_anno)


