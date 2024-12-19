###CSCC_1a_GSE98780_GSE45216 #Affymrtrix, datasets GSE98780 and GSE45216 are from the same Study,so combined together
library(stringi)
library(stringr)
library(FactoMineR)
library(WGCNA); library(lumi); library(ggplot2); library(nlme);library(affy)
library(ggplot2); library(Cairo); library(GEOquery);library(nlme);
library(biomaRt); library(sva);library(limma);library(illuminaio)

##preprocess the GSE98780
setwd("D:\\SCC_DGEP_shared\\SCC_array\\CSCC\\GSE98780_RAW")
dir = getwd()

##
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)

basename(cel.files) 
GSE98780_raw<- ReadAffy(filenames = cel.files) 
sampleNames(GSE98780_raw)
#the length of name maybe 8,9 or 10 
sampleNames(GSE98780_raw)<-stri_sub(sampleNames(GSE98780_raw),1,10)
##rma(Robust Multi-Array Average expression measure)：Background correcting，Normalizing，Calculating Expression
GSE98780_raw_rma <- rma(GSE98780_raw)
#
GSE98780 <- exprs(GSE98780_raw_rma) 
##
GSE98780 = GSE98780[,rownames(pheno_GSE98780_GPL570)]


##annotation
###anno the gene symbol
GSE98780_anno <- as.data.frame(GSE98780)
GSE98780_anno$ID<-rownames(GSE98780_anno)
gpl<-getGEO("GPL570",destdir=".")
iddata<-as.data.frame(Table(gpl)[c('ID','Gene Symbol')])
GSE98780_anno<-merge(x=GSE98780_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE98780_anno <- GSE98780_anno[,-1]
###The expression levels of the same gene names were averaged
GSE98780_anno<- aggregate(GSE98780_anno,by = list(GSE98780_anno$`Gene Symbol`),FUN = mean)
head(GSE98780_anno)

GSE98780_anno <- GSE98780_anno[-c(1),]##blank gene name was dropped
rownames(GSE98780_anno) <- GSE98780_anno$Group.1
GSE98780_anno <- GSE98780_anno[,-c(1,56)]
GSE98780_anno[1:2,1:2]
boxplot(GSE98780_anno)

##prepare the metadata
phe_GSE98780 = read.csv("GSE98780-GPL570_series_matrix.csv", header = F)
phe_GSE98780 = phe_GSE98780[c(29:38),]
phe_GSE98780 = phe_GSE98780 %>% t()
colnames(phe_GSE98780) = phe_GSE98780[1,]
phe_GSE98780 = phe_GSE98780[-1,]
rownames(phe_GSE98780) = phe_GSE98780[,2]
phe_GSE98780 = as.data.frame(phe_GSE98780)
colnames(phe_GSE98780)
pheno_GSE98780 = phe_GSE98780[,c(1,2,10)]
pheno_GSE98780$sample_type = ifelse(pheno_GSE98780$`!Sample_characteristics_ch1`=="tissue: Actinic Keratosis","actinic keratosis","N")
pheno_GSE98780 = pheno_GSE98780[,c(1,2,4)]

colnames(pheno_GSE98780) = c("Sample_title","ID","sample_type")
pheno_GSE98780$organ = "skin"

###match the metadata and the expression data
matchSN = match(colnames(GSE98780_anno), rownames(pheno_GSE98780_GPL570))
GSE98780_anno = GSE98780_anno[,matchSN]
save(file = "GSE98780_match.RData",GSE98780_anno,pheno_GSE98780_GPL570)


###preprocess GSE45216
######
setwd("D:\\SCC_DGEP_shared\\SCC_array\\CSCC\\GSE45216_RAW")
dir = getwd()
##
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)

basename(cel.files) 
GSE45216_raw<- ReadAffy(filenames = cel.files) 
sampleNames(GSE45216_raw) 
#the length of name maybe 8,9 or 10 
sampleNames(GSE45216_raw)<-stri_sub(sampleNames(GSE45216_raw),1,10)

#rma(Robust Multi-Array Average expression measure)：Background correcting，Normalizing，Calculating Expression
GSE45216_raw_rma <- rma(GSE45216_raw)
#
GSE45216 <- exprs(GSE45216_raw_rma) 

##prepare the metadata
phe_GSE45216 = read.csv("GSE45216_series_matrix.csv", header = F)
phe_GSE45216 = phe_GSE45216[c(35:48),]
phe_GSE45216 = phe_GSE45216 %>% t()
colnames(phe_GSE45216) = phe_GSE45216[1,]
phe_GSE45216 = phe_GSE45216[-1,]
rownames(phe_GSE45216) = phe_GSE45216[,2]
phe_GSE45216 = as.data.frame(phe_GSE45216)
colnames(phe_GSE45216)
pheno_GSE45216 = phe_GSE45216[,c(1,2,10,11,12,14)]
pheno_GSE45216$sample_type = ifelse(pheno_GSE45216$`!Sample_characteristics_ch1`=="tissue: Cutaneous SCC","T","actinic keratosis")
pheno_GSE45216$gender = ifelse(pheno_GSE45216$`!Sample_characteristics_ch1.3`=="gender: Male","M","F")
AA = as.data.frame(str_split_fixed(pheno_GSE45216$`!Sample_characteristics_ch1.1`,":",2))
BB = as.data.frame(str_split_fixed(pheno_GSE45216$`!Sample_characteristics_ch1.2`,":",2))
pheno_GSE45216 = cbind(pheno_GSE45216[,c(1,2,7,8)],AA[,2],BB[,2])
pheno_GSE45216$sample_type = c("actinic keratosis","actinic keratosis","actinic keratosis","actinic keratosis","N","N","N","N","N","N","T","T","T","T","T")
colnames(pheno_GSE45216) = c("Sample_title","ID","sample_type","gender","differentiated_state","immune_state")
pheno_GSE45216$organ = "skin"

##match the metadata and the expression data
GSE45216 = GSE45216[,rownames(pheno_GSE45216)]


##annotation
###anno the gene symbol
GSE45216_anno <- as.data.frame(GSE45216)
GSE45216_anno$ID<-rownames(GSE45216_anno)
gpl<-getGEO("GPL570",destdir=".")
iddata<-as.data.frame(Table(gpl)[c('ID','Gene Symbol')])
GSE45216_anno<-merge(x=GSE45216_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE45216_anno <- GSE45216_anno[,-1]
###The expression levels of the same gene names were averaged
GSE45216_anno<- aggregate(GSE45216_anno,by = list(GSE45216_anno$`Gene Symbol`),FUN = mean)
head(GSE45216_anno)

GSE45216_anno <- GSE45216_anno[-c(1),]##blank gene name was dropped
rownames(GSE45216_anno) <- GSE45216_anno$Group.1
GSE45216_anno <- GSE45216_anno[,-c(1,42)]
GSE45216_anno[1:2,1:2]
boxplot(GSE45216_anno)

matchSN = match(colnames(GSE45216_anno), rownames(pheno_GSE45216))
GSE45216_anno = GSE45216_anno[,matchSN]
save(file = "GSE45216_match.RData",GSE45216_anno,pheno_GSE45216)

#clean the datMet
pheno_GSE45216$sample_type = factor(pheno_GSE45216$sample_type,levels = c("N","T"))
pheno_GSE45216$Group = as.numeric(ifelse(pheno_GSE45216$sample_type=="T","1","0"))

pheno_GSE98780_GPL570$sample_type = factor(pheno_GSE98780_GPL570$sample_type,levels = c("N","T"))
pheno_GSE98780_GPL570$Group = as.numeric(ifelse(pheno_GSE98780_GPL570$sample_type=="T","1","0"))
#####
#####
##combined GSE45216 and GSE67610
gene = intersect(rownames(GSE45216_anno),rownames(GSE98780_anno))
sel = c("Sample_title","ID","sample_type","organ","batch")
pheno_GSE98780_com = rbind(pheno_GSE45216[,sel],pheno_GSE98780_GPL570[,sel])
GSE98780_com_anno = cbind(GSE98780_anno[gene,],GSE45216_anno[gene,])
boxplot(GSE98780_com_anno)


##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE98780_com_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE98780_com$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE98780_com_anno)[outliers]); print(table(outliers))
GSE98780_com_anno = GSE98780_com_anno[,!outliers]
pheno_GSE98780_com = pheno_GSE98780_com[!outliers,]

## Batch Correction
mod = model.matrix(~sample_type, data=pheno_GSE98780_com)
batch = factor(pheno_GSE98780_com$batch)
GSE98780_com_anno = ComBat(GSE98780_com_anno, batch=batch, mod=mod)

##remove actinic disease samples
pheno_GSE98780_com = pheno_GSE98780_com[!(pheno_GSE98780_com$sample_type=="actinic keratosis"),]
GSE98780_com_anno = GSE98780_com_anno[,rownames(pheno_GSE98780_com)]

##Regressing Unknown Covariates
pheno_GSE98780_com$sample_type = factor(pheno_GSE98780_com$sample_type,levels = c("N","T"))
pheno_GSE98780_com$Group = as.numeric(ifelse(pheno_GSE98780_com$sample_type=="T","1","0"))

mod = model.matrix(~as.factor(Group),data = pheno_GSE98780_com)
mod0 = model.matrix(~1,data = pheno_GSE98780_com)
n.sv = num.sv(GSE98780_com_anno, mod, method="be")
GSE98780_com_anno = as.matrix(GSE98780_com_anno)

svobj = sva(GSE98780_com_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE98780_com$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 61)

X = svobj$sv
Y = GSE98780_com_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE98780_com_anno = GSE98780_com_anno-t(to_regress)
boxplot(GSE98780_com_anno)
save(file = "GSE98780_com.RData",GSE98780_com_anno,pheno_GSE98780_com)









