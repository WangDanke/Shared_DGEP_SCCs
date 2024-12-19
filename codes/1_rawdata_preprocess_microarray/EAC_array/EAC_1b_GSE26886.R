setwd("D:\\SCC_DGEP_shared\\Adenocarcinoma\\EAC\\GSE26886_RAW")
library(tidyr)
library(dplyr)
library(stringr)
######prepare the metadata
phe_GSE26886 = read.csv("GSE26886_series_matrix.csv", header = F)
phe_GSE26886 = phe_GSE26886[c(29:36),]
phe_GSE26886 = phe_GSE26886 %>% t()
colnames(phe_GSE26886) = phe_GSE26886[1,]
phe_GSE26886 = phe_GSE26886[-1,]
rownames(phe_GSE26886) = phe_GSE26886[,2]
phe_GSE26886 = as.data.frame(phe_GSE26886)

colnames(phe_GSE26886)
pheno_GSE26886 = phe_GSE26886[,c(1,2,8)]
View(pheno_GSE26886)
save(pheno_GSE26886,file = "pheno_GSE26886.RData")

dir = "D:\\SCC_DGEP_shared\\Adenocarcinoma\\EAC\\GSE26886_RAW"
library(affy)
library(stringi)
##
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)

basename(cel.files)
GSE26886_raw<- ReadAffy(filenames = cel.files) 
sampleNames(GSE26886_raw) 
#the length of name maybe 8,9 or 10 
sampleNames(GSE26886_raw)<-stri_sub(sampleNames(GSE26886_raw),1,9)

##rma(Robust Multi-Array Average expression measure)：Background correcting，Normalizing，Calculating Expression
GSE26886_raw_rma <- rma(GSE26886_raw)
#
GSE26886 <- exprs(GSE26886_raw_rma) 
save(GSE26886,file = "GSE26886.RData")
##芯片注释
###anno the gene symbol
library(GEOquery)
GSE26886_anno <- as.data.frame(GSE26886)
GSE26886_anno$ID<-rownames(GSE26886_anno)
gpl<-getGEO("GPL570",destdir=".")
iddata<-as.data.frame(Table(gpl)[c('ID','Gene Symbol')])
GSE26886_anno<-merge(x=GSE26886_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE26886_anno <- GSE26886_anno[,-1]
###The expression levels of the same gene names were averaged
GSE26886_anno<- aggregate(GSE26886_anno,by = list(GSE26886_anno$`Gene Symbol`),FUN = mean)
head(GSE26886_anno)

GSE26886_anno <- GSE26886_anno[-c(1),]##blank gene name was dropped
rownames(GSE26886_anno) <- GSE26886_anno$Group.1
GSE26886_anno <- GSE26886_anno[,-c(1,71)]
GSE26886_anno[1:2,1:2]
boxplot(GSE26886_anno)

matchSN = match(colnames(GSE26886_anno), rownames(pheno_GSE26886))
GSE26886_anno = GSE26886_anno[,matchSN]
save(file = "GSE26886.RData",GSE26886_anno,pheno_GSE26886)

#clean the datMet
#Only esophageal adenocarcinoma tissue and normal tissue were selected
pheno_GSE26886= pheno_GSE26886[c(21:60),]
pheno_GSE26886$sample_type = ifelse(pheno_GSE26886$`!Sample_source_name_ch1`=="esophageal adenocarcinoma","T","N")
pheno_GSE26886 = pheno_GSE26886[,c(1,2,4)]
colnames(pheno_GSE26886) = c("Sample_title","ID","sample_type")
pheno_GSE26886$sample_type = factor(pheno_GSE26886$sample_type,levels = c("N","T"))
pheno_GSE26886$Group = as.numeric(ifelse(pheno_GSE26886$sample_type=="T","1","0"))

GSE26886_anno = GSE26886_anno[,rownames(pheno_GSE26886)]
matchSN = match(colnames(GSE26886_anno), rownames(pheno_GSE26886))
GSE26886_anno = GSE26886_anno[,matchSN]
save(file = "GSE26886_match.RData",GSE26886_anno,pheno_GSE26886)

##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE26886_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE26886$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomright", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE26886_anno)[outliers]); print(table(outliers))
GSE26886_anno = GSE26886_anno[,!outliers]
pheno_GSE26886 = pheno_GSE26886[!outliers,]

####Adjusting for unknown covariates
mod = model.matrix(~as.factor(Group),data = pheno_GSE26886)
mod0 = model.matrix(~1,data = pheno_GSE26886)
n.sv = num.sv(GSE26886_anno, mod, method="be")
GSE26886_anno = as.matrix(GSE26886_anno)

svobj = sva(GSE26886_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE26886$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 38)

X = svobj$sv
Y = GSE26886_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE26886_anno = GSE26886_anno-t(to_regress)
boxplot(GSE26886_anno)

save(file = "GSE26886.RData",GSE26886_anno,pheno_GSE26886)
















