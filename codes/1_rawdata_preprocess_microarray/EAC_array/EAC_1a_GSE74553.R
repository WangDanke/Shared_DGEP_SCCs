###EAC_1a_GSE74553 
library(stringi);library(pd.hugene.2.0.st);library(oligo)
library(stringi)
library(stringr)
library(FactoMineR)
library(WGCNA); library(lumi); library(ggplot2); library(nlme);library(affy)
library(ggplot2); library(Cairo); library(GEOquery);library(nlme);
library(biomaRt); library(sva);library(limma);library(illuminaio);

setwd("D:\\SCC_DGEP_shared\\Adenocarcinoma\\EAC\\GSE74553_RAW")
dir = getwd()

##
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)

basename(cel.files) #



GSE74553_raw<- read.celfiles( cel.files ) #
sampleNames(GSE74553_raw) #
#the length of name maybe 8,9 or 10 
sampleNames(GSE74553_raw)<-stri_sub(sampleNames(GSE74553_raw),1,10)
##rma(Robust Multi-Array Average expression measure)：Background correcting，Normalizing，Calculating Expression
GSE74553_raw_rma <- rma(GSE74553_raw)
GSE74553 <- exprs(GSE74553_raw_rma) 

##Only select esophageal adenocarcinoma tissue samples and normal tissue samples
pheno_GSE74553 = pheno_GSE74553[c(1:8,84:135),]
GSE74553 = GSE74553[,rownames(pheno_GSE74553)]

###anno the gene symbol
GSE74553_anno <- as.data.frame(GSE74553)
GSE74553_anno$ID<-rownames(GSE74553_anno)

###Load the chip annotation file downloaded from ThermorFisher
ANNO=fread(file = "HuGene-2_1-st-v1.na36.hg19.transcript.csv",sep = ",")
probe2gene <- ANNO[,c(2,8)]
library(stringr) 
probe2gene$symbol=trimws(str_split(probe2gene$gene_assignment,'//',simplify = T)[,2])
plot(table(table(probe2gene$symbol)),xlim=c(1,50))
head(probe2gene)

idname = probe2gene[,c(1,3)]
colnames(idname) = c("ID","Gene Symbol")

GSE74553_anno<-merge(x=GSE74553_anno,y=idname,by='ID',all.x=T,all.y=F)
GSE74553_anno <- GSE74553_anno[,-1]
###The expression levels of the same gene names were averaged
GSE74553_anno<- aggregate(GSE74553_anno,by = list(GSE74553_anno$`Gene Symbol`),FUN = mean)
head(GSE74553_anno)

GSE74553_anno <- GSE74553_anno[-c(1),]##blank gene name was dropped
rownames(GSE74553_anno) <- GSE74553_anno$Group.1
GSE74553_anno <- GSE74553_anno[,-c(1,62)]
GSE74553_anno[1:2,1:2]
boxplot(GSE74553_anno)

####prepare the metadata
phe_GSE74553 = read.csv("GSE74553_series_matrix.csv", header = F)
phe_GSE74553 = phe_GSE74553[c(31:45),]
phe_GSE74553 = phe_GSE74553 %>% t()
colnames(phe_GSE74553) = phe_GSE74553[1,]
phe_GSE74553 = phe_GSE74553[-1,]
rownames(phe_GSE74553) = phe_GSE74553[,2]
phe_GSE74553 = as.data.frame(phe_GSE74553)

colnames(phe_GSE74553)
pheno_GSE74553 = phe_GSE74553[,c(1,2,8)]
View(pheno_GSE74553)

#match the metadata and the expression data
matchSN = match(colnames(GSE74553_anno), rownames(pheno_GSE74553))
GSE74553_anno = GSE74553_anno[,matchSN]

#clean the datMet
pheno_GSE74553$sample_type = ifelse(pheno_GSE74553$`!Sample_source_name_ch1`=="normal esophageal squamous","N","T")
pheno_GSE74553 = pheno_GSE74553[,c(1,2,4)]
colnames(pheno_GSE74553) = c("Sample_title","ID","sample_type")
pheno_GSE74553$sample_type = factor(pheno_GSE74553$sample_type,levels = c("N","T"))
pheno_GSE74553$Group = as.numeric(ifelse(pheno_GSE74553$sample_type=="T","1","0"))


##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE74553_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE74553$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomright", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE74553_anno)[outliers]); print(table(outliers))
GSE74553_anno = GSE74553_anno[,!outliers]
pheno_GSE74553 = pheno_GSE74553[!outliers,]

####Adjusting for unknown covariates
mod = model.matrix(~as.factor(Group),data = pheno_GSE74553)
mod0 = model.matrix(~1,data = pheno_GSE74553)
n.sv = num.sv(GSE74553_anno, mod, method="be")
GSE74553_anno = as.matrix(GSE74553_anno)

svobj = sva(GSE74553_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE74553$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 58)

X = svobj$sv
Y = GSE74553_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE74553_anno = GSE74553_anno-t(to_regress)
boxplot(GSE74553_anno)

save(file = "GSE74553.RData",GSE74553_anno,pheno_GSE74553)
