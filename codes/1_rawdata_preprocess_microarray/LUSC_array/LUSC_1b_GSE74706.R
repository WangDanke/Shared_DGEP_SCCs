####LUSC_1b_GSE74706
setwd("D:\\SCC_DGEP_shared\\SCC_array\\LUSC\\GSE74706_RAW")
dir = getwd()

#####prepare metadata
library(tidyr)
library(dplyr)
library(stringr)
###
###GSE74706
phe_GSE74706 = read.csv("GSE74706_series_matrix.csv", header = F)
phe_GSE74706 = phe_GSE74706[c(29:40),]
phe_GSE74706 = phe_GSE74706 %>% t()
colnames(phe_GSE74706) = phe_GSE74706[1,]
phe_GSE74706 = phe_GSE74706[-1,]
rownames(phe_GSE74706) = phe_GSE74706[,2]
phe_GSE74706 = as.data.frame(phe_GSE74706)

colnames(phe_GSE74706)
pheno_GSE74706 = phe_GSE74706[,c(1,2,5)]

pheno_GSE74706$sample_type = ifelse(pheno_GSE74706$X.Sample_characteristics_ch1.2=="tissue type: Squamous Cell Carcinoma","T","N")
rownames(pheno_GSE74706) = pheno_GSE74706$X.Sample_geo_accession
pheno_GSE74706 = pheno_GSE74706[,c(1,2,4)]
colnames(pheno_GSE74706)=c("Sample_title","ID","sample_type")
pheno_GSE74706$organ = "lung"
pheno_GSE74706$batch = "GSE74706"
save(pheno_GSE74706,file = "pheno_GSE74706.RData")

#list .TXT files
txt_files <- list.files(dir, pattern = "\\.txt.gz", full.names = TRUE)

GSE74706_raw <- read.maimages(txt_files,
                              source="agilent", 
                              green.only=TRUE,
                              other.columns = "gIsWellAboveBG")
dim(GSE74706_raw)

##background correction
GSE74706_1 <- backgroundCorrect(GSE74706_raw,method = "normexp") 
##normalization
GSE74706_1 <- normalizeBetweenArrays(GSE74706_1, method="quantile")
class(GSE74706_1)

Control <- GSE74706_1$genes$ControlType==1L
table(Control)

NoSymbol <- is.na(GSE74706_1$genes$ProbeName)
table(NoSymbol)


IsExpr <- rowSums(GSE74706_1$other$gIsWellAboveBG > 0) >= 8
table(IsExpr)

GSE74706_1_filt <- GSE74706_1[!Control & !NoSymbol & IsExpr, ]
dim(GSE74706_1_filt)


GSE74706 = GSE74706_1_filt@.Data[[1]]
boxplot(GSE74706) #


#
colnames(GSE74706) = str_extract(colnames(GSE74706),"GSM\\d*")
GSE74706[1:2,1:2]
#
rownames(GSE74706)= GSE74706_1_filt$genes$ProbeName

###anno the gene symbol
GSE74706_anno <- as.data.frame(GSE74706)
GSE74706_anno$ID<-rownames(GSE74706_anno)
gpl<-getGEO("GPL13497",destdir=".")
iddata<-as.data.frame(Table(gpl)[c('ID','GENE_SYMBOL')])
GSE74706_anno<-merge(x=GSE74706_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE74706_anno <- GSE74706_anno[,-1]
###The expression levels of the same gene names were averaged
GSE74706_anno<- aggregate(GSE74706_anno,by = list(GSE74706_anno$`GENE_SYMBOL`),FUN = mean)
head(GSE74706_anno)

GSE74706_anno <- GSE74706_anno[-c(1),]##blank gene name was dropped
rownames(GSE74706_anno) <- GSE74706_anno$Group.1
GSE74706_anno <- GSE74706_anno[,-c(1,38)]
GSE74706_anno[1:2,1:2]
boxplot(GSE74706_anno)

####match metadata and expression data
GSE74706_anno = GSE74706_anno[,rownames(pheno_GSE74706)]
matchSN = match(colnames(GSE74706_anno), rownames(pheno_GSE74706))
GSE74706_anno = GSE74706_anno[,matchSN]

save(file = "GSE74706_match.RData",GSE74706_anno,pheno_GSE74706)


##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE74706_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE74706$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE74706_anno)[outliers]); print(table(outliers))
GSE74706_anno = GSE74706_anno[,!outliers]
pheno_GSE74706 = pheno_GSE74706[!outliers,]


##Regressing Unknown Covariates
pheno_GSE74706$sample_type = factor(pheno_GSE74706$sample_type,levels = c("N","T"))
pheno_GSE74706$Group = as.numeric(ifelse(pheno_GSE74706$sample_type=="T","1","0"))

mod = model.matrix(~as.factor(Group),data = pheno_GSE74706)
mod0 = model.matrix(~1,data = pheno_GSE74706)
n.sv = num.sv(GSE74706_anno, mod, method="be")
GSE74706_anno = as.matrix(GSE74706_anno)

svobj = sva(GSE74706_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE74706$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 15)

X = svobj$sv
Y = GSE74706_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE74706_anno = GSE74706_anno-t(to_regress)
boxplot(GSE74706_anno)
save(file = "GSE74706.RData",GSE74706_anno,pheno_GSE74706)





