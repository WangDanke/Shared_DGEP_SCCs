###CSCC_1b_GSE108008 #Affymetrix Human Gene 2.0 ST Array
library(stringi);library(pd.hugene.2.0.st);library(oligo)
 
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\CSCC\\GSE108008_RAW")
dir = getwd()
##
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)

basename(cel.files) 

GSE108008_raw<- read.celfiles( cel.files ) 
sampleNames(GSE108008_raw) 
#the length of name maybe 8,9 or 10 
sampleNames(GSE108008_raw)<-stri_sub(sampleNames(GSE108008_raw),1,10)

##rma(Robust Multi-Array Average expression measure)：Background correcting，Normalizing，Calculating Expression
GSE108008_raw_rma <- rma(GSE108008_raw)
#
GSE108008 <- exprs(GSE108008_raw_rma) 

#prepare the metadata
phe_GSE108008 = read.csv("GSE108008_series_matrix.csv", header = F)
phe_GSE108008 = phe_GSE108008[c(36:48),]
phe_GSE108008 = phe_GSE108008 %>% t()
colnames(phe_GSE108008) = phe_GSE108008[1,]
phe_GSE108008 = phe_GSE108008[-1,]
rownames(phe_GSE108008) = phe_GSE108008[,2]
phe_GSE108008 = as.data.frame(phe_GSE108008)
colnames(phe_GSE108008)
pheno_GSE108008 = phe_GSE108008[,c(1,2,8,12,13)]
pheno_GSE108008 = pheno_GSE108008[!(pheno_GSE108008$`!Sample_source_name_ch1`=="Human skin biopsy from actinic keratosis"),]###remove the actinic keratosis samples

pheno_GSE108008$age = substr(pheno_GSE108008$`!Sample_characteristics_ch1.1`,6,8)
pheno_GSE108008$gender = ifelse(pheno_GSE108008$`!Sample_characteristics_ch1`== "gender: male","M","F")
AA = as.data.frame(str_split_fixed(pheno_GSE108008$`!Sample_characteristics_ch1.3`,":",2))

pheno_GSE108008 = pheno_GSE108008[,c(1,2,6,7)]
pheno_GSE108008$sample_type = c("N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T")
pheno_GSE108008$organ = "skin"
colnames(pheno_GSE108008) = c("Sample_title","ID","age","gender","sample_type","organ")

save(pheno_GSE108008,file = "pheno_GSE108008.RData")

##match the expression data and the metadata
GSE108008 = GSE108008[,rownames(pheno_GSE108008)]

save(file = "GSE108008.RData",GSE108008,pheno_GSE108008)

##annotation
###anno the gene symbol
GSE108008_anno <- as.data.frame(GSE108008)
GSE108008_anno$ID<-rownames(GSE108008_anno)

###Load the chip annotation file downloaded from ThermorFisher
ANNO=fread(file = "HuGene-2_1-st-v1.na36.hg19.transcript.csv",sep = ",")
probe2gene <- ANNO[,c(2,8)]

library(stringr) 
probe2gene$symbol=trimws(str_split(probe2gene$gene_assignment,'//',simplify = T)[,2])
plot(table(table(probe2gene$symbol)),xlim=c(1,50))
head(probe2gene)

idname = probe2gene[,c(1,3)]
colnames(idname) = c("ID","Gene Symbol")

GSE108008_anno<-merge(x=GSE108008_anno,y=idname,by='ID',all.x=T,all.y=F)
GSE108008_anno <- GSE108008_anno[,-1]
###The expression levels of the same gene names were averaged
GSE108008_anno<- aggregate(GSE108008_anno,by = list(GSE108008_anno$`Gene Symbol`),FUN = mean)
head(GSE108008_anno)

GSE108008_anno <- GSE108008_anno[-c(1),]##blank gene name was dropped
rownames(GSE108008_anno) <- GSE108008_anno$Group.1
GSE108008_anno <- GSE108008_anno[,-c(1,22)]
GSE108008_anno[1:2,1:2]
boxplot(GSE108008_anno)

matchSN = match(colnames(GSE108008_anno), rownames(pheno_GSE108008))
GSE108008_anno = GSE108008_anno[,matchSN]
#save
save(file = "GSE108008.RData",GSE108008_anno,pheno_GSE108008)

#clean the datMet
pheno_GSE108008$sample_type = factor(pheno_GSE108008$sample_type,levels = c("N","T"))
pheno_GSE108008$Group = as.numeric(ifelse(pheno_GSE108008$sample_type=="T","1","0"))

##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE108008_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE108008$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomright", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE108008_anno)[outliers]); print(table(outliers))
GSE108008_anno = GSE108008_anno[,!outliers]
pheno_GSE108008 = pheno_GSE108008[!outliers,]

####Adjusting for unknown covariates
mod = model.matrix(~as.factor(Group),data = pheno_GSE108008)
mod0 = model.matrix(~1,data = pheno_GSE108008)
n.sv = num.sv(GSE108008_anno, mod, method="be")
GSE108008_anno = as.matrix(GSE108008_anno)

svobj = sva(GSE108008_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE108008$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 4)

X = svobj$sv
Y = GSE108008_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE108008_anno = GSE108008_anno-t(to_regress)
boxplot(GSE108008_anno)

save(file = "GSE108008.RData",GSE108008_anno,pheno_GSE108008)




