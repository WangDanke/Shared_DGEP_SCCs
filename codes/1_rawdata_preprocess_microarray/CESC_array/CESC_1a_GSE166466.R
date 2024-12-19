###CESC_1a_GSE166466 #[Clariom_D_Human] Affymetrix Human Clariom D Assay [transcript (gene) version]
library(stringi);library(pd.hugene.2.0.st);library(oligo)

setwd("D:\\SCC_DGEP_shared\\SCCs_array\\CESC\\GSE166466_RAW")
dir = getwd()

##List files ending with .cel.gz
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)

basename(cel.files) 
GSE166466_raw<- read.celfiles( cel.files ) 
sampleNames(GSE166466_raw) 
#the length of name maybe 8,9 or 10 
sampleNames(GSE166466_raw)<-stri_sub(sampleNames(GSE166466_raw),1,10)

##preprocess the rawdata, rma(Robust Multi-Array Average expression measure)：Background correcting，Normalizing，Calculating Expression
GSE166466_raw_rma <- rma(GSE166466_raw)
##
GSE166466 <- exprs(GSE166466_raw_rma) 

##prepare the metadata
phe_GSE166466 = read.csv("GSE166466_series_matrix.csv", header = F)
phe_GSE166466 = phe_GSE166466[c(26:33),]
phe_GSE166466 = phe_GSE166466 %>% t()
colnames(phe_GSE166466) = phe_GSE166466[1,]
phe_GSE166466 = phe_GSE166466[-1,]
rownames(phe_GSE166466) = phe_GSE166466[,2]
pheno_GSE166466 = phe_GSE166466[,c(1,2)]
pheno_GSE166466 = pheno_GSE166466[-c(8:13),] ##Removal of actinic disease patient samples
pheno_GSE166466 = as.data.frame(pheno_GSE166466)

colnames(pheno_GSE166466) = c("Sample_title","ID")
pheno_GSE166466$sample_type = c("N","N","N","N","N","N","N","T","T","T","T","T","T","T")
pheno_GSE166466$organ = "ovary"
save(pheno_GSE166466,file = "pheno_GSE166466.RData")

##match the metadata and the expression data
GSE166466 = GSE166466[,rownames(pheno_GSE166466)]

save(file = "GSE166466.RData",GSE166466,pheno_GSE166466)

##annotation
###anno the gene symbol
GSE166466_anno <- as.data.frame(GSE166466)
GSE166466_anno$ID<-rownames(GSE166466_anno)

###Load the chip annotation file downloaded from ThermorFisher
library(data.table)
ANNO=fread(file = "Clariom_D_Human.na36.hg38.probeset.csv",sep = ",")

probe2gene <- ANNO[,c(7,11)]

library(stringr) 
probe2gene$symbol=trimws(str_split(probe2gene$gene_assignment,'//',simplify = T)[,2])
plot(table(table(probe2gene$symbol)),xlim=c(1,50))
head(probe2gene)

idname = probe2gene[,c(1,3)]
colnames(idname) = c("ID","Gene Symbol")

GSE166466_anno<-merge(x=GSE166466_anno,y=idname,by='ID',all.x=T,all.y=F)
GSE166466_anno <- GSE166466_anno[,-1]
###The expression levels of the same gene names were averaged
GSE166466_anno<- aggregate(GSE166466_anno[,-15],by = list(GSE166466_anno$`Gene Symbol`),FUN = mean)
head(GSE166466_anno)

GSE166466_anno <- GSE166466_anno[-c(1),]##blank gene name was dropped
rownames(GSE166466_anno) <- GSE166466_anno$Group.1
GSE166466_anno <- GSE166466_anno[,-c(1)]
GSE166466_anno[1:2,1:2]
boxplot(GSE166466_anno)

matchSN = match(colnames(GSE166466_anno), rownames(pheno_GSE166466))
GSE166466_anno = GSE166466_anno[,matchSN]
#
save(file = "GSE166466.RData",GSE166466_anno,pheno_GSE166466)


##remove outliers
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE166466_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE166466$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE166466_anno)[outliers]); print(table(outliers))
GSE166466_anno = GSE166466_anno[,!outliers]
pheno_GSE166466 = pheno_GSE166466[!outliers,]


#clean the datMet
pheno_GSE166466$sample_type = factor(pheno_GSE166466$sample_type,levels = c("N","T"))
pheno_GSE166466$Group = as.numeric(ifelse(pheno_GSE166466$sample_type=="T","1","0"))

###Regressing Unknown Covariates
mod = model.matrix(~as.factor(Group),data = pheno_GSE166466)
mod0 = model.matrix(~1,data = pheno_GSE166466)
n.sv = num.sv(GSE166466_anno, mod, method="be")
GSE166466_anno = as.matrix(GSE166466_anno)

svobj = sva(GSE166466_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE166466$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 13)

X = svobj$sv
Y = GSE166466_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE166466_anno = GSE166466_anno-t(to_regress)
boxplot(GSE166466_anno)
save(file = "GSE166466.RData",pheno_GSE166466,GSE166466_anno)
