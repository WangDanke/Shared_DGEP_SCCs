###batch effect as a covariate in the mixed-effect model.
setwd("D:\\SCC_DGEP_shared\\SCCs_array\\DEG_Meta")
load("D:\\SCC_DGEP_shared\\SCCs_array\\After_preprocess\\CESC\\CESC_EXPR_3_datasets.RData")
load("D:\\SCC_DGEP_shared\\SCCs_array\\After_preprocess\\CSCC\\CSCC_EXPR_2_datasets.RData")
load("D:\\SCC_DGEP_shared\\SCCs_array\\After_preprocess\\HNSC\\HNSC_EXPR_5_datasets.RData")
load("D:\\SCC_DGEP_shared\\SCCs_array\\After_preprocess\\ESCC\\ESCC_EXPR_7_datasets.RData")
load("D:\\SCC_DGEP_shared\\SCCs_array\\After_preprocess\\LUSC\\LUSC_EXPR_4_datasets.RData")
load("D:\\SCC_DGEP_shared\\SCCs_array\\After_preprocess\\EAC\\EAC_EXPR_2_datasets.RData")
load("D:\\SCC_DGEP_shared\\SCCs_array\\After_preprocess\\LUAD\\LUAD_EXPR_4_datasets.RData")

library(nlme)
####CESC
pheno_CESC$sample_type = factor(pheno_CESC$sample_type, levels=c("N", "T"))
pheno_CESC$batch = as.factor(pheno_CESC$batch)

bd_meta = matrix(NA, nrow=nrow(CESC_EXPR), ncol=3)
for(i in 1:nrow(CESC_EXPR)) {
  if(i%%100==0) print(i)
  expr = as.numeric(CESC_EXPR[i,])
  tryCatch({
    bd_meta[i,] = summary(lme(expr~ sample_type + batch, data = pheno_CESC, random = ~1|ID))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

bd_meta=as.data.frame(bd_meta)
colnames(bd_meta) = c("beta", "SE", "p")
rownames(bd_meta) = rownames(CESC_EXPR)
bd_meta$fdr = p.adjust(bd_meta$p, "fdr")
CESC_bd_Meta = bd_meta
save(CESC_bd_Meta,file = "CESC_bd_Meta.RData")


####LUSC
pheno_LUSC$sample_type = factor(pheno_LUSC$sample_type, levels=c("N", "T"))
pheno_LUSC$batch = as.factor(pheno_LUSC$batch)

bd_meta = matrix(NA, nrow=nrow(LUSC_EXPR), ncol=3)
for(i in 1:nrow(LUSC_EXPR)) {
  if(i%%100==0) print(i)
  expr = as.numeric(LUSC_EXPR[i,])
  tryCatch({
    bd_meta[i,] = summary(lme(expr~ sample_type + batch, data = pheno_LUSC, random = ~1|ID))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

bd_meta=as.data.frame(bd_meta)
colnames(bd_meta) = c("beta", "SE", "p")
rownames(bd_meta) = rownames(LUSC_EXPR)
bd_meta$fdr = p.adjust(bd_meta$p, "fdr")
LUSC_bd_Meta = bd_meta
save(LUSC_bd_Meta,file = "LUSC_bd_Meta.RData")

##LUAD
####
pheno_LUAD$sample_type = factor(pheno_LUAD$sample_type, levels=c("N", "T"))
pheno_LUAD$batch = as.factor(pheno_LUAD$batch)

bd_meta = matrix(NA, nrow=nrow(LUAD_EXPR), ncol=3)
for(i in 1:nrow(LUAD_EXPR)) {
  if(i%%100==0) print(i)
  expr = as.numeric(LUAD_EXPR[i,])
  tryCatch({
    bd_meta[i,] = summary(lme(expr~ sample_type + batch, data = pheno_LUAD, random = ~1|ID))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

bd_meta=as.data.frame(bd_meta)
colnames(bd_meta) = c("beta", "SE", "p")
rownames(bd_meta) = rownames(LUAD_EXPR)
bd_meta$fdr = p.adjust(bd_meta$p, "fdr")
LUAD_bd_Meta = bd_meta
save(LUAD_bd_Meta,file = "LUAD_bd_Meta.RData")

####HNSC
pheno_HNSC$sample_type = factor(pheno_HNSC$sample_type, levels=c("N", "T"))
pheno_HNSC$batch = as.factor(pheno_HNSC$batch)

bd_meta = matrix(NA, nrow=nrow(HNSC_EXPR), ncol=3)
for(i in 1:nrow(HNSC_EXPR)) {
  if(i%%100==0) print(i)
  expr = as.numeric(HNSC_EXPR[i,])
  tryCatch({
    bd_meta[i,] = summary(lme(expr~ sample_type + batch, data = pheno_HNSC, random = ~1|ID))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

bd_meta=as.data.frame(bd_meta)
colnames(bd_meta) = c("beta", "SE", "p")
rownames(bd_meta) = rownames(HNSC_EXPR)
bd_meta$fdr = p.adjust(bd_meta$p, "fdr")
HNSC_bd_Meta = bd_meta
save(HNSC_bd_Meta,file = "HNSC_bd_Meta.RData")


####ESCC
pheno_ESCC$sample_type = factor(pheno_ESCC$sample_type, levels=c("N", "T"))
pheno_ESCC$batch = as.factor(pheno_ESCC$batch)

bd_meta = matrix(NA, nrow=nrow(ESCC_EXPR), ncol=3)
for(i in 1:nrow(ESCC_EXPR)) {
  if(i%%100==0) print(i)
  expr = as.numeric(ESCC_EXPR[i,])
  tryCatch({
    bd_meta[i,] = summary(lme(expr~ sample_type + batch, data = pheno_ESCC, random = ~1|ID))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

bd_meta=as.data.frame(bd_meta)
colnames(bd_meta) = c("beta", "SE", "p")
rownames(bd_meta) = rownames(ESCC_EXPR)
bd_meta$fdr = p.adjust(bd_meta$p, "fdr")
ESCC_bd_Meta = bd_meta
save(ESCC_bd_Meta,file = "ESCC_bd_Meta.RData")

####EAC
pheno_EAC$sample_type = factor(pheno_EAC$sample_type, levels=c("N", "T"))
pheno_EAC$batch = as.factor(pheno_EAC$batch)

bd_meta = matrix(NA, nrow=nrow(EAC_EXPR), ncol=3)
for(i in 1:nrow(EAC_EXPR)) {
  if(i%%100==0) print(i)
  expr = as.numeric(EAC_EXPR[i,])
  tryCatch({
    bd_meta[i,] = summary(lme(expr~ sample_type + batch, data = pheno_EAC, random = ~1|ID))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

bd_meta=as.data.frame(bd_meta)
colnames(bd_meta) = c("beta", "SE", "p")
rownames(bd_meta) = rownames(EAC_EXPR)
bd_meta$fdr = p.adjust(bd_meta$p, "fdr")
EAC_bd_Meta = bd_meta
save(EAC_bd_Meta,file = "EAC_bd_Meta.RData")

####CSCC
pheno_CSCC$sample_type = factor(pheno_CSCC$sample_type, levels=c("N", "T"))
pheno_CSCC$batch = as.factor(pheno_CSCC$batch)

bd_meta = matrix(NA, nrow=nrow(CSCC_EXPR), ncol=3)
for(i in 1:nrow(CSCC_EXPR)) {
  if(i%%100==0) print(i)
  expr = as.numeric(CSCC_EXPR[i,])
  tryCatch({
    bd_meta[i,] = summary(lme(expr~ sample_type + batch, data = pheno_CSCC, random = ~1|ID))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

bd_meta=as.data.frame(bd_meta)
colnames(bd_meta) = c("beta", "SE", "p")
rownames(bd_meta) = rownames(CSCC_EXPR)
bd_meta$fdr = p.adjust(bd_meta$p, "fdr")
CSCC_bd_Meta = bd_meta
save(CSCC_bd_Meta,file = "CSCC_bd_Meta.RData")


##Save compiled expression and metadata for permutation testing
multiExpr = vector(mode="list", length=5)
multiExpr[[1]]$datExpr = LUSC_EXPR;  multiExpr[[1]]$datMeta = pheno_LUSC
multiExpr[[2]]$datExpr = CSCC_EXPR; multiExpr[[2]]$datMeta = pheno_CSCC 
multiExpr[[3]]$datExpr = ESCC_EXPR; multiExpr[[3]]$datMeta = pheno_ESCC 
multiExpr[[4]]$datExpr = HNSC_EXPR; multiExpr[[4]]$datMeta = pheno_HNSC 
multiExpr[[5]]$datExpr = CESC_EXPR; multiExpr[[5]]$datMeta = pheno_CESC 
multiExpr[[6]]$datExpr = EAC_EXPR; multiExpr[[6]]$datMeta = pheno_EAC 
multiExpr[[7]]$datExpr = LUAD_EXPR; multiExpr[[7]]$datMeta = pheno_EAC
names(multiExpr) = c("LUSC", "CSCC", "ESCC", "HNSC", "CESC")

genes = intersect(rownames(multiExpr[[1]]$datExpr), rownames(multiExpr[[2]]$datExpr))
for(i in 3:7) genes = intersect(genes, rownames(multiExpr[[i]]$datExpr))

for(i in 1:7) multiExpr[[i]]$datExpr = multiExpr[[i]]$datExpr[match(genes, rownames(multiExpr[[i]]$datExpr)),]
save(file="Microarray_compiledForPermutationTesting.RData", multiExpr)

####Combined the DEG_meta results
CESC_bd_Meta = na.omit(CESC_bd_Meta)
CSCC_bd_Meta = na.omit(CSCC_bd_Meta)
EAC_bd_Meta = na.omit(EAC_bd_Meta)
ESCC_bd_Meta = na.omit(ESCC_bd_Meta)
HNSC_bd_Meta = na.omit(HNSC_bd_Meta)
LUSC_bd_Meta = na.omit(LUSC_bd_Meta)
LUAD_bd_Meta = na.omit(LUAD_bd_Meta)
AA = intersect(rownames(CESC_bd_Meta),rownames(CSCC_bd_Meta))
BB = intersect(rownames(EAC_bd_Meta),rownames(ESCC_bd_Meta))
CC = intersect(rownames(LUSC_bd_Meta),rownames(LUAD_bd_Meta))
DD = intersect(AA,BB)
EE = intersect(CC,rownames(HNSC_bd_Meta))

FF = intersect(DD,EE)

DEG_beta = cbind(LUAD_bd_Meta[FF,]$beta,LUSC_bd_Meta[FF,]$beta,CSCC_bd_Meta[FF,]$beta,ESCC_bd_Meta[FF,]$beta,EAC_bd_Meta[FF,]$beta,HNSC_bd_Meta[FF,]$beta,CESC_bd_Meta[FF,]$beta)
DEG_beta = as.data.frame(DEG_beta)
rownames(DEG_beta)= FF
colnames(DEG_beta) = c("LUAD","LUSC","CSCC","ESCC","EAC","HNSC","CESC")

save(DEG_beta, file =" DEG_beta_na.RData")

#############correlation anlysis
library(corrplot)
tdc = cor(DEG_beta,method = "spearman")
testRes = cor.mtest(DEG_beta, method="spearman",conf.level = 0.95)


















