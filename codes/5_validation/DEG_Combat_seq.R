
setwd("D:\\SCC_DGEP_shared\\SCC_seq\\Combat_DEG_seq")
load("D:\\SCC_DGEP_shared\\SCCs_array\\DEG_combat\\allDiff_5")
#######
###Differential expression analysis of SCC in different organ sites

##prepare the metadata and plot the QC plots, Fig. S4
genes = intersect(rownames(GSE223804_EXPR),rownames(allDiff_5))

CESC_EXPR_seq = cbind(GSE223804_EXPR[genes,],GSE87410_CESC_EXPR[genes,])
                      
pheno_CESC_seq = rbind(pheno_GSE223804[,c("Sample_title","ID","sample_type","organ","ID_SRR","batch")],
                       pheno_GSE87410_CESC[,c("Sample_title","ID","sample_type","organ","ID_SRR","batch")]
                       )
pheno_CESC_seq$organ = "CESC"
pheno_CESC_seq$sample_type = factor(pheno_CESC_seq$sample_type,levels = c("N","T"))
pheno_CESC_seq$Group = as.numeric(ifelse(pheno_CESC_seq$sample_type=="T","1","0"))
rownames(pheno_CESC_seq) = pheno_CESC_seq$ID_SRR

matchSN = match(colnames(CESC_EXPR_seq), rownames(pheno_CESC_seq))
CESC_EXPR_seq = CESC_EXPR_seq[,matchSN]
pheno_CESC_seq = pheno_CESC_seq[matchSN,]

save(pheno_CESC_seq,CESC_EXPR_seq,file = "CESC_EXPR_seq_TPM_log.RData")
###
library("FactoMineR")
library("factoextra")
ddb.pca <- PCA(t(CESC_EXPR_seq.combat), graph = F )
fviz_pca_ind(ddb.pca,
             geom.ind = "point",     # show points only (but not "text")
             col.ind = pheno_CESC_seq$sample_type, # color by groups
             addEllipses = TRUE,     # Concentration ellipses
             legend.title = "batch",
             range=0
)

library(limma)
##
mod = model.matrix(~sample_type, data=pheno_CESC_seq)
batch = factor(pheno_CESC_seq$batch)
CESC_EXPR_seq.combat = ComBat(CESC_EXPR_seq, batch=batch, mod=mod)

save(file = "CESC_EXPR_seq_combat_combined_2datasets.RData",CESC_EXPR_seq.combat,pheno_CESC_seq)

##After QC plot
boxplot(CESC_EXPR_seq.combat, range = 0, col= as.numeric(pheno_CESC_seq$sample_type), main ="Boxplot", ylab = "Intensity")
legend("topright", legend=c("N", "T"), col = 1:2, pch=19)

# Histogram
i = 1; plot(density((CESC_EXPR_seq.combat[,i]), na.rm=T), col = i, main="Hist of Log2 Exp", xlab = "log2 exp")
for(i in 2:dim(CESC_EXPR_seq.combat)[2])
  lines(density((CESC_EXPR_seq.combat[,i]), na.rm=T), col =as.numeric(pheno_CESC_seq$sample_type))
legend("topright", levels(pheno_CESC_seq$sample_type), cex=0.7, text.col = 1:2)



##Differential gene expression analysis.
data = HNSC_EXPR_seq.combat
group <- as.character(pheno_HNSC_seq$sample_type)
group <- factor(group,levels = c("N","T"),ordered = F)
design <- model.matrix(~group)

fit <- lmFit(data,design)
fit2 <- eBayes(fit)
allDiff_HNSC_seq=topTable(fit2,adjust='fdr',coef=2,number=Inf)

save(allDiff_HNSC_seq,file = "allDiff_HNSC_seq.Rdata",quote = F, row.names = T)


###allDiff_5_seq,combined all 8816 genes's logFC in validation datasets



























