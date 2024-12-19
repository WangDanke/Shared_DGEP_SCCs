##Use the Combat function to remove batch effect
library(stringi)
library(stringr)
library(FactoMineR)
library(WGCNA); library(lumi); library(ggplot2); library(nlme);library(affy)
library(ggplot2); library(Cairo); library(GEOquery);library(nlme);
library(biomaRt); library(sva);library(limma);library(illuminaio);library(oligo)

###Only the overlapped genes across all types of SCC are included in the further analysis.
setwd("D:\\SCC_DGEP_shared\\SCCs_array\\DEG_combat")
load("D:\\SCC_DGEP_shared\\SCCs_array\\After_preprocess\\CESC\\CESC_EXPR.RData")
load("D:\\SCC_DGEP_shared\\SCCs_array\\After_preprocess\\HNSC\\HNSC_EXPR.RData")
load("D:\\SCC_DGEP_shared\\SCCs_array\\After_preprocess\\CSCC\\CSCC_EXPR.RData")
load("D:\\SCC_DGEP_shared\\SCCs_array\\After_preprocess\\ESCC\\ESCC_EXPR.RData")
load("D:\\SCC_DGEP_shared\\SCCs_array\\After_preprocess\\LUSC\\LUSC_EXPR.RData")

AA = intersect(rownames(CSCC_EXPR),rownames(CESC_EXPR))
BB = intersect(rownames(HNSC_EXPR),rownames(ESCC_EXPR))
CC = intersect(AA,rownames(LUSC_EXPR))
DD = intersect(CC,BB)

CSCC_EXPR = CSCC_EXPR[DD,]
CESC_EXPR = CESC_EXPR[DD,]
ESCC_EXPR = ESCC_EXPR[DD,]
HNSC_EXPR = HNSC_EXPR[DD,]
LUSC_EXPR = LUSC_EXPR[DD,]

pheno_CESC$Group = as.numeric(ifelse(pheno_CESC$sample_type=="T","1","0"))
pheno_CSCC$Group = as.numeric(ifelse(pheno_CSCC$sample_type=="T","1","0"))
pheno_HNSC$Group = as.numeric(ifelse(pheno_HNSC$sample_type=="T","1","0"))
pheno_LUSC$Group = as.numeric(ifelse(pheno_LUSC$sample_type=="T","1","0"))

sel = c("Sample_title","ID","sample_type","organ","batch","Group")

pheno_CESC = pheno_CESC[,sel]
pheno_CSCC = pheno_CSCC[,sel]
pheno_ESCC = pheno_ESCC[,sel]
pheno_HNSC = pheno_HNSC[,sel]
pheno_LUSC = pheno_LUSC[,sel]

save(pheno_CESC,pheno_CSCC,pheno_ESCC,pheno_HNSC,pheno_LUSC,
     CESC_EXPR,CSCC_EXPR,ESCC_EXPR,HNSC_EXPR,LUSC_EXPR,file = "SCC_combat_combined.RData")

#####LUSC
###use the limma package to find the DEGs
data = LUSC_EXPR
group <- as.character(pheno_LUSC$sample_type)
group <- factor(group,levels = c("N","T"),ordered = F)
design <- model.matrix(~group)

fit <- lmFit(data,design)
fit2 <- eBayes(fit)
allDiff_LUSC=topTable(fit2,adjust='fdr',coef=2,number=Inf)

save(allDiff_LUSC,file = "allDiff_LUSC.Rdata",quote = F, row.names = T)

#####LUAD
data = LUAD_EXPR[DD,]
group <- as.character(pheno_LUAD$sample_type)
group <- factor(group,levels = c("N","T"),ordered = F)
design <- model.matrix(~group)

fit <- lmFit(data,design)
fit2 <- eBayes(fit)
allDiff_LUAD=topTable(fit2,adjust='fdr',coef=2,number=Inf)

save(allDiff_LUAD,file = "allDiff_LUAD.Rdata",quote = F, row.names = T)

#####ESCC
###use the limma package to find the DEGs
#### Differential expression analysis of individual subtype and paired normal tissues
data = ESCC_EXPR
group <- as.character(pheno_ESCC$sample_type)
group <- factor(group,levels = c("N","T"),ordered = F)
design <- model.matrix(~group)

fit <- lmFit(data,design)
fit2 <- eBayes(fit)
allDiff_ESCC=topTable(fit2,adjust='fdr',coef=2,number=Inf)

save(allDiff_ESCC,file = "allDiff_ESCC.Rdata",quote = F, row.names = T)

#####cSCC
###use the limma package to find the DEGs
data = CSCC_EXPR
group <- as.character(pheno_CSCC$sample_type)
group <- factor(group,levels = c("N","T"),ordered = F)
design <- model.matrix(~group)

fit <- lmFit(data,design)
fit2 <- eBayes(fit)
allDiff_CSCC=topTable(fit2,adjust='fdr',coef=2,number=Inf)
save(allDiff_CSCC,file = "allDiff_CSCC.Rdata",quote = F, row.names = T)

#####HNSC
data = HNSC_EXPR
group <- as.character(pheno_HNSC$sample_type)
group <- factor(group,levels = c("N","T"),ordered = F)
design <- model.matrix(~group)

fit <- lmFit(data,design)
fit2 <- eBayes(fit)
allDiff_HNSC=topTable(fit2,adjust='fdr',coef=2,number=Inf)

save(allDiff_HNSC,file = "allDiff_HNSC.Rdata",quote = F, row.names = T)


#####CESC
data = CESC_EXPR
group <- as.character(pheno_CESC$sample_type)
group <- factor(group,levels = c("N","T"),ordered = F)
design <- model.matrix(~group)

fit <- lmFit(data,design)
fit2 <- eBayes(fit)
allDiff_CESC=topTable(fit2,adjust='fdr',coef=2,number=Inf)

save(allDiff_CESC,file = "allDiff_CESC.Rdata",quote = F, row.names = T)

#####EAC
###use the limma package to find the DEGs
data = EAC_EXPR[DD,]
group <- as.character(pheno_EAC$sample_type)
group <- factor(group,levels = c("N","T"),ordered = F)
design <- model.matrix(~group)
fit <- lmFit(data,design)
fit2 <- eBayes(fit)
allDiff_EAC=topTable(fit2,adjust='fdr',coef=2,number=Inf)
save(allDiff_EAC,file = "allDiff_EAC.Rdata",quote = F, row.names = T)

#####LUAD
data = LUAD_EXPR[DD,]
group <- as.character(pheno_LUAD$sample_type)
group <- factor(group,levels = c("N","T"),ordered = F)
design <- model.matrix(~group)
fit <- lmFit(data,design)
fit2 <- eBayes(fit)
allDiff_LUAD=topTable(fit2,adjust='fdr',coef=2,number=Inf)
save(allDiff_LUAD,file = "allDiff_LUAD.Rdata",quote = F, row.names = T)

###########the overlapped genes in all types of carcinoma are included in the correlation analysis######
AA = intersect(intersect(rownames(allDiff_LUSC),rownames(allDiff_EAC)),rownames(allDiff_LUAD))
allDiff_CSCC = allDiff_CSCC[AA,]
allDiff_ESCC = allDiff_ESCC[AA,]
allDiff_EAC = allDiff_EAC[AA,]
allDiff_HNSC = allDiff_HNSC[AA,]
allDiff_CESC = allDiff_CESC[AA,]
allDiff_LUAD = allDiff_LUAD[AA,]
allDiff_LUSC = allDiff_LUSC[AA,]

allDiff = cbind(allDiff_LUAD$logFC,allDiff_LUSC$logFC,allDiff_CSCC$logFC,allDiff_ESCC$logFC,allDiff_EAC$logFC,
                allDiff_HNSC$logFC,allDiff_CESC$logFC)
colnames(allDiff) = c("LUAD","LUSC","CSCC","ESCC","EAC","HNSC","CESC")
rownames(allDiff) = AA
allDiff = as.data.frame(allDiff)
save(allDiff,file = "allDiff_7.RData")


library(corrplot)
library(ggplot2)
library(ggpubr)
class(allDiff)
###correlation analysis
tdc = cor(allDiff,method = "spearman")
testRes = cor.mtest(allDiff, method="spearman",conf.level = 0.95)

corrplot(tdc)
addcol <- colorRampPalette(c("#67BE54", "#FFC125", "#B22222"))
corrplot(tdc, method = "color", col = addcol(100), 
         tl.col = "black", tl.cex = 0.8, tl.srt = 45,tl.pos = "lt",
         p.mat = testRes$p, diag = T, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.2,
         insig = 'label_sig', pch.col = 'grey20', order = 'AOE')


############ Fig. 2B Scatter plot
##
dat = as.data.frame(allDiff_5)
library(data.table)
dat2= reshape2::melt(dat,id=3)
dat2$value = as.numeric(dat2$value)

pcreg = function(ds1, ds2) {
  #Principle components regression to calculate slope 
  r = prcomp(~ds1+ds2)
  slope <- r$rotation[2,1] / r$rotation[1,1]
  intercept <- r$center[2] - slope*r$center[1]
  return(list(slope,intercept))
}

fit_2 = pcreg(dat$ESCC, dat$CSCC)
fit_3 = pcreg(dat$ESCC, dat$LUSC)
fit_4 = pcreg(dat$ESCC, dat$HNSC)
fit_5 = pcreg(dat$ESCC, dat$CESC)

dat2$variable = as.character(dat2$variable)
dat2$variable = gsub("CSCC", paste("CSCC, slope=", signif(fit_2[[1]],2), sep=""), dat2$variable)
dat2$variable = gsub("LUSC", paste("LUSC, slope=", signif(fit_3[[1]],2), sep=""), dat2$variable)
dat2$variable = gsub("HNSC", paste("HNSC, slope=", signif(fit_4[[1]],2), sep=""), dat2$variable)
dat2$variable = gsub("CESC", paste("CESC, slope=", signif(fit_5[[1]],2), sep=""), dat2$variable)

library(ggplot2)
point_colors <- c("CSCC"="#E1C855","LUSC"="#E07B54", "HNSC"="#51B1B7","CESC"='#A5C496')
SCCs_DEGs=ggplot(dat2,aes(x=ESCC,y=value,color=variable)) + 
  scale_color_manual(values = c('#A5C496',"#E1C855","#51B1B7","#E07B54"))+ 
  geom_point(alpha=.7, size = 1.2) + 
  geom_abline(slope=1, lty=2) + xlim(-4,4) + ylim(-4,4) + 
  geom_abline(slope=fit_2[[1]], intercept = fit_2[[2]], color=point_colors["CSCC"],linewidth = 0.7) +
  geom_abline(slope=fit_3[[1]], intercept = fit_3[[2]], color=point_colors["LUSC"],linewidth = 0.7) + 
  geom_abline(slope=fit_4[[1]], intercept = fit_4[[2]], color=point_colors["HNSC"],linewidth = 0.7) +
  geom_abline(slope=fit_5[[1]], intercept = fit_5[[2]], color=point_colors["CESC"],linewidth = 0.7) +
  xlab("ESCC (log2FC)") + ylab("SCCs(log2FC)") +
  coord_fixed(ratio=1)+
  theme(panel.grid=element_blank(), 
        panel.background=element_rect(color="black", fill="transparent"),
        panel.border = element_rect(color = "#363636", fill = NA, size = 1))##theme

SCCs_DEGs


################defined the DEGs，|LogFC| > 1, P < 0.05
allDiff_LUSC$result = ifelse(abs(allDiff_LUSC$logFC)>1&allDiff_LUSC$adj.P.Val < 0.05,"1","0")
allDiff_HNSC$result = ifelse(abs(allDiff_HNSC$logFC)>1&allDiff_HNSC$adj.P.Val < 0.05,"1","0")
allDiff_CESC$result = ifelse(abs(allDiff_CESC$logFC)>1&allDiff_CESC$adj.P.Val < 0.05,"1","0")
allDiff_CSCC$result = ifelse(abs(allDiff_CSCC$logFC)>1&allDiff_CSCC$adj.P.Val < 0.05,"1","0")
allDiff_ESCC$result = ifelse(abs(allDiff_ESCC$logFC)>1&allDiff_ESCC$adj.P.Val < 0.05,"1","0")

#upset plot
library(UpSetR)
allDiff_CESC = allDiff_CESC[rownames(allDiff_ESCC),]
allDiff_CSCC = allDiff_CSCC[rownames(allDiff_ESCC),]
allDiff_HNSC = allDiff_HNSC[rownames(allDiff_ESCC),]
allDiff_LUSC = allDiff_LUSC[rownames(allDiff_ESCC),]

all= cbind(allDiff_CESC$result,allDiff_CSCC$result,allDiff_ESCC$result,allDiff_HNSC$result,allDiff_LUSC$result)
all = as.data.frame(all)
colnames(all)= c("CESC","CSCC","ESCC","HNSC","LUSC")
rownames(all) = rownames(allDiff_ESCC)
all$CESC = as.numeric(all$CESC)
all$CSCC = as.numeric(all$CSCC)
all$ESCC = as.numeric(all$ESCC)
all$HNSC = as.numeric(all$HNSC)
all$LUSC = as.numeric(all$LUSC)
save(all, file="DEGs_5_SCCs.RData")
##
upset(all, nsets = 6, nintersects = 60, mb.ratio = c(0.5, 0.5),
      order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE),
      mainbar.y.label = "Intersection DEGs",point.size = 2.5, 
      line.size = 0.5,number.angles = 30,text.scale = c(1.3, 1.3, 10, 1, 1.5, 1),
      sets.bar.color = c("#E07B54","#E1C855",'#A5C496',"#8891DB","#51B1B7"),
      main.bar.color = c("#E07B54","#E1C855",'#A5C496',"#8891DB","#51B1B7", rep("#363636", 24)),
      queries = list(
        list(
          query = intersects, 
          params = list("LUSC"), 
          color = "#E07B54", 
          active = F), 
        list(
          query = intersects, 
          params = list("CSCC"), 
          color = "#E1C855", 
          active = F), 
        list(
          query = intersects,
          params = list("CESC"), 
          color = '#A5C496',
          active = F), 
        list(
          query = intersects, 
          params = list("ESCC"), 
          color = "#8891DB", 
          active = F), 
        list(
          query = intersects,
          params = list("HNSC"), 
          color = "#51B1B7",
          active = F)
      ),
)

share_DEGs = subset(all,all$CESC=="1"&all$CSCC=="1"&all$ESCC=="1"&all$HNSC=="1"&all$LUSC=="1")
write.csv(all,file = "DEGs.csv")








