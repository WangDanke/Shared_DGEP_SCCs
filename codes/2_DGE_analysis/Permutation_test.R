##创建NullDistribution
set.seed(20240704)
options(stringsAsFactors=F)
#source("http://bioconductor.org/biocLite.R"); biocLite("lme4")
library(lme4)
rootdir = "D:\\SCC_DGEP_shared\\SCCs_array\\DEG_Meta\\null_distribution"
setwd(rootdir)

load("D:\\SCC_DGEP_shared\\SCCs_array\\DEG_Meta\\Microarray_compiledForPermutationTesting.RData") ##Compiled expression & metadata from 2a

genes = rownames(multiExpr[[1]]$datExpr)
allmeta = matrix(NA,nrow=length(genes), length(multiExpr))
colnames(allmeta) = c("LUSC","CSCC","ESCC","HNSC","CESC","EAC","LUAD")
allmeta=  as.data.frame(allmeta)

lmer_apply=function(x, datMeta) {
  if(length(unique(datMeta$ID)) < nrow(datMeta)) {
    return(summary(lmer(x ~ Group + batch + (1 | ID),data=datMeta))$coefficients[2,1])
  } else if(length(unique(datMeta$ID)) > 1) {
    return(summary(lm(x ~ Group + batch,data=datMeta))$coefficients[2,1])
  } else {
    return(summary(lm(x ~ Group + batch,data=datMeta))$coefficients[2,1])
  }
}

n_iterations = 1000
for (iteration in 1:n_iterations) {
  print(paste("Iteration", iteration))
  for(i in 1:length(multiExpr)) {
    tt = matrix(NA, nrow=length(genes), ncol=3)
    subj = unique(as.character(multiExpr[[i]]$datMeta$ID))
    subj_group = data.frame(Subject = subj, Group = multiExpr[[i]]$datMeta$Group[match(subj, multiExpr[[i]]$datMeta$ID)])
    subj_group$Group = subj_group$Group[order(runif(nrow(subj_group)))] ##Randomly shuffle group assignment for each subject
    multiExpr[[i]]$datMeta$Group = subj_group$Group[match(multiExpr[[i]]$datMeta$ID,subj_group$Subject)]
    
    allmeta[,i] = apply(multiExpr[[i]]$datExpr,1,lmer_apply, multiExpr[[i]]$datMeta)
  }
  cor_vec = vector(mode="numeric")
  comparisons = t(combn(seq(1,ncol(allmeta)),2))
  
  for(i in 1:nrow(comparisons)) {
    r = cor(allmeta[,comparisons[i,1]], allmeta[,comparisons[i,2]], method = "spearman", use="pairwise.complete.obs")
    cor_vec = c(cor_vec,r)
  }
  # 
  save(cor_vec, file=paste("D:\\SCC_DGEP_shared\\SCCs_array\\DEG_Meta\\null_distribution\\iteration_", iteration, ".RData", sep=""), row.names=F, col.names=F)
}


##########Fig. 2A
setwd("D:\\SCC_DGEP_shared\\SCCs_array\\DEG_Meta\\null_distribution")

library(ggplot2); library(mada); library(reshape)
library(NMF); library(WGCNA)

#Load null distribution
#
null_cor <- list()

# 
names <- c("LUSC_CSCC","LUSC_ESCC","LUSC_HNSC","LUSC_CESC","LUSC_EAC","LUSC_LUAD","CSCC_ESCC","CSCC_HNSC","CSCC_CESC","CSCC_EAC","CSCC_LUAD",
           "ESCC_HNSC","ESCC_CESC","ESCC_EAC","ESCC_LUAD","HNSC_CESC","HNSC_EAC","HNSC_LUAD","CESC_EAC","CESC_LUAD","EAC_LUAD")

# 
for (name in names) {
  null_cor[[name]] <- data.frame(matrix(nrow = 1000, ncol = 2))
  colnames(null_cor[[name]])=c("null","prob")
}

null_cor$LUSC_CSCC$null = sort(merged_data$LUSC_CSCC)
null_cor$LUSC_CSCC$prob = 1-abs(2*seq(1:length(null_cor$LUSC_CSCC$null))/length(null_cor$LUSC_CSCC$null)-1)

null_cor$LUSC_ESCC$null = sort(merged_data$LUSC_ESCC)
null_cor$LUSC_ESCC$prob = 1-abs(2*seq(1:length(null_cor$LUSC_ESCC$null))/length(null_cor$LUSC_ESCC$null)-1)

null_cor$LUSC_HNSC$null = sort(merged_data$LUSC_HNSC)
null_cor$LUSC_HNSC$prob = 1-abs(2*seq(1:length(null_cor$LUSC_HNSC$null))/length(null_cor$LUSC_HNSC$null)-1)

null_cor$LUSC_CESC$null = sort(merged_data$LUSC_CESC)
null_cor$LUSC_CESC$prob = 1-abs(2*seq(1:length(null_cor$LUSC_CESC$null))/length(null_cor$LUSC_CESC$null)-1)

null_cor$LUSC_EAC$null = sort(merged_data$LUSC_EAC)
null_cor$LUSC_EAC$prob = 1-abs(2*seq(1:length(null_cor$LUSC_EAC$null))/length(null_cor$LUSC_EAC$null)-1)

null_cor$LUSC_LUAD$null = sort(merged_data$LUSC_LUAD)
null_cor$LUSC_LUAD$prob = 1-abs(2*seq(1:length(null_cor$LUSC_LUAD$null))/length(null_cor$LUSC_LUAD$null)-1)

null_cor$CSCC_ESCC$null = sort(merged_data$CSCC_ESCC)
null_cor$CSCC_ESCC$prob = 1-abs(2*seq(1:length(null_cor$CSCC_ESCC$null))/length(null_cor$CSCC_ESCC$null)-1)

null_cor$CSCC_HNSC$null = sort(merged_data$CSCC_HNSC)
null_cor$CSCC_HNSC$prob = 1-abs(2*seq(1:length(null_cor$CSCC_HNSC$null))/length(null_cor$CSCC_HNSC$null)-1)

null_cor$CSCC_CESC$null = sort(merged_data$CSCC_CESC)
null_cor$CSCC_CESC$prob = 1-abs(2*seq(1:length(null_cor$CSCC_CESC$null))/length(null_cor$CSCC_CESC$null)-1)

null_cor$CSCC_EAC$null = sort(merged_data$CSCC_EAC)
null_cor$CSCC_EAC$prob = 1-abs(2*seq(1:length(null_cor$CSCC_EAC$null))/length(null_cor$CSCC_EAC$null)-1)

null_cor$CSCC_LUAD$null = sort(merged_data$CSCC_LUAD)
null_cor$CSCC_LUAD$prob = 1-abs(2*seq(1:length(null_cor$CSCC_LUAD$null))/length(null_cor$CSCC_LUAD$null)-1)

null_cor$ESCC_HNSC$null = sort(merged_data$ESCC_HNSC)
null_cor$ESCC_HNSC$prob = 1-abs(2*seq(1:length(null_cor$ESCC_HNSC$null))/length(null_cor$ESCC_HNSC$null)-1)

null_cor$ESCC_CESC$null = sort(merged_data$ESCC_CSEC)
null_cor$ESCC_CESC$prob = 1-abs(2*seq(1:length(null_cor$ESCC_CESC$null))/length(null_cor$ESCC_CESC$null)-1)

null_cor$ESCC_EAC$null = sort(merged_data$ESCC_EAC)
null_cor$ESCC_EAC$prob = 1-abs(2*seq(1:length(null_cor$ESCC_EAC$null))/length(null_cor$ESCC_EAC$null)-1)

null_cor$ESCC_LUAD$null = sort(merged_data$ESCC_EAC)
null_cor$ESCC_LUAD$prob = 1-abs(2*seq(1:length(null_cor$ESCC_LUAD$null))/length(null_cor$ESCC_LUAD$null)-1)

null_cor$HNSC_CESC$null = sort(merged_data$HNSC_CESC)
null_cor$HNSC_CESC$prob = 1-abs(2*seq(1:length(null_cor$HNSC_CESC$null))/length(null_cor$HNSC_CESC$null)-1)

null_cor$HNSC_EAC$null = sort(merged_data$HNSC_EAC)
null_cor$HNSC_EAC$prob = 1-abs(2*seq(1:length(null_cor$HNSC_EAC$null))/length(null_cor$HNSC_EAC$null)-1)

null_cor$HNSC_LUAD$null = sort(merged_data$HNSC_LUAD)
null_cor$HNSC_LUAD$prob = 1-abs(2*seq(1:length(null_cor$HNSC_LUAD$null))/length(null_cor$HNSC_LUAD$null)-1)

null_cor$CESC_EAC$null = sort(merged_data$CESC_EAC)
null_cor$CESC_EAC$prob = 1-abs(2*seq(1:length(null_cor$CESC_EAC$null))/length(null_cor$CESC_EAC$null)-1)

null_cor$CESC_LUAD$null = sort(merged_data$CESC_LUAD)
null_cor$CESC_LUAD$prob = 1-abs(2*seq(1:length(null_cor$CESC_LUAD$null))/length(null_cor$CESC_LUAD$null)-1)

null_cor$EAC_LUAD$null = sort(merged_data$EAC_LUAD)
null_cor$EAC_LUAD$prob = 1-abs(2*seq(1:length(null_cor$EAC_LUAD$null))/length(null_cor$EAC_LUAD$null)-1)

save(null_cor,file = "null_cor_distribution.RData")

#Make Bargraph (Separate plots of two methods for calculating differential expression) DEG_Combat/DEG_beta
allDiff = allDiff[,c("LUSC","CSCC","ESCC","HNSC","CESC","EAC","LUAD")]
comparisons = t(combn(seq(1,ncol(allDiff)),2))
barplot = data.frame(Mean = NA, SEM=NA, p.fdr=NA)
for (i in 1:dim(comparisons)[1]) {
  x = comparisons[i,1]
  y = comparisons[i,2]
  R = cor.test(allDiff[,x], allDiff[,y])
  rho =cor(allDiff[,x], allDiff[,y], method="spearman", use="pairwise.complete.obs")
  sem = (tanh(atanh(rho + 1.96/sqrt(nrow(allDiff)-3))) - tanh(atanh(rho - 1.96/sqrt(nrow(allDiff)-3))))/3.92
  
  barplot[i,] = c(rho, sem, R$p.value)
  rownames(barplot)[i] = paste(colnames(allDiff)[x],colnames(allDiff)[y],sep="-")
}
barplot$p.fdr = p.adjust(barplot$p.fdr,method="fdr")
barplot$p.bootstrap[1] = null_cor$LUSC_CSCC$prob[findInterval(barplot$Mean, null_cor$LUSC_CSCC$null)][1]
barplot$p.bootstrap[2] = null_cor$LUSC_ESCC$prob[findInterval(barplot$Mean, null_cor$LUSC_ESCC$null)][2]
barplot$p.bootstrap[3] = null_cor$LUSC_HNSC$prob[findInterval(barplot$Mean, null_cor$LUSC_HNSC$null)][3]
barplot$p.bootstrap[4] = null_cor$LUSC_CESC$prob[findInterval(barplot$Mean, null_cor$LUSC_CESC$null)][4]
barplot$p.bootstrap[5] = null_cor$LUSC_EAC$prob[findInterval(barplot$Mean, null_cor$LUSC_EAC$null)][5]
barplot$p.bootstrap[6] = null_cor$LUSC_LUAD$prob[findInterval(barplot$Mean, null_cor$LUSC_LUAD$null)][6]

barplot$p.bootstrap[7] = null_cor$CSCC_ESCC$prob[findInterval(barplot$Mean, null_cor$LUSC_CSCC$null)][7]
barplot$p.bootstrap[8] = null_cor$CSCC_HNSC$prob[findInterval(barplot$Mean, null_cor$LUSC_CSCC$null)][8]
barplot$p.bootstrap[9] = null_cor$CSCC_CESC$prob[findInterval(barplot$Mean, null_cor$LUSC_CSCC$null)][9]
barplot$p.bootstrap[10] = null_cor$CSCC_EAC$prob[findInterval(barplot$Mean, null_cor$CSCC_EAC$null)][10]
barplot$p.bootstrap[11] = null_cor$CSCC_LUAD$prob[findInterval(barplot$Mean, null_cor$CSCC_LUAD$null)][11]

barplot$p.bootstrap[12] = null_cor$ESCC_HNSC$prob[findInterval(barplot$Mean, null_cor$LUSC_CSCC$null)][12]
barplot$p.bootstrap[13] = null_cor$ESCC_CESC$prob[findInterval(barplot$Mean, null_cor$LUSC_CSCC$null)][13]
barplot$p.bootstrap[14] = null_cor$ESCC_EAC$prob[findInterval(barplot$Mean, null_cor$ESCC_EAC$null)][14]
barplot$p.bootstrap[15] = null_cor$ESCC_LUAD$prob[findInterval(barplot$Mean, null_cor$ESCC_LUAD$null)][15]

barplot$p.bootstrap[16] = null_cor$HNSC_CESC$prob[findInterval(barplot$Mean, null_cor$LUSC_CSCC$null)][16]
barplot$p.bootstrap[17] = null_cor$HNSC_EAC$prob[findInterval(barplot$Mean, null_cor$HNSC_EAC$null)][17]
barplot$p.bootstrap[18] = null_cor$HNSC_LUAD$prob[findInterval(barplot$Mean, null_cor$HNSC_LUAD$null)][18]

barplot$p.bootstrap[19] = null_cor$CESC_EAC$prob[findInterval(barplot$Mean, null_cor$CESC_EAC$null)][19]
barplot$p.bootstrap[20] = null_cor$CESC_LUAD$prob[findInterval(barplot$Mean, null_cor$CESC_LUAD$null)][20]

barplot$p.bootstrap[21] = null_cor$EAC_LUAD$prob[findInterval(barplot$Mean, null_cor$EAC_LUAD$null)][21]

barplot$p.symbol = ""
barplot$p.symbol[barplot$p.bootstrap<0.05] = "*"
barplot$p.symbol[barplot$p.bootstrap<0.01] = "**"
barplot$p.symbol[barplot$p.bootstrap<0.001] = "***"
barplot$Comparison = rownames(barplot)
barplot$modality="microarray"

Fig2_A = ggplot(barplot,aes(x = reorder(Comparison, -Mean), y=Mean, label=p.symbol)) + ylim(-0.17,1) +  
  geom_bar(stat="identity",fill="royalblue3",width=0.75) +
  geom_errorbar(aes(ymin=(Mean - SEM), ymax=(Mean + SEM)), position=position_dodge(width=0.8), width=0.25,size=0.25) +   
  #  ggtitle(title) +   
  theme(plot.title = element_text(size=20, face="bold", vjust=2)) +   
  labs(x="", y=expression(paste("Transcriptome correlation (", rho, ")", sep=""))) +     	
  theme(
    axis.text.x=element_text(angle=50, size=12, hjust=1),
    axis.text.y=element_text(size=10, vjust=0.5), 
    legend.title = element_text(size=12), 
    plot.title = element_text(size=12, face="bold"),
    axis.title.x = element_text(size=10, vjust=-0.35, face="bold"),
    axis.title.y = element_text(size=12, vjust=0.5),
    plot.margin=unit(c(2,2,1,2),"mm")
  ) + geom_text(color="red",size=4,aes(y=Mean+ sign(Mean)*SEM + sign(Mean)*.02))
Fig2_A






