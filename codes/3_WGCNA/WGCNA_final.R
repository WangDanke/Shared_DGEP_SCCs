##########WGCNA 
setwd("D:\\SCC_DGEP_shared\\SCCs_array\\WGCNA")
load("D:\\SCC_DGEP_shared\\SCCs_array\\DEG_Meta\\Microarray_compiledForPermutationTesting.RData")

SCC_EXPR = cbind(multiExpr$LUSC$SCC_EXPR,multiExpr$CSCC$SCC_EXPR,multiExpr$ESCC$SCC_EXPR,multiExpr$HNSC$SCC_EXPR,multiExpr$CESC$SCC_EXPR)
pheno_SCC = rbind(multiExpr$LUSC$datMeta,multiExpr$CSCC$datMeta,multiExpr$ESCC$datMeta,multiExpr$HNSC$datMeta,multiExpr$CESC$datMeta)

##QC Pre-Combat,Fig. S6
plot(density(SCC_EXPR[,1]), xlim=c(-5,20), ylim=c(0, 0.5), col = as.numeric(pheno_SCC$batch[1]), xlab="Intensity (log2)", ylab="Density", main="Mega-Analysis: Pre-Combat")
for(i in 2:dim(SCC_EXPR)[[2]])
  lines(density(SCC_EXPR[,i]), xlim=c(0,20), col = as.numeric(pheno_SCC$batch[i]))  
legend("topleft", (levels(pheno_SCC$batch)), col=c(1:8), pch=16, cex=0.5)


##batch correct
library(sva)
mod = model.matrix(~sample_type, data=pheno_SCC)
batch = factor(pheno_SCC$batch)
SCC_EXPR = ComBat(SCC_EXPR, batch=batch, mod=mod)
SCC_EXPR = as.data.frame(SCC_EXPR)
save(pheno_SCC,SCC_EXPR,file = "SCC_combined_combat_for_WGCNA.RData")

##QC Post-Combat Fig. S6
plot(density(SCC_EXPR[,1]), xlim=c(-5,20), ylim=c(0, 0.5), col = as.numeric(pheno_SCC$batch[1]), xlab="Intensity (log2)", ylab="Density", main="Mega-Analysis: Post-Combat")
for(i in 2:dim(SCC_EXPR)[[2]])
  lines(density(SCC_EXPR[,i]), xlim=c(0,20), col = as.numeric(pheno_SCC$batch[i]))  
legend("topleft", (levels(pheno_SCC$batch)), col=c(1:8), pch=16, cex=0.5)

#Dendrogram
par(mfrow=c(1,1))
tree = hclust(dist(t(SCC_EXPR)), method = "average")
plotDendroAndColors(tree, cbind(as.numeric(pheno_SCC$Group), as.numeric(pheno_SCC$batch)), cex.colorLabels=0.6, cex.dendroLabels=0.2,
                    main="Dendrogram\nPost-Combat")

save(file = "./All_datasets_combined_8816x805.RData", SCC_EXPR,pheno_SCC,multiExpr)

## NETWORK ANALYSIS
library(WGCNA)
## ----------------
####select suitable sft, power = 6
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
data <- as.data.frame(t(combat_data_combined))
sft <- pickSoftThreshold(data, powerVector = powers)
str(sft)

# Plot the results:
pdf("D:\\SCC_DGEP_shared\\SCCs_array\\WGCNA\\WGCNA-softthresh.pdf")
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "SFT, signed R^2", type = "n", main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.9, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], type = "n", xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
dev.off()

#######################
#construct "signed"network
TOM = TOMsimilarityFromExpr(t(SCC_EXPR), power = 6,networkType = "signed",TOMType = "signed")
save(TOM,file = "signed_TOM.RData")
geneTree = hclust(1-as.dist(TOM), method="average")

# Iterate WGCNA parameters for robustness -- this takes a while
colors = vector(mode="list")
labels = vector(mode="list")
for (pam in c(FALSE,TRUE)) {
  for (minModSize in c(30,50,100)) {
    for (dthresh in c(0.1, 0.2)) {
      for(ds in c(0:4)) { 
        print(paste("DS=", ds, ",MMS=", minModSize, ",DCOR=",dthresh,",PAM=",pam,sep=""))
        
        tree = cutreeHybrid(dendro = geneTree, minClusterSize= minModSize, pamStage=pam, cutHeight = 0.999, deepSplit=ds, distM=as.matrix(1-as.dist(TOM)))
        merged = mergeCloseModules(exprData= t(SCC_EXPR), colors = tree$labels, cutHeight=dthresh)
        colors = cbind(colors, labels2colors(merged$colors))
        
        labels = c(labels, paste("DS=", ds, ",MMS=", minModSize, ",DCOR=",dthresh,",PAM=",pam,sep=""))
      }
    }
  }
}
save(file="D:\\SCC_DGEP_shared\\SCCs_array\\WGCNA_diffParams_final.RData", geneTree, colors, labels)

colors2 = colors
colors2[,seq(1:30)] = colors[,seq(1,60,by=2)]
colors2[,seq(31:60)] = colors[,seq(2,60,by=2)]

pdf("WGCNA_params_signed.pdf", width = 15, height = 18)
plotDendroAndColors(geneTree, colors, groupLabels=labels, addGuide= TRUE, dendroLabels=FALSE, main="Dendrogram", cex.colorLabels=0.5)
dev.off()

##
#Determine the final Parameters
# --------------------
# Parameters to Use: "DS=4,MMS=100,DCOR=0.1,PAM=FALSE"
wgcna_parameters = list(powers =  6)
wgcna_parameters$minModSize = 50
wgcna_parameters$minHeight = 0.1
wgcna_parameters$bsize = 10000  ##block size needs to be larger than dim(datExpr)[1]
wgcna_parameters$ds = 4  ##deep split parameter contorls number of modules
wgcna_parameters$networkType = "signed"    ## using signed networks
wgcna_parameters$corFnc = "bicor"
wgcna_parameters$pamStage = FALSE

tree = cutreeHybrid(dendro = geneTree, minClusterSize= wgcna_parameters$minModSize, pamStage=wgcna_parameters$pamStage, cutHeight = 0.999, 
                    deepSplit=wgcna_parameters$ds, distM=as.matrix(1-as.dist(TOM)))
merged = mergeCloseModules(exprData= t(SCC_EXPR), colors = tree$labels, cutHeight=wgcna_parameters$minHeight)
colors = labels2colors(merged$colors)
table(colors)
length(table(colors))

#Fig. 3A
pdf("WGCNA_final_signed_4_50_0.1_FALSE_mini.pdf", width = 8, height = 3)
plotDendroAndColors(geneTree,colors,groupLabels = "mod",cex.colorLabels = 0.5,addGuide=T,dendroLabels=F)
dev.off()

MEs = moduleEigengenes(expr = t(SCC_EXPR), colors, softPower = wgcna_parameters$powers)
kMEtable = signedKME(t(SCC_EXPR),MEs$eigengenes)
tableS1 = data.frame(kMEtable[,paste0("kME", labels2colors(1:10))])
colnames(kMEtable) = paste0("kME.CD", 1:10, ".", labels2colors(1:10))

tableS1 = cbind(data.frame(Module.Color=colors, Module.name = paste0("CD",merged$colors)), kMEtable)

write.csv(file="TableS1 - kME table.csv", tableS1)
save(file="finalizedNetwork_20240425.RData", SCC_EXPR, pheno_SCC, geneTree, colors, wgcna_parameters, MEs, kMEtable)

#Fig. S7A
###plot combined WGCNA_params and final network
# Define a matrix of labels for the original and all resampling runs
c_ref = as.character(colors)##这里的reference colors是final 参数的网络的color

# Relabel modules in each of the resampling runs so that full and reampled modules with best overlaps have
# the same labels. This is achieved by the function matchLabels.
for (i in 1:60){
  ci = as.character(colors[,i])
  c_new = matchLabels(ci, c_ref)
  colors[,i] = c_new
}

colors = cbind(colors[,15], colors) ###the best params is the fifteenth
labels = c("Final Modules", labels)

pdf("WGCNA_diffParams_final.pdf",width=5,height=6)
plotDendroAndColors(geneTree,colors,groupLabels = labels,addGuide=T,dendroLabels=F,cex.colorLabels=0.3)
dev.off()

###模块之间的相关性图,Fig 3B
MEs_col = orderMEs(MEs$eigengenes)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)

###确定模块和临床信息的相关性
pheno_SCC_groups = read.csv("pheno_SCC.csv",header = T,row.names = 1)
match_trait <- match(colnames(SCC_EXPR), pheno_SCC_groups$ID)
Traits <- pheno_SCC_groups[match_trait, -c(1,2,3,4,5,6)]
modTraitCor <- cor(MEs$eigengenes[,-5], Traits, use="p")
modTraitP <- corPvalueStudent(modTraitCor, 805)
textMatrix <- paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")",sep = "")
dim(textMatrix) <- dim(modTraitCor)

##并绘制相关性热图Fig 3C
pdf(file = "module_trait_relationship.pdf", wi=4.5, h=7)
labeledHeatmap(
  Matrix = modTraitCor, xLabels = names(Traits), yLabels = names(MEs$eigengenes[,-5]),
  ySymbols = names(MEs$eigengenes[,-5]), colorLabels = FALSE, colors = blueWhiteRed(50),
  textMatrix = textMatrix, setStdMargins = FALSE,
  cex.text = 0.7, zlim = c(-1, 1), main = paste("Module-trait relationships"))
dev.off()

##########对每个模块进行GO富集分析
library(org.Hs.eg.db)
library(clusterProfiler)

net_information = data.frame(matrix(nrow = 8816, ncol = 2))
colnames(net_information) = c("gene","module")
net_information$gene = rownames(SCC_EXPR)
net_information$module = colors
rownames(net_information) = net_information$gene
save(net_information,file = "net_information.RData")

module1 <- subset(net_information, module == "turquoise" )
module2 <- subset(net_information, module == "blue" )
module3 <- subset(net_information, module == "brown" )
module4 <- subset(net_information, module == "yellow" )
module5 <- subset(net_information, module == "green" )
module6 <- subset(net_information, module == "red" )
module7 <- subset(net_information, module == "black" )
module8 <- subset(net_information, module == "pink" )
module9 <- subset(net_information, module == "magenta" )

###GO analysis
GO_module1<-enrichGO(gene=rownames(module1),OrgDb = "org.Hs.eg.db",ont="ALL",
                     pAdjustMethod = "BH",pvalueCutoff=0.05,qvalueCutoff=0.2,keyType="SYMBOL")
GO_module2<-enrichGO(gene=rownames(module2),OrgDb = "org.Hs.eg.db",ont="ALL",
                     pAdjustMethod = "BH",pvalueCutoff=0.05,qvalueCutoff=0.2,keyType="SYMBOL")
GO_module3<-enrichGO(gene=rownames(module3),OrgDb = "org.Hs.eg.db",ont="ALL",
                     pAdjustMethod = "BH",pvalueCutoff=0.05,qvalueCutoff=0.2,keyType="SYMBOL")
GO_module4<-enrichGO(gene=rownames(module4),OrgDb = "org.Hs.eg.db",ont="ALL",
                     pAdjustMethod = "BH",pvalueCutoff=0.05,qvalueCutoff=0.2,keyType="SYMBOL")
GO_module5<-enrichGO(gene=rownames(module5),OrgDb = "org.Hs.eg.db",ont="ALL",
                     pAdjustMethod = "BH",pvalueCutoff=0.05,qvalueCutoff=0.2,keyType="SYMBOL")
GO_module6<-enrichGO(gene=rownames(module6),OrgDb = "org.Hs.eg.db",ont="ALL",
                     pAdjustMethod = "BH",pvalueCutoff=0.05,qvalueCutoff=0.2,keyType="SYMBOL")
GO_module7<-enrichGO(gene=rownames(module7),OrgDb = "org.Hs.eg.db",ont="ALL",
                     pAdjustMethod = "BH",pvalueCutoff=0.05,qvalueCutoff=0.2,keyType="SYMBOL")
GO_module8<-enrichGO(gene=rownames(module8),OrgDb = "org.Hs.eg.db",ont="ALL",
                     pAdjustMethod = "BH",pvalueCutoff=0.05,qvalueCutoff=0.2,keyType="SYMBOL")
GO_module9<-enrichGO(gene=rownames(module9),OrgDb = "org.Hs.eg.db",ont="ALL",
                     pAdjustMethod = "BH",pvalueCutoff=0.05,qvalueCutoff=0.2,keyType="SYMBOL")

##Fig. 3D
#######将各个模块的信息分别展示
GO_results = GO_module9@result[c(1:5),]
GO_results$group = "CD9_magenta"
library(ggplot2)
# 使用排序索引重新排列数据框
GO_results <- GO_results[order(GO_results$group), ]
#terms因子顺序
GO_results$Description <- factor(GO_results$Description, levels = GO_results$Description)

plot = ggplot(GO_results, aes(x = -log10(p.adjust), y = rev(Description), fill = group)) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.55) +
  geom_text(aes(x = 0.1, y = rev(Description), label = Description), size = 3.5, hjust = 0) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(colour = 'black', size = 12),
    axis.line = element_line(colour = 'black', linewidth = 0.5),
    axis.text.x = element_text(colour = 'black', size = 10),
    axis.ticks.x = element_line(colour = 'black'),
    axis.title.x = element_text(colour = 'black', size = 12),
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),  # Adding a black border around the plot
    panel.background = element_blank()  # Removing background color (optional, can keep it white)
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c("#EE00EE")) +
  labs(title = "CD9_magenta", y = " ")

plot





