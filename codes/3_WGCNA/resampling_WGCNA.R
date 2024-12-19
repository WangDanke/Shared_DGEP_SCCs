#######resampling in linux
#!/usr/bin/anaconda3/envs/R4.3/bin/Rscript
library(WGCNA)
load("/mnt/data/wangdanke/EAC_array/SCC_combined_combat_for_WGCNA.RData")
# 设置参数
#resampling
nRuns = 100
power = 6
deepSplit = 4
minModuleSize = 50
networkType = "signed"
TOMType = "signed"
TOMDenom = "mean"
reassignThreshold = 0
mergeCutHeight = 0.1

#running resampling net construction 100 times
tmf0 = system.time ( {
  mods0 = sampledBlockwiseModules(
    nRuns = nRuns,
    replace = TRUE,
    datExpr = t(SCC_EXPR),
    maxBlockSize = 10000,
    networkType = networkType,
    fraction = 0.63,
    TOMType = TOMType,
    TOMDenom = TOMDenom,
    deepSplit = deepSplit,
    randomSeed = 12345,
    mergeCutHeight = mergeCutHeight,
    reassignThreshold = reassignThreshold,
    imputeMissing = FALSE,
    numericLabels = TRUE,
    checkMissingData = FALSE,
    quickCor = 0, verbose = 2 ) } )

print(tmf0)
save(tmf0, mods0, file = "/mnt/data/wangdanke/EAC_array/net_100.RData")

##Fig. S7B
###plot combined WGCNA_params and final network
# Define a matrix of labels for the original and all resampling runs
nGenes = ncol(t(SCC_EXPR))
# Define a matrix of labels for the original and all resampling runs
labels = matrix(0, nGenes, nRuns + 1)
labels[, 1] = colors ##这里的reference colors是final 参数的网络的color，为与上述的图对应起来，所有的网络配色都以最后的网络颜色作为标准

# Relabel modules in each of the resampling runs so that full and reampled modules with best overlaps have
# the same labels. This is achieved by the function matchLabels.
pind = initProgInd()
for (r in 2:(nRuns+1))
{
  labels[, r] = matchLabels(mods0[[r]]$mods$colors, labels[, 1])
  pind = updateProgInd((r-1)/nRuns, pind)
}
# Save the results
save(labels, file = "AA_Unsigned_ModuleExample-matchedLabels.RData")

pdf(file = "AA_signed_sampledModuleExample-dendrogramAndSampledColors_2.pdf", wi=20, h=20)
plotDendroAndColors(mods0[[1]]$mods$dendrograms[[1]],
                    labels,
                    c("Full data set", paste("Resampling", c(1:nRuns))),
                    main = "Gene dendrogram and module labels from resampled data sets",
                    autoColorHeight = FALSE, colorHeight = 0.8,
                    dendroLabels = FALSE, hang = 0.03, guideHang = 0.05,
                    addGuide = TRUE,
                    guideAll = FALSE,
                    cex.main = 2, cex.lab = 1.6, cex.colorLabels = 0.8, marAll = c(0, 5, 3, 0))
dev.off()