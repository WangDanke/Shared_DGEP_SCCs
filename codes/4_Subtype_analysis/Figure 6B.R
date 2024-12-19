###Figure 6B

genes= c("COL1A1","MMP1","SERPINE1","PTHLH","KRT6A","IGF2BP3","FADS3","SPP1","SULF1","FSTL3","TREM1","FAP","KIF26B","BNC1","ABCC1","DUSP6")
library(org.Hs.eg.db)
library(clusterProfiler)

GO_results = GO_results@result[c(1:10),]
GO_results$group = "survival_genes"
GO_results<-enrichGO(gene=genes,OrgDb = "org.Hs.eg.db",ont="ALL",
                                          pAdjustMethod = "BH",pvalueCutoff=0.05,qvalueCutoff=0.2,keyType="SYMBOL")

#terms因子顺序
GO_results$Description <- factor(GO_results$Description, levels = GO_results$Description)

#展示的基因，我们选择每个terms展示5个基因，实际情况可以展示自己关注的基因
GO_results$geneID  <- sapply(strsplit(GO_results$geneID , "/"), function(x) paste(x[1:5], collapse = "/"))

library(ggplot2)
plot = ggplot(GO_results, aes(x = -log10(p.adjust), y = rev(Description), fill = group))+
  geom_bar(stat = "identity", width = 0.5,alpha = 0.55)+
  geom_text(aes(x=0.1,y=rev(Description),label = Description),size=3.5, hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(colour = 'black', size = 12),
        axis.line = element_line(colour = 'black', linewidth =0.5),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black'),
        axis.title.x = element_text(colour = 'black', size = 12),
        legend.position = "none",
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),  # Adding a black border around the plot
panel.background = element_blank()  # Removing background color (optional, can keep it white)
) +
  scale_x_continuous(expand = c(0,0))+
  scale_fill_manual(values = c("#ff4757"))+
  geom_text(data = GO_results,
            aes(x = 0.1, y = rev(Description), label = geneID, color = group),
            size = 3,
            fontface = 'italic', 
            hjust = 0,
            vjust = 2.3)+
  scale_color_manual(values = c("black"))+
  scale_y_discrete(expand = c(0.1,0))+
  labs(title = "Enrichment of genes",
       y=c("survival_related genes enriched"))








