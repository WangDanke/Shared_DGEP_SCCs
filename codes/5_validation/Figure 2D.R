
###The similarity values of SCC pairs obtained from the discovery set and the validation set are recorded in correlations.csv
library(ggplot2)
library(ggpubr)
####
correlation = read.csv("correlation_both.csv",header = T,row.names = 1)
correlation = as.data.frame(t(correlation))
correlation$group = rownames(correlation)
#Fig. 2D
ggplot(correlation, aes(x = Microarray, y = RNA_seq,)) +
  geom_point() +
  geom_smooth(method = "lm", color = "#CD8500",se = T, fill = "#EEC591") +
  stat_cor(method = "pearson", label.x = 0.3, label.y = 0.8) +
  theme(
    axis.text.x=element_text(size=10, hjust=1,color = "#828282"),
    axis.text.y=element_text(size=10, vjust=0.5,color = "#828282"), 
    legend.title = element_text(size=12,color = "grey"), 
    plot.title = element_text(size=12, face="bold",color = "#363636"),
    axis.title.x = element_text(size=12,color = "#363636"),
    axis.title.y = element_text(size=12, vjust=0.5,color = "#363636"),
    plot.margin=unit(c(2,2,16,2),"mm"),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "#363636", fill = NA, size = 1)
  ) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black",size = 0.55)+
  labs(title = "Microarray vs. RNA_seq",
       x = "Transcriptome correlation (Microarray)",
       y = "Transcriptome correlation (RNA_seq)")+
  coord_fixed(ratio = 1, xlim = c(0.2,0.8), ylim = c(0.2,0.8), expand = TRUE, clip = "off")+
  geom_text(aes(label = group), size = 3, vjust = -0.5, hjust = 0.1, color = "#CD8500") 
  
























