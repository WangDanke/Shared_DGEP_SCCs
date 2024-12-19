
allDiff_5 = allDiff_5[rownames(allDiff_5_seq),]
dat = as.data.frame(cbind(allDiff_5[,c("CESC")],allDiff_5_seq[,c("CESC")])) 
colnames(dat) = c("CESC_array","CESC_seq")
library(ggplot2)
library(ggpubr)

ggplot(dat, aes(x = CESC_array, y = CESC_seq,)) +
  geom_point(color =  "#A5C496") +
  geom_smooth(method = "lm", color =  '#A5C496',se = T, fill = "#FFFFFF") +
  stat_cor(method = "pearson", label.x = -1, label.y =6) +
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
  labs(title = "CESC Microarray vs. RNA_seq",
       x = "log2(FC)(Microarray)",
       y = "log2(FC) (RNA_seq)")+
  coord_fixed(ratio = 1, xlim = c(-8,8), ylim = c(-8,8), expand = TRUE, clip = "off")
 

























