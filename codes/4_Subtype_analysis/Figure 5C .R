data = TCGA_SCC_EXPR[,-c(1:3)]
data = data[,rownames(pheno_SCC_TCGA)]
data = log2(data+1)
library(dplyr)
pheno_SCC_TCGA<-pheno_SCC_TCGA %>%
  mutate(Subtype_3SCC_10 = case_when(
    Subtype_3SCC_10 == 1 ~ "Subtype1",
    Subtype_3SCC_10 == 2 ~ "Subtype2",
    Subtype_3SCC_10 == 3 ~ "Subtype3",
    Subtype_3SCC_10 == 4 ~ "Subtype4"
  ))

data = data[,rownames(pheno_SCC_TCGA)]
data <- as.data.frame(t(data))
data$SUBTYPE<-pheno_SCC_TCGA$Subtype_3SCC_10
data$SUBTYPE<-as.factor(data$SUBTYPE)

compaired <- list(c("Subtype1", "Subtype2"),c("Subtype2","Subtype3"),c("Subtype1","Subtype3"),c("Subtype3", "Subtype4"),c("Subtype2","Subtype4"),c("Subtype1","Subtype4"))

library(ggplot2)
library(ggpubr)
##Plot a boxplot of a single gene
ggplot(data, aes(x=SUBTYPE, y=VIM, fill=SUBTYPE)) + 
  geom_violin(trim = F,alpha = 0.8)+
  scale_fill_manual(values = c("Subtype1"="#C6B3D3","Subtype2"="#ED9F9B","Subtype3"="#80BA8A","Subtype4"="#9CD1C8"))+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)+
  theme(panel.grid=element_blank(), 
        panel.background=element_rect(color="black", fill="transparent"))+
  geom_boxplot(width = 0.1, color = "white", outlier.shape = NA)
