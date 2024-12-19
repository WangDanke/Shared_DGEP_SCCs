###correlation of EMT-related genes expression and immune cells infiltration score
data = TCGA_SCC_EXPR[,-c(1:3)]
data = log2(data+1)
data = data[,colnames(nor_data_1_gsva1)]
immune = as.data.frame(t(nor_data_1_gsva1)) 
data = as.data.frame(t(data[c("SNAI1","SNAI2","TWIST1","VIM","ZEB2"),]))

#correlation analysis
# 
cor_results <- data.frame(Column1 = character(), Column2 = character(),
                          Correlation = numeric(), PValue = numeric(), stringsAsFactors = FALSE)

# 
for (col1 in colnames(data)) {
  for (col2 in colnames(immune)) {
    
    # 
    data1 <- data[[col1]]
    data2 <- immune[[col2]]
    
    # 
    test_result <- cor.test(data1, data2, method="spearman")
    
    # 
    cor_results <- rbind(cor_results, data.frame(Column1 = col1, Column2 = col2,
                                                 Correlation = test_result$estimate,
                                                 PValue = test_result$p.value))
  }
}

# 查看结果
print(cor_results)


library(ggplot2)

##Fig. 5E
cor_results$significance <- ifelse(cor_results$PValue < 0.001, "***",
                                   ifelse(cor_results$PValue < 0.01, "**",
                                          ifelse(cor_results$PValue < 0.05, "*", "")))

ggplot(cor_results, aes(Column2,Column1,  fill = Correlation)) +
  geom_tile(color = "black") +  # 绘制矩形块
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +  # 设置颜色渐变
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text.y = element_text(angle = 0, hjust = 1)) +  # 旋转坐标轴标签
  labs(title = "Correlation Heatmap", x = "Immune cell infiltration score", y = "Expression level")+
  geom_text(aes(label = significance), color = "black", size = 3)