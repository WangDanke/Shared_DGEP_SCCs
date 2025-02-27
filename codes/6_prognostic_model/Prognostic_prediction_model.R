##construct and validate the prognostic model
setwd("D:\\SCC_DGEP_shared\\SCCs_array\\Subtype_analysis\\Prognostic_model")
# loading packages
library(survival)
library(survminer)
library(caret)

load("D:\\SCC_DGEP_shared\\SCC_TCGA\\TCGA_4_SCCs_8816x1169.RData")

gene_names <- c("COL1A1", "MMP1", "SERPINE1", "KRT6A", "IGF2BP3", "SPP1")
pheno_SCC_TCGA = pheno_SCC_TCGA[!pheno_SCC_TCGA$project_id=="TCGA-LUSC",]
TCGA_SCC_EXPR = TCGA_SCC_EXPR[,rownames(pheno_SCC_TCGA)]
survival_time = pheno_SCC_TCGA[,c("time","event")]

# 
expr_data <- as.data.frame(t(TCGA_SCC_EXPR[gene_names,]))
expr_data$survival_time <- survival_time$time
expr_data$event <- survival_time$event

expr_data = expr_data[!expr_data$survival_time=="0",]
expr_data$survival_time = expr_data$survival_time/365
expr_data = cbind(log2(expr_data[,1:6]+1),expr_data[,7:8])

# Divide the training set and the validation set (7:3)
set.seed(20201029)
train_index <- createDataPartition(expr_data$event, p = 0.7, list = FALSE)
train_data <- expr_data[train_index, ]
test_data <- expr_data[-train_index, ]

library(glmnet)

# Extract survival data (time and status)
train_y <- train_data[,c('survival_time','event')]
train_y$survival_time <- as.double(train_y$survival_time)
train_y$event <- as.double(train_y$event)
train_y <- as.matrix(survival::Surv(train_y$survival_time, train_y$event))

# Extract the gene expression data in the training set
train_X <- train_data[, 1:6]  

# Lasso regression
fit <- glmnet(train_X, train_y, family = "cox", alpha = 1)
plot(fit)


train_X <- as.matrix(train_X)
# Use cross validation to select the optimal lambda value
cv_fit <- cv.glmnet(train_X, train_y, family = "cox", type.measure = "deviance")

# Draw a lambda path diagram to view the coefficients under different lambda values
plot(cv_fit)

# Get the coefficient corresponding to the optimal lambda (lambda.min)
coef(cv_fit, s = "lambda.min")


selected_genes <- rownames(coef(cv_fit, s = "lambda.min"))[which(coef(cv_fit, s = "lambda.min") != 0)]
selected_genes <- selected_genes[selected_genes != "(Intercept)"]  

# 
formula <- as.formula(paste("Surv(survival_time, event) ~", paste(selected_genes, collapse = "+")))

# 
cox_model <- coxph(formula, data = train_data)
summary(cox_model)  # 查看Cox回归模型的结果

ggforest(model = cox_model,data = train_data, main = 'harzard ratios of candidate genes',fontsize = 1) 


# Calculate risk scores (on both the training and validation sets)）
train_data$risk_score <- predict(cox_model, newdata = train_data, type = "risk")
test_data$risk_score <- predict(cox_model, newdata = test_data, type = "risk")

# Group by quartile and draw survival curves）
q1 <- quantile(train_data$risk_score, 0.25)
q3 <- quantile(train_data$risk_score, 0.75)
train_data$risk_group <- ifelse(
  train_data$risk_score > q3, "High",
  ifelse(train_data$risk_score< q1, "Low", "Medium")
)
train_data_Group= subset(train_data,train_data$risk_group!="Medium")

fit_test <- survfit(Surv(survival_time, event) ~ risk_group, data = train_data_Group)
ggsurvplot(fit_test, data = train_data, pval = TRUE, risk.table = TRUE,
           title = "Training dataset",
           legend.title = "Risk Group", legend.labs = c("High Risk","Low Risk"),
           palette = c( "red","black"))


# Kaplan-Meier生存曲线（测试集）
q1 <- quantile(test_data$risk_score, 0.25)
q3 <- quantile(test_data$risk_score, 0.75)
test_data$risk_group <- ifelse(
  test_data$risk_score > q3, "High",
  ifelse(test_data$risk_score< q1, "Low", "Medium")
)
test_data_Group= subset(test_data,test_data$risk_group!="Medium")


fit_test <- survfit(Surv(survival_time, event) ~ risk_group, data = test_data_Group)
ggsurvplot(fit_test, data = test_data, pval = TRUE, risk.table = TRUE,
           title = "Test dataset",
           legend.title = "Risk Group", legend.labs = c("High Risk","Low Risk"),
           palette = c( "red","black"))



####independent dataset
#####validation in an external dataset
#read the downloaded gene expression matrix of GSE53625
GSE53625 = read.csv("GSE53625_series_matrix.csv")
GSE53625 = GSE53625[-c(1:78),]
colnames(GSE53625) = GSE53625[1,]
GSE53625 = GSE53625[-1,]
rownames(GSE53625) = GSE53625[,1]
iddata = read.csv("ID-Gene_SYMBOL.csv")
colnames(GSE53625)[1] = "probeId"
GSE53625_anno<-merge(x=GSE53625,y=iddata,by='probeId',all.x=F,all.y=F)
GSE53625_anno = GSE53625_anno[,-c(1,361)]

GSE53625 <- GSE53625_anno[!duplicated(GSE53625_anno$geneName), ]
rownames(GSE53625) = GSE53625$geneName
GSE53625 = GSE53625[,-359]
save(GSE53625,file = "GSE53625_EXPR.RData")
GSE53625_pheno = GSE53625_pheno[GSE53625_pheno$Sample_type=="T",]


##validation
load("D:\\SCC_DGEP_shared\\SCCs_array\\Subtype_analysis\\Prognostic_model\\GSE53625_pheno.RData")
load("D:\\SCC_DGEP_shared\\SCCs_array\\Subtype_analysis\\Prognostic_model\\GSE53625_EXPR.RData")

GSE53625_T = GSE53625[,rownames(GSE53625_pheno)]
GSE53625_T_genes <- as.data.frame(t(GSE53625_T[selected_genes,]))
GSE53625_T_genes$survival_time = GSE53625_pheno$survival_time
GSE53625_T_genes$event = GSE53625_pheno$event

GSE53625_T_genes$COL1A1 = as.numeric(GSE53625_T_genes$COL1A1)
GSE53625_T_genes$MMP1 = as.numeric(GSE53625_T_genes$MMP1)
GSE53625_T_genes$SERPINE1 = as.numeric(GSE53625_T_genes$SERPINE1)
GSE53625_T_genes$KRT6A = as.numeric(GSE53625_T_genes$KRT6A)
GSE53625_T_genes$IGF2BP3 = as.numeric(GSE53625_T_genes$IGF2BP3)
GSE53625_T_genes$SPP1 = as.numeric(GSE53625_T_genes$SPP1)

GSE53625_T_genes$event = as.numeric(GSE53625_T_genes$event)
GSE53625_T_genes$survival_time = as.numeric(GSE53625_T_genes$survival_time)

GSE53625_T_genes = as.data.frame(GSE53625_T_genes)

GSE53625_T_genes$risk_score = predict(cox_model, newdata = GSE53625_T_genes, type = "risk")
GSE53625_T_genes = subset(GSE53625_T_genes,GSE53625_T_genes$survival_time < 4.9)

##Group by quartile and draw survival curves
q1 <- quantile(GSE53625_T_genes$risk_score, 0.25)
q3 <- quantile(GSE53625_T_genes$risk_score, 0.75)
GSE53625_T_genes$risk_group <- ifelse(
  GSE53625_T_genes$risk_score > q3, "High",
  ifelse(GSE53625_T_genes$risk_score< q1, "Low", "Medium")
)
GSE53625_T_genes= subset(GSE53625_T_genes,GSE53625_T_genes$risk_group!="Medium")

fit_test <- survfit(Surv(survival_time, event) ~ risk_group, data = GSE53625_T_genes)
ggsurvplot(fit_test, data = GSE53625_T_genes, pval = TRUE, risk.table = TRUE,
           title = "External dataset",
           legend.title = "Risk Group", legend.labs = c("High Risk","Low Risk"),
           palette = c( "red","black"))


##combind all samples （train_test_external datasets)
TCGA_data = rbind(train_data_Group,test_data_Group)
all = rbind(GSE53625_T_genes,TCGA_data[,colnames(GSE53625_T_genes)])


fit_test <- survfit(Surv(survival_time, event) ~ risk_group, data = all)
ggsurvplot(fit_test, data = all, pval = TRUE, risk.table = TRUE,
           title = "Training + test + external datasets",
           legend.title = "Risk Group", legend.labs = c("High Risk","Low Risk"),
           palette = c( "red","black"))
