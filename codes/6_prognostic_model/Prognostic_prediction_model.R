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


########################
##lambda sensitivity analysis
lambda_seq <- cv_fit$lambda

# lambdas that still retain ≥1 non-zero coefficient.
nz <- sapply(lambda_seq, function(l) {
  sum(coef(cv_fit, s = l) != 0)
})

valid_lambda <- lambda_seq[nz >= 2]

lambda_A <- valid_lambda[1]     
lambda_B <- valid_lambda[round(length(valid_lambda)/2)]
lambda_C <- valid_lambda[length(valid_lambda)] 

risk_A <- predict(cv_fit, train_X, s = lambda_A, type = "link")
risk_B <- predict(cv_fit, train_X, s = lambda_B, type = "link")
risk_C <- predict(cv_fit, train_X, s = lambda_C, type = "link")

cor(risk_A, risk_B, method = "spearman")
cor(risk_A, risk_C, method = "spearman")
plot(risk_A, risk_C,
     xlab = "Risk score (λ_min)",
     ylab = "Risk score (λ_1se)",
     main = "Sensitivity analysis of λ choice")
abline(0,1,col="red",lty=2)


# Draw a lambda path diagram to view the coefficients under different lambda values
plot(cv_fit)

# Get the coefficient corresponding to the optimal lambda (lambda.min)
coef(cv_fit, s = "lambda.min")


gene_names<- rownames(coef(cv_fit, s = "lambda.min"))[which(coef(cv_fit, s = "lambda.min") != 0)]
gene_names <- gene_names[gene_names != "(Intercept)"]  

# 
formula <- as.formula(paste("Surv(survival_time, event) ~", paste(gene_names, collapse = "+")))

# 
cox_model <- coxph(formula, data = train_data)
summary(cox_model)  # check the results
ggforest(model = cox_model,data = train_data, main = 'harzard ratios of candidate genes',fontsize = 1) 


# Calculate risk scores (on both the training and validation sets)）
train_data$risk_score <- predict(cox_model, newdata = train_data, type = "risk")
test_data$risk_score <- predict(cox_model, newdata = test_data, type = "risk")

cox_cont <- coxph(Surv(survival_time, event) ~ risk_score, data = test_data)
summary(cox_cont)
library(survival)
library(rms)

dd <- datadist(test_data)
options(datadist = "dd")

fit_spline <- cph(Surv(survival_time, event) ~ rcs(risk_score, 3),
                  data = test_data, x = TRUE, y = TRUE)

plot(Predict(fit_spline, risk_score, fun = exp),
     xlab = "Risk score",
     ylab = "Hazard ratio")



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

#as the continuous variable
train_data_Group$risk_group = factor(train_data_Group$risk_group,levels = c("Low", "High"))
cox_fit <- coxph(Surv(survival_time, event) ~ risk_group, data = train_data_Group)
summary(cox_fit)

# Kaplan-Meier test dataset
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

test_data_Group$risk_group = factor(test_data_Group$risk_group,levels = c("Low", "High"))
cox_fit <- coxph(Surv(survival_time, event) ~ risk_group, data = test_data_Group)
summary(cox_fit)

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
GSE53625_T_genes <- as.data.frame(t(GSE53625_T[gene_names,]))
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


##
mu  <- mean(train_data$risk_score)
sdv <- sd(train_data$risk_score)

GSE53625_T_genes$risk_score_z <- 
  (GSE53625_T_genes$risk_score - mu) / sdv

cox_cont <- coxph(Surv(survival_time, event) ~ risk_score, data = GSE53625_T_genes)
summary(cox_cont)
library(survival)
library(rms)

dd <- datadist(GSE53625_T_genes)
options(datadist = "dd")

fit_spline <- cph(Surv(survival_time, event) ~ rcs(risk_score_z, 3),
                  data = GSE53625_T_genes, x = TRUE, y = TRUE)

plot(Predict(fit_spline, risk_score_z, fun = exp),
     xlab = "Risk score",
     ylab = "Hazard ratio")


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



# Group based on the medium
median_risk <- median(GSE53625_T_genes$risk_score)
GSE53625_T_genes$risk_group_median <- ifelse(
  GSE53625_T_genes$risk_score >= median_risk, "High", "Low"
)

# survival_analysis
fit_test_median <- survfit(Surv(survival_time, event) ~ risk_group_median, 
                           data = GSE53625_T_genes)

# Km
ggsurvplot(fit_test_median, 
           data = GSE53625_T_genes, 
           pval = TRUE, 
           risk.table = TRUE,
           title = "test dataset (Median split)",
           legend.title = "Risk Group", 
           legend.labs = c("High Risk", "Low Risk"),
           palette = c("red", "black"))

GSE53625_T_genes$risk_group_median = factor(GSE53625_T_genes$risk_group_median,levels = c("Low", "High"))
cox_fit <- coxph(Surv(survival_time, event) ~ risk_group_median, data = GSE53625_T_genes)
summary(cox_fit)














