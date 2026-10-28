rm(list=ls())
library(GGally)
library(corrplot)
library(ggplot2)
library(VIM)
library(glmnet)
library(cluster)
library(grpreg)
library(Boruta)
library(ggrepel)
library(randomForest)
library(lme4)
library(lattice)
library(caret)
library(pROC)
library(dcurves)
library(glmnet)
library(Boruta)
library(splines)
library(dplyr)

data <- read.csv('dataclean_final_withNA.csv', header=T)

#preprocessing 
sum(is.na(data))
colnames(data)
colnames(data)[colnames(data) == "Overall.Survival.at.2.years"] = "OS2ys"
colnames(data)[colnames(data) == "mean_dose_vanacava_S"] <- "mean_dose_venacava_S"
colnames(data)[colnames(data) == "eud_venacavaI"] <- "eud_venacava_I"
colnames(data)[colnames(data) == "eud_vanacavaS"] <- "eud_venacava_S"
data[data == ""] <- NA
data$KPS=factor(data$KPS)
data$PDL.1.level=factor(data$PDL.1.level)
data$schedule_ok_chemo=factor(data$schedule_ok_chemo)
data$cycles_chemo_sum = factor(data$cycles_chemo_sum)
data$cycles_immuno = factor(data$cycles_immuno)
data$prescription_dose_radio = factor(data$prescription_dose_radio)
data$number_of_fractions_radio = factor(data$number_of_fractions_radio)
data$number_of_fractions_radio = factor(data$number_of_fractions_radio)
data[] <- lapply(data, function(x) {
  if (!is.numeric(x)) as.factor(x) else x
})
data$OS2ys_num <- ifelse(data$OS2ys == "YES", 1, 0)
data$hospital <- factor(substr(as.character(data$ID), 3, 3), levels = as.character(1:7))
data$KPS <- as.numeric(as.character(data$KPS))
data$aortic_implant_present = as.factor(data$aortic_implant_present)
data$cardiac_implant_present = as.factor(data$cardiac_implant_present)

#put Volume PTV in CRF
cols <- names(data)
idx <- which(cols == "Stage.at.RT")
data <- data[, c(cols[1:idx], "Volume_PTV", setdiff(cols[(idx+1):length(cols)], "Volume_PTV"))]

datacrf = data[, c(1:41, 92)]
datadose = data[, 42:64]
datavol = data[, 65:77]
datametab = data[, 78:91]



#data preparation 
data$Smoking.status <- ifelse(is.na(data$Smoking.status), NA,ifelse(data$Smoking.status %in% c("YES", "EX"), 1, 0))
data$Diabetes.drugs <- ifelse(data$Diabetes.drugs == "YES", 1, 0)
data$Antilipidic.drugs <- ifelse(data$Antilipidic.drugs == "YES", 1, 0)
data$Hypertension <- ifelse(data$Hypertension == "YES", 1, 0)
data$antithrombotic <- ifelse(data$antithrombotic == "YES", 1, 0)
data$Beta.Blockers <- ifelse(data$Beta.Blockers == "YES", 1, 0)
data$Diuretics <- ifelse(data$Diuretics == "YES", 1, 0)
data$risk_sum <- rowSums(data[, c("Smoking.status", "Diabetes.drugs", "Antilipidic.drugs", "Hypertension", "antithrombotic", "Beta.Blockers", "Diuretics")], na.rm = FALSE) 
data$riskfactors <- cut(data$risk_sum, breaks = c(-Inf, 0, 2, 3, 5, Inf), labels = c(0, 1, 2, 3, 4))
data$riskfactors1 <- rowSums(data[, c("Smoking.status","Diabetes.drugs", "Antilipidic.drugs", "Hypertension")], na.rm = FALSE) 
data$risk_sum2 <- rowSums(data[, c("Beta.Blockers", "Diuretics")], na.rm = FALSE)
data$riskfactors2 = ifelse(data$risk_sum2>0,1,0)
data$riskfactors = as.factor(data$riskfactors)
data$riskfactors1 = as.factor(data$riskfactors1)
data$riskfactors2 = as.factor(data$riskfactors2)


# datasets1 
# hole
colnames(data)
data_global <- data[, -c(1,2,4,5,6:12,13,14,16,17,18,20:25,27:41,65:77,93:95,97)] #riskfactors1, riskfactors2
colnames(data_global)
data_global <- na.omit(data_global)
X_glob <- model.matrix(OS2ys_num ~ ., data_global)[, -1]

datacrf$Smoking.status <- ifelse(is.na(datacrf$Smoking.status), NA,ifelse(datacrf$Smoking.status %in% c("YES", "EX"), 1, 0))
datacrf$Diabetes.drugs <- ifelse(datacrf$Diabetes.drugs == "YES", 1, 0)
datacrf$Antilipidic.drugs <- ifelse(datacrf$Antilipidic.drugs == "YES", 1, 0)
datacrf$Hypertension <- ifelse(datacrf$Hypertension == "YES", 1, 0)
datacrf$antithrombotic <- ifelse(datacrf$antithrombotic == "YES", 1, 0)
datacrf$Beta.Blockers <- ifelse(datacrf$Beta.Blockers == "YES", 1, 0)
datacrf$Diuretics <- ifelse(datacrf$Diuretics == "YES", 1, 0)
datacrf$risk_sum <- rowSums(datacrf[, c("Smoking.status", "Diabetes.drugs", "Antilipidic.drugs", "Hypertension", "antithrombotic", "Beta.Blockers", "Diuretics")], na.rm = FALSE) 
datacrf$riskfactors <- cut(datacrf$risk_sum, breaks = c(-Inf, 0, 2, 3, 5, Inf), labels = c(0, 1, 2, 3, 4))
datacrf$riskfactors1 <- rowSums(datacrf[, c("Smoking.status", "Diabetes.drugs", "Antilipidic.drugs", "Hypertension")], na.rm = FALSE) 
datacrf$risk_sum2 <- rowSums(datacrf[, c("Beta.Blockers", "Diuretics")], na.rm = FALSE)
datacrf$riskfactors2 = ifelse(datacrf$risk_sum2>0,1,0)
datacrf$riskfactors = as.factor(datacrf$riskfactors)
datacrf$riskfactors1 = as.factor(datacrf$riskfactors1)
datacrf$riskfactors2 = as.factor(datacrf$riskfactors2)
colnames(datacrf)
Y_glob <- as.factor(data_global$OS2ys_num)

#CRF 
datacrf_model <- datacrf[, -c(1,2,4,5,6:12,13,14,16,17,18,20:25,27:41,43,44,46)]
colnames(datacrf_model)
datacrf_model <- na.omit(datacrf_model)
X_crf <- model.matrix(OS2ys_num ~ ., datacrf_model)[, -1]
Y_crf <- as.factor(datacrf_model$OS2ys_num)

#dose 
datadose_model <- na.omit(datadose)
X_dose <- model.matrix(data$OS2ys_num~ ., datadose_model)[, -1]
Y_dose <- as.factor(data$OS2ys_num)
#metab 
colnames(datametab)
datametab_model <- na.omit(datametab)
X_metab <-model.matrix(data$OS2ys_num~.,datametab)[,-1]
righe_valide <- as.numeric(rownames(X_metab))
ID_validi = datacrf$ID[righe_valide]
Y_metab <- as.factor(data$OS2ys_num[righe_valide])

dim(data_global)
dim(datacrf_model)
dim(datadose_model)
dim(datametab_model)








## hole1 ----
###LASSO
set.seed(123)
cv_glob <- cv.glmnet(X_glob, Y_glob, family = "binomial", alpha = 1)
lasso_glob <- glmnet(X_glob, Y_glob, family = "binomial", alpha = 1, lambda = cv_glob$lambda.min)
print(coef(lasso_glob))
### bootstrapping LASSO 
X <- X_glob
Y <- Y_glob
B <- 1000
n <- nrow(X)
selection_matrix <- matrix(0, nrow = B, ncol = ncol(X))
colnames(selection_matrix) <- colnames(X)
set.seed(123)
for (i in 1:B) {
  boot_idx <- sample(1:n, size = n, replace = TRUE)
  X_boot <- X[boot_idx, ]
  Y_boot <- Y[boot_idx]
  cv_fit <- cv.glmnet(X_boot, Y_boot, family = "binomial", alpha = 1)
  coef_boot <- as.matrix(coef(cv_fit, s = "lambda.min"))
  non_zero_features <- rownames(coef_boot)[coef_boot[, 1] != 0][-1]
  if(length(non_zero_features) > 0) {
    selection_matrix[i, non_zero_features] <- 1
  }
}
selection_freq <- colMeans(selection_matrix)
freq_ordinata <- sort(selection_freq, decreasing = TRUE)
features_di_ferro_lasso <- names(selection_freq[selection_freq >= 0.80])
print(freq_ordinata[freq_ordinata > 0.1])
features_di_ferro_lasso
### Boruta 
set.seed(123)
df_boruta <- as.data.frame(X)
df_boruta$OS2ys_num <- Y
boruta_res <- Boruta(OS2ys_num ~ ., data = df_boruta, pValue = 0.05, mcAdj = TRUE, maxRuns = 2000, doTrace = 1) 
print(boruta_res)
boruta_res$finalDecision
plot(boruta_res, las = 2, cex.axis = 0.7, main = "Boruta Feature Importance")
features_di_ferro_boruta <- getSelectedAttributes(boruta_res, withTentative = FALSE)
features_di_ferro_boruta
feature_super_definitive <- intersect(features_di_ferro_lasso, features_di_ferro_boruta)
print(feature_super_definitive)



## specific1 ----
###CRF----
#### LASSO 
set.seed(123)
cv_crf <- cv.glmnet(X_crf, Y_crf, family = "binomial", alpha = 1)
lasso_crf <- glmnet(X_crf, Y_crf, family = "binomial", alpha = 1, lambda = cv_crf$lambda.min)
print(coef(lasso_crf))
####bootrstrapping LASSO 
X <- X_crf
Y <- Y_crf
B <- 1000
n <- nrow(X)
selection_matrix <- matrix(0, nrow = B, ncol = ncol(X))
colnames(selection_matrix) <- colnames(X)
set.seed(123)
for (i in 1:B) {
  boot_idx <- sample(1:n, size = n, replace = TRUE)
  X_boot <- X[boot_idx, ]
  Y_boot <- Y[boot_idx]
  cv_fit <- cv.glmnet(X_boot, Y_boot, family = "binomial", alpha = 1)
  coef_boot <- as.matrix(coef(cv_fit, s = "lambda.min"))
  non_zero_features <- rownames(coef_boot)[coef_boot[, 1] != 0][-1]
  if(length(non_zero_features) > 0) {
    selection_matrix[i, non_zero_features] <- 1
  }
}
selection_freq <- colMeans(selection_matrix)
freq_ordinata <- sort(selection_freq, decreasing = TRUE)
features_di_ferro_lasso <- names(selection_freq[selection_freq >= 0.80])
print(freq_ordinata[freq_ordinata > 0.1])
features_di_ferro_lasso
####Boruta
set.seed(123)
df_boruta <- as.data.frame(X)
df_boruta$OS2ys_num <- Y
boruta_res <- Boruta(OS2ys_num ~ ., data = df_boruta, pValue = 0.05, mcAdj = TRUE, maxRuns = 2000, doTrace = 1)
print(boruta_res)
boruta_res$finalDecision
plot(boruta_res, las = 2, cex.axis = 0.7, main = "Boruta Feature Importance")
features_di_ferro_boruta <- getSelectedAttributes(boruta_res, withTentative = FALSE)
features_di_ferro_boruta
feature_super_definitive <- intersect(features_di_ferro_lasso, features_di_ferro_boruta)
print(feature_super_definitive)

###dose----
#### LASSO
set.seed(123)
cv_dose <- cv.glmnet(X_dose, Y_dose, family = "binomial", alpha = 1)
lasso_dose <- glmnet(X_dose, Y_dose, family = "binomial", alpha = 1, lambda = cv_dose$lambda.min)
print(coef(lasso_dose))
####bootrstrapping LASSO 
X <- X_dose
Y <- Y_dose
B <- 1000
n <- nrow(X)
selection_matrix <- matrix(0, nrow = B, ncol = ncol(X))
colnames(selection_matrix) <- colnames(X)
set.seed(123)
for (i in 1:B) {
  boot_idx <- sample(1:n, size = n, replace = TRUE)
  X_boot <- X[boot_idx, ]
  Y_boot <- Y[boot_idx]
  cv_fit <- cv.glmnet(X_boot, Y_boot, family = "binomial", alpha = 1)
  coef_boot <- as.matrix(coef(cv_fit, s = "lambda.min"))
  non_zero_features <- rownames(coef_boot)[coef_boot[, 1] != 0][-1]
  if(length(non_zero_features) > 0) {
    selection_matrix[i, non_zero_features] <- 1
  }
}
selection_freq <- colMeans(selection_matrix)
freq_ordinata <- sort(selection_freq, decreasing = TRUE)
features_di_ferro_lasso <- names(selection_freq[selection_freq >= 0.80])
print(freq_ordinata[freq_ordinata > 0.1])
features_di_ferro_lasso
####Boruta 
set.seed(123)
df_boruta <- as.data.frame(X)
df_boruta$OS2ys_num <- Y
boruta_res <- Boruta(OS2ys_num ~ ., data = df_boruta, pValue = 0.05, mcAdj = TRUE, maxRuns = 2000, doTrace = 1) #tried also 2000 but EDIC remains tentative
print(boruta_res)
boruta_res$finalDecision
plot(boruta_res, las = 2, cex.axis = 0.7, main = "Boruta Feature Importance")
features_di_ferro_boruta <- getSelectedAttributes(boruta_res, withTentative = FALSE)
features_di_ferro_boruta
feature_super_definitive <- intersect(features_di_ferro_lasso, features_di_ferro_boruta)
print(feature_super_definitive)

###metab----
#### LASSO
set.seed(123)
cv_metab <- cv.glmnet(X_metab, Y_metab, family = "binomial", alpha = 1)
lasso_metab <- glmnet(X_metab, Y_metab, family = "binomial", alpha = 1, lambda = cv_metab$lambda.min)
print(coef(lasso_metab))
####bootrstrapping LASSO 
X <- X_metab
Y <- Y_metab
B <- 1000
n <- nrow(X)
selection_matrix <- matrix(0, nrow = B, ncol = ncol(X))
colnames(selection_matrix) <- colnames(X)
set.seed(123)
for (i in 1:B) {
  boot_idx <- sample(1:n, size = n, replace = TRUE)
  X_boot <- X[boot_idx, ]
  Y_boot <- Y[boot_idx]
  cv_fit <- cv.glmnet(X_boot, Y_boot, family = "binomial", alpha = 1)
  coef_boot <- as.matrix(coef(cv_fit, s = "lambda.min"))
  non_zero_features <- rownames(coef_boot)[coef_boot[, 1] != 0][-1]
  if(length(non_zero_features) > 0) {
    selection_matrix[i, non_zero_features] <- 1
  }
}
selection_freq <- colMeans(selection_matrix)
freq_ordinata <- sort(selection_freq, decreasing = TRUE)
features_di_ferro_lasso <- names(selection_freq[selection_freq >= 0.80])
print(freq_ordinata[freq_ordinata > 0.1])
features_di_ferro_lasso
####Boruta 
set.seed(123)
df_boruta <- as.data.frame(X)
df_boruta$OS2ys_num <- Y
boruta_res <- Boruta(OS2ys_num ~ ., data = df_boruta, pValue = 0.05, mcAdj = TRUE, maxRuns = 2000, doTrace = 1)
print(boruta_res)
boruta_res$finalDecision
plot(boruta_res, las = 2, cex.axis = 0.7, main = "Boruta Feature Importance")
features_di_ferro_boruta <- getSelectedAttributes(boruta_res, withTentative = FALSE)
features_di_ferro_boruta
feature_super_definitive <- intersect(features_di_ferro_lasso, features_di_ferro_boruta)
print(feature_super_definitive)








#RF features importance -----
df_rf_glob<-as.data.frame(X_glob)
df_rf_glob$OS2ys_num<-Y_glob
set.seed(123)
rf_glob<-randomForest(OS2ys_num~.,data=df_rf_glob,importance=TRUE,ntree=1000)
imp_glob<-importance(rf_glob)
print(head(imp_glob[order(imp_glob[,"MeanDecreaseAccuracy"],decreasing=TRUE),],15))
varImpPlot(rf_glob,main="Random Forest - Variable Importance (Global)",pch=16,col="darkblue")
df_rf_crf<-as.data.frame(X_crf)
df_rf_crf$OS2ys_num<-Y_crf
set.seed(123)
rf_crf<-randomForest(OS2ys_num~.,data=df_rf_crf,importance=TRUE,ntree=1000)
imp_crf<-importance(rf_crf)
print(head(imp_crf[order(imp_crf[,"MeanDecreaseAccuracy"],decreasing=TRUE),],15))
varImpPlot(rf_crf,main="Random Forest - Variable Importance (CRF)",pch=16,col="darkblue")
df_rf_dose<-as.data.frame(X_dose)
df_rf_dose$OS2ys_num<-Y_dose
set.seed(123)
rf_dose<-randomForest(OS2ys_num~.,data=df_rf_dose,importance=TRUE,ntree=1000)
imp_dose<-importance(rf_dose)
print(head(imp_dose[order(imp_dose[,"MeanDecreaseAccuracy"],decreasing=TRUE),],15))
varImpPlot(rf_dose,main="Random Forest - Variable Importance (Dose)",pch=16,col="darkblue")
df_rf_metab<-as.data.frame(X_metab)
df_rf_metab$OS2ys_num<-Y_metab
set.seed(123)
rf_metab<-randomForest(OS2ys_num~.,data=df_rf_metab,importance=TRUE,ntree=1000)
imp_metab<-importance(rf_metab)
print(head(imp_metab[order(imp_metab[,"MeanDecreaseAccuracy"],decreasing=TRUE),],15))
varImpPlot(rf_metab,main="Random Forest - Variable Importance (Metab)",pch=16,col="darkblue")
















#dendogramma metab -----
datametab_clean <- na.omit(datametab)
datametab_clean = datametab_clean[,-c(9,10)]
cor_metab <- cor(datametab_clean, method = "spearman")
dist_metab <- as.dist(1 - abs(cor_metab))
hc_metab <- hclust(dist_metab, method = "average")
plot(hc_metab, main = "Dendrograms metabolic", xlab = "Feature", sub = "", ylab = "dist (1 - |corr(Spearman)|)", cex = 0.8, col = "darkblue")
#correlazione metab -----
var_cardiache <- datametab[, c("aorta_agatston", "aorta_volume_mm3", "heart_agatston", "heart_volume_mm3")]
var_cardiache_clean <- na.omit(var_cardiache)
matrice_correlazione <- cor(var_cardiache_clean, method = "spearman")
print(round(matrice_correlazione, 2))



#t-test/wilcoxon/chi-quadro/fisher - TABLE 
## CRF
##dose 
##metab 



## hospital ----
tab_hospital <- table(data$hospital, data$OS2ys_num)
print(tab_hospital)
test_chi <- chisq.test(tab_hospital)
print(test_chi)
test_fisher <- fisher.test(tab_hospital, simulate.p.value = TRUE, B = 1e5)
print(test_fisher)

#or
data$is_hosp3 <- ifelse(data$hospital == "3", "Ospedale_3", "Altri_Ospedali")
tab_hosp3 <- table(data$is_hosp3, data$OS2ys_num)
print(tab_hosp3)
print(chisq.test(tab_hosp3))
print(fisher.test(tab_hosp3))



#reimport-----
rm(list=ls())
graphics.off()
library(GGally)
library(corrplot)
library(ggplot2)
library(VIM)
library(glmnet)
library(cluster)
library(grpreg)
library(Boruta)
library(ggrepel)
library(randomForest)
library(lme4)
library(lattice)
library(caret)
library(pROC)
library(dcurves)
library(glmnet)
library(Boruta)
library(splines)

data <- read.csv('dataclean_final_withNA.csv', header=T)

#preprocessing 
sum(is.na(data))
colnames(data)
colnames(data)[colnames(data) == "Overall.Survival.at.2.years"] = "OS2ys"
colnames(data)[colnames(data) == "mean_dose_vanacava_S"] <- "mean_dose_venacava_S"
colnames(data)[colnames(data) == "eud_venacavaI"] <- "eud_venacava_I"
colnames(data)[colnames(data) == "eud_vanacavaS"] <- "eud_venacava_S"
data[data == ""] <- NA
data$KPS=factor(data$KPS)
data$PDL.1.level=factor(data$PDL.1.level)
data$schedule_ok_chemo=factor(data$schedule_ok_chemo)
data$cycles_chemo_sum = factor(data$cycles_chemo_sum)
data$cycles_immuno = factor(data$cycles_immuno)
data$prescription_dose_radio = factor(data$prescription_dose_radio)
data$number_of_fractions_radio = factor(data$number_of_fractions_radio)
data$number_of_fractions_radio = factor(data$number_of_fractions_radio)
data[] <- lapply(data, function(x) {
  if (!is.numeric(x)) as.factor(x) else x
})
data$OS2ys_num <- ifelse(data$OS2ys == "YES", 1, 0)
data$hospital <- factor(substr(as.character(data$ID), 3, 3), levels = as.character(1:7))
data$KPS <- as.numeric(as.character(data$KPS))
data$aortic_implant_present = as.factor(data$aortic_implant_present)
data$cardiac_implant_present = as.factor(data$cardiac_implant_present)

#put Volume PTV in CRF
cols <- names(data)
idx <- which(cols == "Stage.at.RT")
data <- data[, c(cols[1:idx], "Volume_PTV", setdiff(cols[(idx+1):length(cols)], "Volume_PTV"))]

datacrf = data[, c(1:41, 92)]
datadose = data[, 42:64]
datavol = data[, 65:77]
datametab = data[, 78:91]



#data preparation 
data$Smoking.status <- ifelse(is.na(data$Smoking.status), NA,ifelse(data$Smoking.status %in% c("YES", "EX"), 1, 0))
data$Diabetes.drugs <- ifelse(data$Diabetes.drugs == "YES", 1, 0)
data$Antilipidic.drugs <- ifelse(data$Antilipidic.drugs == "YES", 1, 0)
data$Hypertension <- ifelse(data$Hypertension == "YES", 1, 0)
data$antithrombotic <- ifelse(data$antithrombotic == "YES", 1, 0)
data$Beta.Blockers <- ifelse(data$Beta.Blockers == "YES", 1, 0)
data$Diuretics <- ifelse(data$Diuretics == "YES", 1, 0)
data$risk_sum <- rowSums(data[, c("Smoking.status", "Diabetes.drugs", "Antilipidic.drugs", "Hypertension", "antithrombotic", "Beta.Blockers", "Diuretics")], na.rm = FALSE) 
data$riskfactors <- cut(data$risk_sum, breaks = c(-Inf, 0, 2, 3, 5, Inf), labels = c(0, 1, 2, 3, 4))
data$riskfactors1 <- rowSums(data[, c("Smoking.status","Diabetes.drugs", "Antilipidic.drugs", "Hypertension")], na.rm = FALSE) 
data$risk_sum2 <- rowSums(data[, c("Beta.Blockers", "Diuretics")], na.rm = FALSE)
data$riskfactors2 = ifelse(data$risk_sum2>0,1,0)
data$riskfactors = as.factor(data$riskfactors)
data$riskfactors1 = as.factor(data$riskfactors1)
data$riskfactors2 = as.factor(data$riskfactors2)



##datasets 2
# hole
colnames(data)
data_global <- data[, -c(1,2,4,5,6:12,14,16,17,18,20:25,27:41,65:77,93:95,96,97)] #riskfactors1, riskfactors2
colnames(data_global)
data_global <- na.omit(data_global)
X_glob <- model.matrix(OS2ys_num ~ ., data_global)[, -1]

datacrf$Smoking.status <- ifelse(is.na(datacrf$Smoking.status), NA,ifelse(datacrf$Smoking.status %in% c("YES", "EX"), 1, 0))
datacrf$Diabetes.drugs <- ifelse(datacrf$Diabetes.drugs == "YES", 1, 0)
datacrf$Antilipidic.drugs <- ifelse(datacrf$Antilipidic.drugs == "YES", 1, 0)
datacrf$Hypertension <- ifelse(datacrf$Hypertension == "YES", 1, 0)
datacrf$antithrombotic <- ifelse(datacrf$antithrombotic == "YES", 1, 0)
datacrf$Beta.Blockers <- ifelse(datacrf$Beta.Blockers == "YES", 1, 0)
datacrf$Diuretics <- ifelse(datacrf$Diuretics == "YES", 1, 0)
datacrf$risk_sum <- rowSums(datacrf[, c("Smoking.status", "Diabetes.drugs", "Antilipidic.drugs", "Hypertension", "antithrombotic", "Beta.Blockers", "Diuretics")], na.rm = FALSE) 
datacrf$riskfactors <- cut(datacrf$risk_sum, breaks = c(-Inf, 0, 2, 3, 5, Inf), labels = c(0, 1, 2, 3, 4))
datacrf$riskfactors1 <- rowSums(datacrf[, c("Smoking.status", "Diabetes.drugs", "Antilipidic.drugs", "Hypertension")], na.rm = FALSE) 
datacrf$risk_sum2 <- rowSums(datacrf[, c("Beta.Blockers", "Diuretics")], na.rm = FALSE)
datacrf$riskfactors2 = ifelse(datacrf$risk_sum2>0,1,0)
datacrf$riskfactors = as.factor(datacrf$riskfactors)
datacrf$riskfactors1 = as.factor(datacrf$riskfactors1)
datacrf$riskfactors2 = as.factor(datacrf$riskfactors2)
colnames(datacrf)
Y_glob <- as.factor(data_global$OS2ys_num)

#CRF 
datacrf_model <- datacrf[, -c(1,2,4,5,6:12,14,16,17,18,20:25,27:41,43,44,45,46)]
colnames(datacrf_model)
datacrf_model <- na.omit(datacrf_model)
X_crf <- model.matrix(OS2ys_num ~ ., datacrf_model)[, -1]
Y_crf <- as.factor(datacrf_model$OS2ys_num)

#dose 
datadose_model <- na.omit(datadose)
X_dose <- model.matrix(data$OS2ys_num~ ., datadose_model)[, -1]
Y_dose <- as.factor(data$OS2ys_num)
#metab 
colnames(datametab)
datametab_model <- na.omit(datametab)
X_metab <-model.matrix(data$OS2ys_num~.,datametab)[,-1]
righe_valide <- as.numeric(rownames(X_metab))
ID_validi = datacrf$ID[righe_valide]
Y_metab <- as.factor(data$OS2ys_num[righe_valide])

dim(data_global)
dim(datacrf_model)
dim(datadose_model)
dim(datametab_model)








## hole2 ----
###LASSO
set.seed(123)
cv_glob <- cv.glmnet(X_glob, Y_glob, family = "binomial", alpha = 1)
lasso_glob <- glmnet(X_glob, Y_glob, family = "binomial", alpha = 1, lambda = cv_glob$lambda.min)
print(coef(lasso_glob))
### bootstrapping LASSO 
X <- X_glob
Y <- Y_glob
B <- 1000
n <- nrow(X)
selection_matrix <- matrix(0, nrow = B, ncol = ncol(X))
colnames(selection_matrix) <- colnames(X)
set.seed(123)
for (i in 1:B) {
  boot_idx <- sample(1:n, size = n, replace = TRUE)
  X_boot <- X[boot_idx, ]
  Y_boot <- Y[boot_idx]
  cv_fit <- cv.glmnet(X_boot, Y_boot, family = "binomial", alpha = 1)
  coef_boot <- as.matrix(coef(cv_fit, s = "lambda.min"))
  non_zero_features <- rownames(coef_boot)[coef_boot[, 1] != 0][-1]
  if(length(non_zero_features) > 0) {
    selection_matrix[i, non_zero_features] <- 1
  }
}
selection_freq <- colMeans(selection_matrix)
freq_ordinata <- sort(selection_freq, decreasing = TRUE)
features_di_ferro_lasso <- names(selection_freq[selection_freq >= 0.80])
print(freq_ordinata[freq_ordinata > 0.1])
features_di_ferro_lasso
### Boruta 
set.seed(123)
df_boruta <- as.data.frame(X)
df_boruta$OS2ys_num <- Y
boruta_res <- Boruta(OS2ys_num ~ ., data = df_boruta, pValue = 0.05, mcAdj = TRUE, maxRuns = 2000, doTrace = 1) #the only one that needed 2000 runs to converge, all the others converged with only 500 runs 
print(boruta_res)
boruta_res$finalDecision
plot(boruta_res, las = 2, cex.axis = 0.7, main = "Boruta Feature Importance")
features_di_ferro_boruta <- getSelectedAttributes(boruta_res, withTentative = FALSE)
features_di_ferro_boruta
feature_super_definitive <- intersect(features_di_ferro_lasso, features_di_ferro_boruta)
print(feature_super_definitive)



## specific2 ----
###CRF----
#### LASSO 
set.seed(123)
cv_crf <- cv.glmnet(X_crf, Y_crf, family = "binomial", alpha = 1)
lasso_crf <- glmnet(X_crf, Y_crf, family = "binomial", alpha = 1, lambda = cv_crf$lambda.min)
print(coef(lasso_crf))
####bootrstrapping LASSO 
X <- X_crf
Y <- Y_crf
B <- 1000
n <- nrow(X)
selection_matrix <- matrix(0, nrow = B, ncol = ncol(X))
colnames(selection_matrix) <- colnames(X)
set.seed(123)
for (i in 1:B) {
  boot_idx <- sample(1:n, size = n, replace = TRUE)
  X_boot <- X[boot_idx, ]
  Y_boot <- Y[boot_idx]
  cv_fit <- cv.glmnet(X_boot, Y_boot, family = "binomial", alpha = 1)
  coef_boot <- as.matrix(coef(cv_fit, s = "lambda.min"))
  non_zero_features <- rownames(coef_boot)[coef_boot[, 1] != 0][-1]
  if(length(non_zero_features) > 0) {
    selection_matrix[i, non_zero_features] <- 1
  }
}
selection_freq <- colMeans(selection_matrix)
freq_ordinata <- sort(selection_freq, decreasing = TRUE)
features_di_ferro_lasso <- names(selection_freq[selection_freq >= 0.80])
print(freq_ordinata[freq_ordinata > 0.1])
features_di_ferro_lasso
####Boruta
set.seed(123)
df_boruta <- as.data.frame(X)
df_boruta$OS2ys_num <- Y
boruta_res <- Boruta(OS2ys_num ~ ., data = df_boruta, pValue = 0.05, mcAdj = TRUE, maxRuns = 2000, doTrace = 1)
print(boruta_res)
boruta_res$finalDecision
plot(boruta_res, las = 2, cex.axis = 0.7, main = "Boruta Feature Importance")
features_di_ferro_boruta <- getSelectedAttributes(boruta_res, withTentative = FALSE)
features_di_ferro_boruta
feature_super_definitive <- intersect(features_di_ferro_lasso, features_di_ferro_boruta)
print(feature_super_definitive)

###dose----
#### LASSO
set.seed(123)
cv_dose <- cv.glmnet(X_dose, Y_dose, family = "binomial", alpha = 1)
lasso_dose <- glmnet(X_dose, Y_dose, family = "binomial", alpha = 1, lambda = cv_dose$lambda.min)
print(coef(lasso_dose))
####bootrstrapping LASSO 
X <- X_dose
Y <- Y_dose
B <- 1000
n <- nrow(X)
selection_matrix <- matrix(0, nrow = B, ncol = ncol(X))
colnames(selection_matrix) <- colnames(X)
set.seed(123)
for (i in 1:B) {
  boot_idx <- sample(1:n, size = n, replace = TRUE)
  X_boot <- X[boot_idx, ]
  Y_boot <- Y[boot_idx]
  cv_fit <- cv.glmnet(X_boot, Y_boot, family = "binomial", alpha = 1)
  coef_boot <- as.matrix(coef(cv_fit, s = "lambda.min"))
  non_zero_features <- rownames(coef_boot)[coef_boot[, 1] != 0][-1]
  if(length(non_zero_features) > 0) {
    selection_matrix[i, non_zero_features] <- 1
  }
}
selection_freq <- colMeans(selection_matrix)
freq_ordinata <- sort(selection_freq, decreasing = TRUE)
features_di_ferro_lasso <- names(selection_freq[selection_freq >= 0.80])
print(freq_ordinata[freq_ordinata > 0.1])
features_di_ferro_lasso
####Boruta 
set.seed(123)
df_boruta <- as.data.frame(X)
df_boruta$OS2ys_num <- Y
boruta_res <- Boruta(OS2ys_num ~ ., data = df_boruta, pValue = 0.05, mcAdj = TRUE, maxRuns = 2000, doTrace = 1)
print(boruta_res)
boruta_res$finalDecision
plot(boruta_res, las = 2, cex.axis = 0.7, main = "Boruta Feature Importance")
features_di_ferro_boruta <- getSelectedAttributes(boruta_res, withTentative = FALSE)
features_di_ferro_boruta
feature_super_definitive <- intersect(features_di_ferro_lasso, features_di_ferro_boruta)
print(feature_super_definitive)

###metab----
#### LASSO
set.seed(123)
cv_metab <- cv.glmnet(X_metab, Y_metab, family = "binomial", alpha = 1)
lasso_metab <- glmnet(X_metab, Y_metab, family = "binomial", alpha = 1, lambda = cv_metab$lambda.min)
print(coef(lasso_metab))
####bootrstrapping LASSO 
X <- X_metab
Y <- Y_metab
B <- 1000
n <- nrow(X)
selection_matrix <- matrix(0, nrow = B, ncol = ncol(X))
colnames(selection_matrix) <- colnames(X)
set.seed(123)
for (i in 1:B) {
  boot_idx <- sample(1:n, size = n, replace = TRUE)
  X_boot <- X[boot_idx, ]
  Y_boot <- Y[boot_idx]
  cv_fit <- cv.glmnet(X_boot, Y_boot, family = "binomial", alpha = 1)
  coef_boot <- as.matrix(coef(cv_fit, s = "lambda.min"))
  non_zero_features <- rownames(coef_boot)[coef_boot[, 1] != 0][-1]
  if(length(non_zero_features) > 0) {
    selection_matrix[i, non_zero_features] <- 1
  }
}
selection_freq <- colMeans(selection_matrix)
freq_ordinata <- sort(selection_freq, decreasing = TRUE)
features_di_ferro_lasso <- names(selection_freq[selection_freq >= 0.80])
print(freq_ordinata[freq_ordinata > 0.1])
features_di_ferro_lasso
####Boruta 
set.seed(123)
df_boruta <- as.data.frame(X)
df_boruta$OS2ys_num <- Y
boruta_res <- Boruta(OS2ys_num ~ ., data = df_boruta, pValue = 0.05, mcAdj = TRUE, maxRuns = 2000, doTrace = 1)
print(boruta_res)
boruta_res$finalDecision
plot(boruta_res, las = 2, cex.axis = 0.7, main = "Boruta Feature Importance")
features_di_ferro_boruta <- getSelectedAttributes(boruta_res, withTentative = FALSE)
features_di_ferro_boruta
feature_super_definitive <- intersect(features_di_ferro_lasso, features_di_ferro_boruta)
print(feature_super_definitive)







