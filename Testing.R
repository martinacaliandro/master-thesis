rm(list=ls())
#graphics.off()
library(randomForest)
library(lme4)
library(lattice)
library(caret)
library(pROC)
library(dcurves)
library(splines)
library(xgboost)
library(lightgbm)
library(grid)
library(survival)
library(survminer)
library(uwot)
library(dbscan)
library(mclust)
library(treeshap)
library(rpart)
library(readxl)

data <- read.csv('dataclean_final_withNA.csv', header=T)
datatest =  read_excel("validazione_martina.xlsx", sheet = 1)
datatestos = read_excel("validazione_martina.xlsx", sheet = 3)


# preprocessing data train -----
sum(is.na(data))
colnames(data)[colnames(data) == "Overall.Survival.at.2.years"] = "OS2ys"
colnames(data)[colnames(data) == "mean_dose_vanacava_S"] <- "mean_dose_venacava_S"
colnames(data)[colnames(data) == "eud_venacavaI"] <- "eud_venacava_I"
colnames(data)[colnames(data) == "eud_vanacavaS"] <- "eud_venacava_S"
data[data == ""] <- NA

data$KPS = factor(data$KPS)
data$PDL.1.level = factor(data$PDL.1.level)
data$schedule_ok_chemo = factor(data$schedule_ok_chemo)
data$cycles_chemo_sum = factor(data$cycles_chemo_sum)
data$cycles_immuno = factor(data$cycles_immuno)
data$prescription_dose_radio = factor(data$prescription_dose_radio)
data$number_of_fractions_radio = factor(data$number_of_fractions_radio)

data[] <- lapply(data, function(x) {
  if (!is.numeric(x)) as.factor(x) else x
})
data$OS2ys_num <- ifelse(data$OS2ys == "YES", 1, 0)
data$hospital <- factor(substr(as.character(data$ID), 3, 3), levels = as.character(1:7))
data$aortic_implant_present = as.factor(data$aortic_implant_present)
data$cardiac_implant_present = as.factor(data$cardiac_implant_present)

cols <- names(data)
idx <- which(cols == "Stage.at.RT")
data <- data[, c(cols[1:idx], "Volume_PTV", setdiff(cols[(idx+1):length(cols)], "Volume_PTV"))]

data$Smoking.status <- ifelse(is.na(data$Smoking.status), NA, ifelse(data$Smoking.status %in% c("YES", "EX"), 1, 0))
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
data$riskfactors3 <- rowSums(data[, c("Diabetes.drugs", "Antilipidic.drugs", "Hypertension")], na.rm = FALSE) 
data$riskfactors2 = ifelse(data$risk_sum2 > 0, 1, 0)

data$riskfactors = as.factor(data$riskfactors)
data$riskfactors1 = as.factor(data$riskfactors1)
data$riskfactors2 = as.factor(data$riskfactors2)
data$riskfactors3 = as.factor(data$riskfactors3)
data$Smoking.status = as.factor(data$Smoking.status)
data$Diabetes.drugs = as.factor(data$Diabetes.drugs)
data$Antilipidic.drugs = factor(data$Antilipidic.drugs)
data$Hypertension = as.factor(data$Hypertension)
data$antithrombotic = as.factor(data$antithrombotic)
data$Beta.Blockers = as.factor(data$Beta.Blockers)
data$Diuretics = as.factor(data$Diuretics)
data$Volume_PTV_log = log1p(data$Volume_PTV) 
data$LAA_log = log1p(data$LAA) 
data$aorta_agatston_log = log1p(data$aorta_agatston) 
data$eud_lungs_log = log1p(data$eud_lungs) 


#short analysis test -----
colnames(datatest)
dim(datatest)
colnames(datatestos)
dim(datatestos)
datatest$risk_factor1=as.factor(datatest$risk_factor1)
summary(datatest)
summary(data$Age)
summary(data$Volume_PTV)
summary(data$LAA)
summary(data$riskfactors1)
table(datatestos$OS2Y_2)
table(datatestos$OS2Y_1)
table(datatestos$OS2Y_1)/249
table(data$OS2ys)
table(data$OS2ys)/456
hist(datatest$Age, breaks=30)
hist(data$Age, breaks=30)
hist(datatest$VOLUME_PTV, breaks=30)
hist(data$Volume_PTV, breaks=30)
hist(datatest$`LAA%`, breaks=30)
hist(data$LAA, breaks=30)
counts <- table(data$riskfactors1)
freq <- prop.table(counts) * 100
bp <- barplot(counts,
              col = "steelblue",
              xlab = "Risk factor 1",
              ylab = "Count",
              ylim = c(0, max(counts) * 1.2))
text(bp, counts,
     labels = paste0(counts, "\n",
                     round(freq, 1), "%"),
     pos = 3, cex = 0.9)
counts <- table(datatest$risk_factor1)
freq <- prop.table(counts) * 100
bp <- barplot(counts,
              col = "steelblue",
              xlab = "Risk factor 1",
              ylab = "Count",
              ylim = c(0, max(counts) * 1.2))
text(bp, counts,
     labels = paste0(counts, "\n",
                     round(freq, 1), "%"),
     pos = 3, cex = 0.9)



#function -----
calc_metrics_vec <- function(p, Y) {
  if(is.list(p)) p <- unlist(p)
  p <- as.numeric(p)
  valid_idx <- !is.na(p)
  metric_names <- c("AUC", "F1", "F1_0", "F1_macro", "Precision", "Precision_0", "Recall", "Recall_0", "Brier", "LogLoss", "Cal_Intercept", "Cal_Slope", "Threshold")
  if(sum(valid_idx) < 10 || length(unique(p[valid_idx])) < 2) {
    res_na <- rep(NA_real_, length(metric_names))
    names(res_na) <- metric_names
    return(res_na)
  }
  p_valid <- p[valid_idx]
  Y_valid <- Y[valid_idx]
  Y_num <- as.numeric(as.character(Y_valid))
  Y_factor <- factor(Y_valid, levels = c(0, 1))
  r <- pROC::roc(Y_valid, p_valid, quiet = TRUE)
  t <- 0.5
  pred_factor <- factor(ifelse(p_valid >= t, 1, 0), levels = c(0, 1))
  cm1 <- caret::confusionMatrix(pred_factor, Y_factor, positive = "1")
  cm0 <- caret::confusionMatrix(pred_factor, Y_factor, positive = "0")
  prec1 <- as.numeric(cm1$byClass["Pos Pred Value"])
  rec1  <- as.numeric(cm1$byClass["Sensitivity"])
  f1_1  <- 2 * (prec1 * rec1) / (prec1 + rec1)
  prec0 <- as.numeric(cm0$byClass["Pos Pred Value"])
  rec0  <- as.numeric(cm0$byClass["Sensitivity"])
  f1_0  <- 2 * (prec0 * rec0) / (prec0 + rec0)
  f1_macro <- mean(c(f1_1, f1_0), na.rm = TRUE)
  eps <- 1e-15
  p_adj <- pmax(pmin(p_valid, 1 - eps), eps)
  brier <- mean((p_valid - Y_num)^2)
  logloss <- -mean(Y_num * log(p_adj) + (1 - Y_num) * log(1 - p_adj))
  logit_p <- log(p_adj / (1 - p_adj))
  c_slope <- NA_real_
  c_int   <- NA_real_
  try({
    fit_calib <- glm(Y_num ~ logit_p, family = binomial())
    c_slope <- as.numeric(coef(fit_calib)[2])
    fit_calib_int <- glm(Y_num ~ 1 + offset(logit_p), family = binomial())
    c_int <- as.numeric(coef(fit_calib_int)[1])
  }, silent = TRUE)
  c(
    AUC           = as.numeric(pROC::auc(r)),
    F1            = f1_1,
    F1_0          = f1_0,
    F1_macro      = f1_macro,
    Precision     = prec1,
    Precision_0   = prec0,
    Recall        = rec1,
    Recall_0      = rec0,
    Brier         = brier,
    LogLoss       = logloss,
    Cal_Intercept = c_int,
    Cal_Slope     = c_slope,
    Threshold     = t
  )
}

#train XGBoost -----
hp_xgb     <- list(max_depth = 2, eta = 0.01, lambda = 15, subsample = 0.7, nrounds = 400)
covariates <- c("Age", "Volume_PTV", "LAA", "riskfactors1")
cont_vars  <- c("Age", "Volume_PTV", "LAA")


minPts_hdb <- 25
nn_umap    <- 30
seed_umap  <- 123 
seed_xgb   <- 123



idx <- complete.cases(data[, covariates]) & !is.na(data$OS2ys_num) & !is.na(data$hospital)
df_sub <- data[idx, ]
X_full <- df_sub[, covariates, drop = FALSE]
Y_full <- df_sub$OS2ys_num
Z_full <- droplevels(df_sub$hospital)
sc_train_params <- scale(X_full[, cont_vars])
center_train <- attr(sc_train_params, "scaled:center")
scale_train  <- attr(sc_train_params, "scaled:scale")

X_full_sc <- X_full
X_full_sc[, cont_vars] <- sc_train_params

cat("\npatients included:", nrow(X_full_sc), "| OS2ys=1 events:", sum(Y_full), "\n")

aggregate_shap <- function(shap_matrix, original_df, original_covs) {
  res_list <- list()
  for (v in original_covs) {
    if (is.numeric(original_df[[v]])) {
      res_list[[v]] <- shap_matrix[, v]
    } else {
      dummy_cols <- grep(paste0("^", v), colnames(shap_matrix), value = TRUE)
      if (length(dummy_cols) == 1) {
        res_list[[v]] <- shap_matrix[, dummy_cols]
      } else if (length(dummy_cols) > 1) {
        res_list[[v]] <- rowSums(shap_matrix[, dummy_cols, drop = FALSE])
      } else {
        res_list[[v]] <- rep(0, nrow(shap_matrix))
      }
    }
  }
  as.matrix(as.data.frame(res_list))
}

X_mat <- model.matrix(~ ., data = X_full_sc)[, -1, drop = FALSE]


set.seed(seed_xgb)
fit_xgb <- xgb.train(
  params  = list(objective = "binary:logistic", 
                 max_depth = hp_xgb$max_depth, 
                 eta = hp_xgb$eta, 
                 lambda = hp_xgb$lambda, 
                 subsample = hp_xgb$subsample),
  data    = xgb.DMatrix(X_mat, label = Y_full),
  nrounds = hp_xgb$nrounds,
  verbose = 0
)

uni_xgb <- xgboost.unify(fit_xgb, X_mat)
shap_train_raw <- as.matrix(treeshap(uni_xgb, X_mat, verbose = FALSE)$shaps)
shap_train_agg <- aggregate_shap(shap_train_raw, X_full, covariates)
set.seed(seed_umap)
um2_shap <- uwot::umap(
  shap_train_agg, 
  n_neighbors = nn_umap, 
  min_dist = 0, 
  n_components = 2, 
  n_threads = 1, 
  n_sgd_threads = 1
)

set.seed(seed_umap) #same UMAP of the original one
umap_model_clB <- uwot::umap(
  shap_train_agg, 
  n_neighbors = nn_umap, 
  min_dist = 0, 
  n_components = 2, 
  n_threads = 1, 
  n_sgd_threads = 1,
  init = um2_shap,     
  n_epochs = 0,        
  ret_model = TRUE     
)

cl_B_train <- dbscan::hdbscan(umap_model_clB$embedding, minPts = minPts_hdb)
cat("\nclusters found in the training set:\n")
print(table(Cluster_Train = cl_B_train$cluster))


# test ----
df_test <- as.data.frame(datatest)

colnames(df_test)[colnames(df_test) == "VOLUME_PTV"]   <- "Volume_PTV"
colnames(df_test)[colnames(df_test) == "LAA%"]         <- "LAA"
colnames(df_test)[colnames(df_test) == "risk_factor1"] <- "riskfactors1"

df_test$Age        <- as.numeric(df_test$Age)
df_test$Volume_PTV <- as.numeric(df_test$Volume_PTV)
df_test$LAA        <- as.numeric(df_test$LAA)
df_test$riskfactors1 <- factor(as.character(df_test$riskfactors1), levels = levels(X_full$riskfactors1))

X_test <- df_test[, covariates, drop = FALSE]

X_test_sc <- X_test
X_test_sc[, cont_vars] <- scale(X_test[, cont_vars], center = center_train, scale = scale_train)

X_test_mat <- model.matrix(~ ., data = X_test_sc)[, -1, drop = FALSE]

p_test <- as.numeric(predict(fit_xgb, xgb.DMatrix(X_test_mat)))
pred_class_test <- ifelse(p_test >= 0.5, 1, 0)

cat("xgboost test set preparation")
cat("tested patients:", length(p_test), "\n")
cat("predicted class distribution (0.5 threshold):\n")
print(table(Pred_Class = pred_class_test))

treeshap_test <- treeshap(uni_xgb, X_test_mat, verbose = FALSE)
shap_test_raw <- as.matrix(treeshap_test$shaps)
shap_test_agg <- aggregate_shap(shap_test_raw, X_test, covariates)

um2_test <- uwot::umap_transform(shap_test_agg, umap_model_clB)

cl_test_pred <- as.integer(predict(cl_B_train, data = umap_model_clB$embedding, newdata = um2_test))

cat("cluster assignment on the test set")
print(table(Cluster_Test = cl_test_pred))

test_results <- data.frame(
  ID             = df_test$ID,
  Prob_OS2ys     = round(p_test, 4),
  Pred_Class     = pred_class_test,
  Cluster_clB    = cl_test_pred,
  UMAP1          = round(um2_test[, 1], 3),
  UMAP2          = round(um2_test[, 2], 3)
)

print(head(test_results, 10))



















# test metrics ------
datatestos <- read_excel("validazione_martina.xlsx", sheet = 3)

datatestos_aligned <- datatestos[match(df_test$ID, datatestos$ID), ]

#Y_test <- ifelse(datatestos_aligned$OS2Y_2 == 0, 1, 0)
Y_test <- ifelse(datatestos_aligned$OS2Y_1 == "Yes", 1, 0)

valid_outcomes_idx <- !is.na(Y_test)
cat("\npatients with available os2ys outcome in test set:", sum(valid_outcomes_idx), "\n")
cat("\n--- CONFUSION MATRIX: PREDICTION VS GROUND TRUTH (threshold = 0.5) ---\n")
cm_table <- table(Predicted = pred_class_test[valid_outcomes_idx], Actual = Y_test[valid_outcomes_idx])
print(cm_table)

if (sum(valid_outcomes_idx) > 0) {
  
  metrics_test <- calc_metrics_vec(
    p = p_test[valid_outcomes_idx], 
    Y = Y_test[valid_outcomes_idx]
  )
  
  cat("\n--- final test set performance metrics ---\n")
  print(round(metrics_test, 4))
  
  cat("\n--- cluster assignments vs real outcome ---\n")
  
  cl_outcome_table <- table(
    Cluster = cl_test_pred[valid_outcomes_idx], 
    Actual_OS2ys = Y_test[valid_outcomes_idx]
  )
  print(cl_outcome_table)
  
  cat("\nobserved survival rate per cluster in the test set:\n")
  cluster_survival_rates <- prop.table(cl_outcome_table, margin = 1)[, "1"]
  print(round(cluster_survival_rates * 100, 1))
  
} else {
  cat("\ncannot calculate metrics: os2y_1 column is completely empty/na.\n")
}


df_prof_train <- X_full
df_prof_train$Cluster <- cl_B_train$cluster
df_prof_train$Set <- "train"

df_prof_test <- X_test
df_prof_test$Cluster <- cl_test_pred
df_prof_test$Set <- "test"

df_prof_all <- rbind(df_prof_train, df_prof_test)
df_prof_valid <- subset(df_prof_all, Cluster > 0)
df_prof_valid$Cluster <- factor(df_prof_valid$Cluster)

if (length(cont_vars) > 0) {
  prof_cont <- aggregate(df_prof_valid[, cont_vars, drop = FALSE], 
                         by = list(Cluster = df_prof_valid$Cluster, Set = df_prof_valid$Set), 
                         FUN = function(x) round(mean(x, na.rm = TRUE), 2))
  prof_cont <- prof_cont[order(prof_cont$Cluster, prof_cont$Set), ]
  print(prof_cont)
}

cat_vars <- setdiff(covariates, cont_vars)
if (length(cat_vars) > 0) {
  for (v in cat_vars) {
    tab <- table(Cluster = df_prof_valid$Cluster, Set = df_prof_valid$Set, Value = df_prof_valid[[v]])
    prop_tab <- prop.table(tab, margin = c(1, 2)) * 100
    print(round(prop_tab, 1))
  }
}

df_melt <- reshape(
  df_prof_valid, 
  direction = "long", 
  varying = cont_vars, 
  v.names = "Value", 
  timevar = "Feature", 
  times = cont_vars
)


cl_test_rules <- rep(0, nrow(df_test))
for(i in seq_len(nrow(df_test))) {
  a <- df_test$Age[i]
  l <- df_test$LAA[i]
  p <- df_test$Volume_PTV[i]
  
  if(a > 71.5 & l <= 0.908 & p <= 288.16) { cl_test_rules[i] <- 1 }
  else if(a <= 72.5 & p <= 289.37) { cl_test_rules[i] <- 2 }
  else if(a <= 75.5 & l <= 1.511 & p > 289.86) { cl_test_rules[i] <- 3 }
  else if(l > 1.511 & p > 289.59) { cl_test_rules[i] <- 4 }
  else if(a > 75.5 & l <= 1.576 & p > 299.42) { cl_test_rules[i] <- 5 }
}

test_results$Cluster_Rules <- cl_test_rules












df_plot_train <- data.frame(
  UMAP1   = umap_model_clB$embedding[, 1],
  UMAP2   = umap_model_clB$embedding[, 2],
  Cluster = factor(ifelse(cl_B_train$cluster == 0, "noise (0)", paste("cluster", cl_B_train$cluster))),
  Set     = "train"
)

df_plot_test_hdb <- data.frame(
  UMAP1   = um2_test[, 1],
  UMAP2   = um2_test[, 2],
  Cluster = factor(ifelse(cl_test_pred == 0, "noise (0)", paste("cluster", cl_test_pred))),
  Set     = "test"
)

df_plot_test_rules <- data.frame(
  UMAP1   = um2_test[, 1],
  UMAP2   = um2_test[, 2],
  Cluster = factor(ifelse(cl_test_rules == 0, "noise (0)", paste("cluster", cl_test_rules))),
  Set     = "test"
)

df_plot_all <- rbind(df_plot_train, df_plot_test_hdb)

p_train_only <- ggplot(df_plot_train, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(alpha = 0.85, size = 2) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom") +
  labs(
    title = "training set (HDBSCAN)",
    x = "UMAP 1 (SHAP)", y = "UMAP 2 (SHAP)"
  )

p_test_hdb <- ggplot(df_plot_test_hdb, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(alpha = 0.95, size = 2.4, shape = 17) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom") +
  labs(
    title = "test set (HDBSCAN)",
    x = "UMAP 1 (SHAP)", y = "UMAP 2 (SHAP)"
  )

p_test_rules <- ggplot(df_plot_test_rules, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(alpha = 0.95, size = 2.4, shape = 17) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom") +
  labs(
    title = "test set (rules)",
    x = "UMAP 1 (SHAP)", y = "UMAP 2 (SHAP)"
  )

p_both <- ggplot(df_plot_all, aes(x = UMAP1, y = UMAP2, color = Cluster, shape = Set)) +
  geom_point(data = subset(df_plot_all, Set == "train"), alpha = 0.35, size = 1.6) +
  geom_point(data = subset(df_plot_all, Set == "test"), alpha = 0.95, size = 2.4) +
  scale_shape_manual(values = c("train" = 16, "test" = 17)) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom") +
  labs(
    title = "combined projection",
    subtitle = "triangles: test, circles: train",
    x = "UMAP 1 (SHAP)", y = "UMAP 2 (SHAP)"
  )

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 4)))
print(p_train_only,  vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_test_hdb,    vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p_test_rules,  vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
print(p_both,        vp = viewport(layout.pos.row = 1, layout.pos.col = 4))
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 3)))
print(p_train_only,  vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_test_hdb,    vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p_test_rules,  vp = viewport(layout.pos.row = 1, layout.pos.col = 3))



make_profile_boxplot <- function(test_cluster_vec, method_name) {
  df_tr <- X_full
  df_tr$Cluster <- cl_B_train$cluster
  df_tr$Set <- "train"
  
  df_te <- X_test
  df_te$Cluster <- test_cluster_vec
  df_te$Set <- paste0("test (", method_name, ")")
  
  df_all <- rbind(df_tr, df_te)
  df_valid <- subset(df_all, Cluster > 0 & !is.na(Cluster))
  df_valid$Cluster <- factor(df_valid$Cluster)
  
  df_m <- reshape(
    df_valid,
    direction = "long",
    varying   = cont_vars,
    v.names   = "Value",
    timevar   = "Feature",
    times     = cont_vars
  )
  
  col_test <- ifelse(method_name == "HDBSCAN", "#D55E00", "#009E73")
  fill_vals <- setNames(c("#0072B2", col_test), c("train", paste0("test (", method_name, ")")))
  
  p <- ggplot(df_m, aes(x = Cluster, y = Value, fill = Set)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.8, position = position_dodge(0.8)) +
    facet_wrap(~ Feature, scales = "free_y") +
    scale_fill_manual(values = fill_vals) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold")) +
    labs(
      title = paste("Clinical profile comparison: train vs test -", method_name),
      x = "Cluster",
      y = "Value",
      fill = "Cohort:"
    )
  return(p)
}
p_prof_hdb   <- make_profile_boxplot(cl_test_pred,  "HDBSCAN")
p_prof_rules <- make_profile_boxplot(cl_test_rules, "Rules")

print(p_prof_hdb)
print(p_prof_rules)

make_cat_profile_plot <- function(test_cluster_vec, method_name) {
  if (length(cat_vars) == 0) return(NULL)
  
  df_tr <- X_full[, cat_vars, drop = FALSE]
  df_tr$Cluster <- cl_B_train$cluster
  df_tr$Set <- "train"
  
  df_te <- X_test[, cat_vars, drop = FALSE]
  df_te$Cluster <- test_cluster_vec
  df_te$Set <- paste0("test (", method_name, ")")
  
  df_all <- rbind(df_tr, df_te)
  df_valid <- subset(df_all, Cluster > 0 & !is.na(Cluster))
  df_valid$Cluster <- factor(df_valid$Cluster)
  
  res_list <- list()
  for (v in cat_vars) {
    tab <- as.data.frame(table(Cluster = df_valid$Cluster, Set = df_valid$Set, Category = df_valid[[v]]))
    tot <- tapply(tab$Freq, list(tab$Cluster, tab$Set), sum)
    tab$Tot <- apply(tab, 1, function(r) tot[as.character(r["Cluster"]), as.character(r["Set"])])
    tab$Perc <- ifelse(tab$Tot > 0, (tab$Freq / tab$Tot) * 100, 0)
    tab$Label <- ifelse(tab$Tot > 0, sprintf("%.1f%%\n(n=%d)", tab$Perc, tab$Freq), "")
    tab$Feature <- v
    res_list[[v]] <- tab
  }
  df_plot_cat <- do.call(rbind, res_list)
  
  col_test <- ifelse(method_name == "HDBSCAN", "#D55E00", "#009E73")
  fill_vals <- setNames(c("#0072B2", col_test), c("train", paste0("test (", method_name, ")")))
  
  p <- ggplot(df_plot_cat, aes(x = Category, y = Perc, fill = Set)) +
    geom_col(position = position_dodge(0.85), width = 0.72, color = "white", linewidth = 0.3) +
    geom_text(
      aes(label = Label),
      position = position_dodge(0.85),
      vjust = -0.25,
      size = 3.0,
      fontface = "bold",
      lineheight = 0.9
    ) +
    facet_grid(Feature ~ Cluster, scales = "free_x", labeller = labeller(Cluster = label_both)) +
    scale_fill_manual(values = fill_vals) +
    scale_y_continuous(limits = c(0, 118), breaks = seq(0, 100, 25)) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      title = paste("Categorical profile: train vs test -", method_name),
      x = "Category / Level",
      y = "Percentage within Cluster (%)",
      fill = "Cohort:"
    )
  return(p)
}

p_cat_hdb   <- make_cat_profile_plot(cl_test_pred,  "HDBSCAN")
p_cat_rules <- make_cat_profile_plot(cl_test_rules, "Rules")

if (!is.null(p_cat_hdb))   print(p_cat_hdb)
if (!is.null(p_cat_rules)) print(p_cat_rules)


















os_summary_train <- aggregate(Y_full, by = list(Cluster = cl_B_train$cluster), FUN = function(x) sprintf("n=%d (%.1f%%)", length(x), mean(x)*100))
os_summary_hdb <- aggregate(Y_test[valid_outcomes_idx], by = list(Cluster = cl_test_pred[valid_outcomes_idx]), FUN = function(x) sprintf("n=%d (%.1f%%)", length(x), mean(x)*100))
os_summary_rules <- aggregate(Y_test[valid_outcomes_idx], by = list(Cluster = cl_test_rules[valid_outcomes_idx]), FUN = function(x) sprintf("n=%d (%.1f%%)", length(x), mean(x)*100))

comp_os <- merge(os_summary_train, os_summary_hdb, by = "Cluster", all = TRUE)
comp_os <- merge(comp_os, os_summary_rules, by = "Cluster", all = TRUE)
colnames(comp_os) <- c("cluster", "train_os2ys", "test_hdbscan_os2ys", "test_rules_os2ys")
print(comp_os)

plot_hm_med <- function(df_d, cl_l, t_text) {
  i_v <- cl_l > 0
  df_d <- df_d[i_v, , drop=FALSE]
  cl_l <- factor(cl_l[i_v])
  os_s <- aggregate(df_d$Y_target, by = list(Cluster = cl_l), FUN = function(x) round(mean(x, na.rm = TRUE) * 100, 1))
  os_m <- setNames(os_s$x, os_s$Cluster)
  n_v <- intersect(cont_vars, colnames(df_d))
  r_n <- df_d[, n_v, drop = FALSE]
  r_s <- function(x) (x - median(x, na.rm = TRUE)) / (sd(x, na.rm = TRUE) + 1e-6)
  r_med <- aggregate(r_n, by = list(Cluster = cl_l), FUN = median)
  z_med <- aggregate(as.data.frame(lapply(r_n, r_s)), by = list(Cluster = cl_l), FUN = median)
  hm_m_l <- reshape(r_med, direction = "long", varying = n_v, v.names = "raw_median", timevar = "feature", times = n_v)
  hm_z_l <- reshape(z_med, direction = "long", varying = n_v, v.names = "z_robust", timevar = "feature", times = n_v)
  hm_df <- merge(hm_m_l[, c("Cluster", "feature", "raw_median")], hm_z_l[, c("Cluster", "feature", "z_robust")], by = c("Cluster", "feature"))
  hm_df$cell_label <- mapply(function(c, f) {
    v_c <- df_d[[f]][cl_l == c]
    sprintf("%.2f\n[%.2f - %.2f]", median(v_c, na.rm=TRUE), quantile(v_c, 0.25, na.rm=TRUE), quantile(v_c, 0.75, na.rm=TRUE))
  }, hm_df$Cluster, hm_df$feature)
  hm_df$feature <- factor(hm_df$feature, levels = n_v)
  hm_df$cluster_label <- paste0("cluster ", hm_df$Cluster, "\n(os: ", os_m[as.character(hm_df$Cluster)], "%)")
  
  hm_df$z_robust <- pmax(pmin(hm_df$z_robust, 2), -2)
  
  p <- ggplot(hm_df, aes(x = feature, y = cluster_label, fill = z_robust)) +
    geom_tile(color = "white", linewidth = 1.2) +
    geom_text(aes(label = cell_label), size = 3.5, fontface = "bold") +
    scale_fill_gradient2(low = "#0072b2", mid = "white", high = "#d55e00", midpoint = 0, limits = c(-2, 2), na.value = "grey94", name = "z-score") +
    theme_minimal(base_size = 12) + labs(title = t_text, x = "variable", y = "cluster")
  return(p)
}

df_tr_hm <- X_full; df_tr_hm$Y_target <- Y_full
p_hm_tr <- plot_hm_med(df_tr_hm, cl_B_train$cluster, "train set median profile")

df_te_hm <- X_test[valid_outcomes_idx, , drop=FALSE]; df_te_hm$Y_target <- Y_test[valid_outcomes_idx]
p_hm_te_hdb <- plot_hm_med(df_te_hm, cl_test_pred[valid_outcomes_idx], "test set hdbscan median profile")
p_hm_te_ru <- plot_hm_med(df_te_hm, cl_test_rules[valid_outcomes_idx], "test set rules median profile")

print(p_hm_tr)
print(p_hm_te_hdb)
print(p_hm_te_ru)

df_sv <- datatestos_aligned[valid_outcomes_idx, ]
y_t_val <- Y_test[valid_outcomes_idx]

print(table(OS_under_24_months = df_sv$OS < 24, Survived_2y_Yes = y_t_val))



data[data$OS2ys_num == 0 & data$OS >= 24, c("ID","OS2ys", "OS2ys_num","OS")]
data[data$OS2ys_num == 1 & data$OS <= 24, c("ID","OS2ys", "OS2ys_num","OS")]
df_sv[y_t_val == 0 & df_sv$OS >= 24,]
df_sv[y_t_val == 1 & df_sv$OS <= 24,]
df_sv$OS <- ifelse(y_t_val == 1 & df_sv$OS < 24, 24, df_sv$OS)
df_sv$os_24m <- pmin(df_sv$OS, 24)
df_sv$ev_24m <- 1 - y_t_val #censurati tutti i pazienti che sopravvivono a 24 mesi. Non censurati quelli deceduti entro i 24 mesi

plot_km <- function(sd, cl_v, t_text) {
  sd$cl <- cl_v
  sd <- sd[sd$cl > 0, ]
  fk <- survfit(Surv(os_24m, ev_24m) ~ cl, data = sd)
  lr <- survdiff(Surv(os_24m, ev_24m) ~ cl, data = sd)
  p_l <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
  kd <- data.frame(time = fk$time, surv = fk$surv * 100, cl = rep(gsub("cl=", "cluster ", names(fk$strata)), fk$strata))
  kd <- rbind(data.frame(time = 0, surv = 100, cl = unique(kd$cl)), kd)
  p <- ggplot(kd, aes(x = time, y = surv, color = cl)) +
    geom_step(linewidth = 1.2) + scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) + scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 3)) +
    theme_minimal(base_size = 12) + theme(legend.position = "right", panel.grid.minor = element_blank()) +
    labs(title = t_text, subtitle = sprintf("log-rank test chisq = %.1f, p = %.3f", lr$chisq, p_l), x = "months", y = "survival estimate (%)")
  return(p)
}

p_km_hdb <- plot_km(df_sv, cl_test_pred[valid_outcomes_idx], "kaplan-meier test set hdbscan")
p_km_ru <- plot_km(df_sv, cl_test_rules[valid_outcomes_idx], "kaplan-meier test set rules")

print(p_km_hdb)
print(p_km_ru)








#cal plot -----
p_val <- p_test[valid_outcomes_idx]
y_val <- Y_test[valid_outcomes_idx]

span_val <- 0.75
grid_n <- 200
x_grid <- seq(0, 1, length.out = grid_n)

loess_fit <- tryCatch(
  loess(y ~ p, data = data.frame(y = y_val, p = p_val), span = span_val, degree = 1, control = loess.control(surface = "direct")),
  error = function(e) NULL
)

if (!is.null(loess_fit)) {
  pred_obj <- tryCatch(
    predict(loess_fit, newdata = data.frame(p = x_grid), se = TRUE),
    error = function(e) NULL
  )
  
  if (!is.null(pred_obj)) {
    y_hat <- as.numeric(pred_obj$fit)
    se_hat <- as.numeric(pred_obj$se.fit)
    
    y_hat[x_grid < min(p_val) | x_grid > max(p_val)] <- NA_real_
    se_hat[x_grid < min(p_val) | x_grid > max(p_val)] <- NA_real_
    
    ci_low <- pmax(0, y_hat - 1.96 * se_hat)
    ci_high <- pmin(1, y_hat + 1.96 * se_hat)
    
    curve_test <- data.frame(x = x_grid, y = y_hat, ymin = ci_low, ymax = ci_high)
  } else {
    curve_test <- data.frame(x = numeric(0), y = numeric(0), ymin = numeric(0), ymax = numeric(0))
  }
} else {
  curve_test <- data.frame(x = numeric(0), y = numeric(0), ymin = numeric(0), ymax = numeric(0))
}

p_cal_test <- ggplot() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_ribbon(data = curve_test, aes(x = x, ymin = ymin, ymax = ymax), fill = "#2166AC", alpha = 0.2, na.rm = TRUE) +
  geom_line(data = curve_test, aes(x = x, y = y), color = "#2166AC", linewidth = 1.3, na.rm = TRUE) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  labs(
    title = "XGBoost calibration - external test set",
    subtitle = "Validation curve (LOESS) with 95% Confidence Interval",
    x = "Predicted probability",
    y = "Observed probability"
  ) +
  theme_classic(base_size = 12)

print(p_cal_test)



