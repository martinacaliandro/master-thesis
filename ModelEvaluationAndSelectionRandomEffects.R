rm(list=ls())
#graphics.off()
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
library(splines)
library(xgboost)
library(lightgbm)

data <- read.csv('dataclean_final_withNA.csv', header=TRUE)

# Preprocessing iniziale
sum(is.na(data))
colnames(data)[colnames(data) == "Overall.Survival.at.2.years"] <- "OS2ys"
colnames(data)[colnames(data) == "mean_dose_vanacava_S"] <- "mean_dose_venacava_S"
colnames(data)[colnames(data) == "eud_venacavaI"] <- "eud_venacava_I"
colnames(data)[colnames(data) == "eud_vanacavaS"] <- "eud_venacava_S"
data[data == ""] <- NA
data$KPS <- factor(data$KPS)
data$PDL.1.level <- factor(data$PDL.1.level)
data$schedule_ok_chemo <- factor(data$schedule_ok_chemo)
data$cycles_chemo_sum <- factor(data$cycles_chemo_sum)
data$cycles_immuno <- factor(data$cycles_immuno)
data$prescription_dose_radio <- factor(data$prescription_dose_radio)
data$number_of_fractions_radio <- factor(data$number_of_fractions_radio)
data[] <- lapply(data, function(x) {
  if (!is.numeric(x)) as.factor(x) else x
})
data$OS2ys_num <- ifelse(data$OS2ys == "YES", 1, 0)
data$hospital <- factor(substr(as.character(data$ID), 3, 3), levels = as.character(1:7))
data$aortic_implant_present <- as.factor(data$aortic_implant_present)
data$cardiac_implant_present <- as.factor(data$cardiac_implant_present)

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
data$riskfactors2 <- ifelse(data$risk_sum2 > 0, 1, 0)
data$riskfactors <- as.factor(data$riskfactors)
data$riskfactors1 <- as.factor(data$riskfactors1)
data$riskfactors2 <- as.factor(data$riskfactors2)
data$riskfactors3 <- as.factor(data$riskfactors3)
data$Smoking.status <- as.factor(data$Smoking.status)
data$Diabetes.drugs <- as.factor(data$Diabetes.drugs)
data$Antilipidic.drugs <- factor(data$Antilipidic.drugs)
data$Hypertension <- as.factor(data$Hypertension)
data$antithrombotic <- as.factor(data$antithrombotic)
data$Beta.Blockers <- as.factor(data$Beta.Blockers)
data$Diuretics <- as.factor(data$Diuretics)
data$Volume_PTV_log <- log1p(data$Volume_PTV) 
data$LAA_log <- log1p(data$LAA) 

# Iperparametri
hp_gmerf  <- list(ntree = 500, nodesize = 20, itmax = 15, toll = 0.05)
hp_mexgb  <- list(max_depth = 2, eta = 0.05, lambda = 5, subsample = 0.7, nrounds = 100, itmax = 15, toll = 0.05)
hp_melgbm <- list(max_depth = 2, learning_rate = 0.05, min_data_in_leaf = 20, nrounds = 100, itmax = 15, toll = 0.05)

# Calcolo Metriche
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

fit_glmer_safe <- function(fmla, data, offset_vec = NULL) {
  tryCatch({
    if (is.null(offset_vec)) {
      suppressWarnings(glmer(fmla, data = data, family = binomial()))
    } else {
      suppressWarnings(glmer(fmla, data = data, family = binomial(), offset = as.numeric(offset_vec)))
    }
  }, error = function(e) NULL)
}

# Funzioni ad effetti misti per 1-split
run_one_split_glmm <- function(X, Y, Z, cont_vars, folds) {
  pooled <- numeric(length(Y)) + 0.5
  covs <- colnames(X)
  fmla <- as.formula(paste("Y ~", paste(covs, collapse = " + "), "+ (1 | H)"))
  for (i in seq_len(length(folds))) {
    test_idx <- folds[[i]]
    X_train <- X[-test_idx, , drop=FALSE]; Y_train <- Y[-test_idx]; Z_train <- Z[-test_idx]
    X_val <- X[test_idx, , drop=FALSE]; Z_val <- Z[test_idx]
    X_train_sc <- X_train; X_val_sc <- X_val
    is_num <- colnames(X) %in% cont_vars
    if (any(is_num)) {
      sc_p <- scale(X_train[, is_num, drop=FALSE])
      X_train_sc[, is_num] <- sc_p
      X_val_sc[, is_num] <- scale(X_val[, is_num, drop=FALSE], center = attr(sc_p, "scaled:center"), scale = attr(sc_p, "scaled:scale"))
    }
    td <- data.frame(X_train_sc, Y = Y_train, H = Z_train)
    vd <- data.frame(X_val_sc, H = Z_val)
    fit <- fit_glmer_safe(fmla, td)
    if (!is.null(fit)) {
      pooled[test_idx] <- as.numeric(predict(fit, newdata = vd, type = "response", allow.new.levels = TRUE))
    } else {
      pooled[test_idx] <- NA
    }
  }
  as.numeric(unlist(pooled))
}

run_one_split_gmerf <- function(X, Y, Z, cont_vars, folds, hp) {
  pooled <- numeric(length(Y)) + 0.5
  for (i in seq_len(length(folds))) {
    test_idx <- folds[[i]]
    X_train <- X[-test_idx, , drop=FALSE]; Y_train <- Y[-test_idx]; Z_train <- Z[-test_idx]
    X_val <- X[test_idx, , drop=FALSE]; Z_val <- Z[test_idx]
    X_train_sc <- X_train; X_val_sc <- X_val
    is_num <- colnames(X) %in% cont_vars
    
    if (any(is_num)) {
      sc_p <- scale(X_train[, is_num, drop=FALSE])
      X_train_sc[, is_num] <- sc_p
      X_val_sc[, is_num] <- scale(X_val[, is_num, drop=FALSE], center = attr(sc_p, "scaled:center"), scale = attr(sc_p, "scaled:scale"))
    }
    
    tmp_df <- data.frame(Y_train = Y_train, X_train_sc)
    eta_static <- as.numeric(glm(Y_train ~ ., data = tmp_df, family = binomial())$linear.predictors)
    
    b <- rep(0, length(unique(Z_train))); names(b) <- unique(Z_train)
    converged <- FALSE; iter <- 0; fitrf <- NULL
    
    while (!converged && iter < hp$itmax) {
      iter <- iter + 1; b_old <- b
      target <- eta_static - as.numeric(b[as.character(Z_train)])
      fitrf <- randomForest(X_train_sc, target, ntree = hp$ntree, nodesize = hp$nodesize)
      fitglmer <- fit_glmer_safe(Y_train ~ 1 + (1 | Z_train), data.frame(Y_train = Y_train, Z_train = Z_train), offset_vec = as.numeric(fitrf$predicted))
      
      if (is.null(fitglmer)) break
      b_new <- as.numeric(ranef(fitglmer)$Z_train[, 1]); names(b_new) <- rownames(ranef(fitglmer)$Z_train)
      b <- b_new
      if (max(abs(b - b_old)) < hp$toll) converged <- TRUE
    }
    
    if (!is.null(fitrf)) {
      b_val <- as.numeric(b[as.character(Z_val)]); b_val[is.na(b_val)] <- 0
      lk <- as.numeric(predict(fitrf, X_val_sc)) + b_val
      pooled[test_idx] <- as.numeric(exp(lk) / (1 + exp(lk)))
    } else {
      pooled[test_idx] <- NA
    }
  }
  as.numeric(unlist(pooled))
}

run_one_split_mexgb <- function(X_mat, Y, Z, cont_vars, folds, hp) {
  pooled <- numeric(length(Y)) + 0.5
  cols_scale <- intersect(colnames(X_mat), cont_vars)
  for (i in seq_len(length(folds))) {
    test_idx <- folds[[i]]
    X_train <- X_mat[-test_idx, , drop = FALSE]; Y_train <- Y[-test_idx]; Z_train <- Z[-test_idx]
    X_val <- X_mat[test_idx, , drop = FALSE]; Z_val <- Z[test_idx]
    
    if (length(cols_scale) > 0) {
      sc_p <- scale(X_train[, cols_scale, drop = FALSE])
      X_train[, cols_scale] <- sc_p
      X_val[, cols_scale] <- scale(X_val[, cols_scale, drop = FALSE], center = attr(sc_p, "scaled:center"), scale = attr(sc_p, "scaled:scale"))
    }
    
    tmp_df <- data.frame(Y_train = Y_train, X_train)
    eta_static <- as.numeric(glm(Y_train ~ ., data = tmp_df, family = binomial())$linear.predictors)
    
    b <- rep(0, length(unique(Z_train))); names(b) <- unique(Z_train)
    converged <- FALSE; iter <- 0; fit <- NULL
    
    while (!converged && iter < hp$itmax) {
      iter <- iter + 1; b_old <- b
      target <- eta_static - as.numeric(b[as.character(Z_train)])
      fit <- xgb.train(list(objective = "reg:squarederror", max_depth = hp$max_depth, eta = hp$eta, lambda = hp$lambda, subsample = hp$subsample), xgb.DMatrix(X_train, label = target), hp$nrounds, verbose = 0)
      f_xgb <- as.numeric(predict(fit, xgb.DMatrix(X_train)))
      
      fitglmer <- fit_glmer_safe(Y_train ~ 1 + (1 | Z_train), data.frame(Y_train = Y_train, Z_train = Z_train), offset_vec = f_xgb)
      if (is.null(fitglmer)) break
      b_new <- as.numeric(ranef(fitglmer)$Z_train[, 1]); names(b_new) <- rownames(ranef(fitglmer)$Z_train)
      b <- b_new
      if (max(abs(b - b_old)) < hp$toll) converged <- TRUE
    }
    
    if (!is.null(fit)) {
      b_val <- as.numeric(b[as.character(Z_val)]); b_val[is.na(b_val)] <- 0
      lk <- as.numeric(predict(fit, xgb.DMatrix(X_val))) + b_val
      pooled[test_idx] <- as.numeric(exp(lk) / (1 + exp(lk)))
    } else {
      pooled[test_idx] <- NA
    }
  }
  as.numeric(unlist(pooled))
}

run_one_split_melgbm <- function(X_mat, Y, Z, cont_vars, folds, hp) {
  pooled <- numeric(length(Y)) + 0.5
  cols_scale <- intersect(colnames(X_mat), cont_vars)
  for (i in seq_len(length(folds))) {
    test_idx <- folds[[i]]
    X_train <- X_mat[-test_idx, , drop = FALSE]; Y_train <- Y[-test_idx]; Z_train <- Z[-test_idx]
    X_val <- X_mat[test_idx, , drop = FALSE]; Z_val <- Z[test_idx]
    
    if (length(cols_scale) > 0) {
      sc_p <- scale(X_train[, cols_scale, drop = FALSE])
      X_train[, cols_scale] <- sc_p
      X_val[, cols_scale] <- scale(X_val[, cols_scale, drop = FALSE], center = attr(sc_p, "scaled:center"), scale = attr(sc_p, "scaled:scale"))
    }
    
    tmp_df <- data.frame(Y_train = Y_train, X_train)
    eta_static <- as.numeric(glm(Y_train ~ ., data = tmp_df, family = binomial())$linear.predictors)
    
    b <- rep(0, length(unique(Z_train))); names(b) <- unique(Z_train)
    converged <- FALSE; iter <- 0; fit <- NULL
    
    while (!converged && iter < hp$itmax) {
      iter <- iter + 1; b_old <- b
      target <- eta_static - as.numeric(b[as.character(Z_train)])
      fit <- lgb.train(list(objective = "regression", metric = "l2", max_depth = hp$max_depth, learning_rate = hp$learning_rate, min_data_in_leaf = hp$min_data_in_leaf, verbose = -1), lgb.Dataset(X_train, label = target), hp$nrounds)
      f_lgb <- as.numeric(predict(fit, X_train))
      
      fitglmer <- fit_glmer_safe(Y_train ~ 1 + (1 | Z_train), data.frame(Y_train = Y_train, Z_train = Z_train), offset_vec = f_lgb)
      if (is.null(fitglmer)) break
      b_new <- as.numeric(ranef(fitglmer)$Z_train[, 1]); names(b_new) <- rownames(ranef(fitglmer)$Z_train)
      b <- b_new
      if (max(abs(b - b_old)) < hp$toll) converged <- TRUE
    }
    
    if (!is.null(fit)) {
      b_val <- as.numeric(b[as.character(Z_val)]); b_val[is.na(b_val)] <- 0
      lk <- as.numeric(predict(fit, X_val)) + b_val
      pooled[test_idx] <- as.numeric(exp(lk) / (1 + exp(lk)))
    } else {
      pooled[test_idx] <- NA
    }
  }
  as.numeric(unlist(pooled))
}

# Funzione di confronto SOLO MODELLI AD EFFETTI MISTI
run_comparison_me_only <- function(covariates_cmp, cont_cmp, n_repeats = 10) {
  
  # Filtro condizionale per l'ospedale 6
  has_risk <- any(c("riskfactors1", "riskfactors2") %in% covariates_cmp)
  
  df_proc <- data
  if (has_risk) {
    df_proc <- subset(df_proc, hospital != "6")
    df_proc$hospital <- droplevels(df_proc$hospital)
  }
  
  X_cmp <- df_proc[, covariates_cmp, drop = FALSE]
  Y_cmp <- df_proc$OS2ys_num
  Z_cmp <- df_proc$hospital
  
  idx_cmp <- complete.cases(X_cmp) & !is.na(Y_cmp) & !is.na(Z_cmp)
  X_cmp <- X_cmp[idx_cmp, , drop = FALSE]
  Y_cmp <- Y_cmp[idx_cmp]
  Z_cmp <- Z_cmp[idx_cmp]
  
  X_mat_cmp <- model.matrix(~ ., data = X_cmp)[, -1, drop = FALSE]
  
  metric_names <- c("AUC", "F1", "F1_0", "F1_macro", "Precision", "Precision_0", "Recall", "Recall_0", "Brier", "LogLoss", "Cal_Intercept", "Cal_Slope", "Threshold")
  model_names <- c("GLMM", "GMERF", "MEXGB", "MELGBM")
  
  results_list <- setNames(vector("list", length(model_names)), model_names)
  for (m in model_names) {
    results_list[[m]] <- matrix(NA_real_, nrow = n_repeats, ncol = length(metric_names), dimnames = list(NULL, metric_names))
  }
  
  for (r in seq_len(n_repeats)) {
    set.seed(r)
    folds_r <- createFolds(as.factor(Y_cmp), k = 5, list = TRUE, returnTrain = FALSE)
    results_list[["GLMM"]][r, ]   <- calc_metrics_vec(run_one_split_glmm(X_cmp, Y_cmp, Z_cmp, cont_cmp, folds_r), Y_cmp)
    results_list[["GMERF"]][r, ]  <- calc_metrics_vec(run_one_split_gmerf(X_cmp, Y_cmp, Z_cmp, cont_cmp, folds_r, hp_gmerf), Y_cmp)
    results_list[["MEXGB"]][r, ]  <- calc_metrics_vec(run_one_split_mexgb(X_mat_cmp, Y_cmp, Z_cmp, cont_cmp, folds_r, hp_mexgb), Y_cmp)
    results_list[["MELGBM"]][r, ] <- calc_metrics_vec(run_one_split_melgbm(X_mat_cmp, Y_cmp, Z_cmp, cont_cmp, folds_r, hp_melgbm), Y_cmp)
  }
  
  do.call(rbind, lapply(model_names, function(m) {
    mat <- results_list[[m]]
    means <- colMeans(mat, na.rm = TRUE)
    ci_low <- apply(mat, 2, quantile, probs = 0.025, na.rm = TRUE)
    ci_high <- apply(mat, 2, quantile, probs = 0.975, na.rm = TRUE)
    data.frame(Model = m, Metric = colnames(mat), Mean = round(means, 4), CI_low = round(ci_low, 4), CI_high = round(ci_high, 4))
  }))
}

# Definizione Lista Subset

subset_list <- list(
  GS5 = list(covs = c("Age","Volume_PTV","schedule_ok_chemo","IMMUNO","riskfactors1","riskfactors2"),
             cont = c("Age", "Volume_PTV")),
  
  GS7 = list(covs = c("Age","KPS", "Volume_PTV","schedule_ok_chemo","IMMUNO","riskfactors2"),
             cont = c("Age", "Volume_PTV")),
  
  Lhole = list(covs = c("Age","Volume_PTV","IMMUNO","eud_lungs","Subcutaneous_Fat", "riskfactors1"),
               cont = c("Age","Volume_PTV","eud_lungs","Subcutaneous_Fat")),
  
  Lspecific1a = list(covs = c("Age","Volume_PTV","schedule_ok_chemo","IMMUNO","riskfactors1","riskfactors2","EDIC","eud_lungs"),
                     cont = c("Age","Volume_PTV","EDIC","eud_lungs")),
  
  Lspecific2a = list(covs = c("Age","Volume_PTV","schedule_ok_chemo","IMMUNO","riskfactors1","riskfactors2","Subcutaneous_Fat","Skeletal_Muscle"),
                     cont = c("Age","Volume_PTV", "Subcutaneous_Fat", "Skeletal_Muscle")),
  
  Lspecific3a = list(covs = c("Age","Volume_PTV","schedule_ok_chemo","IMMUNO","riskfactors1","riskfactors2","eud_lungs","EDIC", "Subcutaneous_Fat","Skeletal_Muscle"),
                     cont = c("Age","Volume_PTV","EDIC","eud_lungs", "Subcutaneous_Fat", "Skeletal_Muscle")),
  
  Lspecific1b = list(covs = c("Age","Volume_PTV","IMMUNO","riskfactors1","EDIC","eud_lungs"),
                     cont = c("Age","Volume_PTV","EDIC","eud_lungs")),
  
  Lspecific2b = list(covs = c("Age","Volume_PTV","IMMUNO","riskfactors1","Subcutaneous_Fat","Skeletal_Muscle"),
                     cont = c("Age","Volume_PTV", "Subcutaneous_Fat", "Skeletal_Muscle")),
  
  Lspecific3b = list(covs = c("Age","Volume_PTV","IMMUNO","riskfactors1","eud_lungs","EDIC", "Subcutaneous_Fat","Skeletal_Muscle"),
                     cont = c("Age","Volume_PTV","EDIC","eud_lungs", "Subcutaneous_Fat", "Skeletal_Muscle")),
  
  Bhole = list(covs = c("Age","Volume_PTV","LAA", "riskfactors1"),
               cont = c("Age","Volume_PTV","LAA")),
  
  Bspecific1a = list(covs = c("Age","Volume_PTV","schedule_ok_chemo","IMMUNO","riskfactors1","riskfactors2"),
                     cont = c("Age","Volume_PTV")),
  
  Bspecific1b = list(covs = c("Age","Volume_PTV","IMMUNO","riskfactors1"),
                     cont = c("Age","Volume_PTV")),
  
  Bspecific2a = list(covs = c("Age","Volume_PTV","schedule_ok_chemo","IMMUNO","riskfactors1","riskfactors2", "LAA","Subcutaneous_Fat", "Skeletal_Muscle", "diameter.kidney.sep..mm."),
                     cont = c("Age","Volume_PTV","Subcutaneous_Fat", "LAA", "Skeletal_Muscle", "diameter.kidney.sep..mm.")),
  
  Bspecific2b = list(covs = c("Age","Volume_PTV","IMMUNO","riskfactors1","LAA","Subcutaneous_Fat","Skeletal_Muscle", "diameter.kidney.sep..mm."),
                     cont = c("Age","Volume_PTV","Subcutaneous_Fat", "LAA", "Skeletal_Muscle", "diameter.kidney.sep..mm.")),
  
  Lassohole = list(covs = c("Age","Volume_PTV","IMMUNO","eud_lungs","Subcutaneous_Fat","aorta_agatston", "riskfactors1", "riskfactors2"),
                   cont = c("Age","Volume_PTV","eud_lungs","Subcutaneous_Fat","aorta_agatston")),
  
  Lassospecific1 = list(covs = c("Age","Volume_PTV","IMMUNO","schedule_ok_chemo","eud_lad", "eud_lungs","EDIC", "riskfactors1", "riskfactors2"),
                        cont = c("Age","Volume_PTV", "eud_lad", "eud_lungs","EDIC")),
  
  Lassospecific2 = list(covs = c("Age","Volume_PTV","IMMUNO","schedule_ok_chemo", "Subcutaneous_Fat","Skeletal_Muscle", "Intermuscular_Fat", "diameter.kidney.sep..mm.", "aortic_implant_present", "aorta_agatston", "riskfactors1", "riskfactors2"),
                        cont = c("Age","Volume_PTV","Subcutaneous_Fat","Skeletal_Muscle", "Intermuscular_Fat","diameter.kidney.sep..mm.",  "aorta_agatston")),
  
  Lassospecific = list(covs = c("Age","Volume_PTV","IMMUNO","schedule_ok_chemo","eud_lad", "eud_lungs","EDIC", "Subcutaneous_Fat", "Skeletal_Muscle", "Intermuscular_Fat", "diameter.kidney.sep..mm.", "aortic_implant_present", "aorta_agatston", "riskfactors1", "riskfactors2"),
                       cont = c("Age","Volume_PTV","eud_lad", "eud_lungs","EDIC", "Subcutaneous_Fat","Skeletal_Muscle", "Intermuscular_Fat","diameter.kidney.sep..mm.",  "aorta_agatston"))
  
)


subset_list <- list(
  GS5_2 = list(covs = c("Age","Volume_PTV","schedule_ok_chemo","IMMUNO","riskfactors1","riskfactors2"),
               cont = c("Age", "Volume_PTV")),
  
  GS7_2 = list(covs = c("Age","KPS", "Volume_PTV","schedule_ok_chemo","IMMUNO","riskfactors2"),
               cont = c("Age", "Volume_PTV")),
  
  Lhole_2 = list(covs = c("Age","KPS", "Volume_PTV","IMMUNO","eud_lungs","mean_dose_atriumL", "Subcutaneous_Fat", "riskfactors2"),
                 cont = c("Age", "Volume_PTV","eud_lungs","mean_dose_atriumL", "Subcutaneous_Fat")),
  
  Lspecific1a_2 = list(covs = c("Age","KPS", "Volume_PTV","schedule_ok_chemo","IMMUNO","riskfactors2","EDIC","eud_lungs"),
                       cont = c("Age", "Volume_PTV","EDIC","eud_lungs")),
  
  Lspecific2a_2 = list(covs = c("Age","KPS", "Volume_PTV","schedule_ok_chemo","IMMUNO","riskfactors2","Subcutaneous_Fat","Skeletal_Muscle"),
                       cont = c("Age", "Volume_PTV", "Subcutaneous_Fat", "Skeletal_Muscle")),
  
  Lspecific3a_2 = list(covs = c("Age","KPS", "Volume_PTV","schedule_ok_chemo","IMMUNO","riskfactors2","eud_lungs","EDIC", "Subcutaneous_Fat","Skeletal_Muscle"),
                       cont = c("Age", "Volume_PTV","EDIC","eud_lungs", "Subcutaneous_Fat", "Skeletal_Muscle")),
  
  Lspecific1b_2 = list(covs = c("Age","KPS", "Volume_PTV","IMMUNO","riskfactors2","EDIC","eud_lungs"),
                       cont = c("Age", "Volume_PTV","EDIC","eud_lungs")),
  
  Lspecific2b_2 = list(covs = c("Age","KPS", "Volume_PTV","IMMUNO","riskfactors2","Subcutaneous_Fat", "Skeletal_Muscle"),
                       cont = c("Age", "Volume_PTV", "Subcutaneous_Fat", "Skeletal_Muscle")),
  
  Lspecific3b_2 = list(covs = c("Age","KPS", "Volume_PTV","IMMUNO","riskfactors2","eud_lungs","EDIC", "Subcutaneous_Fat","Skeletal_Muscle"),
                       cont = c("Age", "Volume_PTV","EDIC","eud_lungs", "Subcutaneous_Fat", "Skeletal_Muscle")),
  
  Bhole_2 = list(covs = c("Age","Volume_PTV","LAA"),
                 cont = c("Age","Volume_PTV","LAA")),
  
  Bspecific1a_2 = list(covs = c("Age","KPS", "Volume_PTV","schedule_ok_chemo","IMMUNO","riskfactors2"),
                       cont = c("Age", "Volume_PTV")),
  
  Bspecific1b_2 = list(covs = c("Age","Volume_PTV","IMMUNO"),
                       cont = c("Age","Volume_PTV")),
  
  Bspecific2a_2 = list(covs = c("Age","KPS", "Volume_PTV","schedule_ok_chemo","IMMUNO","riskfactors2", "LAA","Subcutaneous_Fat", "Skeletal_Muscle", "diameter.kidney.sep..mm."),
                       cont = c("Age", "Volume_PTV","Subcutaneous_Fat", "LAA", "Skeletal_Muscle", "diameter.kidney.sep..mm.")),
  
  Bspecific2b_2 = list(covs = c("Age","Volume_PTV","IMMUNO","LAA","Subcutaneous_Fat","Skeletal_Muscle", "diameter.kidney.sep..mm."),
                       cont = c("Age","Volume_PTV","Subcutaneous_Fat", "LAA", "Skeletal_Muscle", "diameter.kidney.sep..mm.")),
  
  Lassohole_2 = list(covs = c("Age","KPS", "Volume_PTV","IMMUNO","schedule_ok_chemo", "eud_lungs", "mean_dose_atriumL", "Subcutaneous_Fat","aorta_agatston", "aortic_implant_present", "riskfactors2"),
                     cont = c("Age", "Volume_PTV","eud_lungs","mean_dose_atriumL", "Subcutaneous_Fat","aorta_agatston")),
  
  Lassospecific1_2 = list(covs = c("Age","KPS", "Volume_PTV","schedule_ok_chemo","IMMUNO","riskfactors2", "eud_lad", "eud_lungs","EDIC"),
                          cont = c("Age", "Volume_PTV", "eud_lad", "eud_lungs","EDIC")),
  
  Lassospecific2_2 = list(covs = c("Age","KPS", "Volume_PTV","schedule_ok_chemo","IMMUNO","riskfactors2", "Subcutaneous_Fat","Skeletal_Muscle", "Intermuscular_Fat", "diameter.kidney.sep..mm.", "aortic_implant_present", "aorta_agatston"),
                          cont = c("Age", "Volume_PTV","Subcutaneous_Fat","Skeletal_Muscle", "Intermuscular_Fat","diameter.kidney.sep..mm.",  "aorta_agatston")),
  
  Lassospecific_2 = list(covs = c("Age","KPS", "Volume_PTV","schedule_ok_chemo", "IMMUNO","eud_lad", "eud_lungs","EDIC", "Subcutaneous_Fat", "Skeletal_Muscle", "Intermuscular_Fat", "diameter.kidney.sep..mm.", "aortic_implant_present", "aorta_agatston", "riskfactors2"),
                         cont = c("Age", "Volume_PTV","eud_lad", "eud_lungs","EDIC", "Subcutaneous_Fat","Skeletal_Muscle", "Intermuscular_Fat","diameter.kidney.sep..mm.",  "aorta_agatston"))
  
  
)





all_results <- list()
for (subset_name in names(subset_list)) {
  cat("\nIn elaborazione (Solo Effetti Misti):", subset_name, "...\n")
  all_results[[subset_name]] <- run_comparison_me_only(subset_list[[subset_name]]$covs, subset_list[[subset_name]]$cont)
}

final_comparison <- do.call(rbind, Map(function(n, tab) cbind(Subset = n, tab), names(all_results), all_results))
options(max.print = 40000)

print(final_comparison)


res_np <- lapply(names(subset_list), function(s_name) {
  covs <- subset_list[[s_name]]$covs
  cont_vars <- subset_list[[s_name]]$cont
  
  has_risk <- any(c("riskfactors1", "riskfactors2") %in% covs)
  
  df_sub_base <- data
  if (has_risk) {
    df_sub_base <- subset(df_sub_base, hospital != "6")
    df_sub_base$hospital <- droplevels(df_sub_base$hospital)
  }
  
  keep <- complete.cases(df_sub_base[, covs, drop = FALSE]) & 
    !is.na(df_sub_base$OS2ys_num) & 
    !is.na(df_sub_base$hospital)
  
  df_sub <- df_sub_base[keep, covs, drop = FALSE]
  
  for(v in covs) {
    if(!v %in% cont_vars) {
      df_sub[[v]] <- as.factor(df_sub[[v]])
    }
  }
  
  X_mat <- model.matrix(~ ., data = df_sub)[, -1, drop = FALSE]
  
  data.frame(
    Subset      = s_name,
    n           = nrow(df_sub),
    n_hospitals = length(unique(df_sub_base$hospital[keep])),
    p_covs      = length(covs),
    p_dummy     = ncol(X_mat),
    stringsAsFactors = FALSE
  )
})

summary_np <- do.call(rbind, res_np)
summary_np
print(summary_np, row.names = FALSE)