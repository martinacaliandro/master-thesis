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
library(glmnet)
library(Boruta)
library(splines)
library(xgboost)
library(lightgbm)

data <- read.csv('dataclean_final_withNA.csv', header=T)

# Preprocessing
sum(is.na(data))
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
data[] <- lapply(data, function(x) {
  if (!is.numeric(x)) as.factor(x) else x
})
data$OS2ys_num <- ifelse(data$OS2ys == "YES", 1, 0)
data$hospital <- factor(substr(as.character(data$ID), 3, 3), levels = as.character(1:7))
#data$KPS <- as.numeric(as.character(data$KPS))
data$aortic_implant_present = as.factor(data$aortic_implant_present)
data$cardiac_implant_present = as.factor(data$cardiac_implant_present)

cols <- names(data)
idx <- which(cols == "Stage.at.RT")
data <- data[, c(cols[1:idx], "Volume_PTV", setdiff(cols[(idx+1):length(cols)], "Volume_PTV"))]

datacrf = data[, c(1:41, 92)]
datadose = data[, 42:64]
datavol = data[, 65:77]
datametab = data[, 78:91]

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
data$riskfactors3 <- rowSums(data[, c("Diabetes.drugs", "Antilipidic.drugs", "Hypertension")], na.rm = FALSE) 
data$riskfactors2 = ifelse(data$risk_sum2>0,1,0)
data$riskfactors = as.factor(data$riskfactors)
data$riskfactors1 = as.factor(data$riskfactors1)
data$riskfactors2 = as.factor(data$riskfactors2)
data$riskfactors3 = as.factor(data$riskfactors3)
data$Smoking.status = as.factor(data$Smoking.status)
data$Diabetes.drugs = as.factor(data$Diabetes.drugs)
data$Antilipidic.drugs  = factor(data$Antilipidic.drugs)
data$Hypertension = as.factor(data$Hypertension)
data$antithrombotic = as.factor(data$antithrombotic)
data$Beta.Blockers = as.factor(data$Beta.Blockers)
data$Diuretics = as.factor(data$Diuretics)
data$Volume_PTV_log = log1p(data$Volume_PTV) 
data$LAA_log = log1p(data$LAA) 






#analysis of random -----
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

data$KPS = as.factor(data$KPS)

hp_rf <- list(ntree = 500, nodesize = 20)
hp_xgb <- list(max_depth = 2, eta = 0.05, lambda = 5, subsample = 0.7, nrounds = 100)
hp_lgb <- list(max_depth = 2, learning_rate = 0.05, min_data_in_leaf = 20, nrounds = 100)
hp_gmerf <- list(ntree = 500, nodesize = 20, itmax = 15, toll = 0.05)
hp_mexgb <- list(max_depth = 2, eta = 0.05, lambda = 5, subsample = 0.7, nrounds = 100, itmax = 15, toll = 0.05)
hp_melgbm <- list(max_depth = 2, learning_rate = 0.05, min_data_in_leaf = 20, nrounds = 100, itmax = 15, toll = 0.05)


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

evaluate_full_model <- function(model_type, subset_name, subset_info, data_df) {
  covs <- subset_info$covs
  cont_vars <- subset_info$cont
  
  df_model <- data_df[, c(covs, "OS2ys_num", "hospital"), drop = FALSE]
  df_model <- df_model[complete.cases(df_model), ]
  
  X <- df_model[, covs, drop = FALSE]
  Y <- df_model$OS2ys_num
  Z <- df_model$hospital
  
  X_sc <- X
  is_num <- colnames(X) %in% cont_vars
  if (any(is_num)) {
    X_sc[, is_num] <- scale(X[, is_num])
  }
  
  preds <- NULL
  model_obj <- NULL
  
  if (model_type == "GLM") {
    fit <- glm(Y ~ ., data = data.frame(Y = Y, X_sc), family = binomial())
    preds <- predict(fit, type = "response")
    model_obj <- fit
    
  } else if (model_type == "GLMM") {
    fmla <- as.formula(paste("Y ~", paste(colnames(X_sc), collapse = " + "), "+ (1 | Z)"))
    fit <- glmer(fmla, data = data.frame(Y = Y, X_sc, Z = Z), family = binomial())
    preds <- predict(fit, type = "response", allow.new.levels = TRUE)
    model_obj <- fit
    
  } else if (model_type == "RF") {
    fit <- randomForest(X_sc, as.factor(Y), ntree = hp_rf$ntree, nodesize = hp_rf$nodesize, importance = TRUE)
    preds <- predict(fit, X_sc, type = "prob")[, "1"]
    model_obj <- fit
    
  } else if (model_type == "XGB") {
    X_mat <- model.matrix(~ ., data = X_sc)[, -1, drop = FALSE]
    fit <- xgb.train(
      list(objective = "binary:logistic", max_depth = hp_xgb$max_depth, eta = hp_xgb$eta, lambda = hp_xgb$lambda, subsample = hp_xgb$subsample),
      xgb.DMatrix(X_mat, label = Y), hp_xgb$nrounds, verbose = 0
    )
    preds <- as.numeric(predict(fit, xgb.DMatrix(X_mat)))
    model_obj <- fit
    
  } else if (model_type == "LGBM") {
    X_mat <- model.matrix(~ ., data = X_sc)[, -1, drop = FALSE]
    fit <- lgb.train(
      list(objective = "binary", metric = "auc", max_depth = hp_lgb$max_depth, learning_rate = hp_lgb$learning_rate, min_data_in_leaf = hp_lgb$min_data_in_leaf, verbose = -1),
      lgb.Dataset(X_mat, label = Y), hp_lgb$nrounds
    )
    preds <- as.numeric(predict(fit, X_mat))
    model_obj <- fit
    
  } else if (model_type == "GMERF") {
    eta <- as.numeric(glm(Y ~ ., data = data.frame(Y = Y, X_sc), family = binomial())$linear.predictors)
    b <- rep(0, length(unique(Z))); names(b) <- unique(Z)
    converged <- FALSE; iter <- 0; fitrf <- NULL; fitglmer <- NULL
    
    while (!converged && iter < hp_gmerf$itmax) {
      iter <- iter + 1; b_old <- b
      target <- as.numeric(eta) - as.numeric(b[as.character(Z)])
      fitrf <- randomForest(X_sc, target, ntree = hp_gmerf$ntree, nodesize = hp_gmerf$nodesize)
      
      fitglmer <- tryCatch({
        suppressWarnings(glmer(Y ~ 1 + (1 | Z), data = data.frame(Y = Y, Z = Z), family = binomial(), offset = as.numeric(fitrf$predicted)))
      }, error = function(e) NULL)
      
      if (is.null(fitglmer)) break
      b_new <- as.numeric(ranef(fitglmer)$Z[, 1]); names(b_new) <- rownames(ranef(fitglmer)$Z)
      b <- b_new; eta <- as.numeric(predict(fitglmer, type = "link"))
      if (max(abs(b - b_old)) < hp_gmerf$toll) converged <- TRUE
    }
    
    if (!is.null(fitrf)) {
      b_val <- as.numeric(b[as.character(Z)]); b_val[is.na(b_val)] <- 0
      lk <- as.numeric(predict(fitrf, X_sc)) + b_val
      preds <- as.numeric(exp(lk) / (1 + exp(lk)))
      model_obj <- list(rf = fitrf, glmer = fitglmer, random_effects = b, iterations = iter)
    }
    
  } else if (model_type == "MEXGB") {
    X_mat <- model.matrix(~ ., data = X_sc)[, -1, drop = FALSE]
    tmp_df <- data.frame(Y = Y, X_sc)
    eta_link <- as.numeric(glm(Y ~ ., data = tmp_df, family = binomial())$linear.predictors)
    b <- rep(0, length(unique(Z))); names(b) <- unique(Z)
    converged <- FALSE; iter <- 0; fit <- NULL; fitglmer <- NULL
    
    while (!converged && iter < hp_mexgb$itmax) {
      iter <- iter + 1; b_old <- b
      target <- as.numeric(eta_link) - as.numeric(b[as.character(Z)])
      fit <- xgb.train(
        list(objective = "reg:squarederror", max_depth = hp_mexgb$max_depth, eta = hp_mexgb$eta, lambda = hp_mexgb$lambda, subsample = hp_mexgb$subsample),
        xgb.DMatrix(X_mat, label = target), hp_mexgb$nrounds, verbose = 0
      )
      f_xgb <- as.numeric(predict(fit, xgb.DMatrix(X_mat)))
      
      fitglmer <- tryCatch({
        suppressWarnings(glmer(Y ~ 1 + (1 | Z), data = data.frame(Y = Y, Z = Z), family = binomial(), offset = f_xgb))
      }, error = function(e) NULL)
      
      if (is.null(fitglmer)) break
      b_new <- as.numeric(ranef(fitglmer)$Z[, 1]); names(b_new) <- rownames(ranef(fitglmer)$Z)
      b <- b_new; eta_link <- as.numeric(predict(fitglmer, type = "link"))
      if (max(abs(b - b_old)) < hp_mexgb$toll) converged <- TRUE
    }
    
    if (!is.null(fit)) {
      b_val <- as.numeric(b[as.character(Z)]); b_val[is.na(b_val)] <- 0
      lk <- as.numeric(predict(fit, xgb.DMatrix(model.matrix(~ ., data = X_sc)[, -1, drop = FALSE]))) + b_val
      preds <- as.numeric(exp(lk) / (1 + exp(lk)))
      model_obj <- list(xgb = fit, glmer = fitglmer, random_effects = b, iterations = iter)
    }
    
  } else if (model_type == "MELGBM") {
    X_mat <- model.matrix(~ ., data = X_sc)[, -1, drop = FALSE]
    tmp_df <- data.frame(Y = Y, X_sc)
    eta_link <- as.numeric(glm(Y ~ ., data = tmp_df, family = binomial())$linear.predictors)
    b <- rep(0, length(unique(Z))); names(b) <- unique(Z)
    converged <- FALSE; iter <- 0; fit <- NULL; fitglmer <- NULL
    
    while (!converged && iter < hp_melgbm$itmax) {
      iter <- iter + 1; b_old <- b
      target <- as.numeric(eta_link) - as.numeric(b[as.character(Z)])
      fit <- lgb.train(
        list(objective = "regression", metric = "l2", max_depth = hp_melgbm$max_depth, learning_rate = hp_melgbm$learning_rate, min_data_in_leaf = hp_melgbm$min_data_in_leaf, verbose = -1),
        lgb.Dataset(X_mat, label = target), hp_melgbm$nrounds
      )
      f_lgb <- as.numeric(predict(fit, X_mat))
      
      fitglmer <- tryCatch({
        suppressWarnings(glmer(Y ~ 1 + (1 | Z), data = data.frame(Y = Y, Z = Z), family = binomial(), offset = f_lgb))
      }, error = function(e) NULL)
      
      if (is.null(fitglmer)) break
      b_new <- as.numeric(ranef(fitglmer)$Z[, 1]); names(b_new) <- rownames(ranef(fitglmer)$Z)
      b <- b_new; eta_link <- as.numeric(predict(fitglmer, type = "link"))
      if (max(abs(b - b_old)) < hp_melgbm$toll) converged <- TRUE
    }
    
    if (!is.null(fit)) {
      b_val <- as.numeric(b[as.character(Z)]); b_val[is.na(b_val)] <- 0
      lk <- as.numeric(predict(fit, model.matrix(~ ., data = X_sc)[, -1, drop = FALSE])) + b_val
      preds <- as.numeric(exp(lk) / (1 + exp(lk)))
      model_obj <- list(lgb = fit, glmer = fitglmer, random_effects = b, iterations = iter)
    }
  }
  
  metrics <- calc_metrics_vec(preds, Y)
  
  glmer_obj <- if (model_type %in% c("GMERF", "MEXGB", "MELGBM")) model_obj$glmer else if (model_type == "GLMM") model_obj else NULL
  
  sd_int <- NA
  icc_lat <- NA
  if (!is.null(glmer_obj)) {
    vc <- VarCorr(glmer_obj)
    if (length(vc) > 0) {
      varianza_intercept <- as.numeric(vc[[1]][1, 1])
      sd_int <- sqrt(varianza_intercept) # Estrae la vera Deviazione Standard (0.2792)
      icc_lat <- varianza_intercept / (varianza_intercept + (pi^2 / 3)) # ICC formula corretta
    }
  }
  
  dep_effects <- NULL
  try({
    if (model_type %in% c("GLM", "GLMM")) {
      dep_effects <- summary(model_obj)$coefficients
      
    } else if (model_type %in% c("RF", "GMERF")) {
      rf_target_obj <- if(model_type == "GMERF") model_obj$rf else model_obj
      dep_effects <- lapply(cont_vars, function(v) {
        if(v %in% colnames(X_sc)) {
          x_vals <- seq(min(X_sc[, v]), max(X_sc[, v]), length.out = 30)
          pred_grid <- sapply(x_vals, function(val) {
            tmp_data <- X_sc
            tmp_data[, v] <- val
            if (model_type == "RF") {
              mean(predict(rf_target_obj, tmp_data, type = "prob")[, "1"])
            } else {
              mean(predict(rf_target_obj, tmp_data))
            }
          })
          data.frame(x = x_vals, y = pred_grid)
        }
      })
      names(dep_effects) <- cont_vars
      
    } else if (model_type %in% c("XGB", "MEXGB")) {
      xgb_target_obj <- if(model_type == "MEXGB") model_obj$xgb else model_obj
      dep_effects <- lapply(cont_vars, function(v) {
        if(v %in% colnames(X_sc)) {
          x_vals <- seq(min(X_sc[, v]), max(X_sc[, v]), length.out = 30)
          pred_grid <- sapply(x_vals, function(val) {
            tmp_data <- X_sc
            tmp_data[, v] <- val
            tmp_mat <- model.matrix(~ ., data = tmp_data)[, -1, drop = FALSE]
            mean(predict(xgb_target_obj, xgb.DMatrix(tmp_mat)))
          })
          data.frame(x = x_vals, y = pred_grid)
        }
      })
      names(dep_effects) <- cont_vars
      
    } else if (model_type %in% c("LGBM", "MELGBM")) {
      lgb_target_obj <- if(model_type == "MELGBM") model_obj$lgb else model_obj
      dep_effects <- lapply(cont_vars, function(v) {
        if(v %in% colnames(X_sc)) {
          x_vals <- seq(min(X_sc[, v]), max(X_sc[, v]), length.out = 30)
          pred_grid <- sapply(x_vals, function(val) {
            tmp_data <- X_sc
            tmp_data[, v] <- val
            tmp_mat <- model.matrix(~ ., data = tmp_data)[, -1, drop = FALSE]
            mean(predict(lgb_target_obj, tmp_mat))
          })
          data.frame(x = x_vals, y = pred_grid)
        }
      })
      names(dep_effects) <- cont_vars
    }
  }, silent = TRUE)
  
  list(
    model_name = model_type, 
    subset = subset_name, 
    metrics = metrics, 
    fit = model_obj, 
    dependency_analysis = dep_effects,
    sd_intercept = sd_int,
    icc_latent = icc_lat
  )
}

run_loho_glm <- function(X, Y, Z, cont_vars, test_hospital) {
  test_idx <- which(Z == test_hospital)
  train_idx <- which(Z != test_hospital)
  X_train <- X[train_idx, , drop=FALSE]
  Y_train <- Y[train_idx]
  X_val <- X[test_idx, , drop=FALSE]
  X_train_sc <- X_train
  X_val_sc <- X_val
  is_num <- colnames(X) %in% cont_vars
  if (any(is_num)) {
    sc_p <- scale(X_train[, is_num, drop=FALSE])
    X_train_sc[, is_num] <- sc_p
    X_val_sc[, is_num] <- scale(X_val[, is_num, drop=FALSE], center = attr(sc_p, "scaled:center"), scale = attr(sc_p, "scaled:scale"))
  }
  fit <- glm(Y_train ~ ., data = data.frame(Y_train = Y_train, X_train_sc), family = binomial())
  as.numeric(predict(fit, newdata = X_val_sc, type = "response"))
}

run_loho_glmm <- function(X, Y, Z, cont_vars, test_hospital) {
  test_idx <- which(Z == test_hospital)
  train_idx <- which(Z != test_hospital)
  X_train <- X[train_idx, , drop=FALSE]
  Y_train <- Y[train_idx]
  Z_train <- Z[train_idx]
  X_val <- X[test_idx, , drop=FALSE]
  Z_val <- Z[test_idx]
  X_train_sc <- X_train
  X_val_sc <- X_val
  is_num <- colnames(X) %in% cont_vars
  if (any(is_num)) {
    sc_p <- scale(X_train[, is_num, drop=FALSE])
    X_train_sc[, is_num] <- sc_p
    X_val_sc[, is_num] <- scale(X_val[, is_num, drop=FALSE], center = attr(sc_p, "scaled:center"), scale = attr(sc_p, "scaled:scale"))
  }
  covs <- colnames(X)
  fmla <- as.formula(paste("Y ~", paste(covs, collapse = " + "), "+ (1 | H)"))
  td <- data.frame(X_train_sc, Y = Y_train, H = Z_train)
  vd <- data.frame(X_val_sc, H = Z_val)
  fit <- fit_glmer_safe(fmla, td)
  if (!is.null(fit)) {
    as.numeric(predict(fit, newdata = vd, type = "response", allow.new.levels = TRUE))
  } else {
    rep(NA, length(test_idx))
  }
}

run_loho_rf <- function(X, Y, Z, cont_vars, test_hospital, hp) {
  test_idx <- which(Z == test_hospital)
  train_idx <- which(Z != test_hospital)
  X_train <- X[train_idx, , drop=FALSE]
  Y_train <- Y[train_idx]
  X_val <- X[test_idx, , drop=FALSE]
  X_train_sc <- X_train
  X_val_sc <- X_val
  is_num <- colnames(X) %in% cont_vars
  if (any(is_num)) {
    sc_p <- scale(X_train[, is_num, drop=FALSE])
    X_train_sc[, is_num] <- sc_p
    X_val_sc[, is_num] <- scale(X_val[, is_num, drop=FALSE], center = attr(sc_p, "scaled:center"), scale = attr(sc_p, "scaled:scale"))
  }
  fit <- randomForest(X_train_sc, as.factor(Y_train), ntree = hp$ntree, nodesize = hp$nodesize)
  as.numeric(predict(fit, X_val_sc, type = "prob")[, "1"])
}

run_loho_xgb <- function(X_mat, Y, Z, cont_vars, test_hospital, hp) {
  test_idx <- which(Z == test_hospital)
  train_idx <- which(Z != test_hospital)
  X_train <- X_mat[train_idx, , drop = FALSE]
  Y_train <- Y[train_idx]
  X_val <- X_mat[test_idx, , drop = FALSE]
  cols_scale <- intersect(colnames(X_mat), cont_vars)
  if (length(cols_scale) > 0) {
    sc_p <- scale(X_train[, cols_scale, drop = FALSE])
    X_train[, cols_scale] <- sc_p
    X_val[, cols_scale] <- scale(X_val[, cols_scale, drop = FALSE], center = attr(sc_p, "scaled:center"), scale = attr(sc_p, "scaled:scale"))
  }
  fit <- xgb.train(list(objective = "binary:logistic", max_depth = hp$max_depth, eta = hp$eta, lambda = hp$lambda, subsample = hp$subsample), xgb.DMatrix(X_train, label = Y_train), hp$nrounds, verbose = 0)
  as.numeric(predict(fit, xgb.DMatrix(X_val)))
}

run_loho_lgb <- function(X_mat, Y, Z, cont_vars, test_hospital, hp) {
  test_idx <- which(Z == test_hospital)
  train_idx <- which(Z != test_hospital)
  X_train <- X_mat[train_idx, , drop = FALSE]
  Y_train <- Y[train_idx]
  X_val <- X_mat[test_idx, , drop = FALSE]
  cols_scale <- intersect(colnames(X_mat), cont_vars)
  if (length(cols_scale) > 0) {
    sc_p <- scale(X_train[, cols_scale, drop = FALSE])
    X_train[, cols_scale] <- sc_p
    X_val[, cols_scale] <- scale(X_val[, cols_scale, drop = FALSE], center = attr(sc_p, "scaled:center"), scale = attr(sc_p, "scaled:scale"))
  }
  fit <- lgb.train(list(objective = "binary", metric = "auc", max_depth = hp$max_depth, learning_rate = hp$learning_rate, min_data_in_leaf = hp$min_data_in_leaf, verbose = -1), lgb.Dataset(X_train, label = Y_train), hp$nrounds)
  as.numeric(predict(fit, X_val))
}



run_loho_gmerf <- function(X, Y, Z, cont_vars, test_hospital, hp) {
  test_idx <- which(Z == test_hospital)
  train_idx <- which(Z != test_hospital)
  X_train <- X[train_idx, , drop=FALSE]
  Y_train <- Y[train_idx]
  Z_train <- as.character(Z[train_idx])
  X_val <- X[test_idx, , drop=FALSE]
  Z_val <- as.character(Z[test_idx])
  
  X_train_sc <- X_train
  X_val_sc <- X_val
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
    target <- eta_static - as.numeric(b[Z_train])
    
    fitrf <- randomForest(X_train_sc, target, ntree = hp$ntree, nodesize = hp$nodesize)
    fitglmer <- fit_glmer_safe(Y_train ~ 1 + (1 | Z_train), data.frame(Y_train = Y_train, Z_train = Z_train), offset_vec = as.numeric(fitrf$predicted))
    
    if (is.null(fitglmer)) break
    b_new <- as.numeric(ranef(fitglmer)$Z_train[, 1]); names(b_new) <- rownames(ranef(fitglmer)$Z_train)
    b <- b_new
    
    if (max(abs(b - b_old)) < hp$toll) converged <- TRUE
  }
  
  if (is.null(fitrf)) return(rep(NA_real_, length(test_idx)))
  b_val <- as.numeric(b[Z_val]); b_val[is.na(b_val)] <- 0
  lk <- as.numeric(predict(fitrf, X_val_sc)) + b_val
  as.numeric(exp(lk) / (1 + exp(lk)))
}

run_loho_mexgb <- function(X_mat, Y, Z, cont_vars, test_hospital, hp) {
  test_idx <- which(Z == test_hospital)
  train_idx <- which(Z != test_hospital)
  X_train <- X_mat[train_idx, , drop = FALSE]
  Y_train <- Y[train_idx]
  Z_train <- as.character(Z[train_idx])
  X_val <- X_mat[test_idx, , drop = FALSE]
  Z_val <- as.character(Z[test_idx])
  
  cols_scale <- intersect(colnames(X_mat), cont_vars)
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
    target <- eta_static - as.numeric(b[Z_train])
    
    fit <- xgb.train(list(objective = "reg:squarederror", max_depth = hp$max_depth, eta = hp$eta, lambda = hp$lambda, subsample = hp$subsample), xgb.DMatrix(X_train, label = target), hp$nrounds, verbose = 0)
    f_xgb <- as.numeric(predict(fit, xgb.DMatrix(X_train)))
    fitglmer <- fit_glmer_safe(Y_train ~ 1 + (1 | Z_train), data.frame(Y_train = Y_train, Z_train = Z_train), offset_vec = f_xgb)
    
    if (is.null(fitglmer)) break
    b_new <- as.numeric(ranef(fitglmer)$Z_train[, 1]); names(b_new) <- rownames(ranef(fitglmer)$Z_train)
    b <- b_new
    
    if (max(abs(b - b_old)) < hp$toll) converged <- TRUE
  }
  
  if (is.null(fit)) return(rep(NA_real_, length(test_idx)))
  b_val <- as.numeric(b[Z_val]); b_val[is.na(b_val)] <- 0
  lk <- as.numeric(predict(fit, xgb.DMatrix(X_val))) + b_val
  as.numeric(exp(lk) / (1 + exp(lk)))
}

run_loho_melgbm <- function(X_mat, Y, Z, cont_vars, test_hospital, hp) {
  test_idx <- which(Z == test_hospital)
  train_idx <- which(Z != test_hospital)
  X_train <- X_mat[train_idx, , drop = FALSE]
  Y_train <- Y[train_idx]
  Z_train <- as.character(Z[train_idx])
  X_val <- X_mat[test_idx, , drop = FALSE]
  Z_val <- as.character(Z[test_idx])
  
  cols_scale <- intersect(colnames(X_mat), cont_vars)
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
    target <- eta_static - as.numeric(b[Z_train])
    
    fit <- lgb.train(list(objective = "regression", metric = "l2", max_depth = hp$max_depth, learning_rate = hp$learning_rate, min_data_in_leaf = hp$min_data_in_leaf, verbose = -1), lgb.Dataset(X_train, label = target), hp$nrounds)
    f_lgb <- as.numeric(predict(fit, X_train))
    fitglmer <- fit_glmer_safe(Y_train ~ 1 + (1 | Z_train), data.frame(Y_train = Y_train, Z_train = Z_train), offset_vec = f_lgb)
    
    if (is.null(fitglmer)) break
    b_new <- as.numeric(ranef(fitglmer)$Z_train[, 1]); names(b_new) <- rownames(ranef(fitglmer)$Z_train)
    b <- b_new
    
    if (max(abs(b - b_old)) < hp$toll) converged <- TRUE
  }
  
  if (is.null(fit)) return(rep(NA_real_, length(test_idx)))
  b_val <- as.numeric(b[Z_val]); b_val[is.na(b_val)] <- 0
  lk <- as.numeric(predict(fit, X_val)) + b_val
  as.numeric(exp(lk) / (1 + exp(lk)))
}




evaluate_full_model <- function(model_type, subset_name, subset_info, data_df) {
  covs <- subset_info$covs
  cont_vars <- subset_info$cont
  
  df_model <- data_df[, c(covs, "OS2ys_num", "hospital"), drop = FALSE]
  df_model <- df_model[complete.cases(df_model), ]
  
  X <- df_model[, covs, drop = FALSE]
  Y <- df_model$OS2ys_num
  Z <- df_model$hospital
  
  X_sc <- X
  is_num <- colnames(X) %in% cont_vars
  if (any(is_num)) {
    X_sc[, is_num] <- scale(X[, is_num])
  }
  
  preds <- NULL
  model_obj <- NULL
  
  if (model_type == "GLM") {
    fit <- glm(Y ~ ., data = data.frame(Y = Y, X_sc), family = binomial())
    preds <- predict(fit, type = "response")
    model_obj <- fit
    
  } else if (model_type == "GLMM") {
    fmla <- as.formula(paste("Y ~", paste(colnames(X_sc), collapse = " + "), "+ (1 | Z)"))
    fit <- glmer(fmla, data = data.frame(Y = Y, X_sc, Z = Z), family = binomial())
    preds <- predict(fit, type = "response", allow.new.levels = TRUE)
    model_obj <- fit
    
  } else if (model_type == "RF") {
    fit <- randomForest(X_sc, as.factor(Y), ntree = hp_rf$ntree, nodesize = hp_rf$nodesize, importance = TRUE)
    preds <- predict(fit, X_sc, type = "prob")[, "1"]
    model_obj <- fit
    
  } else if (model_type == "XGB") {
    X_mat <- model.matrix(~ ., data = X_sc)[, -1, drop = FALSE]
    fit <- xgb.train(
      list(objective = "binary:logistic", max_depth = hp_xgb$max_depth, eta = hp_xgb$eta, lambda = hp_xgb$lambda, subsample = hp_xgb$subsample),
      xgb.DMatrix(X_mat, label = Y), hp_xgb$nrounds, verbose = 0
    )
    preds <- as.numeric(predict(fit, xgb.DMatrix(X_mat)))
    model_obj <- fit
    
  } else if (model_type == "LGBM") {
    X_mat <- model.matrix(~ ., data = X_sc)[, -1, drop = FALSE]
    fit <- lgb.train(
      list(objective = "binary", metric = "auc", max_depth = hp_lgb$max_depth, learning_rate = hp_lgb$learning_rate, min_data_in_leaf = hp_lgb$min_data_in_leaf, verbose = -1),
      lgb.Dataset(X_mat, label = Y), hp_lgb$nrounds
    )
    preds <- as.numeric(predict(fit, X_mat))
    model_obj <- fit
    
  } else if (model_type == "GMERF") {
    eta <- as.numeric(glm(Y ~ ., data = data.frame(Y = Y, X_sc), family = binomial())$linear.predictors)
    b <- rep(0, length(unique(Z))); names(b) <- unique(Z)
    converged <- FALSE; iter <- 0; fitrf <- NULL; fitglmer <- NULL
    
    while (!converged && iter < hp_gmerf$itmax) {
      iter <- iter + 1; b_old <- b
      target <- as.numeric(eta) - as.numeric(b[as.character(Z)])
      fitrf <- randomForest(X_sc, target, ntree = hp_gmerf$ntree, nodesize = hp_gmerf$nodesize)
      
      fitglmer <- tryCatch({
        suppressWarnings(glmer(Y ~ 1 + (1 | Z), data = data.frame(Y = Y, Z = Z), family = binomial(), offset = as.numeric(fitrf$predicted)))
      }, error = function(e) NULL)
      
      if (is.null(fitglmer)) break
      b_new <- as.numeric(ranef(fitglmer)$Z[, 1]); names(b_new) <- rownames(ranef(fitglmer)$Z)
      b <- b_new
      if (max(abs(b - b_old)) < hp_gmerf$toll) converged <- TRUE
    }
    
    if (!is.null(fitrf)) {
      b_val <- as.numeric(b[as.character(Z)]); b_val[is.na(b_val)] <- 0
      lk <- as.numeric(predict(fitrf, X_sc)) + b_val
      preds <- as.numeric(exp(lk) / (1 + exp(lk)))
      model_obj <- list(rf = fitrf, glmer = fitglmer, random_effects = b, iterations = iter)
    }
    
  } else if (model_type == "MEXGB") {
    X_mat <- model.matrix(~ ., data = X_sc)[, -1, drop = FALSE]
    tmp_df <- data.frame(Y = Y, X_sc)
    eta_link <- as.numeric(glm(Y ~ ., data = tmp_df, family = binomial())$linear.predictors)
    b <- rep(0, length(unique(Z))); names(b) <- unique(Z)
    converged <- FALSE; iter <- 0; fit <- NULL; fitglmer <- NULL
    
    while (!converged && iter < hp_mexgb$itmax) {
      iter <- iter + 1; b_old <- b
      target <- as.numeric(eta_link) - as.numeric(b[as.character(Z)])
      fit <- xgb.train(
        list(objective = "reg:squarederror", max_depth = hp_mexgb$max_depth, eta = hp_mexgb$eta, lambda = hp_mexgb$lambda, subsample = hp_mexgb$subsample),
        xgb.DMatrix(X_mat, label = target), hp_mexgb$nrounds, verbose = 0
      )
      f_xgb <- as.numeric(predict(fit, xgb.DMatrix(X_mat)))
      
      fitglmer <- tryCatch({
        suppressWarnings(glmer(Y ~ 1 + (1 | Z), data = data.frame(Y = Y, Z = Z), family = binomial(), offset = f_xgb))
      }, error = function(e) NULL)
      
      if (is.null(fitglmer)) break
      b_new <- as.numeric(ranef(fitglmer)$Z[, 1]); names(b_new) <- rownames(ranef(fitglmer)$Z)
      b <- b_new
      if (max(abs(b - b_old)) < hp_mexgb$toll) converged <- TRUE
    }
    
    if (!is.null(fit)) {
      b_val <- as.numeric(b[as.character(Z)]); b_val[is.na(b_val)] <- 0
      lk <- as.numeric(predict(fit, xgb.DMatrix(model.matrix(~ ., data = X_sc)[, -1, drop = FALSE]))) + b_val
      preds <- as.numeric(exp(lk) / (1 + exp(lk)))
      model_obj <- list(xgb = fit, glmer = fitglmer, random_effects = b, iterations = iter)
    }
    
  } else if (model_type == "MELGBM") {
    X_mat <- model.matrix(~ ., data = X_sc)[, -1, drop = FALSE]
    tmp_df <- data.frame(Y = Y, X_sc)
    eta_link <- as.numeric(glm(Y ~ ., data = tmp_df, family = binomial())$linear.predictors)
    b <- rep(0, length(unique(Z))); names(b) <- unique(Z)
    converged <- FALSE; iter <- 0; fit <- NULL; fitglmer <- NULL
    
    while (!converged && iter < hp_melgbm$itmax) {
      iter <- iter + 1; b_old <- b
      target <- as.numeric(eta_link) - as.numeric(b[as.character(Z)])
      fit <- lgb.train(
        list(objective = "regression", metric = "l2", max_depth = hp_melgbm$max_depth, learning_rate = hp_melgbm$learning_rate, min_data_in_leaf = hp_melgbm$min_data_in_leaf, verbose = -1),
        lgb.Dataset(X_mat, label = target), hp_melgbm$nrounds
      )
      f_lgb <- as.numeric(predict(fit, X_mat))
      
      fitglmer <- tryCatch({
        suppressWarnings(glmer(Y ~ 1 + (1 | Z), data = data.frame(Y = Y, Z = Z), family = binomial(), offset = f_lgb))
      }, error = function(e) NULL)
      
      if (is.null(fitglmer)) break
      b_new <- as.numeric(ranef(fitglmer)$Z[, 1]); names(b_new) <- rownames(ranef(fitglmer)$Z)
      b <- b_new
      if (max(abs(b - b_old)) < hp_melgbm$toll) converged <- TRUE
    }
    
    if (!is.null(fit)) {
      b_val <- as.numeric(b[as.character(Z)]); b_val[is.na(b_val)] <- 0
      lk <- as.numeric(predict(fit, model.matrix(~ ., data = X_sc)[, -1, drop = FALSE])) + b_val
      preds <- as.numeric(exp(lk) / (1 + exp(lk)))
      model_obj <- list(lgb = fit, glmer = fitglmer, random_effects = b, iterations = iter)
    }
  }
  
  metrics <- calc_metrics_vec(preds, Y)
  
  glmer_obj <- if (model_type %in% c("GMERF", "MEXGB", "MELGBM")) model_obj$glmer else if (model_type == "GLMM") model_obj else NULL
  
  sd_int <- NA
  icc_lat <- NA
  if (!is.null(glmer_obj)) {
    vc <- VarCorr(glmer_obj)
    if (length(vc) > 0) {
      varianza_intercept <- as.numeric(vc[[1]][1, 1])
      sd_int <- sqrt(varianza_intercept)
      icc_lat <- varianza_intercept / (varianza_intercept + (pi^2 / 3))
    }
  }
  
  dep_effects <- NULL
  try({
    if (model_type %in% c("GLM", "GLMM")) {
      dep_effects <- summary(model_obj)$coefficients
      
    } else if (model_type %in% c("RF", "GMERF")) {
      rf_target_obj <- if (model_type == "GMERF") model_obj$rf else model_obj
      dep_effects <- lapply(cont_vars, function(v) {
        if (v %in% colnames(X_sc)) {
          x_vals <- seq(min(X_sc[, v]), max(X_sc[, v]), length.out = 500)
          pred_grid <- sapply(x_vals, function(val) {
            tmp_data <- X_sc
            tmp_data[, v] <- val
            if (model_type == "RF") {
              mean(predict(rf_target_obj, tmp_data, type = "prob")[, "1"])
            } else {
              mean(predict(rf_target_obj, tmp_data))
            }
          })
          data.frame(x = x_vals, y = pred_grid)
        }
      })
      names(dep_effects) <- cont_vars
      
    } else if (model_type %in% c("XGB", "MEXGB")) {
      xgb_target_obj <- if (model_type == "MEXGB") model_obj$xgb else model_obj
      dep_effects <- lapply(cont_vars, function(v) {
        if (v %in% colnames(X_sc)) {
          x_vals <- seq(min(X_sc[, v]), max(X_sc[, v]), length.out = 500)
          pred_grid <- sapply(x_vals, function(val) {
            tmp_data <- X_sc
            tmp_data[, v] <- val
            tmp_mat <- model.matrix(~ ., data = tmp_data)[, -1, drop = FALSE]
            mean(predict(xgb_target_obj, xgb.DMatrix(tmp_mat)))
          })
          data.frame(x = x_vals, y = pred_grid)
        }
      })
      names(dep_effects) <- cont_vars
      
    } else if (model_type %in% c("LGBM", "MELGBM")) {
      lgb_target_obj <- if (model_type == "MELGBM") model_obj$lgb else model_obj
      dep_effects <- lapply(cont_vars, function(v) {
        if (v %in% colnames(X_sc)) {
          x_vals <- seq(min(X_sc[, v]), max(X_sc[, v]), length.out = 500)
          pred_grid <- sapply(x_vals, function(val) {
            tmp_data <- X_sc
            tmp_data[, v] <- val
            tmp_mat <- model.matrix(~ ., data = tmp_data)[, -1, drop = FALSE]
            mean(predict(lgb_target_obj, tmp_mat))
          })
          data.frame(x = x_vals, y = pred_grid)
        }
      })
      names(dep_effects) <- cont_vars
    }
  }, silent = TRUE)
  
  list(
    model_name = model_type, 
    subset = subset_name, 
    metrics = metrics, 
    fit = model_obj, 
    dependency_analysis = dep_effects,
    sd_intercept = sd_int,
    icc_latent = icc_lat
  )
}

analyze_me_trees <- function(subset_name, subset_info, data_df, base_model = "RF", mixed_model = "GMERF") {
  covs <- subset_info$covs
  cont_vars <- subset_info$cont
  
  df_model <- data_df[, c(covs, "OS2ys_num", "hospital"), drop = FALSE]
  df_model <- df_model[complete.cases(df_model), ]
  
  X_full <- df_model[, covs, drop = FALSE]
  Y_full <- df_model$OS2ys_num
  Z_full <- df_model$hospital
  
  cat("\n--- PUNTO 1: Distribuzione Pazienti per Ospedale ---\n")
  print(table(Z_full))
  
  cat(paste("\n--- Training", base_model, "e", mixed_model, "sul dataset completo... ---\n"))
  res_base <- evaluate_full_model(base_model, subset_name, subset_info, data_df)
  res_mixed <- evaluate_full_model(mixed_model, subset_name, subset_info, data_df)
  
  cat("\n--- PUNTO 2: Intervalli di Confidenza Random Effects (estratti da", mixed_model, ") ---\n")
  glmer_obj <- res_mixed$fit$glmer
  if (!is.null(glmer_obj)) {
    p <- lattice::dotplot(ranef(glmer_obj, condVar = TRUE))
    p$main <- "95% CI ran eff"
    print(p)
  } else {
    cat("Nessun random effect stimato (Varianza = 0 o modello non convergente).\n")
  }
  
  cat(paste("\n--- PUNTO 3: Shift delle Curve di Dipendenza (", base_model, "vs", mixed_model, ") ---\n"))
  valid_cont <- intersect(cont_vars, covs)
  if (length(valid_cont) > 0) {
    par(mfrow = c(2, length(valid_cont)))
    for (v in valid_cont) {
      if (!is.null(res_base$dependency_analysis[[v]])) {
        plot(res_base$dependency_analysis[[v]]$x, res_base$dependency_analysis[[v]]$y, type = "l", col = "blue", lwd = 2, main = paste(base_model, ":", v), xlab = v, ylab = "Marginal Effect")
      }
    }
    for (v in valid_cont) {
      if (!is.null(res_mixed$dependency_analysis[[v]])) {
        plot(res_mixed$dependency_analysis[[v]]$x, res_mixed$dependency_analysis[[v]]$y, type = "l", col = "red", lwd = 2, main = paste(mixed_model, ":", v), xlab = v, ylab = "Marginal Effect")
      }
    }
    par(mfrow = c(1, 1))
  }
  
  cat("\n--- PUNTO 4: Stabilità Random Effects in 5-Fold CV ---\n")
  set.seed(120)
  folds <- caret::createFolds(as.factor(Y_full), k = 5, list = TRUE)
  is_num <- colnames(X_full) %in% cont_vars
  
  b_mat <- sapply(folds, function(te) {
    X_train <- X_full[-te, , drop = FALSE]
    Y_train <- Y_full[-te]
    Z_train <- factor(Z_full[-te])
    
    X_train_sc <- X_train
    if (any(is_num)) {
      X_train_sc[, is_num] <- scale(X_train[, is_num])
    }
    
    tmp_df <- data.frame(Y_train = Y_train, X_train_sc)
    eta <- as.numeric(glm(Y_train ~ ., data = tmp_df, family = binomial())$linear.predictors)
    b <- rep(0, length(unique(Z_train)))
    names(b) <- unique(Z_train)
    converged <- FALSE
    iter <- 0
    
    while (!converged && iter < hp_gmerf$itmax) {
      iter <- iter + 1
      b_old <- b
      target <- as.numeric(eta) - as.numeric(b[as.character(Z_train)])
      fitrf <- randomForest(X_train_sc, target, ntree = hp_gmerf$ntree, nodesize = hp_gmerf$nodesize)
      
      fitglmer <- tryCatch({
        suppressWarnings(glmer(Y_train ~ 1 + (1 | Z_train), data = data.frame(Y_train = Y_train, Z_train = Z_train), family = binomial(), offset = as.numeric(fitrf$predicted)))
      }, error = function(e) NULL)
      
      if (is.null(fitglmer)) break
      b_new <- as.numeric(ranef(fitglmer)$Z_train[, 1])
      names(b_new) <- rownames(ranef(fitglmer)$Z_train)
      b <- b_new
      if (max(abs(b - b_old)) < hp_gmerf$toll) converged <- TRUE
    }
    
    out_b <- rep(NA, length(levels(factor(Z_full))))
    names(out_b) <- levels(factor(Z_full))
    out_b[names(b)] <- b
    return(out_b)
  })
  
  colnames(b_mat) <- paste0("Fold", 1:5)
  print(round(b_mat, 3))
  
  cat("\n--- PUNTO 5: Test LOHO (Leave-One-Hospital-Out) ---\n")
  loho_base <- rep(NA_real_, length(Y_full))
  loho_mixed <- rep(NA_real_, length(Y_full))
  
  X_mat_full <- model.matrix(~ ., data = X_full)[, -1, drop = FALSE]
  
  for (h in levels(factor(Z_full))) {
    te <- which(Z_full == h)
    if (length(te) == 0) next
    
    if (base_model == "RF") {
      loho_base[te] <- run_loho_rf(X_full, Y_full, Z_full, cont_vars, h, hp_rf)
    } else if (base_model == "XGB") {
      loho_base[te] <- run_loho_xgb(X_mat_full, Y_full, Z_full, cont_vars, h, hp_xgb)
    } else if (base_model == "LGBM") {
      loho_base[te] <- run_loho_lgb(X_mat_full, Y_full, Z_full, cont_vars, h, hp_lgb)
    }
    
    if (mixed_model == "GMERF") {
      loho_mixed[te] <- run_loho_gmerf(X_full, Y_full, Z_full, cont_vars, h, hp_gmerf)
    } else if (mixed_model == "MEXGB") {
      loho_mixed[te] <- run_loho_mexgb(X_mat_full, Y_full, Z_full, cont_vars, h, hp_mexgb)
    } else if (mixed_model == "MELGBM") {
      loho_mixed[te] <- run_loho_melgbm(X_mat_full, Y_full, Z_full, cont_vars, h, hp_melgbm)
    }
  }
  
  cat(paste("\nMetriche Globali LOHO -", base_model, ":\n"))
  print(round(calc_metrics_vec(loho_base, Y_full), 4))
  
  cat(paste("\nMetriche Globali LOHO -", mixed_model, ":\n"))
  print(round(calc_metrics_vec(loho_mixed, Y_full), 4))
  
  cat(paste("\nMetriche LOHO per singolo ospedale -", mixed_model, ":\n"))
  metrics_per_h <- do.call(rbind, lapply(levels(factor(Z_full)), function(h) {
    i <- which(Z_full == h)
    calc_metrics_vec(loho_mixed[i], Y_full[i])
  }))
  rownames(metrics_per_h) <- levels(factor(Z_full))
  print(round(metrics_per_h, 4))
}

subset_list <- list(
  
  GS5 = list(covs = c("Age","Volume_PTV","schedule_ok_chemo","IMMUNO","riskfactors1","riskfactors2"),
             cont = c("Age", "Volume_PTV")),
  
  GS7 = list(covs = c("Age","KPS", "Volume_PTV","schedule_ok_chemo","IMMUNO","riskfactors2"),
             cont = c("Age", "Volume_PTV")),
  
  Lspecific2a = list(covs = c("Age","Volume_PTV","schedule_ok_chemo","IMMUNO","riskfactors1","riskfactors2","Subcutaneous_Fat","Skeletal_Muscle"),
                     cont = c("Age","Volume_PTV", "Subcutaneous_Fat", "Skeletal_Muscle")),
  Lassospecific1 = list(covs = c("Age","Volume_PTV","IMMUNO","schedule_ok_chemo", "Subcutaneous_Fat","Skeletal_Muscle", "Intermuscular_Fat", "diameter.kidney.sep..mm.", "aortic_implant_present", "aorta_agatston", "riskfactors1", "riskfactors2"),
                        cont = c("Age","Volume_PTV","Subcutaneous_Fat","Skeletal_Muscle", "Intermuscular_Fat","diameter.kidney.sep..mm.",  "aorta_agatston")),
  Bspecific2b_2 = list(covs = c("Age","Volume_PTV","IMMUNO","LAA","Subcutaneous_Fat","Skeletal_Muscle", "diameter.kidney.sep..mm."),
                       cont = c("Age","Volume_PTV","Subcutaneous_Fat", "LAA", "Skeletal_Muscle", "diameter.kidney.sep..mm.")),
  Bhole = list(covs = c("Age","Volume_PTV","LAA", "riskfactors1"),
               cont = c("Age","Volume_PTV","LAA")),
  Bspecific2a = list(covs = c("Age","Volume_PTV","schedule_ok_chemo","IMMUNO","riskfactors1","riskfactors2", "LAA","Subcutaneous_Fat", "Skeletal_Muscle", "diameter.kidney.sep..mm."),
                     cont = c("Age","Volume_PTV","Subcutaneous_Fat", "LAA", "Skeletal_Muscle", "diameter.kidney.sep..mm.")),
  Lhole_2 = list(covs = c("Age","KPS", "Volume_PTV","IMMUNO","eud_lungs","mean_dose_atriumL", "Subcutaneous_Fat", "riskfactors2"),
                 cont = c("Age", "Volume_PTV","eud_lungs","mean_dose_atriumL", "Subcutaneous_Fat")),
  Lassohole_2 = list(covs = c("Age","KPS", "Volume_PTV","IMMUNO","schedule_ok_chemo", "eud_lungs", "mean_dose_atriumL", "Subcutaneous_Fat","aorta_agatston", "aortic_implant_present", "riskfactors2"),
                     cont = c("Age", "Volume_PTV","eud_lungs","mean_dose_atriumL", "Subcutaneous_Fat","aorta_agatston")),
  Lspecific2b = list(covs = c("Age","Volume_PTV","IMMUNO","riskfactors1","Subcutaneous_Fat","Skeletal_Muscle"),
                     cont = c("Age","Volume_PTV", "Subcutaneous_Fat", "Skeletal_Muscle"))
  
)
evaluate_full_model("GMERF", "Lhole_2", subset_list[["Lhole_2"]], data)
analyze_me_trees("Lhole_2", subset_list[["Lhole_2"]], data, base_model = "RF", mixed_model = "GMERF") #metriche globali sono pooled
