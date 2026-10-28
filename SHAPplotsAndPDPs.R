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
library(grid)
library(survival)
library(survminer)
library(uwot)
library(dbscan)
library(mclust)
library(treeshap)
library(rpart)
library(tidyr)
library(dplyr)

data <- read.csv('dataclean_final_withNA.csv', header=T)

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


#corr(fat) ----
plot(data$Subcutaneous_Fat, data$Skeletal_Muscle, 
     pch = 19, col = "steelblue",
     main = paste("Spearman corr =", round(cor(data$Subcutaneous_Fat, data$Skeletal_Muscle, method="spearman", use="complete.obs"), 3)))

plot(data$Subcutaneous_Fat, data$Torso_Fat, 
     pch = 19, col = "steelblue",
     main = paste("Spearman rho =", round(cor(data$Subcutaneous_Fat, data$Torso_Fat, method="spearman", use="complete.obs"), 3)))

#shap plots----
hp_xgb <- list(max_depth = 2, eta = 0.05, lambda = 5, subsample = 0.7, nrounds = 100)
seed_xgb <- 123


covariates <- c("Age", "Volume_PTV", "LAA", "riskfactors1")
cont_vars  <- c("Age", "Volume_PTV", "LAA")

keep <- complete.cases(data[, covariates]) & !is.na(data$OS2ys_num) & !is.na(data$hospital)
df_sub <- data[keep, ]
X_full <- df_sub[, covariates, drop = FALSE]
Y_full <- df_sub$OS2ys_num
Z_full <- droplevels(df_sub$hospital)

X_full_sc <- X_full
X_full_sc[, cont_vars] <- scale(X_full[, cont_vars])

cat("Pazients:", nrow(X_full_sc), "| Events OS2ys=1:", sum(Y_full), "\n")

xgb_params <- function(hp, seed) {
  list(objective = "binary:logistic", max_depth = hp$max_depth, eta = hp$eta,
       lambda = hp$lambda, subsample = hp$subsample, seed = seed, nthread = 1)
}

make_col_map <- function(X_df) {
  mm <- model.matrix(~ ., data = X_df)
  asg <- attr(mm, "assign")[-1]
  data.frame(col = colnames(mm)[-1], var = colnames(X_df)[asg], stringsAsFactors = FALSE)
}

aggregate_shap <- function(shap_matrix, col_map, original_covs) {
  stopifnot(all(colnames(shap_matrix) %in% col_map$col))
  out <- sapply(original_covs, function(v) {
    cc <- col_map$col[col_map$var == v]
    rowSums(shap_matrix[, cc, drop = FALSE])
  })
  out <- as.matrix(out)
  colnames(out) <- original_covs
  out
}

col_map <- make_col_map(X_full_sc)
X_mat <- model.matrix(~ ., data = X_full_sc)[, -1, drop = FALSE]

set.seed(seed_xgb)
fit_xgb <- xgb.train(
  params  = xgb_params(hp_xgb, seed_xgb),
  data    = xgb.DMatrix(X_mat, label = Y_full),
  nrounds = hp_xgb$nrounds,
  verbose = 0
)
p_final <- as.numeric(predict(fit_xgb, xgb.DMatrix(X_mat)))

uni_xgb <- xgboost.unify(fit_xgb, X_mat)
treeshap_res <- treeshap(uni_xgb, X_mat, verbose = FALSE)

shap_raw <- as.matrix(treeshap_res$shaps)
shap_agg <- aggregate_shap(shap_raw, col_map, covariates)
shap_4vars <- shap_agg

n_repeats <- 10
shap_oof_arr <- array(NA_real_, dim = c(nrow(X_full), ncol(shap_raw), n_repeats))
for (r in seq_len(n_repeats)) {
  set.seed(r)
  folds_r <- createFolds(as.factor(Y_full), k = 5, list = TRUE, returnTrain = FALSE)
  for (te in folds_r) {
    Xtr <- X_full[-te, , drop = FALSE]; Xte <- X_full[te, , drop = FALSE]
    sc_k <- scale(Xtr[, cont_vars])
    Xtr[, cont_vars] <- sc_k
    Xte[, cont_vars] <- scale(Xte[, cont_vars], center = attr(sc_k, "scaled:center"), scale = attr(sc_k, "scaled:scale"))
    Xtr_m <- model.matrix(~ ., data = Xtr)[, -1, drop = FALSE]
    Xte_m <- model.matrix(~ ., data = Xte)[, -1, drop = FALSE]
    stopifnot(identical(colnames(Xtr_m), colnames(X_mat)), identical(colnames(Xte_m), colnames(X_mat)))
    
    f_k <- xgb.train(
      params = xgb_params(hp_xgb, r),
      data = xgb.DMatrix(Xtr_m, label = Y_full[-te]), nrounds = hp_xgb$nrounds, verbose = 0
    )
    sh <- as.matrix(treeshap(xgboost.unify(f_k, Xtr_m), Xte_m, verbose = FALSE)$shaps)
    shap_oof_arr[te, , r] <- sh[, colnames(X_mat), drop = FALSE]
  }
}
shap_oof_7cols <- apply(shap_oof_arr, c(1, 2), mean)
colnames(shap_oof_7cols) <- colnames(shap_raw)
shap_oof_agg <- aggregate_shap(shap_oof_7cols, col_map, covariates)
shap_oof_4vars <- shap_oof_agg

to_num <- function(v) {
  if (is.factor(v) || is.character(v)) {
    n <- suppressWarnings(as.numeric(as.character(v)))
    if (all(is.na(n) == is.na(v))) n else as.numeric(factor(v))
  } else {
    as.numeric(v)
  }
}

plot_shap_direction <- function(shap_mat, feat_obj, title_text) {
  shap_mat <- as.matrix(shap_mat)
  feat_df  <- as.data.frame(feat_obj)
  features <- colnames(shap_mat)
  
  if (!all(features %in% colnames(feat_df))) {
    stop("Colonne SHAP non presenti nelle feature: ",
         paste(setdiff(features, colnames(feat_df)), collapse = ", "))
  }
  stopifnot(nrow(shap_mat) == nrow(feat_df))
  
  plot_rows <- lapply(features, function(f) {
    v_num <- to_num(feat_df[[f]])
    rng <- range(v_num, na.rm = TRUE)
    v_scaled <- if (rng[1] == rng[2]) rep(0.5, length(v_num)) else (v_num - rng[1]) / (rng[2] - rng[1])
    data.frame(Feature = f, SHAP = as.numeric(shap_mat[, f]), Value_Scale = v_scaled, stringsAsFactors = FALSE)
  })
  plot_data <- do.call(rbind, plot_rows)
  
  imp <- tapply(abs(plot_data$SHAP), plot_data$Feature, mean, na.rm = TRUE)
  plot_data$Feature <- factor(plot_data$Feature, levels = names(sort(imp, decreasing = FALSE)))
  
  ggplot(plot_data, aes(x = SHAP, y = Feature, color = Value_Scale)) +
    geom_jitter(width = 0, height = 0.22, size = 1.6, alpha = 0.65) +
    scale_color_gradientn(
      colors = c("#0072B2", "#E69F00", "#D55E00"),
      limits = c(0, 1),
      breaks = c(0, 1),
      labels = c("Low", "High"),
      name   = "Feature\nvalue"
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.8) +
    theme_minimal(base_size = 12) +
    labs(title = title_text, x = "SHAP value (log-odds di OS2ys = YES)", y = NULL) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey90"),
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )
}

shap_direction_table <- function(shap_mat, feat_obj) {
  shap_mat <- as.matrix(shap_mat)
  feat_df  <- as.data.frame(feat_obj)
  do.call(rbind, lapply(colnames(shap_mat), function(f) {
    data.frame(
      Feature = f,
      Spearman_value_vs_SHAP = round(suppressWarnings(cor(to_num(feat_df[[f]]), shap_mat[, f], method = "spearman", use = "complete.obs")), 3),
      Mean_abs_SHAP = round(mean(abs(shap_mat[, f])), 4),
      stringsAsFactors = FALSE
    )
  }))
}

p_full_raw <- plot_shap_direction(shap_raw, X_mat, "SHAP summary: Full model Bhole (non aggregated)")
print(p_full_raw)

p_full_agg <- plot_shap_direction(shap_agg, X_full, "SHAP summary: Full model Bhole (aggregated)")
print(p_full_agg)


p_oof_raw <- plot_shap_direction(shap_oof_7cols, X_mat, "SHAP summary: OOF Model (non aggregated)")
print(p_oof_raw)

p_oof_agg <- plot_shap_direction(shap_oof_agg, X_full, "SHAP summary: OOF repeated CV Bhole (aggregated)")
print(p_oof_agg)

print(shap_direction_table(shap_agg, X_full), row.names = FALSE) 
print(shap_direction_table(shap_oof_agg, X_full), row.names = FALSE)
#spearman computed also for ordinated features (and binary too: IMMUNO means from NO to YES what happens to SHAP)
# NOT for schedule ok chemo







#PDP & SHAP dependence plot -----------------
library(patchwork)

hp_rf     <- list(ntree = 500, nodesize = 20)
hp_xgb    <- list(max_depth = 2, eta = 0.05, lambda = 5, subsample = 0.7, nrounds = 100)
hp_lgb    <- list(max_depth = 2, learning_rate = 0.05, min_data_in_leaf = 20, nrounds = 100)
hp_gmerf  <- list(ntree = 500, nodesize = 20, itmax = 15, toll = 0.05)
hp_mexgb  <- list(max_depth = 2, eta = 0.05, lambda = 5, subsample = 0.7, nrounds = 100, itmax = 15, toll = 0.05)
hp_melgbm <- list(max_depth = 2, learning_rate = 0.05, min_data_in_leaf = 20, nrounds = 100, itmax = 15, toll = 0.05)

calc_metrics_vec <- function(p, Y) {
  if(is.list(p)) p <- unlist(p)
  p <- as.numeric(p)
  valid_idx <- !is.na(p)
  if(sum(valid_idx) < 10 || length(unique(p[valid_idx])) < 2) {
    return(c(AUC = NA, Precision = NA, Recall = NA, F1 = NA, Cal_Intercept = NA, Cal_Slope = NA, Threshold = NA))
  }
  p_valid <- p[valid_idx]
  Y_valid <- Y[valid_idx]
  Y_num <- as.numeric(as.character(Y_valid))
  Y_factor <- factor(Y_valid, levels = c(0, 1))
  r <- pROC::roc(Y_valid, p_valid, quiet = TRUE)
  t <- 0.5
  pred_factor <- factor(ifelse(p_valid >= t, 1, 0), levels = c(0, 1))
  cm <- caret::confusionMatrix(pred_factor, Y_factor, positive = "1")
  eps <- 1e-15
  p_adj <- pmax(pmin(p_valid, 1 - eps), eps)
  logit_p <- log(p_adj / (1 - p_adj))
  c_slope <- NA
  c_int <- NA
  try({
    fit_calib <- glm(Y_num ~ logit_p, family = binomial())
    c_slope <- as.numeric(coef(fit_calib)[2])
    fit_calib_int <- glm(Y_num ~ 1 + offset(logit_p), family = binomial())
    c_int <- as.numeric(coef(fit_calib_int)[1])
  }, silent = TRUE)
  c(
    AUC = as.numeric(pROC::auc(r)),
    Precision = as.numeric(cm$byClass["Pos Pred Value"]),
    Recall = as.numeric(cm$byClass["Sensitivity"]),
    F1 = 2 * (cm$byClass["Pos Pred Value"] * cm$byClass["Sensitivity"]) / (cm$byClass["Pos Pred Value"] + cm$byClass["Sensitivity"]),
    Cal_Intercept = c_int,
    Cal_Slope = c_slope,
    Threshold = t
  )
}

fit_glmer_safe <- function(fmla, data, offset_vec = NULL) {
  tryCatch({
    if(is.null(offset_vec)) suppressWarnings(glmer(fmla, data = data, family = binomial())) else suppressWarnings(glmer(fmla, data = data, family = binomial(), offset = as.numeric(offset_vec)))
  }, error = function(e) NULL)
}

evaluate_full_model <- function(model_type, subset_name, subset_info, data_df) {
  covs <- subset_info$covs
  cont_vars <- subset_info$cont
  df_model <- data_df[, c(covs, "OS2ys_num", "hospital"), drop = FALSE]
  df_model <- df_model[complete.cases(df_model), ]
  
  for(v in covs) {
    if(!v %in% cont_vars) {
      df_model[[v]] <- as.factor(df_model[[v]])
    }
  }
  
  X_raw <- df_model[, covs, drop = FALSE]
  X_sc <- X_raw
  Y <- df_model$OS2ys_num
  Z <- df_model$hospital
  
  is_num <- colnames(X_sc) %in% cont_vars
  if(any(is_num)) X_sc[, is_num] <- scale(X_sc[, is_num])
  
  preds <- NULL
  model_obj <- NULL
  
  set.seed(123)
  
  if(model_type == "GLM") {
    fit <- glm(Y ~ ., data = data.frame(Y = Y, X_sc), family = binomial())
    preds <- predict(fit, type = "response")
    model_obj <- fit
  } else if(model_type == "GLMM") {
    fmla <- as.formula(paste("Y ~", paste(colnames(X_sc), collapse = " + "), "+ (1 | Z)"))
    fit <- glmer(fmla, data = data.frame(Y = Y, X_sc, Z = Z), family = binomial())
    preds <- predict(fit, type = "response", allow.new.levels = TRUE)
    model_obj <- fit
  } else if(model_type == "RF") {
    fit <- randomForest(X_sc, as.factor(Y), ntree = hp_rf$ntree, nodesize = hp_rf$nodesize)
    preds <- predict(fit, X_sc, type = "prob")[, "1"]
    model_obj <- fit
  } else if(model_type == "XGB") {
    X_mat <- model.matrix(~ ., data = X_sc)[, -1, drop = FALSE]
    fit <- xgb.train(list(objective = "binary:logistic", max_depth = hp_xgb$max_depth, eta = hp_xgb$eta, lambda = hp_xgb$lambda, subsample = hp_xgb$subsample), xgb.DMatrix(X_mat, label = Y), hp_xgb$nrounds, verbose = 0)
    preds <- as.numeric(predict(fit, xgb.DMatrix(X_mat)))
    model_obj <- fit
  } else if(model_type == "LGBM") {
    X_mat <- model.matrix(~ ., data = X_sc)[, -1, drop = FALSE]
    fit <- lgb.train(list(objective = "binary", metric = "auc", max_depth = hp_lgb$max_depth, learning_rate = hp_lgb$learning_rate, min_data_in_leaf = hp_lgb$min_data_in_leaf, verbose = -1), lgb.Dataset(X_mat, label = Y), hp_lgb$nrounds)
    preds <- as.numeric(predict(fit, X_mat))
    model_obj <- fit
  } else if(model_type == "GMERF") {
    eta_static <- as.numeric(glm(Y ~ ., data = data.frame(Y = Y, X_sc), family = binomial())$linear.predictors)
    b <- rep(0, length(unique(Z)))
    names(b) <- unique(Z)
    converged <- FALSE
    iter <- 0
    fitrf <- NULL
    fitglmer <- NULL
    while(!converged && iter < hp_gmerf$itmax) {
      iter <- iter + 1
      b_old <- b
      target <- as.numeric(eta_static) - as.numeric(b[as.character(Z)])
      fitrf <- randomForest(X_sc, target, ntree = hp_gmerf$ntree, nodesize = hp_gmerf$nodesize)
      fitglmer <- fit_glmer_safe(Y ~ 1 + (1 | Z), data.frame(Y = Y, Z = Z), offset_vec = as.numeric(fitrf$predicted))
      if(is.null(fitglmer)) break
      b_new <- as.numeric(ranef(fitglmer)$Z[, 1])
      names(b_new) <- rownames(ranef(fitglmer)$Z)
      b <- b_new
      if(max(abs(b - b_old)) < hp_gmerf$toll) converged <- TRUE
    }
    if(!is.null(fitrf)) {
      b_val <- as.numeric(b[as.character(Z)])
      b_val[is.na(b_val)] <- 0
      lk <- as.numeric(predict(fitrf, X_sc)) + b_val
      preds <- as.numeric(exp(lk) / (1 + exp(lk)))
      model_obj <- list(rf = fitrf, glmer = fitglmer)
    }
  } else if(model_type == "MEXGB") {
    X_mat <- model.matrix(~ ., data = X_sc)[, -1, drop = FALSE]
    eta_static <- as.numeric(glm(Y ~ ., data = data.frame(Y = Y, X_sc), family = binomial())$linear.predictors)
    b <- rep(0, length(unique(Z)))
    names(b) <- unique(Z)
    converged <- FALSE
    iter <- 0
    fit <- NULL
    fitglmer <- NULL
    while(!converged && iter < hp_mexgb$itmax) {
      iter <- iter + 1
      b_old <- b
      target <- as.numeric(eta_static) - as.numeric(b[as.character(Z)])
      fit <- xgb.train(list(objective = "reg:squarederror", max_depth = hp_mexgb$max_depth, eta = hp_mexgb$eta, lambda = hp_mexgb$lambda, subsample = hp_mexgb$subsample), xgb.DMatrix(X_mat, label = target), hp_mexgb$nrounds, verbose = 0)
      f_xgb <- as.numeric(predict(fit, xgb.DMatrix(X_mat)))
      fitglmer <- fit_glmer_safe(Y ~ 1 + (1 | Z), data.frame(Y = Y, Z = Z), offset_vec = f_xgb)
      if(is.null(fitglmer)) break
      b_new <- as.numeric(ranef(fitglmer)$Z[, 1])
      names(b_new) <- rownames(ranef(fitglmer)$Z)
      b <- b_new
      if(max(abs(b - b_old)) < hp_mexgb$toll) converged <- TRUE
    }
    if(!is.null(fit)) {
      b_val <- as.numeric(b[as.character(Z)])
      b_val[is.na(b_val)] <- 0
      lk <- as.numeric(predict(fit, xgb.DMatrix(X_mat))) + b_val
      preds <- as.numeric(exp(lk) / (1 + exp(lk)))
      model_obj <- list(xgb = fit, glmer = fitglmer)
    }
  } else if(model_type == "MELGBM") {
    X_mat <- model.matrix(~ ., data = X_sc)[, -1, drop = FALSE]
    eta_static <- as.numeric(glm(Y ~ ., data = data.frame(Y = Y, X_sc), family = binomial())$linear.predictors)
    b <- rep(0, length(unique(Z)))
    names(b) <- unique(Z)
    converged <- FALSE
    iter <- 0
    fit <- NULL
    fitglmer <- NULL
    while(!converged && iter < hp_melgbm$itmax) {
      iter <- iter + 1
      b_old <- b
      target <- as.numeric(eta_static) - as.numeric(b[as.character(Z)])
      fit <- lgb.train(list(objective = "regression", metric = "l2", max_depth = hp_melgbm$max_depth, learning_rate = hp_melgbm$learning_rate, min_data_in_leaf = hp_melgbm$min_data_in_leaf, verbose = -1), lgb.Dataset(X_mat, label = target), hp_melgbm$nrounds)
      f_lgb <- as.numeric(predict(fit, X_mat))
      fitglmer <- fit_glmer_safe(Y ~ 1 + (1 | Z), data.frame(Y = Y, Z = Z), offset_vec = f_lgb)
      if(is.null(fitglmer)) break
      b_new <- as.numeric(ranef(fitglmer)$Z[, 1])
      names(b_new) <- rownames(ranef(fitglmer)$Z)
      b <- b_new
      if(max(abs(b - b_old)) < hp_melgbm$toll) converged <- TRUE
    }
    if(!is.null(fit)) {
      b_val <- as.numeric(b[as.character(Z)])
      b_val[is.na(b_val)] <- 0
      lk <- as.numeric(predict(fit, X_mat)) + b_val
      preds <- as.numeric(exp(lk) / (1 + exp(lk)))
      model_obj <- list(lgb = fit, glmer = fitglmer)
    }
  }
  
  metrics <- calc_metrics_vec(preds, Y)
  dep_effects <- list()
  
  for(v in covs) {
    is_categorical <- !v %in% cont_vars
    
    if(is_categorical) {
      x_vals_real <- levels(as.factor(X_raw[[v]]))
      x_vals_numeric <- seq_along(x_vals_real)
      x_vals_model <- x_vals_real
    } else {
      x_vals_real <- seq(min(X_raw[[v]], na.rm = TRUE), max(X_raw[[v]], na.rm = TRUE), length.out = 500)
      x_vals_numeric <- x_vals_real
      x_vals_model <- seq(min(X_sc[[v]], na.rm = TRUE), max(X_sc[[v]], na.rm = TRUE), length.out = 500)
    }
    
    pred_grid <- sapply(1:length(x_vals_numeric), function(i) {
      tmp_data <- X_sc
      if(is_categorical) {
        tmp_data[[v]] <- factor(x_vals_model[i], levels = levels(as.factor(X_sc[[v]])))
      } else {
        tmp_data[[v]] <- as.numeric(x_vals_model[i])
      }
      tmp_mat <- model.matrix(~ ., data = tmp_data)[, -1, drop = FALSE]
      
      pred_arr <- if(model_type == "GLM") {
        predict(model_obj, newdata = data.frame(tmp_data), type = "response")
      } else if(model_type == "GLMM") {
        predict(model_obj, newdata = data.frame(tmp_data), type = "response", re.form = NA)
      } else if(model_type == "RF") {
        predict(model_obj, tmp_data, type = "prob")[, "1"]
      } else if(model_type == "GMERF") {
        predict(model_obj$rf, tmp_data)
      } else if(model_type == "XGB") {
        predict(model_obj, xgb.DMatrix(tmp_mat))
      } else if(model_type == "MEXGB") {
        predict(model_obj$xgb, xgb.DMatrix(tmp_mat))
      } else if(model_type == "LGBM") {
        predict(model_obj, tmp_mat)
      } else if(model_type == "MELGBM") {
        predict(model_obj$lgb, tmp_mat)
      }
      
      if(model_type %in% c("GMERF", "MEXGB", "MELGBM")) pred_arr <- exp(pred_arr) / (1 + exp(pred_arr))
      mean(pred_arr) * 100
    })
    
    dep_effects[[v]] <- data.frame(variable = v, x_val = x_vals_numeric, prob_os = pred_grid, cat_label = if(is_categorical) as.character(x_vals_real) else NA)
  }
  
  list(model_name = model_type, subset = subset_name, metrics = metrics, pdp_data = do.call(rbind, dep_effects), covariates = covs, cont_vars = cont_vars)
}

plot_pdp_dynamic <- function(eval_result) {
  df_pdp <- eval_result$pdp_data
  df_pdp$variable <- factor(df_pdp$variable, levels = eval_result$covariates)
  is_cat <- setdiff(eval_result$covariates, eval_result$cont_vars)
  
  p <- ggplot(df_pdp, aes(x = x_val, y = prob_os)) +
    geom_line(data = subset(df_pdp, !variable %in% is_cat), color = "#0072B2", linewidth = 1.2) +
    geom_point(data = subset(df_pdp, variable %in% is_cat), color = "#0072B2", size = 3) +
    geom_text(data = subset(df_pdp, variable %in% is_cat), aes(label = cat_label), vjust = -1.2, fontface = "bold", color = "black", size = 3.5) +
    facet_wrap(~ variable, scales = "free_x") +
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_text(face = "bold", size = 11), plot.title = element_text(face = "bold")) +
    labs(title = paste("partial dependence plots -", eval_result$model_name), subtitle = paste("Subset:", eval_result$subset), x = "variable value", y = "marginal os2ys pred (%)")
  print(p)
}

subset_list <- list(
  GS5 = list(covs = c("Age","Volume_PTV","schedule_ok_chemo","IMMUNO","riskfactors1","riskfactors2"),
             cont = c("Age", "Volume_PTV")),
  Lspecific2a = list(covs = c("Age","Volume_PTV","schedule_ok_chemo","IMMUNO","riskfactors1","riskfactors2","Subcutaneous_Fat","Skeletal_Muscle"),
                     cont = c("Age","Volume_PTV", "Subcutaneous_Fat", "Skeletal_Muscle")),
  Bhole = list(covs = c("Age","Volume_PTV","LAA", "riskfactors1"),
               cont = c("Age","Volume_PTV","LAA")),
  Bspecific2b_2 = list(covs = c("Age","Volume_PTV","IMMUNO","LAA","Subcutaneous_Fat","Skeletal_Muscle", "diameter.kidney.sep..mm."),
                       cont = c("Age","Volume_PTV","Subcutaneous_Fat", "LAA", "Skeletal_Muscle", "diameter.kidney.sep..mm.")),
  Bspecific1b_2 = list(covs = c("Age","Volume_PTV","IMMUNO"),
                       cont = c("Age","Volume_PTV")),
  Lassohole = list(covs = c("Age","Volume_PTV","IMMUNO","eud_lungs","Subcutaneous_Fat","aorta_agatston", "riskfactors1", "riskfactors2"),
                   cont = c("Age","Volume_PTV","eud_lungs","Subcutaneous_Fat","aorta_agatston")),
  Lhole_2 = list(covs = c("Age","KPS", "Volume_PTV","IMMUNO","eud_lungs","mean_dose_atriumL", "Subcutaneous_Fat", "riskfactors2"),
                 cont = c("Age", "Volume_PTV","eud_lungs","mean_dose_atriumL", "Subcutaneous_Fat"))
)

# GS1 -----
res <- evaluate_full_model("XGB", "GS5", subset_list[["GS5"]], data)
plot_pdp_dynamic(res)

# Bhole ----
res <- evaluate_full_model("XGB", "Bhole", subset_list[["Bhole"]], data)
plot_pdp_dynamic(res)

# Lspecific 2a ----
res <- evaluate_full_model("XGB", "Lspecific2a", subset_list[["Lspecific2a"]], data)
plot_pdp_dynamic(res)

# Bspecific2b_2 ------
res <- evaluate_full_model("XGB", "Bspecific2b_2", subset_list[["Bspecific2b_2"]], data)
plot_pdp_dynamic(res)

# Bspecific1b_2 MELGBM ------
res <- evaluate_full_model("MELGBM", "Bspecific1b_2", subset_list[["Bspecific1b_2"]], data)
plot_pdp_dynamic(res)

# Lassohole GMERF -----
res <- evaluate_full_model("GMERF", "Lassohole", subset_list[["Lassohole"]], data)
plot_pdp_dynamic(res)

# Lhole_2 RF vs GMERF -----
res <- evaluate_full_model("GMERF", "Lhole_2", subset_list[["Lhole_2"]], data)
plot_pdp_dynamic(res)
res <- evaluate_full_model("RF", "Lhole_2", subset_list[["Lhole_2"]], data)
plot_pdp_dynamic(res)


run_shap_pipeline <- function(sub_name, sub_info, with_color = FALSE) {
  covs <- sub_info$covs
  cont_vars <- sub_info$cont
  
  df_model <- data[, c(covs, "OS2ys_num"), drop = FALSE]
  df_model <- na.omit(df_model)
  
  for(v in covs) {
    if(!v %in% cont_vars) {
      df_model[[v]] <- as.factor(df_model[[v]])
    }
  }
  
  X_raw <- df_model[, covs, drop = FALSE]
  X_sc <- X_raw
  is_num <- colnames(X_sc) %in% cont_vars
  if(any(is_num)) {
    X_sc[, is_num] <- scale(X_sc[, is_num])
  }
  
  Y <- df_model$OS2ys_num
  target_factor <- factor(ifelse(Y == 1, "yes", "no"), levels = c("no", "yes"))
  
  formula_str <- paste("~", paste(covs, collapse = " + "))
  X_mat <- model.matrix(as.formula(formula_str), data = X_sc)[, -1, drop = FALSE]
  
  set.seed(123)
  dtrain <- xgb.DMatrix(data = X_mat, label = Y)
  xgb_mod <- xgb.train(
    params = list(objective = "binary:logistic", max_depth = hp_xgb$max_depth, eta = hp_xgb$eta, lambda = hp_xgb$lambda, subsample = hp_xgb$subsample),
    data = dtrain,
    nrounds = hp_xgb$nrounds,
    verbose = 0
  )
  
  unified_mod <- xgboost.unify(xgb_mod, X_mat)
  shap_res <- treeshap(unified_mod, X_mat, verbose = FALSE)
  shap_raw <- as.matrix(shap_res$shaps)
  
  col_map_sub <- make_col_map(X_sc)
  shap_agg_mat <- aggregate_shap(shap_raw, col_map_sub, covs)
  shap_df <- as.data.frame(shap_agg_mat)
  colnames(shap_df) <- covs
  
  plot_list <- list()
  
  for(v in covs) {
    is_cat <- !v %in% cont_vars
    
    df_sub <- data.frame(
      x_val      = X_raw[[v]],
      shap_value = shap_df[[v]],
      os2ys      = target_factor
    )
    
    if(is_cat) {
      df_sub$x_val <- as.factor(df_sub$x_val)
      
      if(with_color) {
        p <- ggplot(df_sub, aes(x = x_val, y = shap_value, color = os2ys)) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
          geom_point(position = position_jitter(width = 0.12, height = 0, seed = 123), alpha = 0.7, size = 1.8) +
          stat_summary(aes(group = 1), fun = mean, geom = "point", shape = 18, size = 3.2, color = "black") +
          stat_summary(aes(group = 1), fun = mean, geom = "line", color = "black", linetype = "dashed", linewidth = 0.8) +
          scale_x_discrete() +
          scale_color_manual(values = c("no" = "#D55E00", "yes" = "#0072B2"))
      } else {
        p <- ggplot(df_sub, aes(x = x_val, y = shap_value)) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
          geom_point(position = position_jitter(width = 0.12, height = 0, seed = 123), color = "#0072B2", alpha = 0.7, size = 1.8) +
          stat_summary(aes(group = 1), fun = mean, geom = "point", shape = 18, size = 3.2, color = "black") +
          stat_summary(aes(group = 1), fun = mean, geom = "line", color = "black", linetype = "dashed", linewidth = 0.8) +
          scale_x_discrete()
      }
      
    } else {
      df_sub$x_val <- as.numeric(df_sub$x_val)
      
      if(with_color) {
        p <- ggplot(df_sub, aes(x = x_val, y = shap_value, color = os2ys)) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
          geom_point(alpha = 0.7, size = 1.8) +
          geom_smooth(aes(group = 1), method = "loess", se = TRUE, color = "black", linewidth = 0.8, span = 0.9) +
          scale_x_continuous(expand = expansion(mult = 0.1)) +
          scale_color_manual(values = c("no" = "#D55E00", "yes" = "#0072B2"))
      } else {
        p <- ggplot(df_sub, aes(x = x_val, y = shap_value)) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
          geom_point(color = "#0072B2", alpha = 0.7, size = 1.8) +
          geom_smooth(aes(group = 1), method = "loess", se = TRUE, color = "black", linewidth = 0.8, span = 0.9) +
          scale_x_continuous(expand = expansion(mult = 0.1))
      }
    }
    
    p <- p +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9)
      ) +
      labs(
        title = v,
        y = "shap value (log-odds)",
        color = if(with_color) "observed os 2ys:" else NULL
      )
    
    plot_list[[v]] <- p
  }
  
  if(with_color) {
    p_grid <- wrap_plots(plot_list, ncol = 3) +
      plot_layout(guides = "collect") +
      plot_annotation(
        title = paste("shap dependence plots -", sub_name),
        theme = theme(plot.title = element_text(face = "bold", size = 12), legend.position = "top")
      ) &
      theme(legend.position = "top")
  } else {
    p_grid <- wrap_plots(plot_list, ncol = 3) +
      plot_annotation(
        title = paste("shap dependence plots -", sub_name),
        theme = theme(plot.title = element_text(face = "bold", size = 12))
      )
  }
  
  return(p_grid)
}

# GS1 with color -----
p_grid_color_GS1 <- run_shap_pipeline("GS5", subset_list[["GS5"]], with_color = TRUE)
print(p_grid_color_GS1)

# GS1 no color -----
p_grid_nocolor_GS1 <- run_shap_pipeline("GS5", subset_list[["GS5"]], with_color = FALSE)
print(p_grid_nocolor_GS1)

# Bhole no color -----
p_grid_nocolor_Bhole <- run_shap_pipeline("Bhole", subset_list[["Bhole"]], with_color = FALSE)
print(p_grid_nocolor_Bhole)

# Lspecific2a no color -----
p_grid_nocolor_Lspecific2a <- run_shap_pipeline("Lspecific2a", subset_list[["Lspecific2a"]], with_color = FALSE)
print(p_grid_nocolor_Lspecific2a)

# Bspecific2b_2 no color -----
p_grid_nocolor_Bspecific2b_2 <- run_shap_pipeline("Bspecific2b_2", subset_list[["Bspecific2b_2"]], with_color = FALSE)
print(p_grid_nocolor_Bspecific2b_2)