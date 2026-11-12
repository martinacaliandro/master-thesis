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









table_global <- table(data$hospital, data$OS2ys)
table_global
chisq.test(table_global)
data$is_center3 <- ifelse(data$hospital == "3", "3", "Other")
table_2x2 <- table(data$is_center3, data$OS2ys)
table_2x2
chisq.test(table_2x2)



#KPS vs R1 -----
library(ggplot2)

data_plot <- data
data_plot$KPS <- addNA(as.factor(data_plot$KPS))
levels(data_plot$KPS)[is.na(levels(data_plot$KPS))] <- "NA"

ggplot(data_plot, aes(x = as.factor(riskfactors1), fill = KPS)) +
  geom_bar(position = "fill", color = "white", linewidth = 0.5, width = 0.8) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(
    values = c(
      "60" = "#FFFFBF", 
      "70" = "#FDAE61", 
      "80" = "#F46D43", 
      "90" = "#D73027", 
      "100" = "#A50026", 
      "NA" = "#CBD5E8"
    ),
    name = "KPS"
  ) +
  labs(
    title = NULL,
    x = "Risk factors 1",
    y = "Proportion"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(1, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



library(ggplot2)
library(dplyr)

data_clean <- data %>%
  filter(!is.na(KPS) & !is.na(riskfactors1))

ggplot(data_clean, aes(x = as.factor(riskfactors1), fill = as.factor(KPS))) +
  geom_bar(position = "fill", color = "white", linewidth = 0.5, width = 0.8) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_brewer(palette = "YlOrRd", name = "KPS") +
  labs(
    title = NULL,
    x = "Risk factors 1",
    y = "Proportion"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(1, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )




#cal plots Bhole, Lhole_2 -----
hp_xgb <- list(max_depth = 2, eta = 0.05, lambda = 5, subsample = 0.7, nrounds = 100)
hp_xgb     <- list(max_depth = 2, eta = 0.01, lambda = 15, subsample = 0.7, nrounds = 400)
run_one_split_xgb <- function(X_mat, Y, cont_vars, folds, hp) {
  pooled <- numeric(length(Y))
  cols_scale <- intersect(colnames(X_mat), cont_vars)
  for (i in seq_len(length(folds))) {
    test_idx <- folds[[i]]
    X_train <- X_mat[-test_idx, , drop = FALSE]; X_val <- X_mat[test_idx, , drop = FALSE]; Y_train <- Y[-test_idx]
    if (length(cols_scale) > 0) {
      sc_p <- scale(X_train[, cols_scale, drop = FALSE])
      X_train[, cols_scale] <- sc_p
      X_val[, cols_scale] <- scale(X_val[, cols_scale, drop = FALSE], center = attr(sc_p, "scaled:center"), scale = attr(sc_p, "scaled:scale"))
    }
    fit <- xgb.train(list(objective = "binary:logistic", max_depth = hp$max_depth, eta = hp$eta, lambda = hp$lambda, subsample = hp$subsample), xgb.DMatrix(X_train, label = Y_train), hp$nrounds, verbose = 0)
    pooled[test_idx] <- as.numeric(predict(fit, xgb.DMatrix(X_val)))
  }
  as.numeric(unlist(pooled))
}

plot_calibration_xgb_10repeats <- function(subset_info, data, n_repeats = 10, k = 5, span = 0.75, grid_n = 200) {
  X <- data[, subset_info$covs, drop = FALSE]; Y <- data$OS2ys_num; cont_vars <- subset_info$cont
  idx <- complete.cases(X) & !is.na(Y); X <- X[idx, , drop = FALSE]; Y <- Y[idx]
  X_mat <- model.matrix(~ ., data = X)[, -1, drop = FALSE]; x_grid <- seq(0, 1, length.out = grid_n)
  curves <- vector("list", n_repeats)
  
  for (r in seq_len(n_repeats)) {
    set.seed(r); folds <- createFolds(factor(Y, levels = c(0, 1)), k = k, list = TRUE, returnTrain = FALSE)
    p <- run_one_split_xgb(X_mat = X_mat, Y = Y, cont_vars = cont_vars, folds = folds, hp = hp_xgb)
    valid <- is.finite(p) & !is.na(Y); p <- p[valid]; y <- Y[valid]
    loess_fit <- tryCatch(loess(y ~ p, data = data.frame(y = y, p = p), span = span, degree = 1, control = loess.control(surface = "direct")), error = function(e) NULL)
    
    if (!is.null(loess_fit)) {
      y_hat <- tryCatch(as.numeric(predict(loess_fit, newdata = data.frame(p = x_grid))), error = function(e) rep(NA_real_, length(x_grid)))
      y_hat[x_grid < min(p) | x_grid > max(p)] <- NA_real_
      curves[[r]] <- data.frame(x = x_grid, y = y_hat, Repeat = paste0("Repeat ", r))
    }
  }
  
  curves_df <- do.call(rbind, curves[!vapply(curves, is.null, logical(1))])
  mean_curve <- aggregate(y ~ x, data = curves_df, FUN = function(z) if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE))
  names(mean_curve) <- c("x", "y_mean")
  
  g <- ggplot() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_line(data = curves_df, aes(x = x, y = y, group = Repeat), color = "grey65", alpha = 0.65, linewidth = 0.6, na.rm = TRUE) +
    geom_line(data = mean_curve, aes(x = x, y = y_mean), color = "#2166AC", linewidth = 1.3, na.rm = TRUE) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2)) + scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    labs(title = NULL, subtitle = NULL, x = "Predicted probability", y = "Observed probability") +
    theme_classic(base_size = 20)
  
  print(g)
  invisible(list(plot = g, curves = curves_df, mean_curve = mean_curve))
}

hp_gmerf <- list(ntree = 500, nodesize = 20, itmax = 15, toll = 0.05)

fit_glmer_safe <- function(formula, data, offset_vec) {
  tryCatch({
    glmer(formula, data = data, family = binomial(), offset = offset_vec, 
          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  }, error = function(e) {
    tryCatch({
      glmer(formula, data = data, family = binomial(), offset = offset_vec)
    }, error = function(e2) NULL)
  })
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

plot_calibration_gmerf_10repeats <- function(subset_info, data, n_repeats = 10, k = 5, span = 0.75, grid_n = 200) {
  X <- data[, subset_info$covs, drop = FALSE]
  Y <- data$OS2ys_num
  Z <- data[[subset_info$z_var]]
  cont_vars <- subset_info$cont
  
  idx <- complete.cases(X) & !is.na(Y) & !is.na(Z)
  X <- X[idx, , drop = FALSE]; Y <- Y[idx]; Z <- Z[idx]
  x_grid <- seq(0, 1, length.out = grid_n)
  curves <- vector("list", n_repeats)
  
  for (r in seq_len(n_repeats)) {
    set.seed(r); folds <- createFolds(factor(Y, levels = c(0, 1)), k = k, list = TRUE, returnTrain = FALSE)
    p <- run_one_split_gmerf(X = X, Y = Y, Z = Z, cont_vars = cont_vars, folds = folds, hp = hp_gmerf)
    valid <- is.finite(p) & !is.na(Y); p <- p[valid]; y <- Y[valid]
    loess_fit <- tryCatch(loess(y ~ p, data = data.frame(y = y, p = p), span = span, degree = 1, control = loess.control(surface = "direct")), error = function(e) NULL)
    
    if (!is.null(loess_fit)) {
      y_hat <- tryCatch(as.numeric(predict(loess_fit, newdata = data.frame(p = x_grid))), error = function(e) rep(NA_real_, length(x_grid)))
      y_hat[x_grid < min(p) | x_grid > max(p)] <- NA_real_
      curves[[r]] <- data.frame(x = x_grid, y = y_hat, Repeat = paste0("Repeat ", r))
    }
  }
  
  curves_df <- do.call(rbind, curves[!vapply(curves, is.null, logical(1))])
  mean_curve <- aggregate(y ~ x, data = curves_df, FUN = function(z) if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE))
  names(mean_curve) <- c("x", "y_mean")
  
  g <- ggplot() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_line(data = curves_df, aes(x = x, y = y, group = Repeat), color = "grey65", alpha = 0.65, linewidth = 0.6, na.rm = TRUE) +
    geom_line(data = mean_curve, aes(x = x, y = y_mean), color = "#2166AC", linewidth = 1.3, na.rm = TRUE) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2)) + scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    labs(title = NULL, subtitle = NULL, x = "Predicted probability", y = "Observed probability") +
    theme_classic(base_size = 20)
  
  print(g)
  invisible(list(plot = g, curves = curves_df, mean_curve = mean_curve))
}

subset_list <- list(
  Bhole = list(
    covs = c("Age","Volume_PTV","LAA", "riskfactors1"), 
    cont = c("Age","Volume_PTV","LAA")
  ),
  Lhole_2 = list(
    covs = c("Age","KPS", "Volume_PTV","IMMUNO","eud_lungs","mean_dose_atriumL", "Subcutaneous_Fat", "riskfactors2"), 
    cont = c("Age", "Volume_PTV","eud_lungs","mean_dose_atriumL", "Subcutaneous_Fat"),
    z_var = "hospital"
  )
)

cal_bhole <- plot_calibration_xgb_10repeats(subset_info = subset_list$Bhole, data = data, n_repeats = 10, k = 5, span = 0.75)
cal_lhole2 <- plot_calibration_gmerf_10repeats(subset_info = subset_list$Lhole_2, data = data, n_repeats = 10, k = 5, span = 0.75)

#raneffCI Lhole_2 ----


#PDP for rand eff -----
library(ggplot2)
library(patchwork)

hp_rf      <- list(ntree = 500, nodesize = 20)
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
  
  hosp_counts <- table(df_model$hospital)
  valid_hospitals <- names(hosp_counts[hosp_counts >= 10])
  df_model <- df_model[df_model$hospital %in% valid_hospitals, ]
  df_model$hospital <- droplevels(factor(df_model$hospital))
  
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
      }
      
      if(model_type %in% c("GMERF", "MEXGB", "MELGBM")) pred_arr <- exp(pred_arr) / (1 + exp(pred_arr))
      mean(pred_arr) * 100
    })
    
    dep_effects[[v]] <- data.frame(variable = v, x_val = x_vals_numeric, prob_os = pred_grid, cat_label = if(is_categorical) as.character(x_vals_real) else NA)
  }
  
  list(model_name = model_type, subset = subset_name, metrics = metrics, pdp_data = do.call(rbind, dep_effects), covariates = covs, cont_vars = cont_vars)
}

plot_pdp_dynamic <- function(eval_result, custom_labels = NULL) {
  df_pdp <- eval_result$pdp_data
  covs <- eval_result$covariates
  cont_vars <- eval_result$cont_vars
  
  plot_list <- list()
  
  for (v in covs) {
    df_v <- subset(df_pdp, variable == v)
    is_cat <- !v %in% cont_vars
    plot_title <- if (!is.null(custom_labels) && v %in% names(custom_labels)) custom_labels[v] else v
    
    p <- ggplot(df_v, aes(x = x_val, y = prob_os))
    
    if (is_cat) {
      p <- p + 
        geom_point(color = "#0072B2", size = 3) +
        geom_text(aes(label = cat_label), vjust = -1.2, fontface = "bold", color = "black", size = 3.5) +
        scale_x_continuous(breaks = df_v$x_val, labels = rep("", nrow(df_v)), expand = ggplot2::expansion(mult = 0.2)) +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        )
    } else {
      p <- p + 
        geom_line(color = "#0072B2", linewidth = 1.2) +
        scale_x_continuous(expand = ggplot2::expansion(mult = 0.15)) +
        theme(
          axis.text.x = element_text(size = 9),
          axis.ticks.x = element_line()
        )
    }
    
    p <- p + 
      coord_cartesian(clip = "off") + 
      scale_y_continuous(expand = ggplot2::expansion(mult = c(0.1, 0.25))) + 
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
        axis.title.x = element_text(face = "plain", size = 10),
        axis.title.y = element_text(face = "plain", size = 10),
        plot.subtitle = element_blank(),
        plot.margin = ggplot2::margin(t = 15, r = 10, b = 10, l = 10)
      ) +
      labs(
        title = plot_title,
        x = "clinical value",
        y = "marginal OS2ys predicted probability"
      )
    
    plot_list[[v]] <- p
  }
  
  final_plot <- wrap_plots(plot_list) + plot_layout(axes = "collect")
  print(final_plot)
}

# Configurazione solo per Lhole_2
subset_list <- list(
  Lhole_2 = list(covs = c("Age","KPS", "Volume_PTV","IMMUNO","eud_lungs","mean_dose_atriumL", "Subcutaneous_Fat", "riskfactors2"),
                 cont = c("Age", "Volume_PTV","eud_lungs","mean_dose_atriumL", "Subcutaneous_Fat"))
)

miei_titoli <- c(
  "Age" = "Age [years]",
  "KPS" = "KPS",
  "Volume_PTV" = "Volume PTV [cm³]",
  "IMMUNO" = "Immunotherapy",
  "eud_lungs" = "Lungs EUD [Gy]",
  "mean_dose_atriumL" = "Mean dose atrium left [Gy]",
  "Subcutaneous_Fat" = "Subcutaneous fat (%)",
  "riskfactors2" = "Risk factors 2"
)

res1 <- evaluate_full_model("GMERF", "Lhole_2", subset_list[["Lhole_2"]], data)
plot_pdp_dynamic(res1, custom_labels = miei_titoli)

res2 <- evaluate_full_model("RF", "Lhole_2", subset_list[["Lhole_2"]], data)
plot_pdp_dynamic(res2, custom_labels = miei_titoli)






#boxplots OS2ys across hospitals -----
library(ggplot2)
library(patchwork)
library(dplyr)
library(scales)

data$hospital <- as.factor(data$hospital)

if(is.numeric(data$OS2ys) || is.integer(data$OS2ys)) {
  data$OS2ys_fac <- factor(data$OS2ys, levels = c(0, 1), labels = c("NO", "YES"))
} else {
  data$OS2ys_fac <- as.factor(data$OS2ys)
}

hp_palette <- c("#ff9999", "#99ccff", "#b6d7a8", "#d5a6bd", "#f9cb9c", "#ffe599", "#d9d2e9")

df_counts <- data %>%
  filter(hospital %in% as.character(1:7)) %>%
  group_by(hospital) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(hospital = factor(hospital, levels = as.character(1:7)))

p_count <- ggplot(df_counts, aes(x = hospital, y = n, fill = hospital)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = -0.5, size = 6, fontface = "bold", color = "black") +
  scale_fill_manual(values = hp_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "plain", size = 20),
    axis.text = element_text(size = 16)
  ) +
  labs(x = "Hospital", y = "Number of patients")

print(p_count)

p1 <- ggplot(data, aes(x = hospital, y = OS, fill = hospital)) +
  geom_boxplot(alpha = 0.85, outlier.shape = NA, color = "black", linewidth = 0.5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7, color = "black") +
  scale_fill_manual(values = hp_palette) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "plain", size = 20, hjust = 0),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "plain", size = 20),
    axis.text = element_text(size = 16)
  ) +
  labs(x = "Hospital", y = "OS")


df_prop <- data %>%
  filter(!is.na(OS2ys_fac) & !is.na(hospital)) %>%
  group_by(hospital, OS2ys_fac) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(hospital) %>%
  mutate(
    prop = count / sum(count),
    label = paste0(round(prop * 100), "%")
  )

p2 <- ggplot(data, aes(x = hospital, fill = OS2ys_fac)) +
  geom_bar(position = "fill", color = "white", linewidth = 0.8) +
  geom_text(data = df_prop, 
            aes(y = prop, label = label, group = OS2ys_fac), 
            position = position_fill(vjust = 0.5), 
            size = 6, fontface = "bold", color = "black") +
  scale_fill_manual(values = c("NO" = "#ff8080", "YES" = "#93c47d")) + 
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "plain", size = 20, hjust = 0),       
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 20),               
    legend.text = element_text(size = 16),                                
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "plain", size = 20),                 
    axis.text = element_text(size = 16)                                   
  ) +
  labs(x = "Hospital", y = "Proportion", fill = "OS2ys")


final_plot <- p1 / p2 + plot_layout(heights = c(1, 1))

print(final_plot)


#boxplot features across hospitals -----
hp_palette <- c("#ff9999", "#99ccff", "#b6d7a8", "#d5a6bd", "#f9cb9c", "#ffe599", "#d9d2e9")
p3 <- ggplot(data, aes(x = hospital, y = Volume_PTV, fill = hospital)) +
  geom_boxplot(alpha = 0.85, outlier.shape = NA, color = "black", linewidth = 0.5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7, color = "black") +
  scale_fill_manual(values = hp_palette) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "plain", size = 14, hjust = 0),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_text(face = "bold")
  ) +
  labs(x = "Hospital", y = "Volume PTV")

p4 <- ggplot(data, aes(x = hospital, y = Age, fill = hospital)) +
  geom_boxplot(alpha = 0.85, outlier.shape = NA, color = "black", linewidth = 0.5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7, color = "black") +
  scale_fill_manual(values = hp_palette) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "plain", size = 14, hjust = 0),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_text(face = "bold")
  ) +
  labs(x = "Hospital", y = "Age")

data_plot <- data %>%
  filter(!is.na(hospital)) %>%
  mutate(
    KPS = as.character(KPS),
    KPS = ifelse(is.na(KPS), "NA", KPS),
    KPS = factor(KPS, levels = c("100", "90", "80", "70", "60", "NA"))
  )

df_prop <- data_plot %>%
  group_by(hospital, KPS) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(hospital) %>%
  mutate(
    prop = count / sum(count),
    label = ifelse(prop > 0.02, paste0(round(prop * 100), "%"), "")
  )

p5 <- ggplot(data_plot, aes(x = hospital, fill = KPS)) + 
  geom_bar(position = "fill", color = "white", linewidth = 0.8) +
  geom_text(data = df_prop, 
            aes(y = prop, label = label, group = KPS), 
            position = position_fill(vjust = 0.5), 
            size = 6, fontface = "bold", color = "black") +
  scale_fill_manual(
    name = "KPS",
    values = c(
      "60"  = "#f26b6b",
      "70"  = "#fc8d59",
      "80"  = "#fee08b",
      "90"  = "#d9ef8b",
      "100" = "#91cf60",
      "NA"  = "#cccccc"
    )
  ) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "plain", size = 20, hjust = 0),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 16),              
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "plain", size = 20),  
    axis.text = element_text(size = 16)                    
  ) +
  guides(fill = guide_legend(nrow = 1, reverse = TRUE)) + 
  labs(x = "Hospital", y = "Proportion", fill = "KPS")

print(p5)

print(p3)
print(p4)
print(p5)


#PDP for base -----------------
library(ggplot2)
library(patchwork)

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

plot_pdp_dynamic <- function(eval_result, custom_labels = NULL) {
  df_pdp <- eval_result$pdp_data
  covs <- eval_result$covariates
  cont_vars <- eval_result$cont_vars
  
  plot_list <- list()
  
  for (v in covs) {
    df_v <- subset(df_pdp, variable == v)
    is_cat <- !v %in% cont_vars
    plot_title <- if (!is.null(custom_labels) && v %in% names(custom_labels)) custom_labels[v] else v
    
    p <- ggplot(df_v, aes(x = x_val, y = prob_os))
    
    if (is_cat) {
      p <- p + 
        geom_point(color = "#0072B2", size = 3) +
        geom_text(aes(label = cat_label), vjust = -1.2, fontface = "plain", color = "black", size = 6) +
        scale_x_continuous(breaks = df_v$x_val, labels = rep("", nrow(df_v)), expand = ggplot2::expansion(mult = 0.2)) +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        )
    } else {
      p <- p + 
        geom_line(color = "#0072B2", linewidth = 1.2) +
        scale_x_continuous(expand = ggplot2::expansion(mult = 0.15)) +
        theme(
          axis.text.x = element_text(size = 20),
          axis.ticks.x = element_line()
        )
    }
    
    p <- p + 
      coord_cartesian(clip = "off") + 
      scale_y_continuous(expand = ggplot2::expansion(mult = c(0.1, 0.25))) + 
      theme_minimal(base_size = 20) +
      theme(
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        axis.title.x = element_text(face = "plain", size = 20),
        axis.title.y = element_text(face = "plain", size = 20),
        axis.text.y = element_text(size = 20),
        plot.subtitle = element_blank(),
        plot.margin = ggplot2::margin(t = 15, r = 10, b = 10, l = 10)
      ) +
      labs(
        title = plot_title,
        x = "clinical value",
        y = "marginal OS2ys predicted probability"
      )
    
    plot_list[[v]] <- p
  }
  
  final_plot <- wrap_plots(plot_list) + plot_layout(axes = "collect")
  print(final_plot)
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

miei_titoli1 <- c(
  "Age" = "Age [years]",
  "Volume_PTV" = "Volume PTV [cm³]",
  "schedule_ok_chemo" = "schedule ok chemo",
  "IMMUNO" = "Immunotherapy",
  "riskfactors1" = "risk factors 1",
  "riskfactors2" = "risk factors 2"
)
miei_titoli2 <- c(
  "Age" = "Age [years]",
  "Volume_PTV" = "Volume PTV [cm³]",
  "LAA" = "LAA (%)",
  "riskfactors1" = "Risk factors 1"
)
miei_titoli3 <- c(
  "Age" = "Age [years]",
  "Volume_PTV" = "Volume PTV [cm³]",
  "IMMUNO" = "Immunotherapy",
  "LAA" = "LAA (%)",
  "Subcutaneous_Fat" = "Subcutaneous Fat (%)",
  "Skeletal_Muscle" = "Skeletal Muscle (%)",
  "diameter.kidney.sep..mm." = "diameter kidney sep [mm]"
)

res1 <- evaluate_full_model("XGB", "GS5", subset_list[["GS5"]], data)
plot_pdp_dynamic(res1, custom_labels = miei_titoli1)

res2 <- evaluate_full_model("XGB", "Bhole", subset_list[["Bhole"]], data)
plot_pdp_dynamic(res2, custom_labels = miei_titoli2)

res3 <- evaluate_full_model("XGB", "Bspecific2b_2", subset_list[["Bspecific2b_2"]], data)
plot_pdp_dynamic(res3, custom_labels = miei_titoli3)






# SHAP dependence plot Bhole -----
library(ggplot2)
library(patchwork)
library(xgboost)
library(treeshap)

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

run_shap_pipeline <- function(sub_name, sub_info, data_df, with_color = FALSE, custom_labels = NULL) {
  covs <- sub_info$covs
  cont_vars <- sub_info$cont
  
  df_model <- data_df[, c(covs, "OS2ys_num"), drop = FALSE]
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
    plot_title <- if (!is.null(custom_labels) && v %in% names(custom_labels)) custom_labels[[v]] else v
    
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
          geom_point(position = position_jitter(width = 0.12, height = 0, seed = 123), alpha = 0.7, size = 2.5) +
          stat_summary(aes(group = 1), fun = mean, geom = "point", shape = 18, size = 4.5, color = "black") +
          stat_summary(aes(group = 1), fun = mean, geom = "line", color = "black", linetype = "dashed", linewidth = 1) +
          scale_x_discrete() +
          scale_color_manual(values = c("no" = "#D55E00", "yes" = "#0072B2"))
      } else {
        p <- ggplot(df_sub, aes(x = x_val, y = shap_value)) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
          geom_point(position = position_jitter(width = 0.12, height = 0, seed = 123), color = "#0072B2", alpha = 0.7, size = 2.5) +
          stat_summary(aes(group = 1), fun = mean, geom = "point", shape = 18, size = 4.5, color = "black") +
          stat_summary(aes(group = 1), fun = mean, geom = "line", color = "black", linetype = "dashed", linewidth = 1) +
          scale_x_discrete()
      }
      
    } else {
      df_sub$x_val <- as.numeric(df_sub$x_val)
      
      if(with_color) {
        p <- ggplot(df_sub, aes(x = x_val, y = shap_value, color = os2ys)) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
          geom_point(alpha = 0.7, size = 2.5) +
          geom_smooth(aes(group = 1), method = "loess", se = TRUE, color = "black", linewidth = 1, span = 0.9) +
          scale_x_continuous(expand = expansion(mult = 0.1)) +
          scale_color_manual(values = c("no" = "#D55E00", "yes" = "#0072B2"))
      } else {
        p <- ggplot(df_sub, aes(x = x_val, y = shap_value)) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
          geom_point(color = "#0072B2", alpha = 0.7, size = 2.5) +
          geom_smooth(aes(group = 1), method = "loess", se = TRUE, color = "black", linewidth = 1, span = 0.9) +
          scale_x_continuous(expand = expansion(mult = 0.1))
      }
    }
    
    p <- p +
      theme_minimal(base_size = 20) +
      theme(
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)
      ) +
      labs(
        title = plot_title,
        y = "shap value (log-odds)",
        color = if(with_color) "observed os 2ys:" else NULL
      )
    
    plot_list[[v]] <- p
  }
  
  if(with_color) {
    p_grid <- wrap_plots(plot_list, ncol = 4) +
      plot_layout(guides = "collect") &
      theme(legend.position = "top", legend.text = element_text(size = 16), legend.title = element_text(size = 18))
  } else {
    p_grid <- wrap_plots(plot_list, ncol = 4)
  }
  
  return(p_grid)
}

hp_xgb     <- list(max_depth = 2, eta = 0.01, lambda = 15, subsample = 0.7, nrounds = 400)

subset_list <- list(
  Bhole = list(covs = c("Age","Volume_PTV","LAA", "riskfactors1"),
               cont = c("Age","Volume_PTV","LAA"))
)

data_copy <- data 

miei_titoli_shap <- c(
  "Age" = "Age [years]",
  "Volume_PTV" = "Volume PTV [cm³]",
  "LAA" = "Low Attenuation Area (%)",
  "riskfactors1" = "Risk factors 1"
)

p_grid_nocolor_Bhole <- run_shap_pipeline(
  sub_name = "B-all-R1", 
  sub_info = subset_list[["Bhole"]], 
  data_df = data_copy, 
  with_color = FALSE, 
  custom_labels = miei_titoli_shap
)

print(p_grid_nocolor_Bhole)












#cluster SHAP vs raw-----
hp_xgb     <- list(max_depth = 2, eta = 0.01, lambda = 15, subsample = 0.7, nrounds = 400)
covariates <- c("Age", "Volume_PTV", "LAA", "riskfactors1")
cont_vars  <- c("Age", "Volume_PTV", "LAA")


minPts_hdb <- 25
nn_umap    <- 30
seed_umap  <- 123 #or 123
seed_xgb   <- 123



idx <- complete.cases(data[, covariates]) & !is.na(data$OS2ys_num) & !is.na(data$hospital)
df_sub <- data[idx, ]
X_full <- df_sub[, covariates, drop = FALSE]
Y_full <- df_sub$OS2ys_num
Z_full <- droplevels(df_sub$hospital)

X_full_sc <- X_full
X_full_sc[, cont_vars] <- scale(X_full[, cont_vars])

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
  params  = list(objective = "binary:logistic", max_depth = hp_xgb$max_depth, eta = hp_xgb$eta, lambda = hp_xgb$lambda, subsample = hp_xgb$subsample),
  data    = xgb.DMatrix(X_mat, label = Y_full),
  nrounds = hp_xgb$nrounds,
  verbose = 0
)
p_final <- as.numeric(predict(fit_xgb, xgb.DMatrix(X_mat)))

uni_xgb <- xgboost.unify(fit_xgb, X_mat)
treeshap_res <- treeshap(uni_xgb, X_mat, verbose = FALSE)

shap_raw <- as.matrix(treeshap_res$shaps)
shap_agg <- aggregate_shap(shap_raw, X_full, covariates)
shap_4vars <- shap_agg 

n_repeats <- 10
shap_oof_arr <- array(NA_real_, dim = c(nrow(X_full), ncol(shap_raw), n_repeats))
for (r in seq_len(n_repeats)) {
  set.seed(r)
  folds_r <- createFolds(as.factor(Y_full), k = 5, list = TRUE, returnTrain = FALSE)
  for (te in folds_r) {
    Xtr <- X_full[-te, , drop = FALSE]; Xte <- X_full[te, , drop = FALSE]
    sc_k <- scale(Xtr[, cont_vars])
    Xtr[, cont_vars] <- sc_k; Xte[, cont_vars] <- scale(Xte[, cont_vars], center = attr(sc_k, "scaled:center"), scale = attr(sc_k, "scaled:scale"))
    Xtr_m <- model.matrix(~ ., data = Xtr)[, -1, drop = FALSE]; Xte_m <- model.matrix(~ ., data = Xte)[, -1, drop = FALSE]
    
    f_k <- xgb.train(
      params = list(objective = "binary:logistic", max_depth = hp_xgb$max_depth, eta = hp_xgb$eta, lambda = hp_xgb$lambda, subsample = hp_xgb$subsample),
      data = xgb.DMatrix(Xtr_m, label = Y_full[-te]), nrounds = hp_xgb$nrounds, verbose = 0
    )
    shap_oof_arr[te, , r] <- as.matrix(treeshap(xgboost.unify(f_k, Xte_m), Xte_m, verbose = FALSE)$shaps)
  }
}
shap_oof_7cols <- apply(shap_oof_arr, c(1, 2), mean)
colnames(shap_oof_7cols) <- colnames(shap_raw)
shap_oof_agg <- aggregate_shap(shap_oof_7cols, X_full, covariates)
shap_oof_4vars <- shap_oof_agg 

cl_A <- dbscan::hdbscan(shap_agg, minPts = minPts_hdb)
print(cl_A)

set.seed(seed_umap)
um2_shap7 <- uwot::umap(shap_raw, n_neighbors = nn_umap, min_dist = 0, n_components = 2, n_threads = 1, n_sgd_threads = 1)
cl_C <- dbscan::hdbscan(um2_shap7, minPts = minPts_hdb)
print(cl_C)

set.seed(seed_umap)
um2_shap <- uwot::umap(shap_agg, n_neighbors = nn_umap, min_dist = 0, n_components = 2, n_threads = 1, n_sgd_threads = 1)
cl_B <- dbscan::hdbscan(um2_shap, minPts = minPts_hdb)
print(cl_B)

##new mapping -----
mapping_manuale <- c("2" = "1", "1" = "2", "3" = "3", "5" = "4", "4" = "5")
map <- factor(mapping_manuale[as.character(cl_B$cluster)], levels = c("1", "2", "3", "4", "5"))
cl_B$cluster <- as.numeric(as.character(map))

set.seed(seed_umap)
um2_oof <- uwot::umap(shap_oof_agg, n_neighbors = nn_umap, min_dist = 0, n_components = 2, n_threads = 1, n_sgd_threads = 1)
cl_oof <- dbscan::hdbscan(um2_oof, minPts = minPts_hdb)
print(cl_oof)

X_raw_mm <- model.matrix(~ ., data = X_full_sc)[, -1, drop = FALSE]
set.seed(seed_umap)
um2_raw <- uwot::umap(X_raw_mm, n_neighbors = nn_umap, min_dist = 0, n_components = 2, n_threads = 1, n_sgd_threads = 1)

Cluster <- factor(cl_B$cluster)
cat("\n--- Valori Clinici VERI (Media) per Cluster (Variabili Continue) ---\n")
if (length(cont_vars) > 0) {
  
  raw_cont_mean <- aggregate(X_full[, cont_vars, drop = FALSE], list(Cluster = Cluster), mean)
  raw_cont_mean <- aggregate(X_full[, cont_vars, drop = FALSE], list(Cluster = Cluster), median)
  raw_cont_mean <- aggregate(X_full[, cont_vars, drop = FALSE], list(Cluster = Cluster), function(x) quantile(x, probs = c(0.25, 0.50, 0.75), na.rm = TRUE))
  
  raw_cont_sd   <- aggregate(X_full[, cont_vars, drop = FALSE], list(Cluster = Cluster), sd)
  tab_clinica_vera <- data.frame(Cluster = raw_cont_mean$Cluster, n = as.numeric(table(Cluster)), OS_2ys = paste0(round(tapply(Y_full, Cluster, mean) * 100, 1), "%"))
  for (cv in cont_vars) tab_clinica_vera[[cv]] <- sprintf("%.3f (±%.3f)", raw_cont_mean[[cv]], raw_cont_sd[[cv]])
  print(tab_clinica_vera, row.names = FALSE)
}


df_plot <- data.frame(
  SHAP_UMAP1 = um2_shap[, 1],
  SHAP_UMAP2 = um2_shap[, 2],
  RAW_UMAP1  = um2_raw[, 1],
  RAW_UMAP2  = um2_raw[, 2],
  Cluster    = Cluster,
  OS2ys      = factor(df_sub$OS2ys, levels = c("NO", "YES")),
  OS_months  = df_sub$OS,
  Pred_Prob  = p_final
)
for (v in covariates) df_plot[[v]] <- df_sub[[v]]
p_shap_2d <- ggplot(df_plot, aes(x = SHAP_UMAP1, y = SHAP_UMAP2, color = Cluster)) +
  geom_point(alpha = 0.85, size = 3) + 
  theme_minimal(base_size = 20) +
  labs(x = "UMAP 1 (SHAP)", y = "UMAP 2 (SHAP)") +
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )

p_raw_2d <- ggplot(df_plot, aes(x = RAW_UMAP1, y = RAW_UMAP2, color = Cluster)) +
  geom_point(alpha = 0.85, size = 3) + 
  theme_minimal(base_size = 20) +
  labs(x = "UMAP 1 (raw data)", y = "UMAP 2 (raw data)") +
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(p_shap_2d, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_raw_2d,  vp = viewport(layout.pos.row = 1, layout.pos.col = 2))

#clusters colored by OS2ys ------
library(ggplot2)
library(grid)

df_plot_copy <- df_plot

p_raw_os <- ggplot(df_plot_copy, aes(x = RAW_UMAP1, y = RAW_UMAP2, color = OS2ys)) +
  geom_point(alpha=0.8, size=3.5) +
  scale_color_manual(
    values = c("NO" = "#D55E00", "YES" = "#0072B2"), 
    name = "OS2ys:",
    labels = c("No", "Yes")
  ) +
  labs(x = "UMAP 1 (raw data)", y = "UMAP 2 (raw data)") +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20, margin = ggplot2::margin(t = 10)),
    axis.title.y = element_text(size = 20, margin = ggplot2::margin(r = 10)),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18),
    legend.position = "bottom"
  )

p_shap_os <- ggplot(df_plot_copy, aes(x = SHAP_UMAP1, y = SHAP_UMAP2, color = OS2ys)) +
  geom_point(alpha=0.8, size=3.5) +
  scale_color_manual(
    values = c("NO" = "#D55E00", "YES" = "#0072B2"), 
    name = "OS2ys:",
    labels = c("No", "Yes")
  ) +
  labs(x = "UMAP 1 (SHAP)", y = "UMAP 2 (SHAP)") +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20, margin = ggplot2::margin(t = 10)),
    axis.title.y = element_text(size = 20, margin = ggplot2::margin(r = 10)),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18),
    legend.position = "bottom"
  )

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))

print(p_raw_os,  vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_shap_os, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))



#SHAP dependence plots Bhole colored by clusters -----
library(ggplot2); library(ggbeeswarm); library(patchwork)

to_num <- function(v) {
  if (is.factor(v) || is.character(v)) {
    n <- suppressWarnings(as.numeric(as.character(v)))
    if (all(is.na(n) == is.na(v))) n else as.numeric(factor(v))
  } else as.numeric(v)
}
shap_dependence <- function(feature, shap_mat, X, color_by = NULL, cluster = NULL,
                            smooth = TRUE, title = NULL, custom_labels = NULL) {
  shap_mat <- as.matrix(shap_mat)
  stopifnot(feature %in% colnames(shap_mat), feature %in% colnames(X))
  d <- data.frame(x = X[[feature]], SHAP = shap_mat[, feature])
  is_cat <- is.factor(d$x) || is.character(d$x) || length(unique(d$x)) <= 5
  col_lab <- NULL
  if (!is.null(color_by)) {
    if (identical(color_by, "cluster")) { d$col <- factor(cluster); col_lab <- "Cluster" }
    else { d$col <- X[[color_by]]; col_lab <- color_by }
  }
  
  plot_title <- if (!is.null(custom_labels) && feature %in% names(custom_labels)) custom_labels[[feature]] else paste("Dependence:", feature)
  
  p <- ggplot(d, aes(x = if (is_cat) factor(x) else x, y = SHAP)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40")
  
  if (is_cat) {
    p <- p + geom_boxplot(outlier.shape = NA, fill = "grey92", color = "grey40", width = 0.5)
    p <- p + if (is.null(color_by)) geom_quasirandom(alpha = 0.6, size = 1.6, color = "#0072B2")
    else geom_quasirandom(aes(color = col), alpha = 0.7, size = 1.6)
    p <- p + stat_summary(aes(group = 1), fun = mean, geom = "point", shape = 18, size = 4.5, color = "black") +
      stat_summary(aes(group = 1), fun = mean, geom = "line", color = "black", linetype = "dashed", linewidth = 1)
  } else {
    p <- p + if (is.null(color_by)) geom_point(alpha = 0.6, size = 1.8, color = "#0072B2")
    else geom_point(aes(color = col), alpha = 0.75, size = 1.8)
    
    if (smooth) {
      p <- p + geom_smooth(aes(group = 1), method = "loess", formula = y ~ x, se = TRUE, color = "black", fill = "grey50", alpha = 0.2, linewidth = 0.8, span = 0.9)
    }
  }
  
  if (!is.null(color_by)) {
    p <- p + if (is.numeric(d$col)) scale_color_viridis_c(name = col_lab) else labs(color = col_lab)
  }
  
  p + theme_minimal(base_size = 20) +
    labs(title = plot_title,
         x = NULL, y = "SHAP value (log-odds)") +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 20),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20)
    )
}
data_copy <- na.omit(data[, colnames(shap_agg), drop = FALSE])

miei_titoli_shap <- c(
  "Age" = "Age [years]",
  "Volume_PTV" = "Volume PTV [cm³]",
  "LAA" = "LAA (%)",
  "riskfactors1" = "Risk factors 1"
)


plots <- lapply(seq_along(colnames(shap_agg)), function(i) {
  f <- colnames(shap_agg)[i]
  p <- shap_dependence(
    feature = f,
    shap_mat = shap_agg,
    X = data_copy,
    color_by = "cluster",
    cluster = Cluster,
    custom_labels = miei_titoli_shap
  )
  if (i < 4) {
    p <- p + theme(legend.position = "none")
  } else {
    p <- p + theme(legend.position = "right", legend.box = "vertical")
  }
  return(p)
})

print(wrap_plots(plots, ncol = 4))


#pred IP (OS2ys) by clusters -------
p_prob_cluster <- ggplot(
  df_plot,
  aes(x = Cluster, y = Pred_Prob, fill = Cluster)
) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1.5) +
  theme_minimal(base_size = 20) +
  labs(
    x = "Cluster",
    y = "Predicted probability for target (OS2ys)"
  ) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18, color = "black"),
    legend.position = "none"
  )

print(p_prob_cluster)

#mean
if (length(cont_vars) > 0) {
  
  raw_cont_mean <- aggregate(X_full[, cont_vars, drop = FALSE], list(Cluster = Cluster), mean)
  raw_cont_sd   <- aggregate(X_full[, cont_vars, drop = FALSE], list(Cluster = Cluster), sd)
  
  tab_clinica_vera <- data.frame(
    Cluster  = raw_cont_mean$Cluster, 
    n        = as.numeric(table(Cluster)), 
    OS_2ys   = paste0(round(tapply(Y_full, Cluster, mean) * 100, 1), "%"),
    Mean_IP  = round(tapply(p_final, Cluster, mean), 3) # <-- Qui calcoli la media delle IP
  )
  
  for (cv in cont_vars) {
    tab_clinica_vera[[cv]] <- sprintf("%.3f (±%.3f)", raw_cont_mean[[cv]], raw_cont_sd[[cv]])
  }
  
  print(tab_clinica_vera, row.names = FALSE)
}


#avg abs SHAP -----
abs_shap <- abs(shap_agg)
mean_abs_by_class <- aggregate(abs_shap, by = list(Class = df_plot$OS2ys), FUN = mean)
df_bar_class <- reshape(mean_abs_by_class, direction = "long", varying = colnames(abs_shap), v.names = "Mean_Abs_SHAP", timevar = "Feature", times = colnames(abs_shap))
tot_feature <- ave(df_bar_class$Mean_Abs_SHAP, df_bar_class$Feature, FUN = sum)
df_bar_class$Perc <- (df_bar_class$Mean_Abs_SHAP / tot_feature) * 100
df_bar_class$Label_Perc <- sprintf("%.1f%%", df_bar_class$Perc)
df_bar_class_custom <- df_bar_class
feature_labels <- c(
  "Age" = "Age",
  "Volume_PTV" = "Volume PTV",
  "LAA" = "LAA",
  "riskfactors1" = "Risk factors 1"
)
df_bar_class_custom$Feature <- ifelse(df_bar_class_custom$Feature %in% names(feature_labels),
                                      feature_labels[df_bar_class_custom$Feature],
                                      df_bar_class_custom$Feature)
p_bar_class <- ggplot(df_bar_class_custom, aes(x = reorder(Feature, Mean_Abs_SHAP, sum), y = Mean_Abs_SHAP, fill = Class)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "white", linewidth = 0.4) +
  geom_text(aes(label = Label_Perc), position = position_stack(vjust = 0.5), color = "white", fontface = "bold", size = 6) +
  coord_flip() + 
  scale_fill_manual(values = c("NO" = "#ff8080", "YES" = "#93c47d")) + 
  theme_minimal(base_size = 20) + 
  labs(
    x = "Variable", 
    y = "Mean(|SHAP value|)", 
    fill = "OS 2ys"
  ) +
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )

print(p_bar_class)
mean_abs_by_cl <- aggregate(abs_shap, by = list(Cluster = Cluster), FUN = mean)
df_bar_cl <- reshape(mean_abs_by_cl, direction = "long", varying = colnames(abs_shap), v.names = "Mean_Abs_SHAP", timevar = "Feature", times = colnames(abs_shap))
tot_feature_cl <- ave(df_bar_cl$Mean_Abs_SHAP, df_bar_cl$Feature, FUN = sum)
df_bar_cl$Perc <- (df_bar_cl$Mean_Abs_SHAP / tot_feature_cl) * 100
df_bar_cl$Label_Perc <- ifelse(df_bar_cl$Perc >= 4, sprintf("%.1f%%", df_bar_cl$Perc), "")
df_bar_cl_custom <- df_bar_cl
df_bar_cl_custom$Feature <- ifelse(df_bar_cl_custom$Feature %in% names(feature_labels),
                                   feature_labels[df_bar_cl_custom$Feature],
                                   df_bar_cl_custom$Feature)
p_bar_cl <- ggplot(df_bar_cl_custom, aes(x = reorder(Feature, Mean_Abs_SHAP, sum), y = Mean_Abs_SHAP, fill = Cluster)) +
  geom_bar(stat = "identity", position = "stack", width = 0.65, color = "white", linewidth = 0.4) +
  geom_text(aes(label = Label_Perc), position = position_stack(vjust = 0.5), color = "white", fontface = "bold", size = 3.5) +
  coord_flip() + theme_minimal(base_size = 12) + labs(x = "Variable", y = "Mean(|SHAP value|)", fill = "Cluster")
print(p_bar_cl)




# waterfall plots -----
library(ggplot2)
library(grid)

custom_names <- c(
  "Volume_PTV"   = "Volume PTV",
  "Age"          = "Age",
  "LAA"          = "LAA",
  "riskfactors1" = "Risk factors 1"
)

shap_df_custom <- as.data.frame(shap_agg)
colnames(shap_df_custom) <- ifelse(
  colnames(shap_df_custom) %in% names(custom_names),
  custom_names[colnames(shap_df_custom)],
  colnames(shap_df_custom)
)

logit_p   <- log(p_final / (1 - p_final))
base_val  <- mean(logit_p - rowSums(shap_df_custom))
prob_base <- round(plogis(base_val) * 100, 1)

make_waterfall_cl <- function(cl_id) {
  idx_cl <- which(Cluster == cl_id)
  cl_shap_means_raw <- colMeans(shap_df_custom[idx_cl, , drop = FALSE])
  
  cl_shap_sorted <- cl_shap_means_raw[order(abs(cl_shap_means_raw), decreasing = TRUE)]
  feat_names     <- names(cl_shap_sorted)
  n_feat         <- length(feat_names)
  
  wf_df <- data.frame(
    Feature      = feat_names,
    Contribution = as.numeric(cl_shap_sorted),
    stringsAsFactors = FALSE
  )
  
  wf_df$CumSum <- base_val + cumsum(wf_df$Contribution)
  wf_df$Start  <- c(base_val, head(wf_df$CumSum, -1))
  wf_df$Pos    <- rev(seq_len(n_feat))
  wf_df$Sign   <- ifelse(wf_df$Contribution >= 0, "protective", "risk")
  
  conn_df <- data.frame(
    x    = head(wf_df$Pos, -1) - 0.28,
    xend = tail(wf_df$Pos, -1) + 0.28,
    y    = head(wf_df$CumSum, -1),
    yend = head(wf_df$CumSum, -1)
  )
  
  final_val  <- tail(wf_df$CumSum, 1)
  prob_final <- round(plogis(final_val) * 100, 1)
  col_final  <- ifelse(final_val >= base_val, "#1e88e5", "#ff0d57")
  
  # Data frame per le linee di riferimento con relative etichette per la legenda
  ref_lines <- data.frame(
    xintercept = c(base_val, final_val),
    Line_Type  = factor(
      c("Population baseline (log-odds)", "Cluster mean prediction (log-odds)"),
      levels = c("Population baseline (log-odds)", "Cluster mean prediction (log-odds)")
    )
  )
  line_colors <- c("Population baseline (log-odds)" = "grey40", "Cluster mean prediction (log-odds)" = col_final)
  
  x_limits <- range(c(base_val, final_val, wf_df$Start, wf_df$CumSum))
  x_span   <- diff(x_limits)
  x_limits <- c(x_limits[1] - 0.20 * x_span, x_limits[2] + 0.20 * x_span)
  
  offset_val <- 0.035 * x_span
  wf_df$Text_X <- ifelse(wf_df$Contribution >= 0,
                         pmax(wf_df$Start, wf_df$CumSum) + offset_val,
                         pmin(wf_df$Start, wf_df$CumSum) - offset_val)
  wf_df$HJust  <- ifelse(wf_df$Contribution >= 0, 0, 1)
  wf_df$Col_Val <- ifelse(wf_df$Contribution >= 0, "#1e88e5", "#ff0d57")
  
  ggplot(wf_df) +
    # Linee tratteggiate con estetica mappata per generare la legenda
    geom_vline(data = ref_lines, 
               aes(xintercept = xintercept, color = Line_Type), 
               linetype = "dashed", linewidth = 0.7) +
    geom_segment(data = conn_df, aes(x = y, xend = yend, y = x, yend = xend), 
                 linetype = "dotted", color = "grey55", linewidth = 0.5) +
    geom_rect(aes(xmin = pmin(Start, CumSum), xmax = pmax(Start, CumSum), 
                  ymin = Pos - 0.28, ymax = Pos + 0.28, fill = Sign), 
              color = "white", linewidth = 0.3) +
    geom_text(aes(x = Text_X, y = Pos, label = sprintf("%+.2f", Contribution), hjust = HJust), 
              color = wf_df$Col_Val, fontface = "bold", size = 9) +
    scale_y_continuous(breaks = wf_df$Pos, labels = wf_df$Feature, expand = ggplot2::expansion(mult = c(0.08, 0.08))) +
    scale_x_continuous(limits = x_limits) +
    scale_fill_manual(values = c("protective" = "#1e88e5", "risk" = "#ff0d57")) +
    scale_color_manual(name = NULL, values = line_colors) +
    guides(
      fill  = "none", # Nasconde la legenda protective/risk
      color = guide_legend(override.aes = list(linetype = "dashed", linewidth = 0.8)) # Mostra solo le linee tratteggiate
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "bottom", # Posizionata sotto per evitare sovrapposizioni con le barre
      legend.text = element_text(size = 20, color = "black"),
      legend.key.width = unit(1.2, "cm"),
      legend.margin = ggplot2::margin(t = -2, b = 4),
      axis.line = element_line(color = "grey40", linewidth = 0.5),
      axis.ticks = element_line(color = "grey40", linewidth = 0.5),
      axis.text.y = element_text(face = "plain", size = 20, color = "black"),
      axis.text.x = element_text(face = "plain", size = 20, color = "black"),
      axis.title.x = element_text(face = "plain", size = 20, color = "black", margin = ggplot2::margin(t = 8)),
      plot.title = element_text(face = "bold", size = 23, hjust = 0),
      plot.subtitle = element_text(face = "plain", size = 20, hjust = 0, color = "black", margin = ggplot2::margin(b = 8)),
      plot.margin = ggplot2::margin(t = 12, r = 18, b = 8, l = 16)
    ) +
    labs(
      title = paste0("Cluster ", cl_id),
      #subtitle = paste0("Cluster mean prediction: ", prob_final, "% (Population baseline: ", prob_base, "%)"),
      x = "SHAP value (log-odds)",
      y = NULL
    )
}
valid_cl_levels <- setdiff(levels(Cluster), "0")
wf_list <- lapply(valid_cl_levels, make_waterfall_cl)

n_wf   <- length(wf_list)
ncol_g <- ifelse(n_wf <= 4, 2, 3)
nrow_g <- ceiling(n_wf / ncol_g)


print(wf_list[[1]])
print(wf_list[[2]])
print(wf_list[[3]])
print(wf_list[[4]])
print(wf_list[[5]])




#heatmaps -----
library(ggplot2)

cat_vars <- c("riskfactors1")
custom_names <- c(
  "Volume_PTV"   = "Volume PTV",
  "Age"          = "Age",
  "LAA"          = "LAA",
  "riskfactors1" = "Risk factors 1"
)

Cluster <- factor(cl_B$cluster)
os_summary <- aggregate(Y_full, by = list(Cluster = Cluster), FUN = function(x) round(mean(x, na.rm = TRUE) * 100, 1))
os_map <- setNames(os_summary$x, os_summary$Cluster)

cl_order_desc <- os_summary$Cluster[order(os_summary$x, decreasing = TRUE)]
cl_levels_y <- rev(cl_order_desc)

label_map <- setNames(
  paste0("Cluster ", os_summary$Cluster, "\n(", os_summary$x, "%)"),
  os_summary$Cluster
)
ordered_labels_y <- label_map[as.character(cl_levels_y)]

num_vars <- intersect(cont_vars, colnames(df_sub))
raw_num  <- df_sub[, num_vars, drop = FALSE]

raw_means <- aggregate(raw_num, by = list(Cluster = Cluster), FUN = mean)
z_means   <- aggregate(scale(raw_num), by = list(Cluster = Cluster), FUN = mean)
hm_raw_long <- reshape(raw_means, direction = "long", varying = num_vars, v.names = "Raw_Mean", timevar = "Feature", times = num_vars)
hm_z_long   <- reshape(z_means, direction = "long", varying = num_vars, v.names = "Z_Score", timevar = "Feature", times = num_vars)
hm_df <- merge(hm_raw_long[, c("Cluster", "Feature", "Raw_Mean")], hm_z_long[, c("Cluster", "Feature", "Z_Score")], by = c("Cluster", "Feature"))
hm_df$Cell_Label <- sprintf("%.3f", hm_df$Raw_Mean)
clusters_vec <- sort(unique(Cluster))

rf_list <- list()
for (cat_v in cat_vars) {
  if (cat_v %in% colnames(df_sub)) {
    rf_list[[cat_v]] <- data.frame(
      Cluster = clusters_vec, Feature = cat_v, Raw_Mean = NA, Z_Score = NA,
      Cell_Label = sapply(clusters_vec, function(cl) {
        tab <- table(df_sub[[cat_v]][Cluster == cl])
        if(length(tab) == 0) return("N/A")
        sprintf("Mode: %s\n(%d%%)", names(which.max(tab)), round(max(tab) / sum(tab) * 100, 0))
      })
    )
  }
}

if (length(rf_list) > 0) {
  hm_df_final <- rbind(hm_df, do.call(rbind, rf_list))
  all_features <- c(num_vars, cat_vars)
} else {
  hm_df_final <- hm_df
  all_features <- num_vars
}

hm_df_final$Cluster_Label <- factor(label_map[as.character(hm_df_final$Cluster)], levels = ordered_labels_y)

hm_df_final$Feature_Label <- ifelse(
  as.character(hm_df_final$Feature) %in% names(custom_names),
  custom_names[as.character(hm_df_final$Feature)],
  as.character(hm_df_final$Feature)
)
ordered_feat_labels <- ifelse(all_features %in% names(custom_names), custom_names[all_features], all_features)
hm_df_final$Feature_Label <- factor(hm_df_final$Feature_Label, levels = ordered_feat_labels)

robust_scale <- function(x) (x - median(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
raw_medians <- aggregate(raw_num, by = list(Cluster = Cluster), FUN = median)
z_medians   <- aggregate(as.data.frame(lapply(raw_num, robust_scale), check.names = FALSE), by = list(Cluster = Cluster), FUN = median)

colnames(raw_medians) <- c("Cluster", num_vars)
colnames(z_medians)   <- c("Cluster", num_vars)

hm_med_long  <- reshape(raw_medians, direction = "long", varying = num_vars, v.names = "Raw_Median", timevar = "Feature", times = num_vars)
hm_zmed_long <- reshape(z_medians, direction = "long", varying = num_vars, v.names = "Z_Robust", timevar = "Feature", times = num_vars)
hm_df_med <- merge(hm_med_long[, c("Cluster", "Feature", "Raw_Median")], hm_zmed_long[, c("Cluster", "Feature", "Z_Robust")], by = c("Cluster", "Feature"))

hm_df_med$Cell_Label <- mapply(function(cl, feat) {
  val_cl <- df_sub[[feat]][Cluster == cl]
  sprintf("%.3f\n[%.3f - %.3f]", median(val_cl, na.rm = TRUE), quantile(val_cl, 0.25, na.rm = TRUE), quantile(val_cl, 0.75, na.rm = TRUE))
}, hm_df_med$Cluster, hm_df_med$Feature)

if (length(rf_list) > 0) {
  rf_df_med_list <- list()
  for (cat_v in cat_vars) {
    if (cat_v %in% names(rf_list)) {
      rf_df_med_list[[cat_v]] <- data.frame(
        Cluster = clusters_vec, Feature = cat_v, Raw_Median = NA, Z_Robust = NA, Cell_Label = rf_list[[cat_v]]$Cell_Label
      )
    }
  }
  hm_df_med_final <- rbind(hm_df_med, do.call(rbind, rf_df_med_list))
} else {
  hm_df_med_final <- hm_df_med
}

hm_df_med_final$Cluster_Label <- factor(label_map[as.character(hm_df_med_final$Cluster)], levels = ordered_labels_y)
hm_df_med_final$Feature_Label <- ifelse(
  as.character(hm_df_med_final$Feature) %in% names(custom_names),
  custom_names[as.character(hm_df_med_final$Feature)],
  as.character(hm_df_med_final$Feature)
)
hm_df_med_final$Feature_Label <- factor(hm_df_med_final$Feature_Label, levels = ordered_feat_labels)

max_z_common <- max(abs(c(hm_df_final$Z_Score, hm_df_med_final$Z_Robust)), na.rm = TRUE)
common_limits <- c(-max_z_common, max_z_common)

# Stampa della Heatmap corretta con Feature_Label (inclusa la moda di riskfactors1) e carattere 20
p_heatmap <- ggplot(hm_df_final, aes(x = Feature_Label, y = Cluster_Label, fill = Z_Score)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = Cell_Label), size = 6, fontface = "bold") +
  scale_fill_gradient2(low = "#1e88e5", mid = "white", high = "#D55E00", 
                       midpoint = 0, limits = common_limits, na.value = "grey94", name = "Z-score") +
  theme_minimal(base_size = 20) + 
  labs(x = "Variable", y = "Cluster") +
  theme(
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )

print(p_heatmap)

hm_df_med_numeric <- subset(hm_df_med_final, Feature %in% num_vars)

p_heatmap_med <- ggplot(hm_df_med_numeric, aes(x = Feature_Label, y = Cluster_Label, fill = Z_Robust)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = Cell_Label), size = 6, fontface = "bold") +
  scale_fill_gradient2(low = "#1e88e5", mid = "white", high = "#D55E00", 
                       midpoint = 0, limits = common_limits, na.value = "grey94", name = "Z-score") +
  theme_minimal(base_size = 20) + 
  labs(x = "Variable", y = "Cluster") +
  theme(
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )

print(p_heatmap)
print(p_heatmap_med)



#plot OS2ys distribution -----
library(ggplot2)
library(patchwork)

ordine_cluster_desiderato <- c(2, 1, 3, 5, 4)

df_bar <- as.data.frame(table(Cluster = Cluster, OS2ys = factor(df_sub$OS2ys, levels = c("NO", "YES"))))
tot_cl <- tapply(df_bar$Freq, df_bar$Cluster, sum)
df_bar$Total_N <- as.numeric(tot_cl[as.character(df_bar$Cluster)])
df_bar$Perc <- (df_bar$Freq / df_bar$Total_N) * 100

df_bar$Cluster <- factor(df_bar$Cluster, levels = ordine_cluster_desiderato)
df_bar <- df_bar[order(df_bar$Cluster), ]

etich_cluster <- paste0("Cluster ", ordine_cluster_desiderato, "\n(n = ", as.numeric(tot_cl[as.character(ordine_cluster_desiderato)]), ")")
df_bar$Cluster_Label <- factor(paste0("Cluster ", df_bar$Cluster, "\n(n = ", df_bar$Total_N, ")"), levels = etich_cluster)

custom_colors <- c("NO" = "#ff8080", "YES" = "#93c47d")

p_bar_perc <- ggplot(df_bar, aes(x = Cluster_Label, y = Perc, fill = OS2ys)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", Perc)), position = position_stack(vjust = 0.5), color = "white", fontface = "bold", size = 4.5) +
  scale_fill_manual(values = custom_colors, labels = c("NO" = "No", "YES" = "Yes")) +
  theme_minimal(base_size = 11) + 
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(), 
    axis.text.x = element_text(face = "bold", size = 9.5)
  ) +
  labs(title = NULL, x = NULL, y = "Percentage (%)")

p_bar_count <- ggplot(df_bar, aes(x = Cluster_Label, y = Freq, fill = OS2ys)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.6) +
  geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), color = "white", fontface = "bold", size = 4.5) +
  scale_fill_manual(values = custom_colors, labels = c("NO" = "No", "YES" = "Yes")) +
  theme_minimal(base_size = 11) + 
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(), 
    axis.text.x = element_text(face = "bold", size = 9.5)
  ) +
  labs(title = NULL, x = NULL, y = "Number of patients")

combined_plot <- (p_bar_perc | p_bar_count) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )

print(combined_plot)

#clustering glm/glmm -----
ref_cluster <- "1"

Cluster <- factor(cl_B$cluster)

if (length(cont_vars) > 0) {
  raw_cont_mean <- aggregate(X_full[, cont_vars, drop = FALSE], list(Cluster = Cluster), mean)
  raw_cont_sd   <- aggregate(X_full[, cont_vars, drop = FALSE], list(Cluster = Cluster), sd)
  tab_clinica_vera <- data.frame(Cluster = raw_cont_mean$Cluster, n = as.numeric(table(Cluster)), OS_2ys = paste0(round(tapply(Y_full, Cluster, mean) * 100, 1), "%"))
  for (cv in cont_vars) tab_clinica_vera[[cv]] <- sprintf("%.3f (±%.3f)", raw_cont_mean[[cv]], raw_cont_sd[[cv]])
  print(tab_clinica_vera, row.names = FALSE)
}


cat_vars <- setdiff(covariates, cont_vars)
if (length(cat_vars) > 0) {
  for (catv in cat_vars) {
    cat(paste("\nVariabile:", catv, "\n"))
    print(round(prop.table(table(Cluster = Cluster, Category = X_full[[catv]]), margin = 1) * 100, 1))
  }
}


print(table(Cluster, Z_full))
res_chi <- chisq.test(table(Cluster, Z_full), simulate.p.value = TRUE, B = 10000)
round(res_chi$stdres, 2)
print(res_chi)


idx_valid_cl <- cl_B$cluster > 0
cl_valid <- droplevels(factor(cl_B$cluster[idx_valid_cl]))
cl_valid <- relevel(cl_valid, ref = ref_cluster)

d_test <- data.frame(Y = Y_full[idx_valid_cl], cl = cl_valid, Z = droplevels(Z_full[idx_valid_cl]))

cat(sprintf("\n--- TEST 1: Regressione Logistica Standard (Senza Ospedale, Reference = Cluster %s) ---\n", ref_cluster))
m_std0 <- glm(Y ~ 1, family = binomial(), data = d_test)
m_std1 <- glm(Y ~ cl, family = binomial(), data = d_test)
print(anova(m_std0, m_std1, test = "Chisq"))

coef_std1 <- as.data.frame(round(summary(m_std1)$coefficients, 4))
coef_std1[["Pr(>|z|)"]] <- sprintf("%.3e", summary(m_std1)$coefficients[, 4])
print(coef_std1)

cat(sprintf("\n--- TEST 2: GLMM a Effetti Misti (Corretto per Ospedale, Reference = Cluster %s) ---\n", ref_cluster))
m_me0 <- glmer(Y ~ 1 + (1 | Z), family = binomial(), data = d_test)
m_me1 <- suppressWarnings(glmer(Y ~ cl + (1 | Z), family = binomial(), data = d_test))
print(anova(m_me0, m_me1))
coef_me1 <- as.data.frame(round(summary(m_me1)$coefficients, 4))
coef_me1[["Pr(>|z|)"]] <- sprintf("%.3e", summary(m_me1)$coefficients[, 4])
print(coef_me1)


#clustering cox -----
df_sub$Cluster <- Cluster
df_sub$OS[df_sub$ID == "LM1001"] <- 24
df_sub$event_death_24m <- 1 - df_sub$OS2ys_num
df_sub$OS_24m <- pmin(df_sub$OS, 24)

df_surv_valid <- df_sub[df_sub$Cluster != "0", ]
df_surv_valid$Cluster <- droplevels(df_surv_valid$Cluster)
df_surv_valid$Cluster <- relevel(df_surv_valid$Cluster, ref = ref_cluster)

cox_mod <- coxph(Surv(OS_24m, event_death_24m) ~ Cluster, data = df_surv_valid)
cox_sum <- summary(cox_mod)
tab_cox <- data.frame(
  Confronto = rownames(cox_sum$conf.int),
  HR = round(cox_sum$conf.int[, "exp(coef)"], 4),
  CI_95 = sprintf("[%.4f - %.4f]", cox_sum$conf.int[, "lower .95"], cox_sum$conf.int[, "upper .95"]),
  p_value = signif(cox_sum$coefficients[, "Pr(>|z|)"], 4)
)
print(tab_cox, row.names = FALSE)


#clustering km ----

fit_km    <- survfit(Surv(OS_24m, event_death_24m) ~ Cluster, data = df_surv_valid)
lr_test   <- survdiff(Surv(OS_24m, event_death_24m) ~ Cluster, data = df_surv_valid)
p_logrank <- 1 - pchisq(lr_test$chisq, length(lr_test$n) - 1)

km_df <- data.frame(time = fit_km$time, surv = fit_km$surv * 100, Cluster = rep(gsub("Cluster=", "Cluster ", names(fit_km$strata)), fit_km$strata))
km_df <- rbind(data.frame(time = 0, surv = 100, Cluster = unique(km_df$Cluster)), km_df)

p_km_24m <- ggplot(km_df, aes(x = time, y = surv, color = Cluster)) +
  geom_step(linewidth = 1.2) + scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) + scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 3)) +
  theme_minimal(base_size = 12) + theme(legend.position = "top", panel.grid.minor = element_blank()) +
  labs(title = "Kaplan-Meier curves in the first 24 months by Cluster", subtitle = sprintf("Log-Rank Test: Chisq = %.1f, p = %.2e", lr_test$chisq, p_logrank), x = "months from RT", y = "Survival estimate (%)")

library(survminer)

library(survminer)

fit_km    <- survfit(Surv(OS_24m, event_death_24m) ~ Cluster, data = df_surv_valid)
lr_test   <- survdiff(Surv(OS_24m, event_death_24m) ~ Cluster, data = df_surv_valid)
p_logrank <- 1 - pchisq(lr_test$chisq, length(lr_test$n) - 1)

nomi_cluster <- sub(".*=", "", names(fit_km$strata)) 

p_km <- ggsurvplot(
  fit_km,
  data = df_surv_valid,
  conf.int = TRUE,
  conf.int.alpha = 0.15,
  xlim = c(0, 24),
  break.x.by = 3,
  ylim = c(0, 1),
  pval = FALSE,         #sprintf("p = %.2e", p_logrank),
  risk.table = TRUE,
  risk.table.fontsize = 6,
  risk.table.height = 0.35,
  legend = "top",                  
  legend.title = "Cluster",
  legend.labs = nomi_cluster,      
  font.legend = c(20, "bold"),     
  font.title = c(20, "bold"),
  font.x = c(20, "plain"),
  font.y = c(20, "plain"),
  font.tickslab = c(20, "plain"),
  xlab = "Months from the beginning of radiotherapy",
  ylab = "Survival probability",
  #title = "Kaplan-Meier curves by Cluster",
  ggtheme = theme_minimal(base_size = 20),
  tables.theme = theme_minimal(base_size = 20)
)

p_km$plot <- p_km$plot + 
  theme(
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 20, face = "bold")
  )

p_km$table <- p_km$table +
  theme(
    text = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title = element_text(size = 20),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

print(p_km)




#clustering OOF comparison ----
library(ggplot2)
library(patchwork)

Cluster_oof <- factor(cl_oof$cluster)
df_oof_cmp <- data.frame(
  Full_UMAP1   = um2_shap[, 1],
  Full_UMAP2   = um2_shap[, 2],
  Cluster_Full = factor(cl_B$cluster),
  OOF_UMAP1    = um2_oof[, 1],
  OOF_UMAP2    = um2_oof[, 2],
  Cluster_OOF  = factor(cl_oof$cluster),
  OS2ys        = factor(df_sub$OS2ys, levels = c("NO", "YES"))
)

p1 <- ggplot(df_oof_cmp, aes(x = Full_UMAP1, y = Full_UMAP2, color = Cluster_Full)) + 
  geom_point(alpha = 0.85, size = 2.2) + 
  theme_minimal(base_size = 15) + 
  labs(x = "UMAP 1", y = "UMAP 2", color = "Full-cohort cluster")

p2 <- ggplot(df_oof_cmp, aes(x = OOF_UMAP1, y = OOF_UMAP2, color = Cluster_OOF)) + 
  geom_point(alpha = 0.85, size = 2.2) + 
  theme_minimal(base_size = 15) + 
  labs(x = "UMAP 1", y = "UMAP 2", color = "Cross-Validated cluster")

p3 <- ggplot(df_oof_cmp, aes(x = OOF_UMAP1, y = OOF_UMAP2, color = Cluster_Full, shape = Cluster_OOF)) + 
  geom_point(alpha = 0.85, size = 2.4) + 
  theme_minimal(base_size = 15) + 
  labs(x = "UMAP 1", y = "UMAP 2", color = "Full-cohort cluster", shape = "Cross-Validated cluster")

p_combined <- (p1 | p2 | p3) & theme(legend.position = "bottom", legend.box = "vertical")

print(p_combined)



#kruskal - tests to assess shifts ----
library(dplyr)
library(purrr)
df_plot <- data.frame(
  SHAP_UMAP1 = um2_shap[, 1],
  SHAP_UMAP2 = um2_shap[, 2],
  RAW_UMAP1  = um2_raw[, 1],
  RAW_UMAP2  = um2_raw[, 2],
  Cluster    = Cluster,
  OS2ys      = factor(df_sub$OS2ys, levels = c("NO", "YES")),
  OS_months  = df_sub$OS,
  Pred_Prob  = p_final
)
for (v in cont_vars) {
  df_plot[[v]] <- df_sub[[v]]
}

num_vars <- intersect(cont_vars, colnames(df_plot))
kruskal_results <- map_dfr(num_vars, function(v) {
  test <- kruskal.test(df_plot[[v]] ~ factor(cl_B$cluster))
  data.frame(
    Variable = v,
    Chi_Squared = round(test$statistic, 3),
    DF = test$parameter,
    P_Value = signif(test$p.value, 3)
  )
})

print(kruskal_results)

library(purrr)



tab_train <- table(df_sub$riskfactors1, cl_B$cluster)
tab_train
chi_test_sim <- chisq.test(tab_train, simulate.p.value = TRUE, B = 10000)

print(chi_test_sim)
print(round(chi_test_sim$residuals, 2))
print(round(chi_test_sim$stdres, 2))



get_mode <- function(v) {
  v <- v[!is.na(v)] 
  if(length(v) == 0) return(NA)
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
get_mode <- function(v) {
  v <- v[!is.na(v)]
  if (length(v) == 0) return(NA)
  
  tb <- table(v)
  if (length(tb) == 0) return(NA)
  names(which.max(tb))
}
tab = tapply(df_sub$riskfactors1, cl_B$cluster, get_mode)
tab
sprintf("Mode: %s\n(%d%%)", names(which.max(tab)), round(max(tab) / sum(tab) * 100, 0))


















#testing -----
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



#test: plots eda -----
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
library(ggplot2)

df_train_plot <- data.frame(
  Age = as.numeric(data$Age),
  Volume_PTV = as.numeric(data$Volume_PTV),
  LAA = as.numeric(data$LAA),
  Risk_Factor_1 = as.factor(data$riskfactors1)
)

patient = df_train_plot[,c("Age", "Volume_PTV", "LAA", "Risk_Factor_1")]
rows= complete.cases(patient) 
df_train_plot = df_train_plot[rows,]


df_test_plot <- data.frame(
  Age = as.numeric(datatest$Age),
  Volume_PTV = as.numeric(datatest$VOLUME_PTV),
  LAA = as.numeric(datatest$`LAA%`),
  Risk_Factor_1 = as.factor(datatest$risk_factor1)
)

p_age_test <- ggplot(df_test_plot, aes(x = Age)) +
  geom_histogram(aes(y = after_stat(density)), fill = "#4682B4", color = "white", bins = 30, alpha = 0.85) +
  theme_minimal(base_size = 20) +
  labs(x = "Age", y = "Density", title = NULL)

p_age_train <- ggplot(df_train_plot, aes(x = Age)) +
  geom_histogram(aes(y = after_stat(density)), fill = "#4682B4", color = "white", bins = 30, alpha = 0.85) +
  theme_minimal(base_size = 20) +
  labs(x = "Age", y = "Density", title = NULL)

p_vol_test <- ggplot(df_test_plot, aes(x = Volume_PTV)) +
  geom_histogram(aes(y = after_stat(density)), fill = "#4682B4", color = "white", bins = 30, alpha = 0.85) +
  theme_minimal(base_size = 20) +
  labs(x = "Volume PTV", y = "Density", title = NULL)

p_vol_train <- ggplot(df_train_plot, aes(x = Volume_PTV)) +
  geom_histogram(aes(y = after_stat(density)), fill = "#4682B4", color = "white", bins = 30, alpha = 0.85) +
  theme_minimal(base_size = 20) +
  labs(x = "Volume PTV", y = "Density", title = NULL)

p_laa_test <- ggplot(df_test_plot, aes(x = LAA)) +
  geom_histogram(aes(y = after_stat(density)), fill = "#4682B4", color = "white", bins = 30, alpha = 0.85) +
  theme_minimal(base_size = 20) +
  labs(x = "LAA", y = "Density", title = NULL)

p_laa_train <- ggplot(df_train_plot, aes(x = LAA)) +
  geom_histogram(aes(y = after_stat(density)), fill = "#4682B4", color = "white", bins = 30, alpha = 0.85) +
  theme_minimal(base_size = 20) +
  labs(x = "LAA", y = "Density", title = NULL)

p_risk_train <- ggplot(df_train_plot[!is.na(df_train_plot$Risk_Factor_1), ], aes(x = Risk_Factor_1)) +
  geom_bar(aes(y = after_stat(prop), group = 1), fill = "#4682B4", color = "white", alpha = 0.85) +
  geom_text(stat = "count", aes(y = after_stat(prop), group = 1, label = sprintf("%.1f%%\n(n=%d)", after_stat(prop) * 100, after_stat(count))), vjust = -0.2, size = 6, lineheight = 0.9) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme_minimal(base_size = 20) +
  labs(x = "Risk factors 1", y = "Proportion", title = NULL)

p_risk_test <- ggplot(df_test_plot[!is.na(df_test_plot$Risk_Factor_1), ], aes(x = Risk_Factor_1)) +
  geom_bar(aes(y = after_stat(prop), group = 1), fill = "#4682B4", color = "white", alpha = 0.85) +
  geom_text(stat = "count", aes(y = after_stat(prop), group = 1, label = sprintf("%.1f%%\n(n=%d)", after_stat(prop) * 100, after_stat(count))), vjust = -0.2, size = 6, lineheight = 0.9) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme_minimal(base_size = 20) +
  labs(x = "Risk factors 1", y = "Proportion", title = NULL)


p_age_train / p_age_test
p_vol_train / p_vol_test
p_laa_train / p_laa_test
p_risk_train / p_risk_test





#test plots eda target -----
library(patchwork)
library(scales)

datatestos$hospital <- factor("Test") 

if("OS2ys" %in% names(datatestos)) {
  if(is.numeric(datatestos$OS2ys) || is.integer(datatestos$OS2ys)) {
    datatestos$OS2ys_fac <- factor(datatestos$OS2ys, levels = c(0, 1), labels = c("NO", "YES"))
  } else {
    datatestos$OS2ys_fac <- as.factor(datatestos$OS2ys)
  }
} else {
  datatestos$OS2ys_fac <- NA
}

if(is.numeric(data$OS2ys) || is.integer(data$OS2ys)) {
  data$OS2ys_fac <- factor(data$OS2ys, levels = c(0, 1), labels = c("NO", "YES"))
} else {
  data$OS2ys_fac <- as.factor(data$OS2ys)
}

data$hospital <- as.factor(data$hospital)
data_combined <- bind_rows(
  data %>% select(any_of(c("hospital", "OS", "OS2ys_fac"))),
  datatestos %>% select(any_of(c("hospital", "OS", "OS2ys_fac")))
)

hospital_levels <- c(as.character(1:7), "Test")
data_combined$hospital <- factor(data_combined$hospital, levels = hospital_levels)
hp_palette <- c("#ff9999", "#99ccff", "#b6d7a8", "#d5a6bd", "#f9cb9c", "#ffe599", "#d9d2e9", "#aaaaaa")

df_counts <- data_combined %>%
  filter(!is.na(hospital)) %>%
  group_by(hospital) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(hospital = factor(hospital, levels = hospital_levels))

p_count <- ggplot(df_counts, aes(x = hospital, y = n, fill = hospital)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = -0.5, size = 6, fontface = "bold", color = "black") +
  scale_fill_manual(values = hp_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "plain", size = 20),
    axis.text = element_text(size = 16)
  ) +
  labs(x = "Hospital", y = "Number of patients")

print(p_count)

p1 <- ggplot(data_combined, aes(x = hospital, y = OS, fill = hospital)) +
  geom_boxplot(alpha = 0.85, outlier.shape = NA, color = "black", linewidth = 0.5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7, color = "black") +
  scale_fill_manual(values = hp_palette) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "plain", size = 20, hjust = 0),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "plain", size = 20),
    axis.text = element_text(size = 16)
  ) +
  labs(x = "Hospital", y = "OS")

print(p1)

df_comparison <- bind_rows(
  data.frame(Cohort = "Train", OS = data$OS),
  data.frame(Cohort = "Test", OS = datatestos$OS)
)

p_comparison <- ggplot(df_comparison, aes(x = Cohort, y = OS, fill = Cohort)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.4, color = "black") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
  scale_fill_manual(values = c("Train" = "#99ccff", "Test" = "#f9cb9c")) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16)
  ) +
  labs(x="", y = "OS")

print(p_comparison)

tab_train <- data %>%
  filter(!is.na(OS2ys_fac)) %>%
  count(OS2ys_fac) %>%
  mutate(Percentage = paste0(round(n / sum(n) * 100, 1), "%"))
print(tab_train)

tab_test <- datatestos %>%
  filter(!is.na(OS2Y_1)) %>%
  count(OS2Y_1) %>%
  mutate(Percentage = paste0(round(n / sum(n) * 100, 1), "%"))
print(tab_test)

library(ggplot2)
library(dplyr)
library(scales)

datatestos$hospital <- factor("Test")
datatestos$OS2ys <- datatestos$OS2Y_1

if("OS2ys" %in% names(datatestos)) {
  if(is.numeric(datatestos$OS2ys) || is.integer(datatestos$OS2ys)) {
    datatestos$OS2ys_fac <- factor(datatestos$OS2ys, levels = c(0, 1), labels = c("NO", "YES"))
  } else {
    datatestos$OS2ys_fac <- toupper(as.character(datatestos$OS2ys))
    datatestos$OS2ys_fac <- factor(datatestos$OS2ys_fac, levels = c("NO", "YES"))
  }
} else {
  datatestos$OS2ys_fac <- NA
}

if(is.numeric(data$OS2ys) || is.integer(data$OS2ys)) {
  data$OS2ys_fac <- factor(data$OS2ys, levels = c(0, 1), labels = c("NO", "YES"))
} else {
  data_vals <- toupper(as.character(data$OS2ys))
  data$OS2ys_fac <- factor(data_vals, levels = c("NO", "YES"))
}

data$hospital <- as.factor(data$hospital)
data_combined <- bind_rows(
  data %>% select(any_of(c("hospital", "OS2ys_fac"))),
  datatestos %>% select(any_of(c("hospital", "OS2ys_fac")))
)

hospital_levels <- c(as.character(1:7), "Test")
data_combined$hospital <- factor(data_combined$hospital, levels = hospital_levels)
hp_palette <- c("#ff9999", "#99ccff", "#b6d7a8", "#d5a6bd", "#f9cb9c", "#ffe599", "#d9d2e9", "#aaaaaa")

df_prop <- data_combined %>%
  filter(!is.na(OS2ys_fac) & !is.na(hospital)) %>%
  group_by(hospital, OS2ys_fac) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(hospital) %>%
  mutate(
    prop = count / sum(count),
    label = paste0(round(prop * 100), "%")
  )

p_prop_combined <- ggplot(data_combined, aes(x = hospital, fill = OS2ys_fac)) +
  geom_bar(position = "fill", color = "white", linewidth = 0.8) +
  geom_text(data = df_prop, 
            aes(y = prop, label = label, group = OS2ys_fac), 
            position = position_fill(vjust = 0.5), 
            size = 6, fontface = "bold", color = "black") +
  scale_fill_manual(values = c("NO" = "#ff8080", "YES" = "#93c47d")) + 
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "plain", size = 20, hjust = 0),        
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 20),                 
    legend.text = element_text(size = 16),                                   
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "plain", size = 20),                  
    axis.text = element_text(size = 16)                                    
  ) +
  labs(x = "Hospital", y = "Proportion", fill = "OS2ys")

print(p_prop_combined)


data$Cohort <- "Train"
datatestos$Cohort <- "Test"

data_combined <- bind_rows(
  data %>% select(any_of(c("Cohort", "OS2ys"))),
  datatestos %>% select(any_of(c("Cohort", "OS2ys")))
)

if(is.numeric(data_combined$OS2ys) || is.integer(data_combined$OS2ys)) {
  data_combined$OS2ys_fac <- factor(data_combined$OS2ys, levels = c(0, 1), labels = c("NO", "YES"))
} else {
  data_vals <- toupper(as.character(data_combined$OS2ys))
  data_combined$OS2ys_fac <- factor(data_vals, levels = c("NO", "YES"))
}

df_prop <- data_combined %>%
  filter(!is.na(OS2ys_fac) & !is.na(Cohort)) %>%
  group_by(Cohort, OS2ys_fac) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Cohort) %>%
  mutate(
    prop = count / sum(count),
    label = paste0(round(prop * 100, 1), "%")
  )

p_prop_global <- ggplot(data_combined, aes(x = Cohort, fill = OS2ys_fac)) +
  geom_bar(position = "fill", color = "white", linewidth = 0.8, width = 0.5) +
  geom_text(data = df_prop, 
            aes(y = prop, label = label, group = OS2ys_fac), 
            position = position_fill(vjust = 0.5), 
            size = 6, fontface = "bold", color = "black") +
  scale_fill_manual(values = c("NO" = "#ff8080", "YES" = "#93c47d")) + 
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "plain", size = 20, hjust = 0),        
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 20),                 
    legend.text = element_text(size = 16),                                   
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "plain", size = 20),                  
    axis.text = element_text(size = 16)                                    
  ) +
  labs(x ="", y = "Proportion", fill = "OS2ys")

print(p_prop_global)

final_plot <- p1 / p_prop_combined + plot_layout(heights = c(1, 1))
final_plot_prop <- p_comparison / p_prop_global + plot_layout(heights = c(1, 1))


print(final_plot)
print(final_plot_prop)


#test eda feature comparison tests ----
colnames(df_train_plot)
colnames(df_test_plot)
dim(df_train_plot)
dim(df_test_plot)
patient = df_train_plot[,c("Age", "Volume_PTV", "LAA", "Risk_Factor_1")]
rows= complete.cases(patient) 
datatrain = df_train_plot[rows,]
dim(datatrain)
colnames(datatrain)
datatrain$Cohort <- "Train"
df_test_plot$Cohort <- "Test"
df_comparison_features <- rbind(datatrain, df_test_plot)
df_comparison_features$Cohort <- factor(df_comparison_features$Cohort, levels = c("Train", "Test"))
shapiro.test(df_comparison_features$Age[df_comparison_features$Cohort == "Train"])
shapiro.test(df_comparison_features$Age[df_comparison_features$Cohort == "Test"])
shapiro.test(df_comparison_features$Volume_PTV[df_comparison_features$Cohort == "Train"])
shapiro.test(df_comparison_features$Volume_PTV[df_comparison_features$Cohort == "Test"])
shapiro.test(df_comparison_features$LAA[df_comparison_features$Cohort == "Train"])
shapiro.test(df_comparison_features$LAA[df_comparison_features$Cohort == "Test"])

wt_age <- wilcox.test(Age ~ Cohort, data = df_comparison_features)
cat("Age:\n"); print(wt_age)
wt_vol <- wilcox.test(Volume_PTV ~ Cohort, data = df_comparison_features)
cat("\nVolume PTV:\n"); print(wt_vol)
wt_vol <- wilcox.test(LAA ~ Cohort, data = df_comparison_features)
cat("\nLAA:\n"); print(wt_vol)
tab_risk <- table(df_comparison_features$Cohort, df_comparison_features$Risk_Factor_1)
cat("Tabella di contingenza:\n")
print(tab_risk)
chi_risk <- chisq.test(tab_risk)
cat("\nChi-squared test:\n"); print(chi_risk)


median(datatrain$Age)
median(datatrain$Volume_PTV)
median(datatrain$LAA)
median(df_test_plot$Age)
median(df_test_plot$Volume_PTV)
median(df_test_plot$LAA)


# test eda plots 2 ------
library(ggplot2)
library(dplyr)
library(patchwork)

datatrain$Cohort <- "Train"
df_test_plot$Cohort  <- "Test"
df_comparison_features <- rbind(datatrain, df_test_plot)
df_comparison_features$Cohort <- factor(df_comparison_features$Cohort, levels = c("Train", "Test"))

cohort_colors <- c("Train" = "#99ccff", "Test" = "#f9cb9c")

p_box_age <- ggplot(df_comparison_features, aes(x = Cohort, y = Age, fill = Cohort)) +
  geom_boxplot(alpha = 0.85, outlier.shape = 16, outlier.alpha = 0.5, color = "black") +
  scale_fill_manual(values = cohort_colors) +
  theme_minimal(base_size = 20) +
  theme(legend.position = "none", plot.title = element_blank()) +
  labs(title = NULL, x = "", y = "Age (years)")

p_box_vol <- ggplot(df_comparison_features, aes(x = Cohort, y = Volume_PTV, fill = Cohort)) +
  geom_boxplot(alpha = 0.85, outlier.shape = 16, outlier.alpha = 0.5, color = "black") +
  scale_fill_manual(values = cohort_colors) +
  theme_minimal(base_size = 20) +
  theme(legend.position = "none", plot.title = element_blank()) +
  labs(title = NULL, x = "", y = "Volume PTV (cm³)")

p_box_laa <- ggplot(df_comparison_features, aes(x = Cohort, y = LAA, fill = Cohort)) +
  geom_boxplot(alpha = 0.85, outlier.shape = 16, outlier.alpha = 0.5, color = "black") +
  scale_fill_manual(values = cohort_colors) +
  theme_minimal(base_size = 20) +
  theme(legend.position = "none", plot.title = element_blank()) +
  labs(title = NULL, x = "", y = "LAA (%)")

df_risk_combined <- df_comparison_features %>%
  filter(!is.na(Risk_Factor_1)) %>%
  group_by(Cohort, Risk_Factor_1) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Cohort) %>%
  mutate(prop = count / sum(count))

p_risk_train <- ggplot(df_train_plot[!is.na(df_train_plot$Risk_Factor_1), ], aes(x = Risk_Factor_1)) +
  geom_bar(aes(y = after_stat(prop), group = 1), fill = "#99ccff", color = "white", alpha = 0.85, linewidth = 0.5) +
  geom_text(stat = "count", aes(y = after_stat(prop), group = 1, label = sprintf("%.1f%%\n(n=%d)", after_stat(prop) * 100, after_stat(count))), vjust = -0.2, size = 5, lineheight = 0.9) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.2))) +
  theme_minimal(base_size = 20) +
  theme(plot.title = element_blank()) +
  labs(title = NULL, x = "Risk factors 1 (Train)", y = "Proportion")


p_risk_test <- ggplot(df_test_plot[!is.na(df_test_plot$Risk_Factor_1), ], aes(x = Risk_Factor_1)) +
  geom_bar(aes(y = after_stat(prop), group = 1), fill = "#f9cb9c", color = "white", alpha = 0.85, linewidth = 0.5) +
  geom_text(stat = "count", aes(y = after_stat(prop), group = 1, label = sprintf("%.1f%%\n(n=%d)", after_stat(prop) * 100, after_stat(count))), vjust = -0.2, size = 5, lineheight = 0.9) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.2))) +
  theme_minimal(base_size = 20) +
  theme(plot.title = element_blank()) +
  labs(title = NULL, x = "Risk factors 1 (Test)", y = "Proportion")

p_box_age
p_box_vol
p_box_laa
p_risk_train
p_risk_test




#test metrics -----
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

hp_xgb     <- list(max_depth = 2, eta = 0.01, lambda = 15, subsample = 0.7, nrounds = 400)
covariates <- c("Age", "Volume_PTV", "LAA", "riskfactors1")
cont_vars  <- c("Age", "Volume_PTV", "LAA")

seed_xgb   <- 123

idx <- complete.cases(data[, covariates]) & !is.na(data$OS2ys_num) & !is.na(data$hospital)
df_sub <- data[idx, ]
X_full <- df_sub[, covariates, drop = FALSE]
Y_full <- df_sub$OS2ys_num
sc_train_params <- scale(X_full[, cont_vars])
center_train <- attr(sc_train_params, "scaled:center")
scale_train  <- attr(sc_train_params, "scaled:scale")

X_full_sc <- X_full
X_full_sc[, cont_vars] <- sc_train_params

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
set.seed(seed_umap)
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
##new mapping ----
cl_B_train <- dbscan::hdbscan(umap_model_clB$embedding, minPts = minPts_hdb)
mapping_manuale <- c("2" = "1", "1" = "2", "3" = "3", "5" = "4", "4" = "5")
map <- factor(mapping_manuale[as.character(cl_B_train$cluster)], levels = c("1", "2", "3", "4", "5"))
cl_B_train$cluster <- as.numeric(as.character(map))

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

datatestos_aligned <- datatestos[match(df_test$ID, datatestos$ID), ]

Y_test <- ifelse(datatestos_aligned$OS2Y_1 == "Yes", 1, 0)

valid_outcomes_idx <- !is.na(Y_test)

metrics_test <- calc_metrics_vec(
  p = p_test[valid_outcomes_idx], 
  Y = Y_test[valid_outcomes_idx]
)

print(round(metrics_test, 4))
treeshap_test <- treeshap(uni_xgb, X_test_mat, verbose = FALSE)
shap_test_raw <- as.matrix(treeshap_test$shaps)
shap_test_agg <- aggregate_shap(shap_test_raw, X_test, covariates)

um2_test <- uwot::umap_transform(shap_test_agg, umap_model_clB)

cl_test_pred <- as.integer(predict(cl_B_train, data = umap_model_clB$embedding, newdata = um2_test))


#test cal plot ----

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
  labs(x = "Predicted probability", y = "Observed probability") +
  theme_classic(base_size = 20) +
  theme(legend.position = "none")

print(p_cal_test)


#test SHAP vs raw plot ------
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
}

cat_vars <- setdiff(covariates, cont_vars)
if (length(cat_vars) > 0) {
  for (v in cat_vars) {
    tab <- table(Cluster = df_prof_valid$Cluster, Set = df_prof_valid$Set, Value = df_prof_valid[[v]])
    prop_tab <- prop.table(tab, margin = c(1, 2)) * 100
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
  
  if(a <= 72.5 & p <= 289.37) { cl_test_rules[i] <- 1 }
  else if(a > 71.5 & l <= 0.908 & p <= 288.16) { cl_test_rules[i] <- 2 }
  else if(a <= 75.5 & l <= 1.511 & p > 289.86) { cl_test_rules[i] <- 3 }
  else if(a > 75.5 & l <= 1.576 & p > 299.42) { cl_test_rules[i] <- 4 }
  else if(l > 1.511 & p > 289.59) { cl_test_rules[i] <- 5 }
}


test_results <- data.frame(
  ID             = df_test$ID,
  Prob_OS2ys     = round(p_test, 4),
  Pred_Class     = pred_class_test,
  Cluster_clB    = cl_test_pred,
  UMAP1          = round(um2_test[, 1], 3),
  UMAP2          = round(um2_test[, 2], 3)
)
test_results$Cluster_Rules <- cl_test_rules

df_plot_train <- data.frame(
  UMAP1   = umap_model_clB$embedding[, 1],
  UMAP2   = umap_model_clB$embedding[, 2],
  Cluster = factor(ifelse(cl_B_train$cluster == 0, "0", as.character(cl_B_train$cluster)), levels = c("0", "1", "2", "3", "4", "5")),
  Set     = "train"
)

df_plot_test_hdb <- data.frame(
  UMAP1   = um2_test[, 1],
  UMAP2   = um2_test[, 2],
  Cluster = factor(ifelse(cl_test_pred == 0, "0", as.character(cl_test_pred)), levels = c("0", "1", "2", "3", "4", "5")),
  Set     = "test"
)

df_plot_test_rules <- data.frame(
  UMAP1   = um2_test[, 1],
  UMAP2   = um2_test[, 2],
  Cluster = factor(ifelse(cl_test_rules == 0, "0", as.character(cl_test_rules)), levels = c("0", "1", "2", "3", "4", "5")),
  Set     = "test"
)

cluster_colors <- c(
  "0" = "gray60", 
  "1" = "#F8766D", 
  "2" = "#A3A500", 
  "3" = "#00BF7D", 
  "4" = "#00B0F6", 
  "5" = "#E76BF3"
)


p_train_only <- ggplot(df_plot_train, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(alpha = 0.85, size = 2) +
  scale_color_manual(values = cluster_colors, breaks = c("1", "2", "3", "4", "5")) + 
  theme_minimal(base_size = 20) +
  theme(legend.position = "bottom", legend.key.width = unit(1.5, "cm"), plot.title = element_text(face = "bold")) +
  guides(color = guide_legend(nrow = 1, title.position = "left", title.vjust = 0.5, override.aes = list(size = 5))) +
  labs(title = "Training clustering", x = "UMAP 1 (SHAP)", y = "UMAP 2 (SHAP)", color = "Cluster:")

p_test_hdb <- ggplot(df_plot_test_hdb, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(alpha = 0.95, size = 2.4, shape = 17) +
  scale_color_manual(values = cluster_colors, drop = FALSE) + # <-- AGGIUNTA
  theme_minimal(base_size = 20) +
  theme(legend.position = "bottom", legend.key.width = unit(1.5, "cm"), plot.title = element_text(face = "bold")) +
  guides(color = guide_legend(nrow = 1, title.position = "left", title.vjust = 0.5, override.aes = list(size = 5))) +
  labs(title = "Test clustering (HDBSCAN)", x = "UMAP 1 (SHAP)", y = "UMAP 2 (SHAP)", color = "Cluster:")

p_test_rules <- ggplot(df_plot_test_rules, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(alpha = 0.95, size = 2.4, shape = 17) +
  scale_color_manual(values = cluster_colors, drop = FALSE) + # <-- AGGIUNTA
  theme_minimal(base_size = 20) +
  theme(legend.position = "bottom", legend.key.width = unit(1.5, "cm"), plot.title = element_text(face = "bold")) +
  guides(color = guide_legend(nrow = 1, title.position = "left", title.vjust = 0.5, override.aes = list(size = 5))) +
  labs(title = "Test clustering (Skope-Rules)", x = "UMAP 1 (SHAP)", y = "UMAP 2 (SHAP)", color = "Cluster:")

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 3)))
print(p_train_only,  vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_test_hdb,    vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p_test_rules,  vp = viewport(layout.pos.row = 1, layout.pos.col = 3))




#test shap vs raw - train/test sovrapposti ----
library(ggplot2)
library(grid)

p_train_only <- ggplot(df_plot_train, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(alpha = 0.85, size = 2) +
  scale_color_manual(values = cluster_colors, breaks = c("1", "2", "3", "4", "5")) + 
  theme_minimal(base_size = 16) +
  theme(legend.position = "bottom", legend.key.width = unit(1, "cm"), plot.title = element_text(face = "bold")) +
  guides(color = guide_legend(nrow = 1, title.position = "left", title.vjust = 0.5, override.aes = list(size = 4))) +
  labs(title = "Training clustering", x = "UMAP 1", y = "UMAP 2", color = "Cluster:")

p_test_hdb_bg <- ggplot() +
  geom_point(data = df_plot_train, aes(x = UMAP1, y = UMAP2, color = Cluster), alpha = 0.2, size = 1.2) +
  geom_point(data = df_plot_test_hdb, aes(x = UMAP1, y = UMAP2, color = Cluster), alpha = 0.9, size = 2.2, shape = 17) +
  scale_color_manual(values = cluster_colors, drop = FALSE) + 
  theme_minimal(base_size = 16) +
  theme(legend.position = "bottom", legend.key.width = unit(1, "cm"), plot.title = element_text(face = "bold")) +
  guides(color = guide_legend(nrow = 1, title.position = "left", title.vjust = 0.5, override.aes = list(size = 4, alpha = 1))) +
  labs(title = "Test HDBSCAN on Train", x = "UMAP 1", y = "UMAP 2", color = "Cluster:")

p_test_rules_bg <- ggplot() +
  geom_point(data = df_plot_train, aes(x = UMAP1, y = UMAP2, color = Cluster), alpha = 0.2, size = 1.2) +
  geom_point(data = df_plot_test_rules, aes(x = UMAP1, y = UMAP2, color = Cluster), alpha = 0.9, size = 2.2, shape = 17) +
  scale_color_manual(values = cluster_colors, drop = FALSE) + 
  theme_minimal(base_size = 16) +
  theme(legend.position = "bottom", legend.key.width = unit(1, "cm"), plot.title = element_text(face = "bold")) +
  guides(color = guide_legend(nrow = 1, title.position = "left", title.vjust = 0.5, override.aes = list(size = 4, alpha = 1))) +
  labs(title = "Test Skope-Rules on Train", x = "UMAP 1", y = "UMAP 2", color = "Cluster:")

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 1)))
print(p_test_hdb_bg,   vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_test_rules_bg, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))


#test SHAP vs raw - solo HDBSCAN e Skope-Rules per summary -----
p_test_hdb <- ggplot(df_plot_test_hdb, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(alpha = 0.85, size = 3, shape = 16) +
  scale_color_manual(values = cluster_colors, drop = FALSE) +
  theme_minimal(base_size = 20) +
  labs(title = "", x = "UMAP 1 (SHAP)", y = "UMAP 2 (SHAP)", color = "Cluster:") +
  theme(
    plot.title = element_text(face = "bold", size = 22),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.position = "right",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))

p_test_rules <- ggplot(df_plot_test_rules, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(alpha = 0.85, size = 3, shape = 16) +
  scale_color_manual(values = cluster_colors, drop = FALSE) +
  theme_minimal(base_size = 20) +
  labs(title = "", x = "UMAP 1 (SHAP)", y = "UMAP 2 (SHAP)", color = "Cluster:") +
  theme(
    plot.title = element_text(face = "bold", size = 22),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.position = "right",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))


print(p_test_hdb)
print(p_test_rules)





#test R1 clusters distr plots -----
make_cat_stacked_plot <- function(test_cluster_vec, method_name) {
  if (length(cat_vars) == 0) return(NULL)
  
  df_tr <- X_full[, cat_vars, drop = FALSE]
  df_tr$Cluster <- cl_B_train$cluster
  df_tr$Set <- "Train"
  
  df_te <- X_test[, cat_vars, drop = FALSE]
  df_te$Cluster <- test_cluster_vec
  df_te$Set <- "Test"
  
  df_all <- rbind(df_tr, df_te)
  df_valid <- subset(df_all, Cluster > 0 & !is.na(Cluster))
  df_valid$Cluster <- factor(df_valid$Cluster)
  df_valid$Set <- factor(df_valid$Set, levels = c("Train", "Test"))
  
  res_list <- list()
  for (v in cat_vars) {
    tab <- as.data.frame(table(Cluster = df_valid$Cluster, Set = df_valid$Set, Category = df_valid[[v]]))
    tot <- tapply(tab$Freq, list(tab$Cluster, tab$Set), sum)
    tab$Tot <- apply(tab, 1, function(r) tot[as.character(r["Cluster"]), as.character(r["Set"])])
    tab$Perc <- ifelse(tab$Tot > 0, (tab$Freq / tab$Tot), 0)
    tab$Feature <- v
    res_list[[v]] <- tab
  }
  df_plot_cat <- do.call(rbind, res_list)
  
  p <- ggplot(df_plot_cat, aes(x = Set, y = Perc, fill = Category)) +
    geom_col(position = "fill", color = "white", linewidth = 0.5, width = 0.8) +
    facet_wrap(~ Cluster, nrow = 1, labeller = labeller(Cluster = function(x) paste("Cluster", x))) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_brewer(palette = "YlOrRd", name = "Risk factors 1:") +
    theme_minimal(base_size = 20) +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(1, "cm"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(title = NULL, x = NULL, y = "Proportion")
  
  return(p)
}

p_cat_stack_hdb <- make_cat_stacked_plot(cl_test_pred, "HDBSCAN")
p_cat_stack_rules <- make_cat_stacked_plot(cl_test_rules, "Skope-Rules")

grid.newpage()
print(p_cat_stack_hdb)
print(p_cat_stack_rules)



#test clinical boxplots with cl 0-----
make_profile_boxplot <- function(test_cluster_vec, method_name) {
  df_tr <- X_full
  df_tr$Cluster <- cl_B_train$cluster
  df_tr$Set <- "train"
  
  df_te <- X_test
  df_te$Cluster <- test_cluster_vec
  df_te$Set <- paste0("test (", method_name, ")")
  
  df_all <- rbind(df_tr, df_te)
  df_valid <- subset(df_all, Cluster >= 0 & !is.na(Cluster))
  df_valid$Cluster <- factor(ifelse(df_valid$Cluster == 0, "0", as.character(df_valid$Cluster)), levels = c("0", "1", "2", "3", "4", "5"))
  
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
    theme_minimal(base_size = 20) +
    theme(
      legend.position = "bottom", 
      strip.text = element_text(face = "bold"),
      legend.key.size = unit(1.5, "cm"),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 22, face = "bold")
    ) +
    labs(
      title = NULL,
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



#test clinical boxplots without cl 0-----
make_profile_boxplot <- function(test_cluster_vec, method_name) {
  df_tr <- X_full
  df_tr$Cluster <- cl_B_train$cluster
  df_tr$Set <- "train"
  
  df_te <- X_test
  df_te$Cluster <- test_cluster_vec
  df_te$Set <- paste0("test (", method_name, ")")
  
  df_all <- rbind(df_tr, df_te)
  
  df_valid <- subset(df_all, (Set == "train" | Cluster > 0) & !is.na(Cluster))
  df_valid$Cluster <- factor(df_valid$Cluster)
  
  df_m <- reshape(
    df_valid,
    direction = "long",
    varying   = cont_vars,
    v.names   = "Value",
    timevar   = "Feature",
    times     = cont_vars
  )
  
  df_m$Feature <- factor(df_m$Feature, 
                         levels = c("Age", "Volume_PTV", "LAA"),
                         labels = c("Age", "Volume PTV", "LAA"))
  
  col_test <- ifelse(method_name == "HDBSCAN", "#D55E00", "#009E73")
  fill_vals <- setNames(c("#0072B2", col_test), c("train", paste0("test (", method_name, ")")))
  
  p <- ggplot(df_m, aes(x = Cluster, y = Value, fill = Set)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.8, position = position_dodge(0.8)) +
    facet_wrap(~ Feature, scales = "free_y") +
    scale_fill_manual(values = fill_vals) +
    theme_minimal(base_size = 25) +
    theme(
      legend.position = "bottom", 
      strip.text = element_text(face = "bold", size = 25),
      axis.title.x = element_text(size = 25, margin = ggplot2::margin(t = 15)),
      axis.title.y = element_text(size = 25, margin = ggplot2::margin(r = 15)),
      axis.text = element_text(size = 20),
      legend.key.size = unit(2, "cm"),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 25, face = "bold"),
      legend.margin = ggplot2::margin(t = 20)
    ) +
    labs(
      title = NULL,
      x = "Cluster",
      y = "Clinical value",
      fill = "Cohort:"
    )
  return(p)
}

p_prof_hdb   <- make_profile_boxplot(cl_test_pred,  "HDBSCAN")
p_prof_rules <- make_profile_boxplot(cl_test_rules, "Skope-Rules")

print(p_prof_hdb)
print(p_prof_rules)


#test clinical boxplots zoomed-----
make_profile_boxplot <- function(test_cluster_vec, method_name) {
  df_tr <- X_full
  df_tr$Cluster <- cl_B_train$cluster
  df_tr$Set <- "train"
  
  df_te <- X_test
  df_te$Cluster <- test_cluster_vec
  df_te$Set <- paste0("test (", method_name, ")")
  
  df_all <- rbind(df_tr, df_te)
  
  df_valid <- subset(df_all, (Set == "train" | Cluster > 0) & !is.na(Cluster))
  df_valid$Cluster <- factor(df_valid$Cluster)
  
  df_m <- reshape(
    df_valid,
    direction = "long",
    varying   = cont_vars,
    v.names   = "Value",
    timevar   = "Feature",
    times     = cont_vars
  )
  
  df_m$Feature <- factor(df_m$Feature, 
                         levels = c("Age", "Volume_PTV", "LAA"),
                         labels = c("Age", "Volume PTV", "LAA"))
  
  df_m_plot <- subset(df_m, !(Feature == "LAA" & Value > 5))
  
  col_test <- ifelse(method_name == "HDBSCAN", "#D55E00", "#009E73")
  fill_vals <- setNames(c("#0072B2", col_test), c("train", paste0("test (", method_name, ")")))
  
  p <- ggplot(df_m_plot, aes(x = Cluster, y = Value, fill = Set)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.8, position = position_dodge(0.8)) +
    facet_wrap(~ Feature, scales = "free_y") +
    scale_fill_manual(values = fill_vals) +
    theme_minimal(base_size = 25) +
    theme(
      legend.position = "bottom", 
      strip.text = element_text(face = "bold", size = 25),
      axis.title.x = element_text(size = 25, margin = ggplot2::margin(t = 15)),
      axis.title.y = element_text(size = 25, margin = ggplot2::margin(r = 15)),
      axis.text = element_text(size = 20),
      legend.key.size = unit(2, "cm"),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 25, face = "bold"),
      legend.margin = ggplot2::margin(t = 20)
    ) +
    labs(
      title = NULL,
      x = "Cluster",
      y = "Clinical value",
      fill = "Cohort:"
    )
  return(p)
}

p_prof_hdb   <- make_profile_boxplot(cl_test_pred,  "HDBSCAN")
p_prof_rules <- make_profile_boxplot(cl_test_rules, "Skope-Rules")

print(p_prof_hdb)
print(p_prof_rules)

#test table n, os ----
os_summary_train <- aggregate(Y_full, by = list(Cluster = cl_B_train$cluster), FUN = function(x) sprintf("n=%d (%.1f%%)", length(x), mean(x)*100))
os_summary_hdb <- aggregate(Y_test[valid_outcomes_idx], by = list(Cluster = cl_test_pred[valid_outcomes_idx]), FUN = function(x) sprintf("n=%d (%.1f%%)", length(x), mean(x)*100))
os_summary_rules <- aggregate(Y_test[valid_outcomes_idx], by = list(Cluster = cl_test_rules[valid_outcomes_idx]), FUN = function(x) sprintf("n=%d (%.1f%%)", length(x), mean(x)*100))

comp_os <- merge(os_summary_train, os_summary_hdb, by = "Cluster", all = TRUE)
comp_os <- merge(comp_os, os_summary_rules, by = "Cluster", all = TRUE)
colnames(comp_os) <- c("cluster", "train_os2ys", "test_hdbscan_os2ys", "test_rules_os2ys")
print(comp_os)






#test clustering km ----
df_tr_hm <- X_full; df_tr_hm$Y_target <- Y_full
df_te_hm <- X_test[valid_outcomes_idx, , drop=FALSE]; df_te_hm$Y_target <- Y_test[valid_outcomes_idx]

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
  
  fit_km  <- survfit(Surv(os_24m, ev_24m) ~ cl, data = sd)
  lr_test <- survdiff(Surv(os_24m, ev_24m) ~ cl, data = sd)
  p_logrank <- 1 - pchisq(lr_test$chisq, length(lr_test$n) - 1)
  
  nomi_cluster <- sub(".*=", "", names(fit_km$strata)) 
  
  p_km <- ggsurvplot(
    fit_km,
    data = sd,
    conf.int = TRUE,
    conf.int.alpha = 0.15,
    xlim = c(0, 24),
    break.x.by = 3,
    ylim = c(0, 1),
    pval = sprintf("p = %.2e", p_logrank),
    risk.table = TRUE,
    legend = "top",                  
    legend.title = "Cluster",
    legend.labs = nomi_cluster,      
    font.legend = c(13, "bold"),     
    xlab = "Months from the beginning of radiotherapy",
    ylab = "Survival probability",
    title = t_text,
    ggtheme = theme_minimal()
  )
  
  p_km$plot <- p_km$plot + 
    theme(
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 13)
    )
  
  return(p_km)
}

p_km_hdb <- plot_km(df_sv, cl_test_pred[valid_outcomes_idx], "Kaplan-Meier Test set (HDBSCAN)")
p_km_ru <- plot_km(df_sv, cl_test_rules[valid_outcomes_idx], "Kaplan-Meier Test set (Skope-Rules)")

print(p_km_hdb)
print(p_km_ru)



plot_km <- function(sd, cl_v, t_text) {
  sd$cl <- cl_v
  sd <- sd[sd$cl > 0, ]
  
  fit_km <- survfit(Surv(os_24m, ev_24m) ~ cl, data = sd)
  
  nomi_cluster <- sub(".*=", "", names(fit_km$strata))
  
  p_km <- ggsurvplot(
    fit_km,
    data = sd,
    conf.int = TRUE,
    conf.int.alpha = 0.15,
    xlim = c(0, 24),
    break.x.by = 3,
    ylim = c(0, 1),
    risk.table = TRUE,
    risk.table.fontsize = 6,
    risk.table.height = 0.35,       
    legend = "top",
    legend.title = "Cluster",
    legend.labs = nomi_cluster,
    font.legend = c(20, "bold"),
    font.title = c(20, "bold"),
    font.x = c(20, "plain"),
    font.y = c(20, "plain"),
    font.tickslab = c(20, "plain"),
    xlab = "Months from the beginning of radiotherapy",
    ylab = "Survival probability",
    title = t_text,
    ggtheme = theme_minimal(base_size = 20),
    tables.theme = theme_minimal(base_size = 20)
  )
  
  p_km$plot <- p_km$plot +
    theme(
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 20)
    )
  
  p_km$table <- p_km$table +
    theme(
      text = element_text(size = 20),
      axis.text.x = element_text(size = 20),
      axis.title = element_text(size = 20),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    )
  
  return(p_km)
}

p_km_hdb <- plot_km(
  df_sv, cl_test_pred[valid_outcomes_idx],
  "Kaplan-Meier Test set (HDBSCAN)"
)

p_km_ru <- plot_km(
  df_sv, cl_test_rules[valid_outcomes_idx],
  "Kaplan-Meier Test set (Skope-Rules)"
)

print(p_km_hdb)
print(p_km_ru)




#test SHAP dependence plot -----
data_copy_test <- na.omit(X_test[, colnames(shap_test_agg), drop = FALSE])

test_plots <- lapply(seq_along(colnames(shap_test_agg)), function(i) {
  f <- colnames(shap_test_agg)[i]
  p <- shap_dependence(
    feature = f,
    shap_mat = shap_test_agg,
    X = X_test,
    color_by = "cluster",
    cluster = cl_test_pred,
    custom_labels = miei_titoli_shap
  )
  if (i < 4) {
    p <- p + theme(legend.position = "none")
  } else {
    p <- p + theme(legend.position = "right", legend.box = "vertical")
  }
  return(p)
})

print(wrap_plots(test_plots, ncol = 4))



#test mean abs SHAP plots ----
feature_labels <- c(
  "Age" = "Age",
  "Volume_PTV" = "Volume PTV",
  "LAA" = "LAA",
  "riskfactors1" = "Risk factors 1"
)

abs_shap_test <- abs(shap_test_agg)
test_class_labels <- ifelse(datatestos_aligned$OS2Y_1 == "Yes", "YES", "NO")
mean_abs_by_class_test <- aggregate(abs_shap_test, by = list(Class = test_class_labels), FUN = mean)
df_bar_class_test <- reshape(mean_abs_by_class_test, direction = "long", varying = colnames(abs_shap_test), v.names = "Mean_Abs_SHAP", timevar = "Feature", times = colnames(abs_shap_test))
tot_feature_test <- ave(df_bar_class_test$Mean_Abs_SHAP, df_bar_class_test$Feature, FUN = sum)
df_bar_class_test$Perc <- (df_bar_class_test$Mean_Abs_SHAP / tot_feature_test) * 100
df_bar_class_test$Label_Perc <- sprintf("%.1f%%", df_bar_class_test$Perc)
df_bar_class_test_custom <- df_bar_class_test
df_bar_class_test_custom$Feature <- ifelse(df_bar_class_test_custom$Feature %in% names(feature_labels),
                                           feature_labels[df_bar_class_test_custom$Feature],
                                           df_bar_class_test_custom$Feature)
p_bar_class_test <- ggplot(df_bar_class_test_custom, aes(x = reorder(Feature, Mean_Abs_SHAP, sum), y = Mean_Abs_SHAP, fill = Class)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "white", linewidth = 0.4) +
  geom_text(aes(label = Label_Perc), position = position_stack(vjust = 0.5), color = "white", fontface = "bold", size = 6) +
  coord_flip() + 
  scale_fill_manual(values = c("NO" = "#ff8080", "YES" = "#93c47d")) + 
  theme_minimal(base_size = 20) + 
  labs(
    x = "Variable", 
    y = "Mean(|SHAP value|)", 
    fill = "OS 2ys"
  ) +
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )

print(p_bar_class_test)


#test pred IP (OS2ys) by cluster ----
df_test_plot <- data.frame(
  Cluster = factor(cl_test_pred),
  Pred_Prob = p_test
)

p_prob_cluster_test <- ggplot(
  df_test_plot,
  aes(x = Cluster, y = Pred_Prob, fill = Cluster)
) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1.5) +
  theme_minimal(base_size = 20) +
  labs(
    x = "Cluster (HDBSCAN)",
    y = "Predicted probability for target (OS2ys)"
  ) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18, color = "black"),
    legend.position = "none"
  )

print(p_prob_cluster_test)

df_test_plot_rules <- data.frame(
  Cluster = factor(cl_test_rules),
  Pred_Prob = p_test
)

p_prob_cluster_rules_test <- ggplot(
  df_test_plot_rules,
  aes(x = Cluster, y = Pred_Prob, fill = Cluster)
) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1.5) +
  theme_minimal(base_size = 20) +
  labs(
    x = "Cluster (Skope-Rules)",
    y = "Predicted probability for target (OS2ys)"
  ) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18, color = "black"),
    legend.position = "none"
  )

print(p_prob_cluster_rules_test)





df_test_plot <- data.frame(
  Cluster = factor(cl_test_rules),
  Pred_Prob = p_test,
  OS2ys = Y_test
)

tapply(df_test_plot$Pred_Prob, df_test_plot$Cluster, mean)
prop.table(table(df_test_plot$Cluster, df_test_plot$OS2ys), margin = 1)

df_test_plot_pred <- data.frame(
  Cluster = factor(cl_test_pred),
  Pred_Prob = p_test,
  OS2ys = Y_test
)

tapply(df_test_plot_pred$Pred_Prob, df_test_plot_pred$Cluster, mean)
prop.table(table(df_test_plot_pred$Cluster, df_test_plot_pred$OS2ys), margin = 1)

ip_train <- predict(fit_xgb, newdata = X_mat)

df_train_plot <- data.frame(
  Cluster = factor(cl_B_train$cluster),
  Pred_Prob = ip_train,
  OS2ys = Y_full
)

tapply(df_train_plot$Pred_Prob, df_train_plot$Cluster, mean)
prop.table(table(df_train_plot$Cluster, df_train_plot$OS2ys), margin = 1)




#KM estimate of OS2ys ----
library(survival)
df_sub = data
df_sub$OS[df_sub$ID == "LM1001"] <- 24
df_sub$OS_24m <- pmin(df_sub$OS, 24)

df_sub$event_death_24m <- 1 - df_sub$OS2ys_num # 1 = evento, 0 = censurato
km_fit_overall <- survfit(Surv(OS_24m, event_death_24m) ~ 1, data = df_sub)
plot(km_fit_overall, 
     xlab = "Mesi", 
     ylab = "Probabilità di Sopravvivenza",
     main = "Curva di Kaplan-Meier Globale (OS a 24 mesi)")
km_summary_24m <- summary(km_fit_overall, times = 24)
data.frame(
  Time = km_summary_24m$time,
  n.risk = km_summary_24m$n.risk,
  n.event = km_summary_24m$n.event,
  Survival_Estimate = km_summary_24m$surv,
  Std_Error = km_summary_24m$std.err,
  Lower_CI = km_summary_24m$lower,
  Upper_CI = km_summary_24m$upper
)

library(survival)
df_sub = datatestos
df_sub$OS_24m <- pmin(df_sub$OS, 24)
df_sub$event_death_24m <-df_sub$OS2Y # 1 = evento, 0 = censurato
km_fit_overall <- survfit(Surv(OS_24m, event_death_24m) ~ 1, data = df_sub)
plot(km_fit_overall, 
     xlab = "Mesi", 
     ylab = "Probabilità di Sopravvivenza",
     main = "Curva di Kaplan-Meier Globale (OS a 24 mesi)")
km_summary_24m <- summary(km_fit_overall, times = 24)
data.frame(
  Time = km_summary_24m$time,
  n.risk = km_summary_24m$n.risk,
  n.event = km_summary_24m$n.event,
  Survival_Estimate = km_summary_24m$surv,
  Std_Error = km_summary_24m$std.err,
  Lower_CI = km_summary_24m$lower,
  Upper_CI = km_summary_24m$upper
)


#test: calibration who are patients with where the model makes more mistakes ----
flat_patients_idx <- which(p_test < 0.5)

df_flat_investigation <- data.frame(
  ID = df_test$ID[flat_patients_idx],
  Predicted_Prob = p_test[flat_patients_idx],
  Observed_OS2ys = Y_test[flat_patients_idx],
  Cluster_HDBSCAN = cl_test_pred[flat_patients_idx],
  Cluster_Rules = cl_test_rules[flat_patients_idx],
  Age = df_test$Age[flat_patients_idx],
  Volume_PTV = df_test$Volume_PTV[flat_patients_idx],
  LAA = df_test$LAA[flat_patients_idx]
)

print(head(df_flat_investigation, 20))
dim(df_flat_investigation) 
table(df_flat_investigation$Cluster_Rules)
table(cl_test_rules)
table(cl_test_pred)
table(cl_B_train$cluster)

summary_flat_cluster <- aggregate(
  cbind(Age, Volume_PTV, LAA, Observed_OS2ys, Predicted_Prob) ~ Cluster_HDBSCAN, 
  data = df_flat_investigation, 
  FUN = function(x) c(Mean = mean(x, na.rm = TRUE), Median = median(x, na.rm = TRUE))
)
print(summary_flat_cluster)

table_cluster_flat <- table(df_flat_investigation$Cluster_HDBSCAN, df_flat_investigation$Observed_OS2ys)
print(table_cluster_flat)

#test: median in clusters -----
tapply(df_train_plot$Volume_PTV, cl_B_train$cluster, median, na.rm = TRUE)
tapply(df_test_plot$Volume_PTV, cl_test_pred, median, na.rm = TRUE)
tapply(df_test_plot$Volume_PTV, cl_test_rules, median, na.rm = TRUE)

tapply(df_train_plot$LAA, cl_B_train$cluster, median, na.rm = TRUE)
tapply(df_test_plot$LAA, cl_test_pred, median, na.rm = TRUE)
tapply(df_test_plot$LAA, cl_test_rules, median, na.rm = TRUE)



#test: test age, volume PTV, LAA, R1 difference between clusters test ----
library(dplyr)
df_train_plot <- data.frame(
  Age = as.numeric(data$Age),
  Volume_PTV = as.numeric(data$Volume_PTV),
  LAA = as.numeric(data$LAA),
  Risk_Factor_1 = as.factor(data$riskfactors1)
)

patient = df_train_plot[,c("Age", "Volume_PTV", "LAA", "Risk_Factor_1")]
rows= complete.cases(patient) 
df_train_plot = df_train_plot[rows,]


df_test_plot <- data.frame(
  Age = as.numeric(datatest$Age),
  Volume_PTV = as.numeric(datatest$VOLUME_PTV),
  LAA = as.numeric(datatest$`LAA%`),
  Risk_Factor_1 = as.factor(datatest$risk_factor1)
)

patient = df_train_plot[,c("Age", "Volume_PTV", "LAA", "Risk_Factor_1")]
rows= complete.cases(patient) 
datatrain = df_train_plot[rows,]
table(cl_test_rules, datatest$risk_factor1)
table(cl_test_pred, datatest$risk_factor1)
table(cl_B_train$cluster, datatrain$Risk_Factor_1)
prop.table(table(cl_test_rules, datatest$risk_factor1), margin = 1) * 100
prop.table(table(cl_test_pred, datatest$risk_factor1), margin = 1) * 100
prop.table(table(cl_B_train$cluster, datatrain$Risk_Factor_1), margin = 1) * 100

dim(datatrain)
colnames(datatrain)
datatrain$Cohort <- "Train"
df_test_plot$Cohort <- "Test"
df_comparison_features <- rbind(datatrain, df_test_plot)
df_comparison_features$Cohort <- factor(df_comparison_features$Cohort, levels = c("Train", "Test"))



#R1
datatrain$Cluster_HDBSCAN <- cl_B_train$cluster
df_test_plot$Cluster_HDBSCAN <- cl_test_pred
df_train_plot$Cluster_Rules <- cl_B_train$cluster
df_test_plot$Cluster_Rules <- cl_test_rules
datatrain$Cluster_Rules <- cl_B_train$cluster
df_test_plot$Cluster_Rules <- cl_test_rules
df_comparison_features <- rbind(datatrain, df_test_plot)
df_comparison_features$Cohort <- factor(df_comparison_features$Cohort, levels = c("Train", "Test"))

#standard chi squared
sub_c3_hdb <- subset(df_comparison_features, Cluster_HDBSCAN == 3)
tab_c3_hdb <- table(sub_c3_hdb$Cohort, sub_c3_hdb$Risk_Factor_1)
res_c3_hdb <- chisq.test(tab_c3_hdb)


sub_c3_rules <- subset(df_comparison_features, Cluster_Rules == 3)
tab_c3_rules <- table(sub_c3_rules$Cohort, sub_c3_rules$Risk_Factor_1)
res_c3_rules <- chisq.test(tab_c3_rules)




#age, volume PTV, LAA - pairwise wilcoxon Bonf adjusted
vars <- c("Age", "Volume_PTV", "LAA")
alpha <- 0.05
k <- 10
alpha_bonferroni <- alpha / k



for (v in vars) {
  cat(sprintf("\nVariable: %s\n", v))
  for (c in 1:5) {
    sub_c <- subset(df_comparison_features, Cluster_HDBSCAN == c)
    dt_c <- sub_c$Cohort == "Train"
    dt_test <- sub_c$Cohort == "Test"
    if(sum(dt_c) > 1 && sum(dt_test) > 1) {
      form <- as.formula(paste(v, "~ Cohort"))
      res <- wilcox.test(form, data = sub_c)
      cat(sprintf("  Cluster %d (p-value: %e, Bonferroni alpha < %f): %s\n", 
                  c, res$p.value, alpha_bonferroni, ifelse(res$p.value < alpha_bonferroni, "Significant", "Non significant")))
    }
  }
}


for (v in vars) {
  cat(sprintf("\nVariable: %s\n", v))
  for (c in 1:5) {
    sub_c <- subset(df_comparison_features, Cluster_Rules == c)
    dt_c <- sub_c$Cohort == "Train"
    dt_test <- sub_c$Cohort == "Test"
    if(sum(dt_c) > 1 && sum(dt_test) > 1) {
      form <- as.formula(paste(v, "~ Cohort"))
      res <- wilcox.test(form, data = sub_c)
      cat(sprintf("  Cluster %d (p-value: %e, Bonferroni alpha < %f): %s\n", 
                  c, res$p.value, alpha_bonferroni, ifelse(res$p.value < alpha_bonferroni, "Significant", "Non significant")))
    }
  }
}


wilcox.test(LAA ~ Cohort, data = subset(df_comparison_features, Cluster_HDBSCAN == 5))
wilcox.test(LAA ~ Cohort, data = subset(df_comparison_features, Cluster_Rules == 5))



#test: HDBSCAN vs SR contigency table ---
contingency_table <- table(HDBSCAN = cl_test_pred, SkopeRules = cl_test_rules)
print(contingency_table)


#test: kruskal tests to assess shifts ----
patient = data[,c("Age", "Volume_PTV", "LAA", "riskfactors1")]
rows= complete.cases(patient) 
prova = data[rows,]

df_test_plot <- data.frame(
  Age = as.numeric(datatest$Age),
  Volume_PTV = as.numeric(datatest$VOLUME_PTV),
  LAA = as.numeric(datatest$`LAA%`),
  Risk_Factor_1 = as.factor(datatest$risk_factor1),
  Cluster = factor(cl_test_pred)
)
df_test_plot2 <- data.frame(
  Age = as.numeric(datatest$Age),
  Volume_PTV = as.numeric(datatest$VOLUME_PTV),
  LAA = as.numeric(datatest$`LAA%`),
  Risk_Factor_1 = as.factor(datatest$risk_factor1),
  Cluster = factor(cl_test_rules)
)
df_test_plot <- df_test_plot[complete.cases(df_test_plot) & df_test_plot$Cluster != 0, ]
df_test_plot2 <- df_test_plot2[complete.cases(df_test_plot2) & df_test_plot2$Cluster != 0, ]

kruskal.test(Age ~ Cluster, data = df_test_plot)
kruskal.test(Age ~ Cluster, data = df_test_plot2)
kruskal.test(Volume_PTV ~ Cluster, data = df_test_plot)
kruskal.test(Volume_PTV ~ Cluster, data = df_test_plot2)
kruskal.test(LAA ~ Cluster, data = df_test_plot)
kruskal.test(LAA ~ Cluster, data = df_test_plot2)


tab_risk_test <- table(df_test_plot$Risk_Factor_1, df_test_plot$Cluster)
tab_clean <- tab_risk_test[, colSums(tab_risk_test) > 0]
tab_clean <- tab_clean[rowSums(tab_clean) > 0, ]
tab_risk_test2 <- table(df_test_plot2$Risk_Factor_1, df_test_plot2$Cluster)
tab_clean2 <- tab_risk_test2[, colSums(tab_risk_test2) > 0]
tab_clean2 <- tab_clean2[rowSums(tab_clean2) > 0, ]
# chi_risk_test <- chisq.test(tab_clean)
# print(tab_clean)
# print(chi_risk_test)
# print(round(chi_risk_test$residuals, 2))
tab_clean
tab_clean2
chi_sim <- chisq.test(tab_clean, simulate.p.value = TRUE, B = 10000)
print(chi_sim)
print(round(chi_sim$residuals, 2))
chi_sim2 <- chisq.test(tab_clean2, simulate.p.value = TRUE, B = 10000)
print(chi_sim2)
print(round(chi_sim2$residuals, 2))


#spearnman tests 
cor.test(as.numeric(as.character(df_test_plot$Cluster)), as.numeric(as.character(df_test_plot$Risk_Factor_1)), method = "spearman", alternative = "two.sided", exact = FALSE)
cor.test(as.numeric(as.character(df_test_plot2$Cluster)), as.numeric(as.character(df_test_plot2$Risk_Factor_1)), method = "spearman", alternative = "two.sided", exact = FALSE)
#no monotonic association 

#test: cal investigation ----

p_val <- p_test[valid_outcomes_idx]
y_val <- Y_test[valid_outcomes_idx]

library(ggplot2)

df_cal <- data.frame(
  p_val = p_val,
  y_val = as.numeric(as.character(y_val))
)

ggplot(df_cal, aes(x = p_val, y = y_val)) +
  geom_jitter(aes(color = factor(y_val)), height = 0.05, alpha = 0.5, size = 2) +
  
  geom_smooth(method = "loess", color = "black", se = FALSE, linetype = "dashed") +
  scale_y_continuous(breaks = c(0, 1), labels = c("0", "1")) +
  scale_color_manual(values = c("0" = "#E69F00", "1" = "#0072B2"), name = "Esito Reale") +
  theme_minimal() +
  labs(
    title = "Distribuzione 0 vs 1 lungo le Probabilità Predette",
    x = "Probabilità Predetta (p_val)",
    y = "Esito Reale (Y)"
  )


#test: clusters os2ys rate - overall os2ysrate ----
ip_train <- predict(fit_xgb, newdata = X_mat)
df_train_plot <- data.frame(
  Cluster = factor(cl_B_train$cluster),
  OS2ys = Y_full,
  prob = ip_train
)

df_test_plot <- data.frame(
  Cluster = factor(cl_test_pred),
  OS2ys = Y_test,
  prob = p_val
)
df_test_plot2 <- data.frame(
  Cluster = factor(cl_test_rules),
  OS2ys = Y_test,
  prob = p_val
)


mean(Y_full)
mean(Y_test)
prop.table(table(df_train_plot$Cluster, df_train_plot$OS2ys), margin = 1) * 100
prop.table(table(df_test_plot$Cluster, df_test_plot$OS2ys), margin = 1) * 100
prop.table(table(df_test_plot2$Cluster, df_test_plot2$OS2ys), margin = 1) * 100

tapply(df_train_plot$prob, df_train_plot$Cluster, mean, na.rm = TRUE) * 100
tapply(df_test_plot$prob, df_test_plot$Cluster, mean, na.rm = TRUE) * 100
tapply(df_test_plot2$prob, df_test_plot2$Cluster, mean, na.rm = TRUE) * 100

#test: R1 investigation -----
p_train <- predict(fit_xgb, newdata = X_mat)
p_val <- p_test[valid_outcomes_idx]

patient = data[,c("Age", "Volume_PTV", "LAA", "riskfactors1")]
rows= complete.cases(patient) 
prova = data[rows,]

df_train_plot <- data.frame(
  cluster = factor(cl_B_train$cluster),
  risk1 = prova$riskfactors1,
  prob = p_train
)

df_test_plot <- data.frame(
  cluster = factor(cl_test_pred),
  risk1 = datatest$risk_factor1,
  prob = p_val
)

df_test_plot2 <- data.frame(
  cluster = factor(cl_test_rules),
  risk1 = datatest$risk_factor1,
  prob = p_val
)


tapply(df_train_plot$prob, df_train_plot$risk1, median, na.rm = TRUE)
tapply(df_test_plot$prob, df_test_plot$risk1, median, na.rm = TRUE)
tapply(df_test_plot2$prob, df_test_plot2$risk1, median, na.rm = TRUE)


tapply(df_train_plot$prob, list(df_train_plot$cluster, df_train_plot$risk1), function(x) median(x, na.rm = TRUE))
tapply(df_test_plot$prob, list(df_test_plot$cluster, df_test_plot$risk1), function(x) median(x, na.rm = TRUE))
tapply(df_test_plot2$prob, list(df_test_plot2$cluster, df_test_plot2$risk1), function(x) median(x, na.rm = TRUE))


table(df_train_plot$cluster, df_train_plot$risk1)
table(df_test_plot$cluster, df_test_plot$risk1)
table(df_test_plot2$cluster, df_test_plot$risk1)





