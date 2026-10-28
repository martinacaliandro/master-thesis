#legenda:
#shap_raw non aggregati --> cl A
#shap_raw non aggregati --> UMAP --> cl C
#shap_agg o shap_4vars (questa non utilizzata) --> UMAP --> cl B
#shap_oof_7cols OOF non aggregati
#shap_oof_agg --> UMAP --> cl_oof
#nota: per plot con UMAP dati raw, non aggrego dummies (non sono SHAP)
#SHAP su TR set non spiegano rumore: confronto con OOF e modalità di validazione adottate a monte
#OOF: fluttuazioni di SHAP annullate, rimane solo il segnale forte e robusto (Volume), riflettendone la OS

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

data <- read.csv('dataclean_final_withNA.csv', header=T)

# preprocessing -----
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

# patient = data[,c("ID", "Age", "Volume_PTV", "LAA", "riskfactors1")]
# patient= na.omit(patient) 

# 1) hyperparameters -----
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

X_full_sc <- X_full
X_full_sc[, cont_vars] <- scale(X_full[, cont_vars])

cat("Patients:", nrow(X_full_sc), "Events OS2ys=1:", sum(Y_full), "\n")


# 2) full XGBoost fit + shap computation -----
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


# 3) UMAP + HDBSCAN -----
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

set.seed(seed_umap)
um2_oof <- uwot::umap(shap_oof_agg, n_neighbors = nn_umap, min_dist = 0, n_components = 2, n_threads = 1, n_sgd_threads = 1)
cl_oof <- dbscan::hdbscan(um2_oof, minPts = minPts_hdb)
print(cl_oof)

X_raw_mm <- model.matrix(~ ., data = X_full_sc)[, -1, drop = FALSE]
set.seed(seed_umap)
um2_raw <- uwot::umap(X_raw_mm, n_neighbors = nn_umap, min_dist = 0.1, n_components = 2, n_threads = 1, n_sgd_threads = 1)




# 4) clusters charcaterization -----
Cluster <- factor(cl_B$cluster)

n_OS <- cbind(n = table(Cluster), OS_rate = round(tapply(Y_full, Cluster, mean), 3))
print(n_OS)

print(round(aggregate(as.data.frame(shap_agg), list(Cluster = Cluster), mean)[, -1], 3))

print(round(aggregate(as.data.frame(shap_raw), list(Cluster = Cluster), mean)[, -1], 3))

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

#chi-squared test
print(table(Cluster, Z_full))
res_chi <- chisq.test(table(Cluster, Z_full), simulate.p.value = TRUE, B = 10000)
round(res_chi$stdres, 2)
print(res_chi)

idx_valid_cl <- cl_B$cluster > 0
d_test <- data.frame(Y = Y_full[idx_valid_cl], cl = droplevels(factor(cl_B$cluster[idx_valid_cl])), Z = droplevels(Z_full[idx_valid_cl]))

#glm
m_std0 <- glm(Y ~ 1, family = binomial(), data = d_test)
m_std1 <- glm(Y ~ cl, family = binomial(), data = d_test)
print(anova(m_std0, m_std1, test = "Chisq"))
print(summary(m_std1)$coefficients, 4)

#glmm
m_me0 <- glmer(Y ~ 1 + (1 | Z), family = binomial(), data = d_test)
m_me1 <- suppressWarnings(glmer(Y ~ cl + (1 | Z), family = binomial(), data = d_test))
print(anova(m_me0, m_me1))
print(summary(m_me1)$coefficients, 4)


# 5) plots -----

## 5.1) shap vs raw ------
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
  geom_point(alpha = 0.85, size = 2.2) + theme_minimal(base_size = 12) +
  labs(title = "UMAP 2D space (SHAP)", x = "UMAP 1 (SHAP)", y = "UMAP 2 (SHAP)")

p_raw_2d <- ggplot(df_plot, aes(x = RAW_UMAP1, y = RAW_UMAP2, color = Cluster)) +
  geom_point(alpha = 0.85, size = 2.2) + theme_minimal(base_size = 12) +
  labs(title = "B) UMAP 2D space (raw data)", x = "UMAP 1 (raw data)", y = "UMAP 2 (raw data)")

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(p_shap_2d, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_raw_2d,  vp = viewport(layout.pos.row = 1, layout.pos.col = 2))

## 5.2) shap vs raw colored by factors -----
roc_obj <- pROC::roc(Y_full, p_final, quiet = TRUE)
thresh_youden <- as.numeric(pROC::coords(roc_obj, "best", ret = "threshold", best.method = "youden")[1])
df_plot$Pred_Class <- factor(ifelse(p_final >= thresh_youden, "YES", "NO"), levels = c("NO", "YES"))

var_a <- covariates[1]
var_b <- ifelse(length(covariates) > 1, covariates[2], covariates[1])

p_raw_cl <- ggplot(df_plot, aes(x = RAW_UMAP1, y = RAW_UMAP2, color = Cluster)) + geom_point(alpha=0.8, size=1.8) + theme_minimal(base_size=11) + labs(title="1. Raw UMAP: Cluster")
p_raw_os <- ggplot(df_plot, aes(x = RAW_UMAP1, y = RAW_UMAP2, color = OS2ys)) + geom_point(alpha=0.8, size=1.8) + scale_color_manual(values=c("NO"="#D55E00", "YES"="#0072B2")) + theme_minimal(base_size=11) + labs(title="2. Raw UMAP: OS2ys")
p_raw_prob <- ggplot(df_plot, aes(x = RAW_UMAP1, y = RAW_UMAP2, color = Pred_Prob)) + geom_point(alpha=0.85, size=1.8) + scale_color_gradient2(low="#D55E00", mid="#F0E442", high="#0072B2", midpoint=0.5) + theme_minimal(base_size=11) + labs(title="3. Raw UMAP: Pred prob")
p_raw_pred <- ggplot(df_plot, aes(x = RAW_UMAP1, y = RAW_UMAP2, color = Pred_Class)) + geom_point(alpha=0.8, size=1.8) + scale_color_manual(values=c("NO"="#D55E00", "YES"="#0072B2")) + theme_minimal(base_size=11) + labs(title=paste0("4. Raw UMAP: Pred class (",round(thresh_youden,2),")"))
p_raw_va <- ggplot(df_plot, aes(x = RAW_UMAP1, y = RAW_UMAP2, color = .data[[var_a]])) + geom_point(alpha=0.8, size=1.8) + theme_minimal(base_size=11) + labs(title=paste("5. Raw UMAP:", var_a)); if(is.numeric(df_plot[[var_a]])) p_raw_va <- p_raw_va + scale_color_viridis_c()
p_raw_vb <- ggplot(df_plot, aes(x = RAW_UMAP1, y = RAW_UMAP2, color = .data[[var_b]])) + geom_point(alpha=0.8, size=1.8) + theme_minimal(base_size=11) + labs(title=paste("6. Raw UMAP:", var_b)); if(is.numeric(df_plot[[var_b]])) p_raw_vb <- p_raw_vb + scale_color_viridis_c()
p_raw_OSm <- ggplot(df_plot, aes(x = RAW_UMAP1, y = RAW_UMAP2, color = OS_months)) + geom_point(alpha=0.8, size=1.8) + scale_color_viridis_c(option="viridis") + theme_minimal(base_size=11) + labs(title="7. Raw UMAP: OS (months)")

p_shap_cl <- ggplot(df_plot, aes(x = SHAP_UMAP1, y = SHAP_UMAP2, color = Cluster)) + geom_point(alpha=0.8, size=1.8) + theme_minimal(base_size=11) + labs(title="1. SHAP UMAP: Cluster")
p_shap_os <- ggplot(df_plot, aes(x = SHAP_UMAP1, y = SHAP_UMAP2, color = OS2ys)) + geom_point(alpha=0.8, size=1.8) + scale_color_manual(values=c("NO"="#D55E00", "YES"="#0072B2")) + theme_minimal(base_size=11) + labs(title="2. SHAP UMAP: OS2ys")
p_shap_prob <- ggplot(df_plot, aes(x = SHAP_UMAP1, y = SHAP_UMAP2, color = Pred_Prob)) + geom_point(alpha=0.85, size=1.8) + scale_color_gradient2(low="#D55E00", mid="#F0E442", high="#0072B2", midpoint=0.5) + theme_minimal(base_size=11) + labs(title="3. SHAP UMAP: Pred prob")
p_shap_pred <- ggplot(df_plot, aes(x = SHAP_UMAP1, y = SHAP_UMAP2, color = Pred_Class)) + geom_point(alpha=0.8, size=1.8) + scale_color_manual(values=c("NO"="#D55E00", "YES"="#0072B2")) + theme_minimal(base_size=11) + labs(title=paste0("4. SHAP UMAP: Pred class (",round(thresh_youden,2),")"))
p_shap_va <- ggplot(df_plot, aes(x = SHAP_UMAP1, y = SHAP_UMAP2, color = .data[[var_a]])) + geom_point(alpha=0.8, size=1.8) + theme_minimal(base_size=11) + labs(title=paste("5. SHAP UMAP:", var_a)); if(is.numeric(df_plot[[var_a]])) p_shap_va <- p_shap_va + scale_color_viridis_c()
p_shap_vb <- ggplot(df_plot, aes(x = SHAP_UMAP1, y = SHAP_UMAP2, color = .data[[var_b]])) + geom_point(alpha=0.8, size=1.8) + theme_minimal(base_size=11) + labs(title=paste("6. SHAP UMAP:", var_b)); if(is.numeric(df_plot[[var_b]])) p_shap_vb <- p_shap_vb + scale_color_viridis_c()
p_shap_OSm <- ggplot(df_plot, aes(x = SHAP_UMAP1, y = SHAP_UMAP2, color = OS_months)) + geom_point(alpha=0.8, size=1.8) + scale_color_viridis_c(option="viridis") + theme_minimal(base_size=11) + labs(title="7. SHAP UMAP: OS (months)")

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 7)))
print(p_raw_cl,   vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_raw_os,   vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p_raw_prob, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
print(p_raw_pred, vp = viewport(layout.pos.row = 1, layout.pos.col = 4))
print(p_raw_va,   vp = viewport(layout.pos.row = 1, layout.pos.col = 5))
print(p_raw_vb,   vp = viewport(layout.pos.row = 1, layout.pos.col = 6))
print(p_raw_OSm,  vp = viewport(layout.pos.row = 1, layout.pos.col = 7))
print(p_shap_cl,  vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(p_shap_os,  vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(p_shap_prob,vp = viewport(layout.pos.row = 2, layout.pos.col = 3))
print(p_shap_pred,vp = viewport(layout.pos.row = 2, layout.pos.col = 4))
print(p_shap_va,  vp = viewport(layout.pos.row = 2, layout.pos.col = 5))
print(p_shap_vb,  vp = viewport(layout.pos.row = 2, layout.pos.col = 6))
print(p_shap_OSm, vp = viewport(layout.pos.row = 2, layout.pos.col = 7))

pushViewport(viewport(layout = grid.layout(2, 5)))
print(p_raw_os,   vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_raw_prob, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p_raw_pred, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
print(p_raw_va,   vp = viewport(layout.pos.row = 1, layout.pos.col = 4))
print(p_raw_OSm,  vp = viewport(layout.pos.row = 1, layout.pos.col = 5))
print(p_shap_os,  vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(p_shap_prob,vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(p_shap_pred,vp = viewport(layout.pos.row = 2, layout.pos.col = 3))
print(p_shap_va,  vp = viewport(layout.pos.row = 2, layout.pos.col = 4))
print(p_shap_OSm, vp = viewport(layout.pos.row = 2, layout.pos.col = 5))

## 5.3) mean abs shap ------
abs_shap <- abs(shap_agg)
mean_abs_by_class <- aggregate(abs_shap, by = list(Class = df_plot$OS2ys), FUN = mean)
df_bar_class <- reshape(mean_abs_by_class, direction = "long", varying = colnames(abs_shap), v.names = "Mean_Abs_SHAP", timevar = "Feature", times = colnames(abs_shap))
tot_feature <- ave(df_bar_class$Mean_Abs_SHAP, df_bar_class$Feature, FUN = sum)
df_bar_class$Perc <- (df_bar_class$Mean_Abs_SHAP / tot_feature) * 100
df_bar_class$Label_Perc <- sprintf("%.1f%%", df_bar_class$Perc)

p_bar_class <- ggplot(df_bar_class, aes(x = reorder(Feature, Mean_Abs_SHAP, sum), y = Mean_Abs_SHAP, fill = Class)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "white", linewidth = 0.4) +
  geom_text(aes(label = Label_Perc), position = position_stack(vjust = 0.5), color = "white", fontface = "bold", size = 3.8) +
  coord_flip() + scale_fill_manual(values = c("NO" = "#D55E00", "YES" = "#0072B2")) +
  theme_minimal(base_size = 12) + labs(title = "Mean |SHAP| per Target class OS2ys", x = "Variable", y = "Mean(|SHAP value|)", fill = "OS 2ys")
print(p_bar_class)

mean_abs_by_cl <- aggregate(abs_shap, by = list(Cluster = Cluster), FUN = mean)
df_bar_cl <- reshape(mean_abs_by_cl, direction = "long", varying = colnames(abs_shap), v.names = "Mean_Abs_SHAP", timevar = "Feature", times = colnames(abs_shap))
tot_feature_cl <- ave(df_bar_cl$Mean_Abs_SHAP, df_bar_cl$Feature, FUN = sum)
df_bar_cl$Perc <- (df_bar_cl$Mean_Abs_SHAP / tot_feature_cl) * 100
df_bar_cl$Label_Perc <- ifelse(df_bar_cl$Perc >= 4, sprintf("%.1f%%", df_bar_cl$Perc), "")

p_bar_cl <- ggplot(df_bar_cl, aes(x = reorder(Feature, Mean_Abs_SHAP, sum), y = Mean_Abs_SHAP, fill = Cluster)) +
  geom_bar(stat = "identity", position = "stack", width = 0.65, color = "white", linewidth = 0.4) +
  geom_text(aes(label = Label_Perc), position = position_stack(vjust = 0.5), color = "white", fontface = "bold", size = 3.5) +
  coord_flip() + theme_minimal(base_size = 12) + labs(title = "Mean |SHAP| per Cluster", x = "Variable", y = "Mean(|SHAP value|)", fill = "Cluster")
print(p_bar_cl)

## 5.4) waterfall plots ------
shap_df_all <- as.data.frame(shap_agg)
mean_signed_by_cl <- aggregate(shap_df_all, by = list(Cluster = Cluster), FUN = mean)
df_profile <- reshape(mean_signed_by_cl, direction = "long", varying = colnames(shap_df_all), v.names = "Mean_SHAP", timevar = "Feature", times = colnames(shap_df_all))
df_profile$Sign <- ifelse(df_profile$Mean_SHAP >= 0, "Aumenta OS (+)", "Riduce OS (-)")
df_profile$Cluster_Lab <- paste("Cluster", df_profile$Cluster)

p_zero_bars <- ggplot(df_profile, aes(x = reorder(Feature, abs(Mean_SHAP), mean), y = Mean_SHAP, fill = Sign)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.4) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.7) +
  geom_text(aes(y = Mean_SHAP + ifelse(Mean_SHAP >= 0, 0.06, -0.06), label = sprintf("%+.2f", Mean_SHAP)), fontface = "bold", size = 3.4) +
  coord_flip() + facet_wrap(~ Cluster_Lab, nrow = 1) +
  scale_fill_manual(values = c("Aumenta OS (+)" = "#0072B2", "Riduce OS (-)" = "#D55E00")) +
  scale_y_continuous(expand = expansion(mult = c(0.22, 0.22))) + theme_minimal(base_size = 12) +
  theme(legend.position = "top", strip.text = element_text(face = "bold", size = 11), panel.grid.minor = element_blank()) +
  labs(title = "SHAP path per cluster", x = "Variable", y = "SHAP Contribution", fill = "Direction:")
print(p_zero_bars)

logit_p   <- log(p_final / (1 - p_final))
base_val  <- mean(logit_p - rowSums(shap_df_all))
prob_base <- round(plogis(base_val) * 100, 1)

make_waterfall_cl <- function(cl_id) {
  idx_cl <- which(Cluster == cl_id)
  cl_shap_means_raw <- colMeans(shap_df_all[idx_cl, , drop = FALSE])
  cl_shap_sorted <- cl_shap_means_raw[order(abs(cl_shap_means_raw), decreasing = TRUE)]
  feat_names <- names(cl_shap_sorted)
  
  feat_labels <- sapply(feat_names, function(f) {
    val_cl <- df_sub[[f]][idx_cl]
    if (is.factor(val_cl) || is.character(val_cl)) {
      tab <- table(val_cl)
      if(length(tab) == 0) return(f)
      paste0(f, "\n[Mode: ", names(which.max(tab)), " (", round(max(tab)/sum(tab)*100,0), "%)]")
    } else {
      paste0(f, "\n[Mean: ", round(mean(val_cl, na.rm = TRUE), 1), "]")
    }
  })
  
  wf_df <- data.frame(Pos = rev(seq_along(cl_shap_sorted)), Feature = feat_names, Feature_Lab = feat_labels, Contribution = as.numeric(cl_shap_sorted), stringsAsFactors = FALSE)
  wf_df$CumSum <- base_val + cumsum(wf_df$Contribution)
  wf_df$Start  <- c(base_val, head(wf_df$CumSum, -1))
  wf_df$Sign   <- ifelse(wf_df$Contribution >= 0, "Positive (OS increases)", "Negative (reduces OS)")
  
  conn_df <- data.frame(x = head(wf_df$Pos, -1) - 0.28, xend = tail(wf_df$Pos, -1) + 0.28, y = head(wf_df$CumSum, -1), yend = head(wf_df$CumSum, -1))
  final_val  <- tail(wf_df$CumSum, 1)
  prob_final <- round(plogis(final_val) * 100, 1)
  col_line   <- ifelse(final_val >= base_val, "#0072B2", "#D55E00")
  
  ggplot(wf_df) +
    geom_hline(yintercept = base_val, linetype = "dashed", color = "grey40", linewidth = 0.7) +
    geom_hline(yintercept = final_val, linetype = "dotted", color = col_line, linewidth = 0.8) +
    geom_segment(data = conn_df, aes(x = x, xend = xend, y = y, yend = yend), linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_rect(aes(xmin = Pos - 0.28, xmax = Pos + 0.28, ymin = pmin(Start, CumSum), ymax = pmax(Start, CumSum), fill = Sign), color = "white", linewidth = 0.3) +
    geom_text(aes(x = Pos, y = (Start + CumSum) / 2, label = sprintf("%+.2f", Contribution)), color = "white", fontface = "bold", size = 3.1) +
    scale_x_continuous(breaks = wf_df$Pos, labels = wf_df$Feature_Lab) +
    scale_fill_manual(values = c("Positive (OS increases)" = "#0072B2", "Negative (reduces OS)" = "#D55E00"), drop = FALSE) +
    coord_flip() + theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), legend.position = "top", axis.text.y = element_text(face = "bold", size = 8.5)) +
    labs(title = paste0("SHAP Waterfall - Cluster ", cl_id, " (n=", length(idx_cl), ")"),
         subtitle = paste0("Log-odds: from ", round(base_val, 2), " to ", round(final_val, 2), " | OS2ys: form ", prob_base, "% to ", prob_final, "%"),
         x = NULL, y = "Cumulative SHAP contribution", fill = "Direction:")
}

valid_cl_levels <- setdiff(levels(Cluster), "0")
wf_list <- lapply(valid_cl_levels, make_waterfall_cl)

n_wf   <- length(wf_list)
ncol_g <- ifelse(n_wf <= 4, 2, 3)
nrow_g <- ceiling(n_wf / ncol_g)

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow_g, ncol_g)))
for (i in seq_along(wf_list)) {
  print(wf_list[[i]], vp = viewport(layout.pos.row = ceiling(i / ncol_g), layout.pos.col = (i - 1) %% ncol_g + 1))
}

## 5.5) heatmaps -----
Cluster <- factor(cl_B$cluster)
os_summary <- aggregate(Y_full, by = list(Cluster = Cluster), FUN = function(x) round(mean(x, na.rm = TRUE) * 100, 1))
os_map <- setNames(os_summary$x, os_summary$Cluster)
num_vars <- intersect(cont_vars, colnames(df_sub))
raw_num  <- df_sub[, num_vars, drop = FALSE]

raw_means <- aggregate(raw_num, by = list(Cluster = Cluster), FUN = mean)
z_means   <- aggregate(scale(raw_num), by = list(Cluster = Cluster), FUN = mean)
hm_raw_long <- reshape(raw_means, direction = "long", varying = num_vars, v.names = "Raw_Mean", timevar = "Feature", times = num_vars)
hm_z_long   <- reshape(z_means, direction = "long", varying = num_vars, v.names = "Z_Score", timevar = "Feature", times = num_vars)
hm_df <- merge(hm_raw_long[, c("Cluster", "Feature", "Raw_Mean")], hm_z_long[, c("Cluster", "Feature", "Z_Score")], by = c("Cluster", "Feature"))
hm_df$Cell_Label <- sprintf("%.3f", hm_df$Raw_Mean)
clusters_vec <- sort(unique(Cluster))

if (length(cat_vars) > 0) {
  rf_list <- list()
  for (cat_v in cat_vars) {
    rf_list[[cat_v]] <- data.frame(
      Cluster = clusters_vec, Feature = cat_v, Raw_Mean = NA, Z_Score = NA,
      Cell_Label = sapply(clusters_vec, function(cl) {
        tab <- table(df_sub[[cat_v]][Cluster == cl])
        if(length(tab) == 0) return("N/A")
        sprintf("Mode: %s\n(%d%%)", names(which.max(tab)), round(max(tab) / sum(tab) * 100, 0))
      })
    )
  }
  hm_df_final <- rbind(hm_df, do.call(rbind, rf_list))
  hm_df_final$Feature <- factor(hm_df_final$Feature, levels = c(num_vars, cat_vars))
} else {
  hm_df_final <- hm_df
  hm_df_final$Feature <- factor(hm_df_final$Feature, levels = num_vars)
}

hm_df_final$Cluster_Label <- paste0("Cluster ", hm_df_final$Cluster, "\n(OS: ", os_map[as.character(hm_df_final$Cluster)], "%)")
p_heatmap <- ggplot(hm_df_final, aes(x = Feature, y = Cluster_Label, fill = Z_Score)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = Cell_Label), size = 4.2, fontface = "bold") +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0, na.value = "grey94", name = "Z-score") +
  theme_minimal(base_size = 12) + labs(title = "Heatmap of mean clinical profile per Cluster", x = "Variable", y = "Cluster")
print(p_heatmap)

robust_scale <- function(x) (x - median(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
raw_medians <- aggregate(raw_num, by = list(Cluster = Cluster), FUN = median)
z_medians   <- aggregate(as.data.frame(lapply(raw_num, robust_scale)), by = list(Cluster = Cluster), FUN = median)
hm_med_long  <- reshape(raw_medians, direction = "long", varying = num_vars, v.names = "Raw_Median", timevar = "Feature", times = num_vars)
hm_zmed_long <- reshape(z_medians, direction = "long", varying = num_vars, v.names = "Z_Robust", timevar = "Feature", times = num_vars)
hm_df_med <- merge(hm_med_long[, c("Cluster", "Feature", "Raw_Median")], hm_zmed_long[, c("Cluster", "Feature", "Z_Robust")], by = c("Cluster", "Feature"))

hm_df_med$Cell_Label <- mapply(function(cl, feat) {
  val_cl <- df_sub[[feat]][Cluster == cl]
  sprintf("%.3f\n[%.3f - %.3f]", median(val_cl, na.rm=T), quantile(val_cl, 0.25, na.rm=T), quantile(val_cl, 0.75, na.rm=T))
}, hm_df_med$Cluster, hm_df_med$Feature)

if (length(cat_vars) > 0) {
  rf_df_med_list <- list()
  for (cat_v in cat_vars) {
    rf_df_med_list[[cat_v]] <- data.frame(
      Cluster = clusters_vec, Feature = cat_v, Raw_Median = NA, Z_Robust = NA, Cell_Label = rf_list[[cat_v]]$Cell_Label
    )
  }
  hm_df_med_final <- rbind(hm_df_med, do.call(rbind, rf_df_med_list))
  hm_df_med_final$Feature <- factor(hm_df_med_final$Feature, levels = c(num_vars, cat_vars))
} else {
  hm_df_med_final <- hm_df_med
  hm_df_med_final$Feature <- factor(hm_df_med_final$Feature, levels = num_vars)
}

hm_df_med_final$Cluster_Label <- paste0("Cluster ", hm_df_med_final$Cluster, "\n(OS: ", os_map[as.character(hm_df_med_final$Cluster)], "%)")
p_heatmap_med <- ggplot(hm_df_med_final, aes(x = Feature, y = Cluster_Label, fill = Z_Robust)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = Cell_Label), size = 3.8, fontface = "bold") +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0, na.value = "grey94", name = "Z-score") +
  theme_minimal(base_size = 12) + labs(title = "Heatmap of median [IQR] clinical profile", x = "Variable", y = "Cluster")
print(p_heatmap_med)



## 5.6) target -----
df_bar <- as.data.frame(table(Cluster = Cluster, OS2ys = factor(df_sub$OS2ys, levels = c("NO", "YES"))))
tot_cl <- tapply(df_bar$Freq, df_bar$Cluster, sum)
df_bar$Total_N <- as.numeric(tot_cl[as.character(df_bar$Cluster)])
df_bar$Perc <- (df_bar$Freq / df_bar$Total_N) * 100
df_bar$Cluster_Label <- factor(paste0("Cluster ", df_bar$Cluster, "\n(N = ", df_bar$Total_N, ")"), levels = paste0("Cluster ", levels(Cluster), "\n(N = ", as.numeric(tot_cl[levels(Cluster)]), ")"))
os_mean_global <- mean(Y_full) * 100

p_bar_perc <- ggplot(df_bar, aes(x = Cluster_Label, y = Perc, fill = OS2ys)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.6) +
  geom_hline(yintercept = os_mean_global, linetype = "dashed", color = "grey20", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", Perc)), position = position_stack(vjust = 0.5), color = "white", fontface = "bold", size = 3.5) +
  scale_fill_manual(values = c("NO" = "#D55E00", "YES" = "#2B83BA"), labels = c("NO" = "No", "YES" = "Yes")) +
  theme_minimal(base_size = 11) + theme(legend.position = "top", panel.grid.major.x = element_blank(), axis.text.x = element_text(face = "bold", size = 9.5)) +
  labs(title = "A) Proportions of OS2ys by cluster", x = NULL, y = "Percentage (%)", fill = "2-Year OS:")

p_bar_count <- ggplot(df_bar, aes(x = Cluster_Label, y = Freq, fill = OS2ys)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.6) +
  geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), color = "white", fontface = "bold", size = 3.5) +
  scale_fill_manual(values = c("NO" = "#D55E00", "YES" = "#2B83BA"), labels = c("NO" = "No", "YES" = "Yes")) +
  theme_minimal(base_size = 11) + theme(legend.position = "top", panel.grid.major.x = element_blank(), axis.text.x = element_text(face = "bold", size = 9.5)) +
  labs(title = "B) Patient counts by cluster", x = NULL, y = "Number of Patients (n)", fill = "2-Year OS:")

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(p_bar_perc,   vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_bar_count,  vp = viewport(layout.pos.row = 1, layout.pos.col = 2))


# 6) Skope-Rules -----
Cluster_valid <- droplevels(factor(cl_B$cluster[cl_B$cluster > 0]))
X_full_valid  <- X_full[cl_B$cluster > 0, , drop = FALSE]
Y_full_valid  <- Y_full[cl_B$cluster > 0]

df_for_python <- X_full_valid
df_for_python$Cluster <- as.character(Cluster_valid)
df_for_python$OS2ys   <- as.numeric(Y_full_valid)
write.csv(df_for_python, "dataset_clusters_for_skope.csv", row.names = FALSE) # dataset for Python


# 7) survival analysis -----
df_sub$Cluster <- Cluste
df_sub$OS[df_sub$ID == "LM1001"] <- 24
df_sub$event_death_24m <- 1 - df_sub$OS2ys_num
df_sub$OS_24m <- pmin(df_sub$OS, 24)


os_stats <- aggregate(df_sub$OS, list(Cluster = Cluster), function(x) c(Median=median(x,na.rm=T), IQR_low=quantile(x,0.25,na.rm=T), IQR_high=quantile(x,0.75,na.rm=T), Mean=mean(x,na.rm=T), SD=sd(x,na.rm=T)))
tab_os_res <- data.frame(Cluster = os_stats$Cluster, N = as.numeric(table(Cluster)), Eventi_24m = as.numeric(tapply(df_sub$event_death_24m, Cluster, sum)), OS2ys_Vivi = paste0(round(tapply(df_sub$OS2ys_num, Cluster, mean) * 100, 1), "%"), OS_Mediana_Mesi = sprintf("%.1f [%.1f - %.1f]", os_stats$x[, 1], os_stats$x[, 2], os_stats$x[, 3]), OS_Media_Mesi = sprintf("%.1f (±%.1f)", os_stats$x[, 4], os_stats$x[, 5]))
print(tab_os_res, row.names = FALSE)


df_surv_valid <- df_sub[df_sub$Cluster != "0", ]
df_surv_valid$Cluster <- droplevels(df_surv_valid$Cluster)

cox_mod <- coxph(Surv(OS_24m, event_death_24m) ~ Cluster, data = df_surv_valid)
cox_sum <- summary(cox_mod)
tab_cox <- data.frame(Confronto = rownames(cox_sum$conf.int), HR = round(cox_sum$conf.int[, "exp(coef)"], 2), CI_95 = sprintf("[%.2f - %.2f]", cox_sum$conf.int[, "lower .95"], cox_sum$conf.int[, "upper .95"]), p_value = signif(cox_sum$coefficients[, "Pr(>|z|)"], 3))
print(tab_cox, row.names = FALSE)

fit_km    <- survfit(Surv(OS_24m, event_death_24m) ~ Cluster, data = df_surv_valid)
lr_test   <- survdiff(Surv(OS_24m, event_death_24m) ~ Cluster, data = df_surv_valid)
p_logrank <- 1 - pchisq(lr_test$chisq, length(lr_test$n) - 1)

km_df <- data.frame(time = fit_km$time, surv = fit_km$surv * 100, Cluster = rep(gsub("Cluster=", "Cluster ", names(fit_km$strata)), fit_km$strata))
km_df <- rbind(data.frame(time = 0, surv = 100, Cluster = unique(km_df$Cluster)), km_df)

p_km_24m <- ggplot(km_df, aes(x = time, y = surv, color = Cluster)) +
  geom_step(linewidth = 1.2) + scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) + scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 3)) +
  theme_minimal(base_size = 12) + theme(legend.position = "top", panel.grid.minor = element_blank()) +
  labs(title = "Kaplan-Meier curves in the first 24 months by Cluster", subtitle = sprintf("Log-Rank Test: Chisq = %.1f, p = %.2e", lr_test$chisq, p_logrank), x = "months from RT", y = "Survival estimate (%)")

grid.newpage()
print(p_km_24m, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))


##new KM plot
library(survival)
library(survminer)
fit_km <- survfit(Surv(OS_24m, event_death_24m) ~ Cluster, data = df_surv_valid)
lr_test <- survdiff(Surv(OS_24m, event_death_24m) ~ Cluster, data = df_surv_valid)
p_val_num <- 1 - pchisq(lr_test$chisq, length(lr_test$n) - 1)
sub_label <- sprintf("Log-rank test: p = %.2e", p_val_num)
p_km_survminer <- ggsurvplot(
  fit_km,
  data = df_surv_valid,
  pval = FALSE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.25,
  risk.table.y.text = FALSE,
  risk.table.title = "Number at risk (every 3 months)",
  tables.theme = theme_cleantable(),
  xlim = c(0, 24),
  break.time.by = 3,
  ylim = c(0, 1),
  palette = "npg",
  censor = TRUE,
  censor.shape = 124,
  censor.size = 3,
  size = 1.1,
  legend.title = "Cluster",
  legend.labs = levels(df_surv_valid$Cluster),
  xlab = "Time (months)",
  ylab = "Survival probability",
  title = "Overall survival at 24 months",
  ggtheme = theme_classic(base_size = 12)
)
p_km_survminer$plot <- p_km_survminer$plot + 
  labs(subtitle = sub_label)
print(p_km_survminer)



# 8) stability and sensitivity analysis ------
ari_AB    <- mclust::adjustedRandIndex(cl_A$cluster, cl_B$cluster)
ami_AB    <- aricode::AMI(cl_A$cluster, cl_B$cluster)
ari_B_oof <- mclust::adjustedRandIndex(cl_B$cluster, cl_oof$cluster)
ami_B_oof <- aricode::AMI(cl_B$cluster, cl_oof$cluster)
ari_BC    <- mclust::adjustedRandIndex(cl_B$cluster, cl_C$cluster)
ami_BC    <- aricode::AMI(cl_B$cluster, cl_C$cluster)

n_seed_test <- 30
lab_seeds <- sapply(1:n_seed_test, function(s) {
  set.seed(s)
  u_s <- uwot::umap(shap_agg, n_neighbors = nn_umap, min_dist = 0, n_components = 2, n_threads = 1, n_sgd_threads = 1)
  dbscan::hdbscan(u_s, minPts = minPts_hdb)$cluster
})
ari_mat <- outer(1:n_seed_test, 1:n_seed_test, Vectorize(function(i, j) mclust::adjustedRandIndex(lab_seeds[, i], lab_seeds[, j])))
ami_mat <- outer(1:n_seed_test, 1:n_seed_test, Vectorize(function(i, j) aricode::AMI(lab_seeds[, i], lab_seeds[, j])))
stab_2d <- mean(ari_mat[upper.tri(ari_mat)])
stab_2d_ami <- mean(ami_mat[upper.tri(ami_mat)])

lab_seeds_xgb <- sapply(1:n_seed_test, function(s_xgb) {
  set.seed(s_xgb)
  fit_s <- xgb.train(
    params  = list(objective = "binary:logistic", max_depth = hp_xgb$max_depth, eta = hp_xgb$eta, lambda = hp_xgb$lambda, subsample = hp_xgb$subsample),
    data    = xgb.DMatrix(X_mat, label = Y_full), nrounds = hp_xgb$nrounds, verbose = 0
  )
  s_raw <- as.matrix(treeshap(xgboost.unify(fit_s, X_mat), X_mat, verbose = FALSE)$shaps)
  s_agg <- aggregate_shap(s_raw, X_full, covariates)
  set.seed(seed_umap)
  u_s <- uwot::umap(s_agg, n_neighbors = nn_umap, min_dist = 0, n_components = 2, n_threads = 1, n_sgd_threads = 1)
  dbscan::hdbscan(u_s, minPts = minPts_hdb)$cluster
})
ari_mat_xgb <- outer(1:n_seed_test, 1:n_seed_test, Vectorize(function(i, j) mclust::adjustedRandIndex(lab_seeds_xgb[, i], lab_seeds_xgb[, j])))
ami_mat_xgb <- outer(1:n_seed_test, 1:n_seed_test, Vectorize(function(i, j) aricode::AMI(lab_seeds_xgb[, i], lab_seeds_xgb[, j])))
stab_xgb <- mean(ari_mat_xgb[upper.tri(ari_mat_xgb)])
stab_xgb_ami <- mean(ami_mat_xgb[upper.tri(ami_mat_xgb)])

cat("Cluster found on 30 seed XGBoost (on SHAP aggregated):", apply(lab_seeds_xgb, 2, function(cl) length(unique(cl[cl > 0]))), "\n")
cat("\nARI (No-UMAP vs UMAP):", round(ari_AB, 3), " AMI:", round(ami_AB, 3),
    "\nARI (UMAP vs OOF-UMAP):", round(ari_B_oof, 3), " AMI:", round(ami_B_oof, 3), 
    "\nStability UMAP-2D (30 seed) ARI:", round(stab_2d, 3), " AMI:", round(stab_2d_ami, 3),
    "\nStability XGBoost (30 seed) ARI:", round(stab_xgb, 3), " AMI:", round(stab_xgb_ami, 3), "\n")


## Bagged xgb -----
n_blocks <- 10
n_bag    <- 10

lab_bag_xgb <- sapply(1:n_blocks, function(b) {
  shap_list <- lapply(((b - 1) * n_bag + 1):(b * n_bag), function(s_xgb) {
    set.seed(s_xgb)
    fit_s <- xgb.train(
      params  = list(objective = "binary:logistic", max_depth = hp_xgb$max_depth, eta = hp_xgb$eta, lambda = hp_xgb$lambda, subsample = hp_xgb$subsample),
      data    = xgb.DMatrix(X_mat, label = Y_full), nrounds = hp_xgb$nrounds, verbose = 0
    )
    s_r <- as.matrix(treeshap(xgboost.unify(fit_s, X_mat), X_mat, verbose = FALSE)$shaps)
    aggregate_shap(s_r, X_full, covariates)
  })
  shap_mean_b <- Reduce("+", shap_list) / n_bag
  set.seed(seed_umap)
  u_s <- uwot::umap(shap_mean_b, n_neighbors = nn_umap, min_dist = 0, n_components = 2, n_threads = 1, n_sgd_threads = 1)
  dbscan::hdbscan(u_s, minPts = minPts_hdb)$cluster
})
ari_mat_bag <- outer(1:n_blocks, 1:n_blocks, Vectorize(function(i, j) mclust::adjustedRandIndex(lab_bag_xgb[, i], lab_bag_xgb[, j])))
ami_mat_bag <- outer(1:n_blocks, 1:n_blocks, Vectorize(function(i, j) aricode::AMI(lab_bag_xgb[, i], lab_bag_xgb[, j])))
stab_bag_xgb <- mean(ari_mat_bag[upper.tri(ari_mat_bag)])
stab_bag_xgb_ami <- mean(ami_mat_bag[upper.tri(ami_mat_bag)])

cat("Stability XGBoost (Ensemble 10 seed on SHAP aggregated) ARI:", round(stab_bag_xgb, 3), " AMI:", round(stab_bag_xgb_ami, 3), "\n")
cat("cluster per run :", apply(lab_bag_xgb, 2, function(cl) length(unique(cl[cl > 0]))), "\n")
cat("noise pts per run :", apply(lab_bag_xgb, 2, function(cl) sum(cl == 0)), "\n")


ari_single_vs_bag <- sapply(1:n_blocks, function(b) mclust::adjustedRandIndex(cl_B$cluster, lab_bag_xgb[, b]))
ami_single_vs_bag <- sapply(1:n_blocks, function(b) aricode::AMI(cl_B$cluster, lab_bag_xgb[, b]))
cat("\nmean ARI (one vs bagged):", round(mean(ari_single_vs_bag), 3), " mean AMI:", round(mean(ami_single_vs_bag), 3), "\n")

cl_bag1 <- lab_bag_xgb[, 1]
print(table(Modello_Singolo_clB = cl_B$cluster, Bagged_XGB = cl_bag1))#contingency table

profilo_clinico_cmp <- function(cl_vec, nome_metodo) {
  cl_f <- factor(cl_vec)
  res_df <- data.frame(Metodo = nome_metodo, Cluster = levels(cl_f), n = as.numeric(table(cl_f)), OS_2ys_perc = round(as.numeric(tapply(Y_full, cl_f, mean)) * 100, 1))
  if (length(cont_vars) > 0) {
    m_cont <- aggregate(X_full[, cont_vars, drop = FALSE], list(Cluster = cl_f), mean)
    res_df <- cbind(res_df, round(m_cont[, -1, drop = FALSE], 2))
  }
  return(res_df)
}

cmp_tab <- rbind(profilo_clinico_cmp(cl_B$cluster, "1. One (seed 123)"), profilo_clinico_cmp(cl_bag1, "2. Bagged (10 seed)"))
print(cmp_tab, row.names = FALSE)

## plots full vs OOF SHAP-----
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
p_full <- ggplot(df_oof_cmp, aes(x = Full_UMAP1, y = Full_UMAP2, color = Cluster_Full)) + geom_point(alpha = 0.85, size = 2.2) + theme_minimal(base_size = 12) + labs(title = "UMAP 2D (SHAP Full Model )", x = "UMAP 1", y = "UMAP 2", color = "Cluster Full")
p_oof <- ggplot(df_oof_cmp, aes(x = OOF_UMAP1, y = OOF_UMAP2, color = Cluster_OOF)) + geom_point(alpha = 0.85, size = 2.2) + theme_minimal(base_size = 12) + labs(title = "UMAP 2D (SHAP Out-Of-Fold)", subtitle = paste("ARI =", round(ari_B_oof, 3)), x = "UMAP 1", y = "UMAP 2", color = "Cluster OOF")

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(p_full, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_oof,  vp = viewport(layout.pos.row = 1, layout.pos.col = 2))

print(table(Cluster_Full = df_oof_cmp$Cluster_Full, Cluster_OOF = df_oof_cmp$Cluster_OOF))

p_full_col <- ggplot(df_oof_cmp, aes(x = Full_UMAP1, y = Full_UMAP2, color = Cluster_Full)) + geom_point(alpha = 0.85, size = 2.2) + theme_minimal(base_size = 12) + labs(title = "UMAP 2D (SHAP Full Model)", x = "UMAP 1", y = "UMAP 2")
p_oof_col <- ggplot(df_oof_cmp, aes(x = OOF_UMAP1, y = OOF_UMAP2, color = Cluster_Full, shape = Cluster_OOF)) + geom_point(alpha = 0.85, size = 2.4) + theme_minimal(base_size = 12) + labs(title = "UMAP 2D (SHAP Out-Of-Fold)", subtitle = "colored by Cluster Full")

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(p_full_col, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_oof_col, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))





# additional material ------
df_oof_cmp <- data.frame(
  Full_UMAP1   = um2_shap[, 1],
  Full_UMAP2   = um2_shap[, 2],
  OOF_UMAP1    = um2_oof[, 1],
  OOF_UMAP2    = um2_oof[, 2],
  RAW_UMAP1    = um2_raw[, 1],
  RAW_UMAP2    = um2_raw[, 2],
  Cluster_Full = factor(ifelse(cl_B$cluster == 0, "0 (Noise)", as.character(cl_B$cluster))),
  Cluster_OOF  = factor(ifelse(cl_oof$cluster == 0, "0 (Noise)", as.character(cl_oof$cluster))),
  OS2ys        = factor(df_sub$OS2ys, levels = c("NO", "YES"))
)

## plot 1 ------
df_plot$Cluster_clA <- factor(ifelse(cl_A$cluster == 0, "0 (Noise)", as.character(cl_A$cluster)))
livelli_clA <- levels(df_plot$Cluster_clA)
colori_clA  <- setNames(scales::hue_pal()(length(livelli_clA)), livelli_clA)
if ("0 (Noise)" %in% livelli_clA) colori_clA["0 (Noise)"] <- "grey70"

p_shap_clB <- ggplot(df_plot, aes(x = SHAP_UMAP1, y = SHAP_UMAP2, color = Cluster)) + geom_point(alpha = 0.85, size = 2) + theme_minimal(base_size = 11) + labs(title = "SHAP UMAP space", x = "UMAP 1", y = "UMAP 2", color = "cluster full")
p_raw_clB <- ggplot(df_plot, aes(x = RAW_UMAP1, y = RAW_UMAP2, color = Cluster)) + geom_point(alpha = 0.85, size = 2) + theme_minimal(base_size = 11) + labs(title = "Raw data UMAP space", x = "UMAP 1", y = "UMAP 2", color = "cluster full")
p_shap_clA <- ggplot(df_plot, aes(x = SHAP_UMAP1, y = SHAP_UMAP2, color = Cluster_clA)) + geom_point(alpha = 0.85, size = 2) + scale_color_manual(values = colori_clA) + theme_minimal(base_size = 11) + labs(title = "SHAP UMAP space", x = "UMAP 1", y = "UMAP 2", color = "cluster raw")
p_raw_clA <- ggplot(df_plot, aes(x = RAW_UMAP1, y = RAW_UMAP2, color = Cluster_clA)) + geom_point(alpha = 0.85, size = 2) + scale_color_manual(values = colori_clA) + theme_minimal(base_size = 11) + labs(title = "Raw data UMAP space", x = "UMAP 1", y = "UMAP 2", color = "cluster raw")

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
print(p_shap_clB, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_raw_clB,  vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p_shap_clA, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(p_raw_clA,  vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

## plot 2 ------
p_shap_full <- ggplot(df_oof_cmp, aes(x = Full_UMAP1, y = Full_UMAP2, color = Cluster_Full)) + geom_point(alpha = 0.85, size = 2) + theme_minimal(base_size = 11) + labs(title = "SHAP UMAP space", x = "UMAP 1 (Full)", y = "UMAP 2 (Full)", color = "cluster full")
p_shap_oof <- ggplot(df_oof_cmp, aes(x = OOF_UMAP1, y = OOF_UMAP2, color = Cluster_OOF)) + geom_point(alpha = 0.85, size = 2) + theme_minimal(base_size = 11) + labs(title = "SHAP OOF UMAP space", x = "UMAP 1 (OOF)", y = "UMAP 2 (OOF)", color = "cluster OOF")
p_raw_full <- ggplot(df_oof_cmp, aes(x = RAW_UMAP1, y = RAW_UMAP2, color = Cluster_Full)) + geom_point(alpha = 0.85, size = 2) + theme_minimal(base_size = 11) + labs(title = "Raw data UMAP space", x = "UMAP 1", y = "UMAP 2", color = "cluster full")
p_raw_oof <- ggplot(df_oof_cmp, aes(x = RAW_UMAP1, y = RAW_UMAP2, color = Cluster_OOF)) + geom_point(alpha = 0.85, size = 2) + theme_minimal(base_size = 11) + labs(title = "Raw data UMAP space", x = "UMAP 1", y = "UMAP 2", color = "cluster OOF")

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
print(p_shap_full, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_shap_oof,  vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p_raw_full,  vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(p_raw_oof,   vp = viewport(layout.pos.row = 2, layout.pos.col = 2))







