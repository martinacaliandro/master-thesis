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

data <- read.csv('dataclean_final_withNA.csv', header=T)
head(data)
dim(data)
View(data)

#create 3 dataset ----
sum(is.na(data))
colnames(data)
colnames(data)[colnames(data) == "Overall.Survival.at.2.years"] = "OS2ys"
colnames(data)[colnames(data) == "mean_dose_vanacava_S"] <- "mean_dose_venacava_S"
colnames(data)[colnames(data) == "eud_venacavaI"] <- "eud_venacava_I"
colnames(data)[colnames(data) == "eud_vanacavaS"] <- "eud_venacava_S"
data[data == ""] <- NA
data[data$OS..yes.no. == "", "OS..yes.no."] = NA

# fix LM6001 -----
#data[data$CHEMO == "NO" & data$schedule_ok_chemo != 99, c("ID", "CHEMO", "schedule_ok_chemo")]
#data$CHEMO[data$CHEMO == "NO" & data$schedule_ok_chemo != 99] <- "YES"
#----------------
#make factors cat num cols 
data$KPS=factor(data$KPS)
data$PDL.1.level=factor(data$PDL.1.level)
data$schedule_ok_chemo=factor(data$schedule_ok_chemo)
data$cycles_chemo_sum = factor(data$cycles_chemo_sum)
data$cycles_immuno = factor(data$cycles_immuno)
data$prescription_dose_radio = factor(data$prescription_dose_radio)
data$number_of_fractions_radio = factor(data$number_of_fractions_radio)
#make factors all cat cols
data[] <- lapply(data, function(x) {
  if (!is.numeric(x)) as.factor(x) else x
})
#----------------
datacrf = data[,1:40]
datadose = data[,41:63]
datavol = data[,64:76]
datametab = data[,77:89]

#identify cat cols in data CRF ----
datacrf_num <- datacrf[, sapply(datacrf, is.numeric)]
colnames(datacrf_num)
datacrf_cat <- datacrf[, !sapply(datacrf, is.numeric)]
colnames(datacrf_cat)
dim(datacrf)[2]
dim(datacrf_num)[2]
dim(datacrf_cat)[2]








#look at all correlations ----
cor = cor(cbind(datacrf_num,datadose,datametab), use='complete.obs')
corrplot(cor,
         method = "color",
         addCoef.col = "black",
         number.cex = 0.3,
         tl.cex = 0.5,
         tl.col = "grey30",
         col = colorRampPalette(c("darkblue", "white", "darkred"))(200))
#zooming 
cor = cor(datadose)
corrplot(cor, method='col',addCoef.col = 'black', number.cex = 0.6, tl.cex = 0.5,tl.col='grey', col = colorRampPalette(c("darkblue", "white", "darkred"))(200))

cor = cor(datametab)
corrplot(cor, method='col',addCoef.col = 'black', number.cex = 0.6, col = colorRampPalette(c("darkblue", "white", "darkred"))(200))







# some scatter plots for VOLUMETRIC and DOSIMETRIC outliers ---------
ggplot(data, aes(x = seq_along(data$A_Aorta), y = data$A_Aorta, color = OS2ys)) +
  geom_point(size = 2)

ggplot(data, aes(x = seq_along(data$A_Lad), y = data$A_Lad, color = OS2ys)) +
  geom_point(size = 2)
data[data$A_Lad>12,"ID"]
ggplot(data, aes(x = seq_along(data$mean_dose_lad), y = data$mean_dose_lad, color = OS2ys)) +
  geom_point(size = 2)
data[data$mean_dose_lad>35,"ID"]
data[data$ID == "LM7031", c("A_Lad","mean_dose_lad")]
data[data$ID == "LM1001", c("A_Lad","mean_dose_lad")] #NO
data[data$ID == "LM1012", c("A_Lad","mean_dose_lad")]
data[data$ID == "LM1016", c("A_Lad","mean_dose_lad")] #NO
data[data$ID == "LM3027", c("A_Lad","mean_dose_lad")]
data[data$ID == "LM7053", c("A_Lad","mean_dose_lad")] #NO
data[data$ID == "LM7096", c("A_Lad","mean_dose_lad")]
ggplot(data, aes(x = seq_along(data$eud_lad), y = data$eud_lad, color = OS2ys)) +
  geom_point(size = 2)

ggplot(data, aes(x = seq_along(data$A_Pulm), y = data$A_Pulm, color = OS2ys)) +
  geom_point(size = 2)

ggplot(data, aes(x = seq_along(data$Atrium_L), y = data$Atrium_L, color = OS2ys)) +
  geom_point(size = 2)
data[data$Atrium_L>200,"ID"]
ggplot(data, aes(x = seq_along(data$mean_dose_atriumL), y = data$mean_dose_atriumL, color = OS2ys)) +
  geom_point(size = 2)
ggplot(data, aes(x = seq_along(data$eud_atriumL), y = data$eud_atriumL, color = OS2ys)) +
  geom_point(size = 2)

ggplot(data, aes(x = seq_along(data$Atrium_R), y = data$Atrium_R, color = OS2ys)) +
  geom_point(size = 2)
data[data$Atrium_R>200,"ID"]
ggplot(data, aes(x = seq_along(data$mean_dose_atriumR), y = data$mean_dose_atriumR, color = OS2ys)) +
  geom_point(size = 2)
data[data$mean_dose_atriumR>30,"ID"]
data[data$ID == "LM1073", c("Atrium_R","mean_dose_atriumR")]
ggplot(data, aes(x = seq_along(data$eud_atriumR), y = data$eud_atriumR, color = OS2ys)) +
  geom_point(size = 2)

ggplot(data, aes(x = seq_along(data$Heart), y = data$Heart, color = OS2ys)) +
  geom_point(size = 2)
data[data$Heart>1500,"ID"]
ggplot(data, aes(x = seq_along(data$mean_dose_Heart), y = data$mean_dose_Heart, color = OS2ys)) +
  geom_point(size = 2)
data[data$mean_dose_Heart>20,"ID"]
data[data$ID == "LM1001", c("Heart","mean_dose_Heart")]
data[data$ID == "LM1023", c("Heart","mean_dose_Heart")]
data[data$ID == "LM1077", c("Heart","mean_dose_Heart")]
data[data$ID == "LM1087", c("Heart","mean_dose_Heart")]
data[data$ID == "LM2006", c("Heart","mean_dose_Heart")]
data[data$ID == "LM2071", c("Heart","mean_dose_Heart")]
data[data$ID == "LM6069", c("Heart","mean_dose_Heart")]
data[data$ID == "LM7053", c("Heart","mean_dose_Heart")]
data[data$ID == "LM7057", c("Heart","mean_dose_Heart")]
data[data$ID == "LM7093", c("Heart","mean_dose_Heart")]
ggplot(data, aes(x = seq_along(data$eud_Heart), y = data$eud_Heart, color = OS2ys)) +
  geom_point(size = 2)



ggplot(data, aes(x = seq_along(data$Lung_L), y = data$Lung_L, color = OS2ys)) +
  geom_point(size = 2)
data[data$Lung_L>3000,"ID"]

ggplot(data, aes(x = seq_along(data$Lung_R), y = data$Lung_R, color = OS2ys)) +
  geom_point(size = 2)
data[data$Lung_R>3000,"ID"] # biggest lung R also biggest lung L for LM1056, LM4058, LM7056

ggplot(data, aes(x = seq_along(data$Lungs), y = data$Lungs, color = OS2ys)) +
  geom_point(size = 2)
data[data$Lungs>6000,"ID"] # LM1056, LM4058, LM7056 included
ggplot(data, aes(x = seq_along(data$mean_dose_lungs), y = data$mean_dose_lungs, color = OS2ys)) +
  geom_point(size = 2)
data[data$mean_dose_lungs>20,"ID"]
data[data$ID == "LM1013", c("Lungs","mean_dose_lungs")] #low
data[data$ID == "LM3012", c("Lungs","mean_dose_lungs")]
data[data$ID == "LM3044", c("Lungs","mean_dose_lungs")]
data[data$ID == "LM3046", c("Lungs","mean_dose_lungs")]
data[data$ID == "LM3050", c("Lungs","mean_dose_lungs")]
data[data$ID == "LM3052", c("Lungs","mean_dose_lungs")]
ggplot(data, aes(x = seq_along(data$eud_lungs), y = data$eud_lungs, color = OS2ys)) +
  geom_point(size = 2)
data[data$eud_lungs<20,"ID"]
data[data$ID == "LM3028", c("Lungs","eud_lungs")]
data[data$ID == "LM3031", c("Lungs","eud_lungs")]
data[data$ID == "LM3033", c("Lungs","eud_lungs")]
data[data$ID == "LM6040", c("Lungs","eud_lungs")]

ggplot(data, aes(x = seq_along(data$Ventricle_L), y = data$Ventricle_L, color = OS2ys)) +
  geom_point(size = 2)
data[data$Ventricle_L>350,"ID"]
data[data$Ventricle_L>400,"ID"]
ggplot(data, aes(x = seq_along(data$mean_dose_ventricleL), y = data$mean_dose_ventricleL, color = OS2ys)) +
  geom_point(size = 2)
data[data$mean_dose_ventricleL>20,"ID"]
data[data$ID == "LM1012", c("Ventricle_L","mean_dose_ventricleL")]
data[data$ID == "LM3058", c("Ventricle_L","mean_dose_ventricleL")]
data[data$ID == "LM7053", c("Ventricle_L","mean_dose_ventricleL")]
ggplot(data, aes(x = seq_along(data$eud_ventricleL), y = data$eud_ventricleL, color = OS2ys)) +
  geom_point(size = 2)

ggplot(data, aes(x = seq_along(data$Ventricle_R), y = data$Ventricle_R, color = OS2ys)) +
  geom_point(size = 2)
data[data$Ventricle_R>200,"ID"]
ggplot(data, aes(x = seq_along(data$mean_dose_ventricleR), y = data$mean_dose_ventricleR, color = OS2ys)) +
  geom_point(size = 2)
data[data$mean_dose_ventricleR>20,"ID"]
data[data$ID == "LM1001", c("Ventricle_R","mean_dose_ventricleR")]
data[data$ID == "LM1012", c("Ventricle_R","mean_dose_ventricleR")]
data[data$ID == "LM7057", c("Ventricle_R","mean_dose_ventricleR")]
data[data$ID == "LM7068", c("Ventricle_R","mean_dose_ventricleR")]
ggplot(data, aes(x = seq_along(data$eud_ventricleR), y = data$eud_ventricleR, color = OS2ys)) +
  geom_point(size = 2)
data[data$eud_ventricleR>40,"ID"]

ggplot(data, aes(x = seq_along(data$Venacava_I), y = data$Venacava_I, color = OS2ys)) +
  geom_point(size = 2)
data[data$Venacava_I>75,"ID"]
ggplot(data, aes(x = seq_along(data$mean_dose_venacava_I), y = data$mean_dose_venacava_I, color = OS2ys)) +
  geom_point(size = 2)
data[data$mean_dose_venacava_I>20,"ID"]
data[data$ID == "LM1058", c("Venacava_I","mean_dose_venacava_I")]
data[data$ID == "LM1077", c("Venacava_I","mean_dose_venacava_I")]
data[data$ID == "LM1082", c("Venacava_I","mean_dose_venacava_I")]
ggplot(data, aes(x = seq_along(data$eud_venacava_I), y = data$eud_venacava_I, color = OS2ys)) +
  geom_point(size = 2)
data[data$eud_venacava_I>40,"ID"]
data[data$ID == "LM1028", c("Venacava_I","eud_venacava_I")]
data[data$ID == "LM1077", c("Venacava_I","eud_venacava_I")]
data[data$ID == "LM1082", c("Venacava_I","eud_venacava_I")]
data[data$ID == "LM3049", c("Venacava_I","eud_venacava_I")]
data[data$ID == "LM5025", c("Venacava_I","eud_venacava_I")]

ggplot(data, aes(x = seq_along(data$Venacava_S), y = data$Venacava_S, color = OS2ys)) +
  geom_point(size = 2)
## some plots ----------
### LAD ----
data$label_ID <- ifelse(data$A_Lad > 12 | data$mean_dose_lad > 35, 
                        as.character(data$ID), 
                        "")
ggplot(data, aes(x = A_Lad, y = mean_dose_lad, color = OS2ys)) +
  geom_point(size = 3, alpha = 0.8) + 
  geom_vline(xintercept = 12, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_hline(yintercept = 35, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_label_repel(aes(label = label_ID), 
                   size = 4,
                   fontface = "bold",
                   box.padding = 0.8,
                   point.padding = 0.5,
                   max.overlaps = Inf,
                   show.legend = FALSE,
                   fill = alpha("white", 0.8),
                   color = "black") +
  theme_minimal(base_size = 14) +
  labs(x = "A_Lad",
       y = "Mean Dose LAD",
       color = "OS2ys")

### Atrium R ----
col1 <- "Atrium_R"
col2 <- "mean_dose_atriumR"
val1 <- 200
val2 <- 30
col_color <- "OS2ys"
data$label_ID <- ifelse(data[[col1]] > val1 | data[[col2]] > val2, 
                        as.character(data$ID), 
                        "")
ggplot(data, aes(x = .data[[col1]], y = .data[[col2]], color = .data[[col_color]])) +
  geom_point(size = 3, alpha = 0.8) + 
  geom_vline(xintercept = val1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_hline(yintercept = val2, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_label_repel(aes(label = label_ID), 
                   size = 4,
                   fontface = "bold",
                   box.padding = 0.8,
                   point.padding = 0.5,
                   max.overlaps = Inf,
                   show.legend = FALSE,
                   fill = alpha("white", 0.8),
                   color = "black") +
  theme_minimal(base_size = 14) +
  labs(x = col1,
       y = col2,
       color = col_color)
### Heart ----
col1 <- "Heart"
col2 <- "mean_dose_Heart"
val1 <- 1500
val2 <- 20
col_color <- "OS2ys"
data$label_ID <- ifelse(data[[col1]] > val1 | data[[col2]] > val2, 
                        as.character(data$ID), 
                        "")
ggplot(data, aes(x = .data[[col1]], y = .data[[col2]], color = .data[[col_color]])) +
  geom_point(size = 3, alpha = 0.8) + 
  geom_vline(xintercept = val1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_hline(yintercept = val2, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_label_repel(aes(label = label_ID), 
                   size = 4,
                   fontface = "bold",
                   box.padding = 0.8,
                   point.padding = 0.5,
                   max.overlaps = Inf,
                   show.legend = FALSE,
                   fill = alpha("white", 0.8),
                   color = "black") +
  theme_minimal(base_size = 14) +
  labs(x = col1,
       y = col2,
       color = col_color)
### Lungs ----
col1 <- "Lungs"
col2 <- "mean_dose_lungs"
val1 <- 6000
val2 <- 20
col_color <- "OS2ys"
data$label_ID <- ifelse(data[[col1]] > val1 | data[[col2]] > val2, 
                        as.character(data$ID), 
                        "")
ggplot(data, aes(x = .data[[col1]], y = .data[[col2]], color = .data[[col_color]])) +
  geom_point(size = 3, alpha = 0.8) + 
  geom_vline(xintercept = val1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_hline(yintercept = val2, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_label_repel(aes(label = label_ID), 
                   size = 4,
                   fontface = "bold",
                   box.padding = 0.8,
                   point.padding = 0.5,
                   max.overlaps = Inf,
                   show.legend = FALSE,
                   fill = alpha("white", 0.8),
                   color = "black") +
  theme_minimal(base_size = 14) +
  labs(x = col1,
       y = col2,
       color = col_color)
### Ventricle L ----
col1 <- "Ventricle_L"
col2 <- "mean_dose_ventricleL"
val1 <- 350
val2 <-20
col_color <- "OS2ys"
data$label_ID <- ifelse(data[[col1]] > val1 | data[[col2]] > val2, 
                        as.character(data$ID), 
                        "")
ggplot(data, aes(x = .data[[col1]], y = .data[[col2]], color = .data[[col_color]])) +
  geom_point(size = 3, alpha = 0.8) + 
  geom_vline(xintercept = val1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_hline(yintercept = val2, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_label_repel(aes(label = label_ID), 
                   size = 4,
                   fontface = "bold",
                   box.padding = 0.8,
                   point.padding = 0.5,
                   max.overlaps = Inf,
                   show.legend = FALSE,
                   fill = alpha("white", 0.8),
                   color = "black") +
  theme_minimal(base_size = 14) +
  labs(x = col1,
       y = col2,
       color = col_color)
### Ventricle L ----
col1 <- "Ventricle_R"
col2 <- "mean_dose_ventricleR"
val1 <- 200
val2 <- 20
col_color <- "OS2ys"
data$label_ID <- ifelse(data[[col1]] > val1 | data[[col2]] > val2, 
                        as.character(data$ID), 
                        "")
ggplot(data, aes(x = .data[[col1]], y = .data[[col2]], color = .data[[col_color]])) +
  geom_point(size = 3, alpha = 0.8) + 
  geom_vline(xintercept = val1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_hline(yintercept = val2, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_label_repel(aes(label = label_ID), 
                   size = 4,
                   fontface = "bold",
                   box.padding = 0.8,
                   point.padding = 0.5,
                   max.overlaps = Inf,
                   show.legend = FALSE,
                   fill = alpha("white", 0.8),
                   color = "black") +
  theme_minimal(base_size = 14) +
  labs(x = col1,
       y = col2,
       color = col_color)
### Venacava I ----
col1 <- "Venacava_I"
col2 <- "mean_dose_venacava_I"
val1 <- 75
val2 <- 20
col_color <- "OS2ys"
data$label_ID <- ifelse(data[[col1]] > val1 | data[[col2]] > val2, 
                        as.character(data$ID), 
                        "")
ggplot(data, aes(x = .data[[col1]], y = .data[[col2]], color = .data[[col_color]])) +
  geom_point(size = 3, alpha = 0.8) + 
  geom_vline(xintercept = val1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_hline(yintercept = val2, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_label_repel(aes(label = label_ID), 
                   size = 4,
                   fontface = "bold",
                   box.padding = 0.8,
                   point.padding = 0.5,
                   max.overlaps = Inf,
                   show.legend = FALSE,
                   fill = alpha("white", 0.8),
                   color = "black") +
  theme_minimal(base_size = 14) +
  labs(x = col1,
       y = col2,
       color = col_color)













#outliers IMPUTATION (for plots)----
dataclean = data
#aorta agatston
length(data$aorta_agatston[data$aorta_agatston >= 1000]) #168
length(data$aorta_agatston[data$aorta_agatston < 0]) #0
length(data$aorta_agatston[data$aorta_agatston == 0]) #102
data$aorta_agatston[data$aorta_agatston >= 1000] <- NA #/!\ PUT NA then impute

#fat 5033
length(data$Subcutaneous_Fat[data$Subcutaneous_Fat >= 1])#1
print(data[which(data$Subcutaneous_Fat >= 1),1]) #ID: 5033
data$Subcutaneous_Fat[data$Subcutaneous_Fat >= 1] <- NA #/!\ PUT NA then impute
data$Torso_Fat[data$Torso_Fat >= 0.9] <- NA
data$Skeletal_Muscle[data$Skeletal_Muscle >= 1] <- NA
data$Intermuscular_Fat[data$Intermuscular_Fat >= 0.1] <- NA

#agatston
length(data[which(data$heart_agatston>5000),1]) #18
length(data[which(data$heart_agatston>10000),1])
data[which(data$heart_agatston>20000),1]
data[which(data$heart_agatston>50000),1] #LM1039
data$heart_agatston[data$heart_agatston >= 5000] <- NA
data$aorta_agatston[data$aorta_agatston >= 5000] <- NA
data$heart_volume_mm3[data$heart_volume_mm3 >= 5000] <- NA
data$aorta_volume_mm3[data$aorta_volume_mm3 >= 5000] <- NA


#volume
data[which(data$heart_volume_mm3>5000),1]

dim(data)
#data = kNN(data, k=5, imp_var=FALSE)
dim(data)

#investigate linear corr wihtin DOSIMETRIC features ----
##Age ----
ggplot(data, aes(x = Age, y = diameter.aortic.arc.1..mm., color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) 

ggplot(data, aes(x = Age, y = diameter.aortic.arc.2..mm., color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) 

ggplot(data, aes(x = Age, y = diameter.kidney.sep..mm., color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) 

ggplot(data, aes(x = Age, y = aorta_agatston, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) # no relevant discoveries

ggplot(data, aes(x = Age, y = Volume_PTV, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) 

##BMI---
hist(datacrfnew$BMI,  breaks = 50, col = "lightblue") #obesity moderate from BMI>35
ggplot(data, aes(x = BMI, y = EDIC, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) #weak lin
ggplot(data, aes(x = BMI, y = Subcutaneous_Fat, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) #lin
data[which(data$Subcutaneous_Fat>1),c("ID","Subcutaneous_Fat")]
ggplot(data, aes(x = BMI, y = Torso_Fat, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) 
ggplot(data, aes(x = BMI, y = Skeletal_Muscle, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2)
ggplot(data, aes(x = BMI, y = Intermuscular_Fat, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2)
cor(data$BMI, data$Subcutaneous_Fat)
cor(data$BMI, data$Torso_Fat)
cor(data$BMI, data$Skeletal_Muscle)
cor(data$BMI, data$Intermuscular_Fat)
##KPS ----
ggplot(data, aes(x = KPS, y = OS, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) 


## dose aorta----
### subset: AORTA - PULM
ggplot(data, aes(x = data$mean_dose_aorta, y = data$mean_dose_pulm, color = data$discontinued_radio)) +
  geom_point(alpha = 0.7, size = 2)  #strong lin




ggplot(data, aes(x = data$mean_dose_aorta, y = data$eud_pulm, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) # s shape
ggplot(data, aes(x = data$eud_aorta, y = data$eud_pulm, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2)
ggplot(data, aes(x = data$eud_aorta, y = data$mean_dose_pulm, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2)




## lad ---- 
### NO RELEVANT FINDINGS
ggplot(data, aes(x = data$mean_dose_lad, y = data$mean_dose_pulm, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2)
ggplot(data, aes(x = sqrt(data$mean_dose_lad), y = data$mean_dose_pulm, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2)
cor(sqrt(data$mean_dose_lad), data$mean_dose_pulm)
ggplot(data, aes(x = data$mean_dose_lad, y = data$eud_ventricleL, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) #ventaglio
ggplot(data, aes(x = data$mean_dose_lad, y = data$eud_ventricleR, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) #ventaglio

## pulm ----
### NO RELEVANT FINDINGS
ggplot(data, aes(x = data$mean_dose_pulm, y = data$eud_atriumL, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) #s shape 

## atriumL ----
### NO RELEVANT FINDINGS
ggplot(data, aes(x = data$eud_atriumL, y = data$eud_aorta, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) 

ggplot(data, aes(x = data$eud_atriumL^2, y = data$eud_lad, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) 
cor(data$eud_atriumL,data$eud_lad)
cor(data$eud_atriumL^2,data$eud_lad)

ggplot(data, aes(x = data$eud_atriumL, y = data$mean_dose_pulm, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) 

##atriumR ----
### NO RELEVANT FINDINGS
ggplot(data, aes(x = data$mean_dose_atriumR, y = data$mean_dose_ventricleR, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) #ventaglio 
ggplot(data, aes(x = data$mean_dose_atriumR, y = data$mean_dose_venacava_I, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2)
ggplot(data, aes(x = data$mean_dose_atriumR, y = data$eud_venacavaI, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) #lin but many outliers

## ventricleL ----
ggplot(data, aes(x = data$mean_dose_ventricleL, y = data$mean_dose_lad, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) #ventaglio
ggplot(data, aes(x = data$mean_dose_ventricleL, y = data$mean_dose_ventricleR, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) #ventaglio but more linear
ggplot(data, aes(x = sqrt(data$mean_dose_ventricleL), y = data$mean_dose_Heart, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) #strong lin
cor(sqrt(data$mean_dose_ventricleL), data$mean_dose_Heart)

## ventricleR ----
ggplot(data, aes(x = sqrt(data$mean_dose_ventricleR), y = data$mean_dose_Heart, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2)
cor(data$mean_dose_ventricleR,data$mean_dose_Heart)
cor(sqrt(data$mean_dose_ventricleR),data$mean_dose_Heart) #higher 

## lungs ----
### NO RELEVANT FINDINGS (just related to all other features)

ggplot(data, aes(x = data$mean_dose_lungs, y = data$mean_dose_aorta, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) #weak lin
ggplot(data, aes(x = data$mean_dose_lungs, y = data$mean_dose_atriumR, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) #ventaglio
ggplot(data, aes(x = data$mean_dose_lungs, y = data$mean_dose_Heart, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) # weak lin

## fat ----
### NO RELEVANT FINDINGS
ggplot(data, aes(x = data$Subcutaneous_Fat, y = data$Torso_Fat, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2)
ggplot(data, aes(x = data$Subcutaneous_Fat, y = data$Skeletal_Muscle, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2)
ggplot(data, aes(x = data$Subcutaneous_Fat, y = data$Intermuscular_Fat, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2)
ggplot(data, aes(x = data$Skeletal_Muscle, y = data$Intermuscular_Fat, color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2)

## agatston-volume ----
### NO RELEVANT FINDINGS
ggplot(data, aes(x = data$heart_agatston, y = data$aorta_agatston )) +
  geom_point(alpha = 0.7, size = 2, color = "darkred") +   labs(x = "Heart agatston",
                                             y = "Aorta agatston",
                                             color = "") 
ggplot(data, aes(x = data$heart_volume_mm3, y = data$aorta_volume_mm3, color = Stage.at.RT)) +
  geom_point(alpha = 0.7, size = 2) +   labs(x = "Heart volume mm3",
                                                y = "Aorta volume mm3",
                                                color = "") 
ggplot(data, aes(x = data$aorta_agatston, y = data$aorta_volume_mm3, color = Stage.at.RT)) +
  geom_point(alpha = 0.7, size = 2) +   labs(x = "Aorta agatston",
                                             y = "Aorta volume mm3",
                                             color = "") 
ggplot(data, aes(x = data$heart_agatston, y = data$heart_volume_mm3, color = Stage.at.RT)) +
  geom_point(alpha = 0.7, size = 2) +   labs(x = "Heart agatston",
                                             y = "Heart volume mm3",
                                             color = "") 
par(mfrow=c(1,3))
ggplot(data, aes(x = data$aorta_agatston, y = data$diameter.aortic.arc.1..mm., color = Stage.at.RT)) +
  geom_point(alpha = 0.7, size = 2, color="darkred") +   labs(x = "Aorta agatston",
                                             y = "Aortic diameter arc 1 [mm]",
                                             color = "") 
ggplot(data, aes(x = data$aorta_agatston, y = data$diameter.aortic.arc.2..mm., color = Stage.at.RT)) +
  geom_point(alpha = 0.7, size = 2, color="darkred") +   labs(x = "Aorta agatston",
                                             y = "Aortic diameter arc 2 [mm]",
                                             color = "") 
ggplot(data, aes(x = data$aorta_agatston, y = data$diameter.kidney.sep..mm., color = Stage.at.RT)) +
  geom_point(alpha = 0.7, size = 2, color="darkred") +   labs(x = "Aorta agatston",
                                             y = "Aortic diameter kidney sep [mm]",
                                             color = "")

## diameters ----
ggplot(data, aes(x = data$diameter.aortic.arc.1..mm., y = data$diameter.aortic.arc.2..mm., color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) #lin
ggplot(data, aes(x = data$diameter.aortic.arc.1..mm., y = data$diameter.kidney.sep..mm., color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) #lin
ggplot(data, aes(x = data$diameter.aortic.arc.2..mm., y = data$diameter.kidney.sep..mm., color = OS2ys)) +
  geom_point(alpha = 0.7, size = 2) #lin 



## hypertension and diuretics  categorical plots ----
barplot_perc <- function(x, col = "lightblue") {
  tab <- table(x)
  
  # Calcoliamo un tetto massimo per l'asse Y (il massimo della tabella + 20% di spazio)
  tetto_y <- max(tab) * 1.2 
  
  # Aggiungiamo ylim al barplot
  bp <- barplot(tab, col = col, ylim = c(0, tetto_y))
  
  # Aggiungiamo il testo (con as.numeric per pulire l'output)
  text(x = bp, y = tab,
       labels = paste0(as.numeric(tab), "\n", round(100 * as.numeric(prop.table(tab)), 1), "%"),
       pos = 3)
}

barplot_perc(datacrfnew$Sex)
barplot_perc(data$Hypertension)
library(ggplot2)
ggplot(data, aes(x = Hypertension, fill = Diuretics)) +
  geom_bar(position = "fill", color = "black") +
  theme_minimal() +
  labs(title = "Uso di Diuretici nei pazienti Ipertesi vs Non Ipertesi",
       x = "Ipertensione (Sì/No)",
       y = "Proporzione",
       fill = "Uso Diuretici") +
  scale_y_continuous(labels = scales::percent)


ggplot(data, aes(x = Hypertension, fill = Beta.Blockers)) +
  geom_bar(position = "fill", color = "black") +
  theme_minimal() +
  labs(title = "Uso di Diuretici nei pazienti Ipertesi vs Non Ipertesi",
       x = "Ipertensione (Sì/No)",
       y = "Proporzione",
       fill = "Uso Diuretici") +
  scale_y_continuous(labels = scales::percent)

ggplot(data, aes(x = Diuretics, y = OS, fill = Diuretics)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.2, fill = "white", color = "black") +
  theme_minimal() +
  labs(title = "Distribuzione dei Mesi di Sopravvivenza per Uso di Diuretici",
       x = "Uso Diuretici",
       y = "Sopravvivenza (Mesi)") +
  theme(legend.position = "none")

ggplot(data, aes(x = Hypertension, y = Torso_Fat, fill = Hypertension)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Grasso del Torso in pazienti Ipertesi vs Non Ipertesi",
       x = "Ipertensione",
       y = "Volume Torso Fat (mm³)") +
  theme(legend.position = "none")











# TARGET DEPENDENCIES -----
boxplot(OS~data$Sex,data)
boxplot(OS_new~data$Sex,data)
boxplot(OS~data$Diabetes.drugs,data)
boxplot(OS_new~data$Diabetes.drugs,data)
boxplot(OS~data$Antilipidic.drugs,data)
boxplot(OS~Hypertension,data)
boxplot(OS~data$antithrombotic,data)
boxplot(OS~data$Beta.Blockers,data)
boxplot(OS~data$Diuretics,data)
boxplot(OS~data$Smoking.status,data)
boxplot(OS~data$KPS,data)
boxplot(OS~data$Stage.at.RT,data)
boxplot(OS~data$PDL.1.level,data)
boxplot(OS~data$CHEMO,data)
boxplot(OS~data$schedule_ok_chemo,data)
boxplot(OS~data$cycles_chemo_sum,data)
boxplot(OS~data$discontinued_chemo,data)
boxplot(OS~data$IMMUNO,data)
boxplot(OS~data$cycles_immuno,data)
boxplot(OS~data$discontinued_immuno,data)
boxplot(OS~data$prescription_dose_radio,data)
boxplot(OS~data$number_of_fractions_radio,data)
boxplot(OS~data$eqd2_dose_radio,data)
boxplot(OS~data$discontinued_radio,data)

boxplot(data$Age~data$OS..yes.no.,data)
boxplot(data$BMI~data$OS..yes.no.,data)

boxplot(data$mean_dose_aorta~data$OS..yes.no.,data)
boxplot(data$mean_dose_lad~data$OS..yes.no.,data)
boxplot(data$mean_dose_pulm~data$OS..yes.no.,data)
boxplot(data$mean_dose_atriumL~data$OS2ys,data)
boxplot(data$mean_dose_atriumR~data$OS..yes.no.,data)
boxplot(data$mean_dose_ventricleL~data$OS..yes.no.,data)
boxplot(data$mean_dose_ventricleR~data$OS..yes.no.,data)
boxplot(data$mean_dose_Heart~data$OS..yes.no.,data)
boxplot(data$mean_dose_venacava_I~data$OS..yes.no.,data)
boxplot(data$mean_dose_venacava_S~data$OS..yes.no.,data)
boxplot(data$mean_dose_lungs~data$OS..yes.no.,data)

boxplot(data$A_Aorta~data$OS..yes.no.,data)
boxplot(data$A_Lad~data$OS..yes.no.,data)
boxplot(data$A_Pulm~data$OS..yes.no.,data)
boxplot(data$Atrium_L~data$OS..yes.no.,data)
boxplot(data$Atrium_R~data$OS..yes.no.,data)
boxplot(data$Heart~data$OS..yes.no.,data)
boxplot(data$Lung_L~data$OS..yes.no.,data)
boxplot(data$Lungs~data$OS..yes.no.,data)
boxplot(data$Ventricle_L~data$OS..yes.no.,data)
boxplot(data$Ventricle_R~data$OS..yes.no.,data)
boxplot(data$Venacava_I~data$OS..yes.no.,data)
boxplot(data$Venacava_S~data$OS..yes.no.,data)

boxplot(data$Subcutaneous_Fat~data$OS..yes.no.,data)
boxplot(data$Torso_Fat~data$OS..yes.no.,data)
boxplot(data$Skeletal_Muscle~data$OS..yes.no.,data)
boxplot(data$Intermuscular_Fat~data$OS..yes.no.,data)
boxplot(data$heart_agatston~data$OS..yes.no.,data)
boxplot(data$heart_volume_mm3~data$OS..yes.no.,data)
boxplot(data$aorta_agatston~data$OS..yes.no.,data)
boxplot(data$aorta_volume_mm3~data$OS..yes.no.,data)
boxplot(data$LAA~data$OS..yes.no.,data)
boxplot(data$Volume_PTV~data$OS..yes.no.,data)
boxplot(data$diameter.aortic.arc.1..mm.~data$OS..yes.no.,data)
boxplot(data$diameter.aortic.arc.2..mm.~data$OS..yes.no.,data)
boxplot(data$diameter.kidney.sep..mm.~data$OS..yes.no.,data)






# HOSPITAL ANALYSIS ----
data$hospital <- factor(substr(as.character(data$ID), 3, 3),
                        levels = as.character(1:7))
barplot(table(data$Sex, data$hospital),col = terrain.colors(length(unique(data$Sex))), legend.text = TRUE)
barplot(table(data$Diabetes.drugs, data$hospital),col = terrain.colors(length(unique(data$Diabetes.drugs))), legend.text = TRUE)
barplot(table(data$Hypertension, data$hospital),col = terrain.colors(length(unique(data$Hypertension))), legend.text = TRUE)
barplot(table(data$Diuretics, data$hospital),col = terrain.colors(length(unique(data$Diuretics))), legend.text = TRUE)
barplot(table(data$Smoking.status, data$hospital),col = terrain.colors(length(unique(data$Smoking.status))), legend.text = TRUE)
barplot(table(data$KPS, data$hospital),col = terrain.colors(length(unique(data$KPS))), legend.text = TRUE) #PUT
barplot(table(data$Stage.at.RT, data$hospital),col = terrain.colors(length(unique(data$Stage.at.RT))), legend.text = TRUE) #PUT
barplot(table(data$PDL.1.level, data$hospital),col = terrain.colors(length(unique(data$PDL.1.level))), legend.text = TRUE) #PUT WITH IMMNUO
barplot(table(data$CHEMO, data$hospital),col = terrain.colors(length(unique(data$CHEMO))), legend.text = TRUE)
barplot(table(data$schedule_ok_chemo, data$hospital),col = terrain.colors(length(unique(data$schedule_ok_chemo))), legend.text = TRUE) #PUT CONTAINS ALSO INFO OF CHEMO
barplot(table(data$cycles_chemo_sum, data$hospital),col = terrain.colors(length(unique(data$cycles_chemo_sum))), legend.text = TRUE)
barplot(table(data$discontinued_chemo, data$hospital),col = terrain.colors(length(unique(data$discontinued_chemo))), legend.text = TRUE)
barplot(table(data$IMMUNO, data$hospital),col = terrain.colors(length(unique(data$IMMUNO))), legend.text = TRUE) #3 DIDN'T DO IMMUNO
barplot(table(data$discontinued_immuno, data$hospital),col = terrain.colors(length(unique(data$discontinued_immuno))), legend.text = TRUE)
barplot(table(data$discontinued_radio, data$hospital),col = terrain.colors(length(unique(data$discontinued_radio))), legend.text = TRUE)
barplot(table(data$OS2ys, data$hospital),col = terrain.colors(length(unique(data$OS2ys))), legend.text = TRUE)


boxplot(data$Age~data$hospital,data)
boxplot(data$BMI~data$hospital,data)
boxplot(data$Volume_PTV~data$hospital,data)

boxplot(data$mean_dose_aorta~data$hospital,data)
boxplot(data$mean_dose_lad~data$hospital,data)
boxplot(data$mean_dose_pulm~data$hospital,data)
boxplot(data$mean_dose_atriumL~data$hospital,data)
boxplot(data$mean_dose_atriumR~data$hospital,data)
boxplot(data$mean_dose_ventricleL~data$hospital,data)
boxplot(data$mean_dose_ventricleR~data$hospital,data)
boxplot(data$mean_dose_Heart~data$hospital,data)
boxplot(data$mean_dose_venacava_I~data$hospital,data)
boxplot(data$mean_dose_venacava_S~data$hospital,data)
boxplot(data$mean_dose_lungs~data$hospital,data)

boxplot(data$A_Aorta~data$hospital,data)
boxplot(data$A_Lad~data$hospital,data)
boxplot(data$A_Pulm~data$hospital,data)
boxplot(data$Atrium_L~data$hospital,data)
boxplot(data$Atrium_R~data$hospital,data)
boxplot(data$Heart~data$hospital,data)
boxplot(data$Lung_L~data$hospital,data)
boxplot(data$Lungs~data$hospital,data)
boxplot(data$Ventricle_L~data$hospital,data)
boxplot(data$Ventricle_R~data$hospital,data)
boxplot(data$Venacava_I~data$hospital,data)
boxplot(data$Venacava_S~data$hospital,data)

boxplot(data$Subcutaneous_Fat~data$hospital,data)
boxplot(data$Torso_Fat~data$hospital,data)
boxplot(data$Skeletal_Muscle~data$hospital,data)
boxplot(data$Intermuscular_Fat~data$hospital,data)
boxplot(data$heart_agatston~data$hospital,data)
boxplot(data$heart_volume_mm3~data$hospital,data)
boxplot(data$aorta_agatston~data$hospital,data)
boxplot(data$aorta_volume_mm3~data$hospital,data)
boxplot(data$LAA~data$hospital,data)
boxplot(data$Volume_PTV~data$hospital,data) #only this
boxplot(data$diameter.aortic.arc.1..mm.~data$hospital,data)
boxplot(data$diameter.aortic.arc.2..mm.~data$hospital,data)
boxplot(data$diameter.kidney.sep..mm.~data$hospital,data)


boxplot(data$OS~data$hospital,data)
boxplot(data$OS_new~data$hospital,data)
barplot(table(data$OS2ys, data$hospital),col = terrain.colors(length(unique(data$OS2ys))), legend.text = TRUE) #3 SHOWS BETTER OS
barplot(table(data$OS..yes.no., data$hospital),col = terrain.colors(length(unique(data$OS..yes.no.))), legend.text = TRUE)

barplot(table(data$KPS, data$hospital),col = terrain.colors(length(unique(data$KPS))), legend.text = TRUE) #PUT

barplot(table(data$Stage.at.RT, data$hospital),col = terrain.colors(length(unique(data$Stage.at.RT))), legend.text = TRUE) 

barplot(table(data$PDL.1.level, data$hospital),col = terrain.colors(length(unique(data$PDL.1.level))), legend.text = TRUE)
barplot(table(data$IMMUNO, data$hospital),col = terrain.colors(length(unique(data$IMMUNO))), legend.text = TRUE)
barplot(table(data$schedule_ok_chemo, data$hospital),col = terrain.colors(length(unique(data$schedule_ok_chemo))), legend.text = TRUE) #PUT CONTAINS ALSO INFO OF CHEMO

boxplot(data$Volume_PTV~data$hospital,data) 

#PLOT REPORT ----
ggplot(data, aes(x = eud_aorta,
                 y = eud_pulm,
                 color = discontinued_radio)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(x = "Eud to aorta",
       y = "Eud to pulm",
       color = "Radiotherapy discontinued") 
ggplot(data, aes(x = eud_Heart,
                 y = eud_lungs,
                 color = discontinued_radio)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(x = "Eud to heart",
       y = "Eud to lungs",
       color = "Radiotherapy discontinued") 
ggplot(data, aes(x = data$heart_agatston, y = data$aorta_agatston )) +
  geom_point(alpha = 0.7, size = 2, color = "darkred") +   labs(x = "Heart agatston",
                                                                y = "Aorta agatston",
                                                                color = "") 
ggplot(data, aes(x = data$aorta_agatston, y = data$diameter.aortic.arc.1..mm., color = Stage.at.RT)) +
  geom_point(alpha = 0.7, size = 2, color="darkred") +   labs(x = "Aorta agatston",
                                                              y = "Aortic diameter arc 1 [mm]",
                                                              color = "") 
ggplot(data, aes(x = data$aorta_agatston, y = data$diameter.aortic.arc.2..mm., color = Stage.at.RT)) +
  geom_point(alpha = 0.7, size = 2, color="darkred") +   labs(x = "Aorta agatston",
                                                              y = "Aortic diameter arc 2 [mm]",
                                                              color = "") 
ggplot(data, aes(x = data$aorta_agatston, y = data$diameter.kidney.sep..mm., color = Stage.at.RT)) +
  geom_point(alpha = 0.7, size = 2, color="darkred") +   labs(x = "Aorta agatston",
                                                              y = "Aortic diameter kidney sep [mm]",
                                                              color = "")
ggpairs(datametab[, 1:4], lower = list(continuous = wrap("points", alpha = 0.5, color = "#2C7FB8")), diag = list(continuous = wrap("densityDiag", fill = "#2C7FB8", alpha = 0.5)))

## target - osdep ---------
par(mfrow=c(1,3))
boxplot(OS_new ~ Sex,
        data = data,
        col = c("pink", "#2C7FB8"),
        border = "grey30",
        xlab = "Sex",
        ylab = "Overall Survival new (months)")
tmp <- subset(data, Smoking.status != "UNKNOWN")
tmp$Smoking.status <- droplevels(tmp$Smoking.status)
boxplot(OS_new ~ Smoking.status,
        data = tmp,
        col = c("#FDAE61","#1A9850", "#D73027"),
        xlab = "Smoking status",
        ylab = "Overall Survival new (months)")

boxplot(OS_new ~ data$KPS,
        data = data,
        col = c(
          "#D73027",  # rosso
          "#FDAE61",  # arancione
          "#FEE08B",  # giallo
          "#66BD63",  # verde chiaro
          "#1A9850"   # verde
        ),
        border = "grey30",
        xlab = "KPS",
        ylab = "Overall Survival new (months)")

boxplot(OS_new ~ data$Diabetes.drugs,
        data = data,
        col = c("#2C7FB8", "#D73027"),
        border = "grey30",
        xlab = "Diabetes drugs",
        ylab = "Overall Survival new (months)")
boxplot(OS_new ~ data$Hypertension,
        data = data,
        col = c("#2C7FB8", "#D73027"),
        border = "grey30",
        xlab = "Hypertension",
        ylab = "Overall Survival new (months)")

tmp <- subset(data, Diuretics != "UNKNOWN")
tmp$Diuretics <- droplevels(tmp$Diuretics)
boxplot(OS_new ~ Diuretics,
        data = tmp,
        col = c("#2C7FB8", "#D73027"),
        xlab = "Diuretics",
        ylab = "Overall Survival new (months)")

par(mfrow=c(1,2))
boxplot(OS_new ~ data$Stage.at.RT,
        data = data,
        col = c(
          "#FEE08B",  # giallo
          "#FDAE61",  # arancione
          "#D73027"   # rosso
        ),
        border = "grey30",
        xlab = "Stage at radiotherapy",
        ylab = "Overall Survival new (months)")
#levels(data$PDL.1.level) <- c("-1","0-1%", "1-49%", "50-100%")
tmp <- subset(data, PDL.1.level != -1)
tmp$PDL.1.level <- droplevels(tmp$PDL.1.level)
boxplot(OS_new ~ tmp$PDL.1.level,
        data = tmp,
        col = c(
          "#FDAE61",  # arancione
          "#FEE08B",  # giallo
          "#66BD63"),
        border = "grey30",
        xlab = "PDL-1 level",
        ylab = "Overall Survival new (months)")

par(mfrow=c(1,2))
boxplot(OS_new ~ data$CHEMO,
        data = data,
        col = c(
          "#FEE08B",  # giallo
          "#FDAE61"),
        border = "grey30",
        xlab = "CHEMO",
        ylab = "Overall Survival new (months)")
boxplot(OS_new ~ data$IMMUNO,
        data = data,
        col = c(
          "#FEE08B",  # giallo
          "#FDAE61"),
        border = "grey30",
        xlab = "IMMUNO",
        ylab = "Overall Survival new (months)")
par(mfrow=c(1,3))
tmp <- subset(data, discontinued_chemo != "-1" & discontinued_chemo != "0.0")
tmp$discontinued_chemo <- droplevels(tmp$discontinued_chemo)
boxplot(OS_new ~ tmp$discontinued_chemo,
        data = tmp,
        col = c("#2C7FB8", "#D73027"),
        border = "grey30",
        xlab = "Discontinued CHEMO",
        ylab = "Overall Survival new (months)")
tmp <- subset(data, discontinued_immuno != "-1" & discontinued_immuno != "0.0")
tmp$discontinued_immuno <- droplevels(tmp$discontinued_immuno)
boxplot(OS_new ~ tmp$discontinued_immuno,
        data = tmp,
        col = c("#2C7FB8", "#D73027"),
        border = "grey30",
        xlab = "Discontinued IMMUNO",
        ylab = "Overall Survival new (months)")
tmp <- subset(data, discontinued_radio != "-1" & discontinued_radio != "0.0")
tmp$discontinued_radio <- droplevels(tmp$discontinued_radio)
boxplot(OS_new ~ tmp$discontinued_radio,
        data = tmp,
        col = c("#2C7FB8", "#D73027"),
        border = "grey30",
        xlab = "Discontinued RADIO",
        ylab = "Overall Survival new (months)")


par(mfrow=c(1,2))
boxplot(Age ~ data$OS..yes.no.,
        data = data,
        col = c(
          "#D73027",  
          "#66BD63"  # verde chiaro
        ),
        border = "grey30",
        xlab = "OS (yes/no)",
        ylab = "Age")
boxplot(BMI ~ data$OS..yes.no.,
        data = data,
        col = c(
          "#D73027",  
          "#66BD63"  # verde chiaro
        ),
        border = "grey30",
        xlab = "OS (yes/no)",
        ylab = "BMI")

par(mfrow=c(1,4))
boxplot(data$mean_dose_aorta ~ data$OS..yes.no.,
        data = data,
        col = c(
          "#D73027",  
          "#66BD63"  # verde chiaro
        ),
        border = "grey30",
        xlab = "OS (yes/no)",
        ylab = "Mean dose to aorta")
boxplot(data$mean_dose_lad ~ data$OS..yes.no.,
        data = data,
        col = c(
          "#D73027",  
          "#66BD63"  # verde chiaro
        ),
        border = "grey30",
        xlab = "OS (yes/no)",
        ylab = "Mean dose to lad")
boxplot(data$mean_dose_pulm ~ data$OS..yes.no.,
        data = data,
        col = c(
          "#D73027",  
          "#66BD63"  # verde chiaro
        ),
        border = "grey30",
        xlab = "OS (yes/no)",
        ylab = "Mean dose to pulm")
boxplot(data$mean_dose_lungs ~ data$OS..yes.no.,
        data = data,
        col = c(
          "#D73027",  
          "#66BD63"  # verde chiaro
        ),
        border = "grey30",
        xlab = "OS (yes/no)",
        ylab = "Mean dose to lungs")
par(mfrow=c(1,3))
boxplot(data$A_Lad ~ data$OS..yes.no.,
        data = data,
        col = c(
          "#D73027",  
          "#66BD63"  # verde chiaro
        ),
        border = "grey30",
        xlab = "OS (yes/no)",
        ylab = "Volume Lad [cm3]")
boxplot(data$Ventricle_R ~ data$OS..yes.no.,
        data = data,
        col = c(
          "#D73027",  
          "#66BD63"  # verde chiaro
        ),
        border = "grey30",
        xlab = "OS (yes/no)",
        ylab = "Volume right Ventricle [cm3]")
boxplot(data$Heart ~ data$OS..yes.no.,
        data = data,
        col = c(
          "#D73027",  
          "#66BD63"  # verde chiaro
        ),
        border = "grey30",
        xlab = "OS (yes/no)",
        ylab = "Volume Heart [cm3]")


#IQR is used for visualization
par(mfrow=c(2,3)) 
Q1 <- quantile(data$Subcutaneous_Fat, 0.25, na.rm = TRUE)
Q3 <- quantile(data$Subcutaneous_Fat, 0.75, na.rm = TRUE)
IQRv <- Q3 - Q1

tmp <- data[
  data$Subcutaneous_Fat >= (Q1 - 1.5 * IQRv) &
    data$Subcutaneous_Fat <= (Q3 + 1.5 * IQRv),
]
boxplot(tmp$Subcutaneous_Fat ~ tmp$OS..yes.no.,
        data = tmp,
        col = c(
          "#D73027",  
          "#66BD63"  # verde chiaro
        ),
        border = "grey30",
        xlab = "OS (yes/no)",
        ylab = "Subcutaneous Fat")


Q1 <- quantile(data$Torso_Fat, 0.25, na.rm = TRUE)
Q3 <- quantile(data$Torso_Fat, 0.75, na.rm = TRUE)
IQRv <- Q3 - Q1

tmp <- data[
  data$Torso_Fat >= (Q1 - 1.5 * IQRv) &
    data$Torso_Fat <= (Q3 + 1.5 * IQRv),
]
boxplot(tmp$Torso_Fat ~ tmp$OS..yes.no.,
        data = tmp,
        col = c(
          "#D73027",  
          "#66BD63"  # verde chiaro
        ),
        border = "grey30",
        xlab = "OS (yes/no)",
        ylab = "Torso Fat")
Q1 <- quantile(data$heart_agatston, 0.25, na.rm = TRUE)
Q3 <- quantile(data$heart_agatston, 0.75, na.rm = TRUE)
IQRv <- Q3 - Q1

tmp <- data[
  data$heart_agatston >= (Q1 - 1.5 * IQRv) &
    data$heart_agatston <= (Q3 + 1.5 * IQRv),
]

boxplot(tmp$heart_agatston ~ tmp$OS..yes.no.,
        data = tmp,
        col = c(
          "#D73027",  
          "#66BD63"  # verde chiaro
        ),
        border = "grey30",
        xlab = "OS (yes/no)",
        ylab = "Heart agatston")

Q1 <- quantile(data$aorta_agatston, 0.25, na.rm = TRUE)
Q3 <- quantile(data$aorta_agatston, 0.75, na.rm = TRUE)
IQRv <- Q3 - Q1

tmp <- data[
  data$aorta_agatston >= (Q1 - 1.5 * IQRv) &
    data$aorta_agatston <= (Q3 + 1.5 * IQRv),
]

boxplot(tmp$aorta_agatston ~ tmp$OS..yes.no.,
        data = tmp,
        col = c(
          "#D73027",  
          "#66BD63"  # verde chiaro
        ),
        border = "grey30",
        xlab = "OS (yes/no)",
        ylab = "Aorta agatston")

Q1 <- quantile(data$LAA, 0.25, na.rm = TRUE)
Q3 <- quantile(data$LAA, 0.75, na.rm = TRUE)
IQRv <- Q3 - Q1

tmp <- data[
  data$LAA >= (Q1 - 1.5 * IQRv) &
    data$LAA <= (Q3 + 1.5 * IQRv),
]

boxplot(tmp$LAA ~ tmp$OS..yes.no.,
        data = tmp,
        col = c(
          "#D73027",  
          "#66BD63"  # verde chiaro
        ),
        border = "grey30",
        xlab = "OS (yes/no)",
        ylab = "LAA")

Q1 <- quantile(data$Volume_PTV, 0.25, na.rm = TRUE)
Q3 <- quantile(data$Volume_PTV, 0.75, na.rm = TRUE)
IQRv <- Q3 - Q1

tmp <- data[
  data$Volume_PTV >= (Q1 - 1.5 * IQRv) &
    data$Volume_PTV <= (Q3 + 1.5 * IQRv),
]

boxplot(tmp$Volume_PTV ~ tmp$OS..yes.no.,
        data = tmp,
        col = c(
          "#D73027",  
          "#66BD63"  # verde chiaro
        ),
        border = "grey30",
        xlab = "OS (yes/no)",
        ylab = "Volume PTV")


## structure an -----
par(mfrow=c(1,1))

data$hospital <- factor(substr(as.character(data$ID), 3, 3),
                        levels = as.character(1:7))

boxplot(Age ~ hospital,
        data = data,
        col = "#2C7FB8",
        border = "grey30",
        outline = FALSE,
        main = "Age distribution",
        xlab = "Hospital",
        ylab = "Age")

stage_cols <- c("IIIA" = "#FEE08B",
                "IIIB" = "#FDAE61",
                "IIIC" = "#D73027")

tab_stage <- table(data$hospital, data$Stage.at.RT)
prop_stage <- prop.table(tab_stage, margin = 1)

max_y_stage <- max(prop_stage, na.rm = TRUE) * 1.35

bp_stage <- barplot(t(prop_stage),
                    beside = TRUE,
                    col = stage_cols[colnames(prop_stage)],
                    border = "grey30",
                    ylim = c(0, max_y_stage),
                    legend.text = colnames(prop_stage),
                    args.legend = list(x = "topright",
                                       cex = 0.9,
                                       bty = "o",
                                       box.col = "grey30",
                                       bg = rgb(1, 1, 1, alpha = 0.7)),
                    xlab = "Hospital",
                    ylab = "Proportion",
                    main = "Stage distribution")

box()

labels_stage <- ifelse(t(prop_stage) > 0,
                       paste0(round(t(prop_stage) * 100, 0), "%"),
                       "")

text(x = bp_stage,
     y = t(prop_stage),
     labels = labels_stage,
     srt=90,
     adj = c(-0.5, 0.5),
     cex = 0.7,
     font = 2,
     col = "black")

kps_levels <- sort(unique(data$KPS))

kps_cols <- colorRampPalette(
  c("#D73027", "#F46D43", "#FDAE61", "#A6D96A", "#1A9850")
)(length(kps_levels))

names(kps_cols) <- kps_levels

tab_kps <- table(data$hospital, data$KPS)
prop_kps <- prop.table(tab_kps, margin = 1)

max_y_kps <- max(prop_kps, na.rm = TRUE) * 1.35

bp_kps <- barplot(t(prop_kps),
                  beside = TRUE,
                  col = kps_cols[colnames(prop_kps)],
                  border = "grey30",
                  ylim = c(0, max_y_kps),
                  legend.text = colnames(prop_kps),
                  args.legend = list(x = "topright",
                                     cex = 0.9,
                                     bty = "o",
                                     box.col = "grey30",
                                     bg = rgb(1, 1, 1, alpha = 0.7)),
                  xlab = "Hospital",
                  ylab = "Proportion",
                  main = "KPS distribution")

box()

labels_kps <- ifelse(t(prop_kps) > 0,
                     paste0(round(t(prop_kps) * 100, 0), "%"),
                     "")

offset_y <- max_y_kps * 0.02

text(x = bp_kps,
     y = t(prop_kps) + offset_y,
     labels = labels_kps,
     srt = 90,
     adj = c(0, 0.5),
     cex = 0.7,
     font = 2,
     col = "black")











