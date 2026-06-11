# =============================
# Gut Microbial Metabolites and Lipidomics - MCI_to_AD Pipeline 
# Author: Alejandro I. Trejo-Castro 
# Last version 11.11.25
# Goal: 
# 1. Make the cohort construction explicit and reproducible
# 2. Prevent data leakage
# =============================

# 0) Reproducibility & Setup
set.seed(0)
options(stringsAsFactors = FALSE)

# Helper to load or stop with a clear error
load_pkg <- function(pkgs) {
  for (p in pkgs) {
    if (!require(p, character.only = TRUE)) {
      stop(sprintf("Package '%s' is required but not installed.", p))
    }
  }
}

load_pkg(c(
  "dplyr","stringi","tidyverse","DataExplorer",
  "FRESA.CAD","Rcpp","stringr","miscTools","Hmisc","pROC",
  "nlme","rpart","gplots","RColorBrewer","class","cvTools",
  "glmnet","randomForest","survival","e1071","MASS","naivebayes",
  "mRMRe","epiR","DescTools","irr","survminer","BeSS","ggplot2",
  "robustbase","mda","twosamples","Rfast","whitening","corrplot","caret",
  "pheatmap","RColorBrewer"
))

# 1) Paths & Output directories
# Adjust 'root' to your environment 
root <- "G:/Mi unidad/Investigación Trejo-Castro/E31 - DBC Metabolic Health and Disease/E31.1 Alzheimer_GutMetabolites_Lipidomics"
setwd(root)

# Create output dirs if not present
out_db   <- file.path("3. Results Databases/1. Exp MCI_to_AD")
out_figs <- file.path("4. Results Plots/1. Exp MCI_to_AD")
dir.create(out_db, recursive = TRUE, showWarnings = FALSE)
dir.create(out_figs, recursive = TRUE, showWarnings = FALSE)

pdf(file = file.path(out_figs, "Exp_MCI5_to_AD2.pdf"), width = 8, height = 6)

# Utility
safe_log2 <- function(x) {
  # Handles zeros/negatives by returning NA; caller may impute if desired
  y <- suppressWarnings(log2(as.numeric(x)))
  y[!is.finite(y)] <- NA
  return(y)
}
make_flag <- function(RID, VISCODE2) paste(RID, VISCODE2)

# 2) Diagnostic & Demographic data 
dx <- read.csv(file.path("1. ADNI/Assessments","DXSUM_PDXCONV_18Mar2024.csv"))
dx[dx == -4] <- NA
dx <- replace(dx, dx == "", NA)
dx <- dplyr::select(dx, c("RID","VISCODE2","EXAMDATE","DIAGNOSIS"))
dx <- dx[order(dx$RID, dx$EXAMDATE),]

dem <- read.csv(file.path("1. ADNI/Subject Characteristics","PTDEMOG_19Mar2024.csv"))
dem[dem == -4] <- NA
dem <- replace(dem, dem == "", NA)
dem <- dplyr::select(dem, c("RID","PTGENDER","PTDOB","PTDOBYY","PTEDUCAT"))
dem$PTDOBMONTH <- substr(dem$PTDOB, 1, 2)
dem$DOB <- as.Date(paste(dem$PTDOBYY, dem$PTDOBMONTH, 15, sep = "-"))
dem <- dplyr::select(dem, c("RID","DOB","PTGENDER","PTEDUCAT"))

dx$DOB      <- dem$DOB[match(dx$RID, dem$RID)]
dx$PTGENDER <- dem$PTGENDER[match(dx$RID, dem$RID)]
dx$PTEDUCAT <- dem$PTEDUCAT[match(dx$RID, dem$RID)]

# Keep only subjects with MCI at baseline (VISCODE2 == "bl", DIAGNOSIS == 2)
mciBl <- dx[dx$RID %in% dx$RID[dx$VISCODE2=="bl" & dx$DIAGNOSIS==2],]
# Remove subjects with only one visit
mciBl <- mciBl[!(mciBl$RID %in% as.numeric(names(which(table(mciBl$RID)==1)))),]
# Age at each visit
mciBl$AGE <- round(as.numeric(as.Date(mciBl$EXAMDATE) - mciBl$DOB)/365, 1)
# Remove cases who revert to DX=1 after BL DX=2
mciBl <- mciBl[!(mciBl$RID %in% mciBl$RID[mciBl$DIAGNOSIS==1]),]
data <- mciBl

# 3) Define study population (MCI-stable vs MCI-progressors)
data$Classification <- ifelse(data$RID %in% data$RID[data$DIAGNOSIS==3], "MCI-p","MCI-s")

# Stable
mcidata <- data[data$Classification=="MCI-s",]
mcidata$LV <- mcidata$VISCODE2
mcidata$LV <- ifelse(mcidata$LV=="bl", 0, ifelse(mcidata$LV %in% c("m102","m108","m114","m120","m126","m132",
                                                                   "m138","m144","m150","m156","m174","m180","m186",
                                                                   "m192"),
                                                 as.numeric(stringi::stri_sub(mcidata$LV,-3,-1)),
                                                 as.numeric(stringi::stri_sub(mcidata$LV,-2,-1))))
mcidata$LV <- suppressWarnings(as.numeric(mcidata$LV))
mcidata$CLASS2 <- ifelse(mcidata$RID %in% mcidata$RID[mcidata$LV>=60], 1, 0)
mcidata <- mcidata[mcidata$CLASS2==1,]
write.csv(mcidata,
          file.path(out_db, "1. MCI_S_Population.csv"),
          row.names = FALSE)
mcidata <- mcidata[mcidata$VISCODE2=="bl",]
mcidata <- subset(mcidata, select = -c(LV,CLASS2))
write.csv(mcidata,
          file.path(out_db, "2. MCI_S_Population_BL_Info.csv"),
          row.names = FALSE)
# Progressors
addata <- data[data$Classification=="MCI-p",]
addata$LV <- addata$VISCODE2
addata$LV <- ifelse(addata$LV=="bl", 0,
                    ifelse(addata$LV %in% c("m102","m108","m114","m120","m126","m132",
                                            "m138","m144","m150","m156","m174","m180","m186",
                                            "m192"),
                           as.numeric(stringi::stri_sub(addata$LV,-3,-1)),
                           as.numeric(stringi::stri_sub(addata$LV,-2,-1))))
DX3 <- addata[addata$DIAGNOSIS==3,] 
DX3<-group_by(DX3,RID)
DX3<-slice(DX3,1)
DX3 <- DX3[DX3$LV<=24,]
addata$classDX3 <- ifelse(addata$RID %in% DX3$RID, 1, 0)
addata <- addata[addata$classDX3==1,]
write.csv(addata,
          file.path(out_db, "3. MCI_P_Population.csv"),
          row.names = FALSE)
addata <- addata[addata$VISCODE2=="bl",]
addata <- subset(addata, select = -c(LV,classDX3))
write.csv(addata,
          file.path(out_db, "4. MCI_P_Population_BL_Info.csv"),
          row.names = FALSE)
population <- rbind(addata, mcidata)
# NOTE: Age filter according to LOAD
population <- population[population$AGE>=64.5,]
population <- population[order(population$RID),]
write.csv(population,
          file.path(out_db, "5. Population_Age_Check.csv"),
          row.names = FALSE)


#  4) Neuropsychology & formatting 
GDS  <- read.csv(file.path("1. ADNI/Assessments","GDSCALE_19Mar2024.csv"))
GDS[GDS == -4] <- NA; GDS <- replace(GDS, GDS=="", NA); GDS <- replace(GDS, GDS=="sc", "bl")
GDS  <- dplyr::select(GDS, c("RID","VISCODE2","GDTOTAL"))

MMSE <- read.csv(file.path("1. ADNI/Assessments","MMSE_19Mar2024.csv"))
MMSE[MMSE == -4] <- NA; MMSE <- replace(MMSE, MMSE=="", NA); MMSE <- replace(MMSE, MMSE=="sc", "bl")
MMSE <- dplyr::select(MMSE, c("RID","VISCODE2","MMSCORE"))

UWN  <- read.csv(file.path("1. ADNI/Assessments","UWNPSYCHSUM_19Mar2024.csv"))
UWN[UWN == -4] <- NA; UWN <- replace(UWN, UWN=="", NA)
UWN  <- dplyr::select(UWN, c("RID","VISCODE2","ADNI_MEM","ADNI_EF","ADNI_LAN","ADNI_VS"))

ADAS1 <- read.csv(file.path("1. ADNI/Assessments","ADASSCORES_28Mar2024.csv"))
ADAS1[ADAS1 == -4] <- NA; ADAS1 <- replace(ADAS1, ADAS1=="", NA)
ADAS1 <- dplyr::select(ADAS1, c("RID","VISCODE","TOTAL11"))
ADAS2 <- read.csv(file.path("1. ADNI/Assessments","ADAS_ADNIGO23_28Mar2024.csv"))
ADAS2[ADAS2 == -4] <- NA; ADAS2 <- replace(ADAS2, ADAS2=="", NA)
ADAS2 <- dplyr::select(ADAS2, c("RID","VISCODE2","TOTSCORE"))
colnames(ADAS2) <- c("RID","VISCODE","TOTAL11")
ADAS <- rbind(ADAS1, ADAS2)

population$flag <- make_flag(population$RID, population$VISCODE2)
GDS$flag  <- make_flag(GDS$RID, GDS$VISCODE2)
MMSE$flag <- make_flag(MMSE$RID, MMSE$VISCODE2)
UWN$flag  <- make_flag(UWN$RID, UWN$VISCODE2)
ADAS$flag <- make_flag(ADAS$RID, ADAS$VISCODE)

population$MMSCORE  <- MMSE$MMSCORE[match(population$flag, MMSE$flag)]
population$GDTOTAL  <- GDS$GDTOTAL[match(population$flag, GDS$flag)]
population$ADAS11   <- ADAS$TOTAL11[match(population$flag, ADAS$flag)]
population$ADNI_MEM <- UWN$ADNI_MEM[match(population$flag, UWN$flag)]
population$ADNI_EF  <- UWN$ADNI_EF[match(population$flag, UWN$flag)]
population$ADNI_LAN <- UWN$ADNI_LAN[match(population$flag, UWN$flag)]
population$ADNI_VS  <- UWN$ADNI_VS[match(population$flag, UWN$flag)]

population <- subset(population, select = -c(EXAMDATE,DIAGNOSIS,DOB))
population <- population %>% dplyr::select(Classification,RID,VISCODE2,PTGENDER,PTEDUCAT,AGE, dplyr::everything())
write.csv(population,
          file.path(out_db, "6. Population_Neuropsychological.csv"),
          row.names = FALSE)


# 5) Join BMI, APOE, CSF, Lipidomics, Gut metabolites 
vitalsigns <- read.csv("1. ADNI/Vital Signs/VITALS_19Mar2024.csv", stringsAsFactors = FALSE)
vitalsigns <- subset(vitalsigns, select=c(RID, VISCODE2, VSWEIGHT, VSWTUNIT, VSHEIGHT, VSHTUNIT, VSBPSYS, VSBPDIA))
vitalsigns[vitalsigns == -4] <- NA
vitalsigns <- replace(vitalsigns, vitalsigns == "", NA)
# VSWTUNIT: 1 = lb, 2 = kg
# VSHTUNIT: 1 = in, 2 = cm
vitalsigns <- vitalsigns[vitalsigns$VSWTUNIT %in% c(1,2),]
vitalsigns[vitalsigns == -1] <- NA

# Convert to kg/cm
vitalsigns$VSWEIGHT <- ifelse(vitalsigns$VSWTUNIT==1, vitalsigns$VSWEIGHT*0.45359237, vitalsigns$VSWEIGHT)
vitalsigns$VSHEIGHT <- as.numeric(vitalsigns$VSHEIGHT)
vitalsigns$VSHEIGHT <- ifelse(vitalsigns$VSHTUNIT==1, vitalsigns$VSHEIGHT*2.54, vitalsigns$VSHEIGHT)
vitalsigns$VSHEIGHT <- vitalsigns$VSHEIGHT/100
vitalsigns <- subset(vitalsigns, select = -c(VSWTUNIT,VSHTUNIT))

# Fix height per RID using first available
vitalsigns <- group_by(vitalsigns,RID)
height<-slice(vitalsigns,1)
vitalsigns$VSHEIGHT<-height$VSHEIGHT[match(vitalsigns$RID,height$RID)]
vitalsigns$BMI<-vitalsigns$VSWEIGHT/(vitalsigns$VSHEIGHT*vitalsigns$VSHEIGHT)  

vitalsigns$BMIClass <- dplyr::case_when(
  vitalsigns$BMI < 18.5 ~ "Underweight",
  vitalsigns$BMI >= 18.5 & vitalsigns$BMI < 25 ~ "Normal",
  vitalsigns$BMI >= 25   & vitalsigns$BMI < 30 ~ "Overweight",
  TRUE ~ "Obesity"
)

vitalsigns$ObesityClass <- dplyr::case_when(
  vitalsigns$BMIClass=="Obesity" & vitalsigns$BMI>=30 & vitalsigns$BMI<35 ~ "Class 1",
  vitalsigns$BMIClass=="Obesity" & vitalsigns$BMI>=35 & vitalsigns$BMI<40 ~ "Class 2",
  vitalsigns$BMIClass=="Obesity" & vitalsigns$BMI>=40 ~ "Class 3",
  TRUE ~ "No"
)  

vitalsigns$flag         <- make_flag(vitalsigns$RID, vitalsigns$VISCODE2)
population$VSBPSYS      <- vitalsigns$VSBPSYS[match(population$flag, vitalsigns$flag)]
population$VSBPDIA      <- vitalsigns$VSBPDIA[match(population$flag, vitalsigns$flag)]
population$BMI          <- vitalsigns$BMI[match(population$flag, vitalsigns$flag)]
population$BMIClass     <- vitalsigns$BMIClass[match(population$flag, vitalsigns$flag)]
population$ObesityClass <- vitalsigns$ObesityClass[match(population$flag, vitalsigns$flag)]

write.csv(population,
          file.path(out_db, "7. Population_VitalSigns.csv"),
          row.names = FALSE)

# APOE
APOE <- read.csv(file.path("1. ADNI/APOE Results","APOERES_27Mar2024.csv")) %>%
  subset(select = c(RID,APGEN1,APGEN2))
APOE$APOE4 <- ifelse(APOE$APGEN1==4 | APOE$APGEN2==4, 1, 0)
population$APOE4 <- APOE$APOE4[match(population$RID, APOE$RID)]

write.csv(population,
          file.path(out_db, "8. Population_APOE.csv"),
          row.names = FALSE)
# CSF
CSF <- read.csv(file.path("1. ADNI/CSF Metabolites","UPENNBIOMK_ROCHE_ELECSYS_27Mar2024.csv"))
CSF <-subset(CSF,select=c(RID,VISCODE2,ABETA42,TAU,PTAU))
CSF <-CSF[which(CSF$VISCODE2=="bl"),]
population$ABETA42 <- CSF$ABETA42[match(population$RID, CSF$RID)]
population$TAU     <- CSF$TAU[match(population$RID, CSF$RID)]
population$PTAU    <- CSF$PTAU[match(population$RID, CSF$RID)]
#Keep patients that have all the information
population<-population[complete.cases(population),]
write.csv(population,
          file.path(out_db, "9. Population_CSF.csv"),
          row.names = FALSE)

# Lipidomics
lipidomics <- read.csv(file.path("1. ADNI/Lipidomics","ADMCLIPIDOMICSMEIKLELABLONG_08_13_21_19Mar2024.csv"))
lipidomics <- subset(lipidomics, select = -c(update_stamp,VISCODE,GUSPECID,COHORT,SAMPLE.ID,EXAMDATE))
lipidomics[ ,3:783] <- lapply(lipidomics[ ,3:783], safe_log2)
lipidomics$flag <- make_flag(lipidomics$RID, lipidomics$VISCODE2)

# Join
j <- ncol(population)
for (i in 3:(ncol(lipidomics)-1)) {
  population[ , i + j - 2] <- lipidomics[ , i][match(population$flag, lipidomics$flag)]
  colnames(population)[i + j - 2] <- colnames(lipidomics)[i]
}
#Keep patients that have all the information
population<-population[complete.cases(population),]
write.csv(population,
          file.path(out_db, "10. Population_Lipidomics.csv"),
          row.names = FALSE)

# Gut metabolites
gut <- read.csv(file.path("1. ADNI/Gut Metabolites","ADMCGUTMETABOLITESLONG_12_13_21_27Mar2024.csv"))
gut <- subset(gut, select = -c(EXAMDATE,VISCODE,DRAWDATE,SAMPLE_ID,VID,DRAWTIME,update_stamp))
gut[ ,3:106] <- lapply(gut[ ,3:106], safe_log2)
gut$flag <- make_flag(gut$RID, gut$VISCODE2)

j <- ncol(population)
for (i in 3:(ncol(gut)-1)) {
  population[ , i + j - 2] <- gut[ , i][match(population$flag, gut$flag)]
  colnames(population)[i + j - 2] <- colnames(gut)[i]
}

#Keep patients that have all the information
population<-population[complete.cases(population),]
write.csv(population,
          file.path(out_db, "11. Population_GutMicrobialMetabolites.csv"),
          row.names = FALSE)

# 6) Types & Save cohort snapshot 
population$Classification <- factor(population$Classification, levels=c("MCI-s","MCI-p"))
population$PTGENDER       <- factor(population$PTGENDER)
population$BMIClass       <- factor(population$BMIClass)
population$ObesityClass   <- factor(population$ObesityClass)
population$APOE4          <- factor(population$APOE4, levels=c(0,1))


write.csv(population, file.path(out_db, "12.DatasetMCI5AD2_rawCompleteCases.csv"), row.names = FALSE)

population <- subset(population, select = -c(flag, RID, VISCODE2))

# Exploratory plots
DataExplorer::plot_intro(population)
DataExplorer::plot_bar(population, by="Classification")

# 7) Train/Test split BEFORE feature selection (prevent leakage)
# We keep a hold-out (20%) for final unbiased evaluation
set.seed(0)
idx_val <- caret::createDataPartition(population$Classification, p=0.80, list=FALSE)
train <- population[idx_val, ]
test  <- population[-idx_val, ]
write.csv(train, file.path(out_db, "13.DatasetMCI5AD2_TRAIN.csv"), row.names = FALSE)
write.csv(test, file.path(out_db, "14.DatasetMCI5AD2_TEST.csv"), row.names = FALSE)

# 8) Feature filtering & mRMR on TRAIN ONLY 
# Prepare numeric frame for selection
train_sel <- train
train_sel$Classification <- ifelse(train_sel$Classification=="MCI-p", 1, 0)
train_sel$APOE4 <- ifelse(train_sel$APOE4==1, 1, 0)
train_sel <- data.frame(lapply(train_sel, function(x) if(is.factor(x)) as.numeric(x) else x))

q_values_idi <- univariate_Logit(data=train_sel, Outcome="Classification",
                                 pvalue=0.05, adjustMethod="BH", uniTest="zIDI", thr=0.999)
q_values_nri <- univariate_Logit(data=train_sel, Outcome="Classification",
                                 pvalue=0.05, adjustMethod="BH", uniTest="zNRI", thr=0.999)
q_values_ftest <- univariate_residual(data=train_sel, Outcome="Classification",
                                      pvalue=0.05, adjustMethod="BH", uniTest="Ftest", type="LOGIT", thr=0.999)
q_values_w <- univariate_Wilcoxon(data=train_sel, Outcome="Classification",
                                  pvalue=0.05, adjustMethod="BH", thr=0.999)
q_values_k <- univariate_correlation(data=train_sel, Outcome="Classification",
                                     pvalue=0.05, adjustMethod="BH", method="kendall", thr=0.999)

# Join q-values across filters
vecs <- list(
  IDI = q_values_idi,
  NRI = q_values_nri,
  F_Test = q_values_ftest,
  Wilcoxon   = q_values_w,
  Kendall   = q_values_k
)

features <- Reduce(union, lapply(vecs, names))

qmat <- data.frame(
  IDI = vecs$IDI[features],
  NRI = vecs$NRI[features],
  F_Test   = vecs$F_Test[features],
  Wilcoxon   = vecs$Wilcoxon[features],
  Kendall   = vecs$Kendall[features],
  row.names = features,
  check.names = FALSE
)

obs <- colSums(!is.na(qmat))

write.csv(qmat, file.path(out_db, "15.AllSignificantFeatures.csv"), row.names = TRUE)

features <- rownames(qmat)
other <- c("PTGENDER","PTEDUCAT","AGE","MMSCORE","GDTOTAL","ADAS11","ADNI_MEM",
             "ADNI_EF","ADNI_LAN","ADNI_VS","VSBPSYS","VSBPDIA","BMI","BMIClass",
           "ObesityClass","APOE4","ABETA42","TAU","PTAU")
other <-data.frame(other)
other$present<-other$other%in%features

alpha <- 0.05
min_hits <- 3

sig_mat <- qmat[, c("IDI","NRI","F_Test","Wilcoxon","Kendall")] < alpha



# Heatmap (features x tests)
if (nrow(qmat) > 2) {
  gplots::heatmap.2(as.matrix(qmat), Rowv=FALSE, dendrogram="col", main="Features q.values", cexRow=0.2, trace="none")
}

rn <- rownames(qmat)

qmat <- data.frame(lapply(qmat, function(x) {
  if(is.numeric(x)) sprintf("%.3e", x) else x
}))

rownames(qmat)<-rn

hits <- rowSums(sig_mat, na.rm = TRUE)
qmat$hits <- hits

keep <- hits >= min_hits
selected_features <- rownames(qmat)[keep]

cat(sprintf("Seleccionadas %d de %d features (criterio: ≥ %d tests con q < %.3f)\n",
            sum(keep), nrow(qmat), min_hits, alpha))
print(table(hits))
features <- rownames(qmat)
other$present2<-other$other%in%features
# Save Results
write.csv(qmat[order(-qmat$hits), ],
          file.path(out_db, "16.HitsSignificantFeatures.csv"))
write.csv(data.frame(feature = selected_features),
          file.path(out_db, "17.Selected_features_ge3of5.csv"),
          row.names = FALSE)

#Clinical <- c("ADNI_MEM","ADAS11","MMSCORE","AGE","ABETA42",
#              "ADNI_EF","APOE4","PTAU","TAU","ADNI_LAN","ADNI_VS","GDTOTAL")
#outcome_col <- "Classification"

#metabolomic <- setdiff(selected_features,c(outcome_col, Clinical))


# mRMR on filtered set
train_mrmr <- train_sel[c("Classification", selected_features)]
mrmr_res <- mRMR.classic_FRESA(data=train_mrmr, Outcome="Classification")
mrmr_df  <- as.data.frame(mrmr_res)
mrmr_df$Feature <- rownames(mrmr_df)
mrmr_df <- mrmr_df %>% dplyr::select(Feature, dplyr::everything())
# Avoid ambiguous column name "mrmr_results"
colnames(mrmr_df)[2] <- "mRMR_score"
mrmr_df <- mrmr_df[order(mrmr_df$mRMR_score, decreasing=TRUE),]

theme_set(theme_bw(base_size=10))
# Bar plot of mRMR ranking
ggplot(mrmr_df, aes(x=reorder(Feature, mRMR_score), y=mRMR_score)) +
  geom_bar(stat="identity") + coord_flip() +
  ggtitle("mRMR features") + xlab("Feature") + ylab("Mutual Information score")

mrmr_df <- mrmr_df %>%
  mutate(cumFrac = cumsum(mRMR_score) / sum(mRMR_score))

N = which(mrmr_df$cumFrac >= 0.9)[1]
# # Select top N features (tune N if needed)
#topN <- min(, nrow(mrmr_df_met))  # TODO: decide per-validation
selected_features <- c("Classification", 
                       mrmr_df$Feature[1:N])

mrmr_df2<-mrmr_df[1:N,]

theme_set(theme_bw(base_size=10))
# Bar plot of mRMR ranking
ggplot(mrmr_df2, aes(x=reorder(Feature, mRMR_score), y=mRMR_score)) +
  geom_bar(stat="identity") + coord_flip() +
  ggtitle("mRMR features") + xlab("Feature") + ylab("Mutual Information score")

features <- rownames(mrmr_df2)
other$present3<-other$other%in%features
lip<-rownames(mrmr_df2)%in%colnames(lipidomics)
table(lip)
gut<-rownames(mrmr_df2)%in%colnames(gut)
table(gut)
# 
# Reduce TRAIN and TEST to selected features
train_final <- train_sel[ selected_features]
test_final  <- test[ selected_features]
# Apply same transforms to TEST
test_final$Classification <- ifelse(test_final$Classification=="MCI-p", 1, 0)
test_final$APOE4 <- ifelse(test_final$APOE4==1, 1, 0)
#test_final <- data.frame(lapply(test_final, function(x) if(is.factor(x)) as.numeric(x) else x))
#test_final <- test_final[ , intersect(colnames(test_final), colnames(train_final))]
# 
# # 9) Table for selected features 
summary_tbl <- data.frame(Feature=colnames(train_final), stringsAsFactors = FALSE)
summary_tbl <- summary_tbl[summary_tbl$Feature!="Classification",,drop=FALSE]
# 
# # Attach mRMR & q-values
summary_tbl$mRMR <- mrmr_df$mRMR_score[match(summary_tbl$Feature, mrmr_df$Feature)]
summary_tbl$IDI  <- qmat$IDI[match(summary_tbl$Feature, rownames(qmat))]
summary_tbl$NRI  <- qmat$NRI[match(summary_tbl$Feature, rownames(qmat))]
summary_tbl$F_Test  <- qmat$F_Test[match(summary_tbl$Feature, rownames(qmat))]
summary_tbl$Wilcoxon    <- qmat$Wilcoxon[match(summary_tbl$Feature, rownames(qmat))]
summary_tbl$Kendall    <- qmat$Kendall[match(summary_tbl$Feature, rownames(qmat))]
summary_tbl$hits <- qmat$hits[match(summary_tbl$Feature, rownames(qmat))]
# 
write.csv(summary_tbl, file.path(out_db, "18.Summary_Table.csv"), row.names = FALSE)
write.csv(train_final, file.path(out_db, "19.Train_final.csv"), row.names = FALSE)
write.csv(test_final, file.path(out_db, "20.Test_final.csv"), row.names = FALSE)
# 

dev.off()


