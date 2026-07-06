# =============================
# Age and Sex Analysis
# Author: Alejandro I. Trejo-Castro 
# Last version 30.11.25
# Goal: 
# 1. Analysis of Sex and Age Features
# =============================

root <- "G:/Mi unidad/Investigación Trejo-Castro/E31 - DBC Metabolic Health and Disease/E31.1 Alzheimer_GutMetabolites_Lipidomics"
setwd(root)

data<-read.csv("3. Results Databases/1. Exp MCI_to_AD/13.DatasetMCI5AD2_TRAIN.csv")
wilcox.test(data[which(data$Classification=="MCI-p"),]$AGE,
            data[which(data$Classification=="MCI-s"),]$AGE,)

quantile(data[which(data$Classification=="MCI-p"),]$AGE)
quantile(data[which(data$Classification=="MCI-s"),]$AGE)


table(data$PTGENDER)
# 1 male, 2 female 
# Create a contingency table
data$PTGENDER<-ifelse(data$PTGENDER==1,"Male","Female")
contingency_table <- table(data$PTGENDER, data$Classification)
print(contingency_table)

chi_square_result <- chisq.test(contingency_table)
print(chi_square_result)


library(ggplot2)

ggplot(data, aes(x = PTGENDER, fill = Classification)) +
  geom_bar(position = "fill") +
  labs(title = "Classification Distribution by PTGENDER",
       y = "Proportion",
       fill = "Classification") +
  theme_minimal()

ggplot(data, aes(x = Classification, fill = PTGENDER)) +
  geom_bar(position = "fill") +
  labs(title = "Classification Distribution by PTGENDER",
       y = "Proportion",
       fill = "PTGENDER") +
  theme_minimal()
