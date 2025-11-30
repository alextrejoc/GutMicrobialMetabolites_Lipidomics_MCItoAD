# =============================
# Plot significant features
# Author: Alejandro I. Trejo-Castro 
# Last version 11.11.25
# Goal: 
# 1. Make the plot of significant features 
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

load_pkg(c("ggplot2","ggrain","ggpubr"))

# 1) Paths & Output directories
# Adjust 'root' to your environment
root <- "G:/Mi unidad/Investigación Trejo-Castro/E31 - DBC Metabolic Health and Disease/E31.1 Alzheimer_GutMetabolites_Lipidomics"
setwd(root)

# 2) Train data 
train<-read.csv("3. Results Databases/1. Exp MCI_to_AD/19.Train_final.csv")
train <- dplyr::select(train, c("Classification","ADAS11", "ABETA42", "LPC.O.20.1.","LPI.20.4...SN2.", "C4_0","C8_0"))

pval<-read.csv("3. Results Databases/1. Exp MCI_to_AD/18.Summary_Table.csv")
pval<-pval[which(pval$Feature%in%colnames(train)),]
pval<-t(pval)
colnames(pval) <- pval[1, ]
pval <- pval[-c(1,2,8), ]   

pval[is.na(pval)] <- "n.s."


# Classification as Factor
train$Classification<-ifelse(train$Classification==1,"MCI-p", "MCI-s")
train$Classification <- factor(train$Classification,
                               levels = c("MCI-s","MCI-p"),
                               labels = c("MCI-Stable","MCI-Progressive"))

# 3) Plot

## ADAS

pvals_ADAS <- pval[, "ADAS11"]
pvals_text <- paste0(names(pvals_ADAS), ": ", format(pvals_ADAS, scientific = TRUE, digits = 3), collapse = "\n")


theme_set(theme_bw())
plot_ADAS <-ggplot(train, aes(Classification, ADAS11, fill = Classification, colour = Classification)) +
  geom_rain() +
  scale_fill_manual(values = c("#e6e6fa","#6e6ec5")) +
  scale_color_manual(values = c("black","black")) + ylab("Alzheimer's Disease Assessment Scale 11 items score") + 
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.caption = element_text(size = 14)
) + guides(fill = 'none', color = 'none')+ guides(fill = 'none', color = 'none') +
  labs(caption = paste("Adjusted p-values:\n", pvals_text))




## ABETA42
pvals_AB42 <- pval[, "ABETA42"]

pvals_text <- paste0(names(pvals_AB42), ": ", format(pvals_AB42, scientific = TRUE, digits = 3), collapse = "\n")


theme_set(theme_bw())
plot_AB42 <-ggplot(train, aes(Classification, ABETA42, fill = Classification, colour = Classification)) +
  geom_rain() +
  scale_fill_manual(values = c("#e6e6fa","#6e6ec5")) +
  scale_color_manual(values = c("black","black"))+ ylab("CSF AB42 [pg/ml]") + 
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.caption = element_text(size = 14)
  ) + guides(fill = 'none', color = 'none')+ guides(fill = 'none', color = 'none') +
  labs(caption = paste("Adjusted p-values:\n", pvals_text))


## LPC.O.20.1.
pvals_LPC<- pval[, "LPC.O.20.1."]
pvals_text <- paste0(names(pvals_LPC), ": ", format(pvals_LPC, scientific = TRUE, digits = 3), collapse = "\n")

theme_set(theme_bw())
plot_LPC <-ggplot(train, aes(Classification, LPC.O.20.1., fill = Classification, colour = Classification)) +
  geom_rain() +
  scale_fill_manual(values = c("#e6e6fa","#6e6ec5")) +
  scale_color_manual(values = c("black","black"))+ ylab("LPC(O-20:1). Lipid metabolite log2[nM]") + 
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.caption = element_text(size = 14)
  ) + guides(fill = 'none', color = 'none')+ guides(fill = 'none', color = 'none') +
  labs(caption = paste("Adjusted p-values:\n", pvals_text))


## LPI.20.4...SN2.
pvals_LPI<- pval[, "LPI.20.4...SN2."]
pvals_text <- paste0(names(pvals_LPI), ": ", format(pvals_LPI, scientific = TRUE, digits = 3), collapse = "\n")


theme_set(theme_bw())
plot_LPI <-ggplot(train, aes(Classification, LPI.20.4...SN2., fill = Classification, colour = Classification)) +
  geom_rain() +
  scale_fill_manual(values = c("#e6e6fa","#6e6ec5")) +
  scale_color_manual(values = c("black","black")) + ylab("LPI(20:4) [sn2]. Lipid metabolite log2[nM]") + 
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.caption = element_text(size = 14)
  ) + guides(fill = 'none', color = 'none')+ guides(fill = 'none', color = 'none') +
  labs(caption = paste("Adjusted p-values:\n", pvals_text))

## C4_0
pvals_C4<- pval[, "C4_0"]
pvals_text <- paste0(names(pvals_C4), ": ", format(pvals_C4, scientific = TRUE, digits = 3), collapse = "\n")

theme_set(theme_bw())
plot_C4_0 <-ggplot(train, aes(Classification, C4_0, fill = Classification, colour = Classification)) +
  geom_rain() +
  scale_fill_manual(values = c("#e6e6fa","#6e6ec5")) +
  scale_color_manual(values = c("black","black"))+ ylab("Butyric Acid log2[M]") + 
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.caption = element_text(size = 14)
  ) + guides(fill = 'none', color = 'none')+ guides(fill = 'none', color = 'none') +
  labs(caption = paste("Adjusted p-values:\n", pvals_text))


## C8_0
pvals_C8<- pval[, "C8_0"]
pvals_text <- paste0(names(pvals_C8), ": ", format(pvals_C8, scientific = TRUE, digits = 3), collapse = "\n")

theme_set(theme_bw())
plot_C8_0 <-ggplot(train, aes(Classification, C8_0, fill = Classification, colour = Classification)) +
  geom_rain() +
  scale_fill_manual(values = c("#e6e6fa","#6e6ec5")) +
  scale_color_manual(values = c("black","black"))+ ylab("Caprylic Acid log2[M]") + 
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.caption = element_text(size = 14)
  ) + guides(fill = 'none', color = 'none')+ guides(fill = 'none', color = 'none') +
  labs(caption = paste("Adjusted p-values:\n", pvals_text))


figure <- ggarrange(plot_ADAS, plot_AB42, plot_LPC,
                    plot_LPI,plot_C4_0,plot_C8_0,
                    ncol = 3, nrow = 2)

# 4) Save
ggsave("4. Results Plots/1. Exp MCI_to_AD/Plot_ADAS.tiff", plot = plot_ADAS, width = 10, height = 10, dpi = 600)
ggsave("4. Results Plots/1. Exp MCI_to_AD/Plot_AB42.tiff", plot = plot_AB42, width = 10, height = 10, dpi = 600)
ggsave("4. Results Plots/1. Exp MCI_to_AD/Plot_LPI.tiff", plot = plot_LPI, width = 10, height = 10, dpi = 600)
ggsave("4. Results Plots/1. Exp MCI_to_AD/Plot_LPC.tiff", plot = plot_LPC, width = 10, height = 10, dpi = 600)
ggsave("4. Results Plots/1. Exp MCI_to_AD/Plot_C4_0.tiff", plot = plot_C4_0, width = 10, height = 10, dpi = 600)
ggsave("4. Results Plots/1. Exp MCI_to_AD/Plot_C8.tiff", plot = plot_C8_0, width = 10, height = 10, dpi = 600)

ggsave("4. Results Plots/1. Exp MCI_to_AD/Plot_Features.tiff", plot = figure, width = 20, height = 17.5, dpi = 600)



