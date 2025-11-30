# =============================
# Plot mRMR score
# Author: Alejandro I. Trejo-Castro 
# Last version 11.11.25
# Goal: 
# 1. Make the plot of mRMR features with score
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

load_pkg(c("ggplot2"))

# 1) Paths & Output directories
# Adjust 'root' to your environment
root <- "G:/Mi unidad/Investigación Trejo-Castro/E31 - DBC Metabolic Health and Disease/E31.1 Alzheimer_GutMetabolites_Lipidomics"
setwd(root)

# 2) mRMR data and features description
summary<-read.csv("3. Results Databases/1. Exp MCI_to_AD/18.Summary_Table.csv")
features<-read.csv("3. Results Databases/1. Exp MCI_to_AD/Names_Final_Features.csv")

features$mRMR_score<-summary$mRMR[match(features$Feature,summary$Feature)]
features <- features[order(features$mRMR_score, decreasing=TRUE),]
features<-features[which(features$mRMR_score>0),]

# 3) Plot
theme_set(theme_bw())
# Bar plot of mRMR ranking
plot_mRMR<- ggplot(features, aes(x=reorder(Description, mRMR_score), y=mRMR_score)) +
            geom_bar(stat="identity", colour = "Black", fill="#9e9ee2") + coord_flip() +
            ggtitle("mRMR Features") + xlab("Feature") + ylab("Mutual Information score") +
            theme(
                   axis.text = element_text(size = 14),
                   axis.title = element_text(size = 16, face = "bold"),
                  plot.title = element_text(size = 18, face = "bold")
                  )

# 4) Save
ggsave("4. Results Plots/1. Exp MCI_to_AD/Plot_mRMR.tiff", plot = plot_mRMR, width = 10, height = 7, dpi = 600)