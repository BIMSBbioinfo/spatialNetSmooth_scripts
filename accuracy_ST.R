library(tidyr)
library(reshape2)
library(Seurat)

accuracy <- NULL
accuracy_raw <- NULL
letters <- c("A1", "B1", "C1", "D1", "E1")
for (letter in letters) {
  path <- paste("~/BA/images ST/", letter, "/", sep="")
  setwd(path)
  truth <- read.csv("truth.csv", row.names=1)
  raw <- read.csv("gsea_raw.csv", row.names=1)
  
  nn_spatial <- read.csv(paste("nn_spatial_scores_", letter, ".csv", sep=""), row.names=1)[,16]
  d <-roc_quant(truth$V1, nn_spatial)
  f <- roc_quant(truth$V1, raw$V1)
  accuracy_raw <- append(accuracy_raw, f$accuracy)
  accuracy <- append(accuracy, d$accuracy)
}

letter <- "H1"
path <- paste("~/BA/images ST/", letter, "/", sep="")
setwd(path)
truth <- read.csv("truth.csv", row.names=1)
raw <- read.csv("gsea_raw.csv", row.names=1)

nn_spatial <- read.csv(paste("nn_spatial_scores_", letter, ".csv", sep=""), row.names=1)[,16]
d <-roc_quant(truth$V1, nn_spatial)
f <- roc_quant(truth$V1, raw$V1)
accuracy_raw <- append(accuracy_raw, f$accuracy)
accuracy <- append(accuracy, d$accuracy)


accuracy <- rbind(accuracy_raw, accuracy)
rownames(accuracy) <- c("unsmoothed", "smoothed")
colanmes <-  c(letters, "H1")

setwd("~/BA/images ST/")
write.csv(accuracy, "accuracy_comp_ST.csv")



