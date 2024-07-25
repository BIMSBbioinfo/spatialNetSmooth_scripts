library(pscl)

truth <- readRDS("~/score_smoothing/truth.Rds")

gsea_raw <-read.csv("~/BA/images Visium/raw.csv", row.names=1)
gsea_raw <- as.data.frame(gsea_raw)

gsea_raw <- cbind(truth, gsea_raw)
colnames(gsea_raw) <- c("truth","gsea")
log_model <- glm(truth ~ gsea, data = gsea_raw, family = "binomial")
r_raw <- pR2(log_model)
write.csv(r_raw, "~/BA/images Visium/R2/raw_R2.csv")
l <- c("nn", "spatial", "union", "snn", "inter_snn")
for (method in l) {
  temp <- read.csv(paste("~/BA/images Visium/", method, "_scores.csv", sep=""), row.names = 1)
  temp <- cbind(truth, temp)
  alphas <- c("alpha2", "alpha4", "alpha6", "alpha8")
  colnames(temp) <-c("truth", alphas)
  model_temp <- list()
  for(i in alphas){
    formula <- paste("truth~",i, sep="")
    model<- glm(formula, data = temp, family = "binomial")
    model_temp[[paste(i, collapse = "_")]] <- pR2(model)
  }
  write.csv(model_temp, paste("~/BA/images Visium/R2/", method, "_R2.csv", sep=""))
}
l <- c("alpha", "nn_spatial", "snn_spatial")
alphas <- c("alpha2alpha2","alpha2alpha4", "alpha2alpha6", "alpha2alpha8", "alpha4alpha2", "alpha4alpha4", "alpha4alpha6", "alpha4alpha8", "alpha6alpha2", "alpha6alpha4", "alpha6alpha6", "alpha6alpha8", "alpha8alpha2", "alpha8alpha4", "alpha8alpha6", "alpha8alpha8")

for (method in l) {
  temp <- read.csv(paste("~/BA/images Visium/", method, "_scores.csv", sep=""), row.names = 1)
  temp <- cbind(truth, temp)
  colnames(temp) <-c("truth", alphas)
  model_temp <- list()
  for(i in alphas){
    formula <- paste("truth~",i, sep="")
    model<- glm(formula, data = temp, family = "binomial")
    model_temp[[paste(i, collapse = "_")]] <- pR2(model)
  }
  write.csv(model_temp, paste("~/BA/images Visium/R2/", method, "_R2.csv", sep=""))
}
l <- c("nn", "spatial", "union", "snn", "inter_snn", "raw")
rocs <- NULL
for (method in l) {
  score <- read.csv(paste("~/BA/images Visium/", method, "_scores.csv", sep=""), row.names = 1)
  score <- as.vector(score[,4])
  rocs <- append(rocs, roc_quant(truth, score))
  pdf(paste("~/BA/images Visium/", method, "_roc_plot.pdf", sep=""))
  plot_quant(score, coordinates, truth)
  dev.off()
      
}
m <- c("nn_spatial", "snn_spatial", "alpha")
rocs2 <- NULL
for (method in m) {
  score <- read.csv(paste("~/BA/images Visium/", method, "_scores.csv", sep=""), row.names = 1)
  score <- as.vector(score[,16])
  rocs2 <- append(rocs2, roc_quant(truth, score))
  pdf(paste("~/BA/images Visium/", method, "_roc_plot.pdf", sep=""))
  plot_quant(score, coordinates, truth)
  dev.off()
  
}
best_threshs <- as.matrix(c(rocs[1], rocs[3], rocs[5], rocs[7], rocs[9], rocs[11], rocs2[1], rocs2[3], rocs2[5]), ncol=1)
rownames(best_threshs) <- c(l, m)
colnames(best_threshs) <-"threshold"
accuracy <- as.matrix(c(rocs[2], rocs[4], rocs[6], rocs[8], rocs[10], rocs[12], rocs2[2], rocs2[4], rocs2[6]), ncol=1)
rownames(accuracy) <- c(l, m)
colnames(accuracy) <-"accuracy"
write.csv(accuracy, "~/BA/images Visium/accuracy.csv", row.names = T)
write.csv(best_threshs, "~/BA/images Visium/thresholds_roc.csv", row.names = T)
