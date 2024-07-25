library(pscl)
#letter <- "C1"

gsea_raw <-read.csv(paste("~/BA/images ST/", letter, "/gsea_raw.csv", sep=""), row.names=1)
gsea_raw <- as.data.frame(gsea_raw)
colnames(gsea_raw) <- "gsea"
gsea_raw <- cbind(truth, gsea_raw)
log_model <- glm(truth ~ gsea, data = gsea_raw, family = "binomial")
r_raw <- pR2(log_model)
write.csv(r_raw, paste("~/BA/images ST/", letter, "/R2/raw_R2.csv", sep=""))
l <- c("nn", "spatial", "union")
for (method in l) {
  temp <- read.csv(paste("~/BA/images ST/", letter, "/", method, "_scores_", letter, ".csv", sep=""), row.names = 1)
  temp <- cbind(truth, temp)
  alphas <- c("alpha2", "alpha4", "alpha6", "alpha8")
  colnames(temp) <-c("truth", "alpha2", "alpha4", "alpha6", "alpha8")
  write.csv(temp, paste("~/BA/images ST/", letter,"/", method,"_scores_", letter, ".csv", sep=""))
  model_temp <- list()
  for(i in alphas){
    formula <- paste("truth~",i, sep="")
    model<- glm(formula, data = temp, family = "binomial")
    model_temp[[paste(i, collapse = "_")]] <- pR2(model)
  }
  write.csv(model_temp, paste("~/BA/images ST/", letter, "/R2/", method, "_R2.csv", sep=""))
}
l <- c("alpha", "nn_spatial")
alphas <- c("alpha2alpha2","alpha2alpha4", "alpha2alpha6", "alpha2alpha8", "alpha4alpha2", "alpha4alpha4", "alpha4alpha6", "alpha4alpha8", "alpha6alpha2", "alpha6alpha4", "alpha6alpha6", "alpha6alpha8", "alpha8alpha2", "alpha8alpha4", "alpha8alpha6", "alpha8alpha8")

for (method in l) {
  temp <- read.csv(paste("~/BA/images ST/", letter, "/", method, "_scores_", letter, ".csv", sep=""), row.names = 1)
  temp <- cbind(truth, temp)
  colnames(temp) <-c("truth", alphas)
  write.csv(temp, paste("~/BA/images ST/", letter,"/", method,"_scores_", letter,".csv", sep=""))
  model_temp <- list()
  for(i in alphas){
    formula <- paste("truth~",i, sep="")
    model<- glm(formula, data = temp, family = "binomial")
    model_temp[[paste(i, collapse = "_")]] <- pR2(model)
  }
  write.csv(model_temp, paste("~/BA/images ST/", letter, "/R2/", method, "_R2.csv", sep=""))
}
library(pROC)
#get precision recall curve
truth <- as.vector(truth)
scores<-read.csv(paste("~/BA/images ST/", letter, "/nn_spatial_scores_", letter, ".csv", sep=""), row.names=1)
scores <- as.vector(scores[,17])
roc_c <-roc(response= truth, predictor=scores)
plot(roc_c)
coords <- coords(roc_c, "best", best.method="closest.topleft")
best_threshold <- coords$threshold
pred <- ifelse(scores >=best_threshold, 1, 0)
conf_mat <- table(truth, pred)
accuracy <- sum(diag(conf_mat))/sum(conf_mat)
#thresholds_H1 <- NULL
thresholds_D1$nn_spatial$threshold <- best_threshold
thresholds_D1$nn_spatial$accuracy <- accuracy
L1 <- cbind(thresholds_A1$nn_spatial, thresholds_B1$nn_spatial, thresholds_C1$nn_spatial, thresholds_D1$nn_spatial, thresholds_E1$nn_spatial, thresholds_H1$nn_spatial)
colnames(L1) <- c("A1", "B1", "C1", "D1", "E1","H1")
write.csv(L1, "~/BA/best_thresholds_nn_spatial.csv", row.names=T)
L2 <- cbind(thresholds_A1$raw, thresholds_B1$raw, thresholds_C1$raw, thresholds_D1$raw, thresholds_E1$raw, thresholds_H1$raw)
colnames(L2) <- c("A1", "B1", "C1", "D1", "E1","H1")
write.csv(L2, "~/BA/best_thresholds_raw.csv", row.names=T)

