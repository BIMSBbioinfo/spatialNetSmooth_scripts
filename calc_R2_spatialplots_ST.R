library(pscl)


letter <- "A1"
filepath_output = "~/BA/images ST/"
setwd(paste(filepath_output, letter, sep=""))
truth <- read.csv("truth.csv", row.names=1)
gsea_nn <- read.csv(paste("nn_scores_", letter, ".csv", sep=""), row.names=1)
gsea_nn_spatial <- read.csv(paste("nn_spatial_scores_", letter, ".csv", sep=""), row.names=1)
gsea_spatial <- read.csv(paste("spatial_scores_", letter, ".csv", sep=""), row.names=1)
gsea_union<- read.csv(paste("union_scores_", letter, ".csv", sep=""), row.names=1)
gsea_alpha <- read.csv(paste("alpha_scores_", letter, ".csv", sep=""), row.names=1)
gsea_raw <- read.csv("gsea_raw.csv", row.names=1)
colnames(gsea_raw) <- "gsea"
truth <- as.vector(truth$V1)

#making log-regression-models
log_model <- glm(truth ~ gsea, data = cbind(truth,gsea_raw), family = "binomial")
r_raw <- pR2(log_model)
write.csv(r_raw, "raw_R2.csv")

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
  write.csv(model_temp, paste(method, "_R2.csv", sep=""))
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
  write.csv(model_temp, paste(method, "_R2.csv", sep=""))
}

#select best parameters according to R^2 (look at csv's for each method by selecting the column in dataframe
#make spatial plots
pdf(paste("plot_snn_", letter, ".pdf", sep=""))
plot_quant(as.vector(gsea_nn[[4]]), coordinates, truth)
dev.off()
pdf(paste("plot_spatial_", letter, ".pdf", sep=""))
plot_quant(as.vector(gsea_spatial[[4]]), coordinates, truth)
dev.off()
pdf(paste("plot_union_", letter, ".pdf", sep=""))
plot_quant(as.vector(gsea_union[[4]]), coordinates, truth)
dev.off()

pdf(paste("plot_alpha_", letter, ".pdf", sep=""))
plot_quant(as.vector(alpha[,16]), coordinates, truth)
dev.off()
pdf(paste("plot_snn_spatial", letter, ".pdf", sep=""))
plot_quant(as.vector(nn_spatial[,16]), coordinates, truth)
dev.off()
