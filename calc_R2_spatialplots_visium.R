library(pscl)


filepath_output = "~/BA/images Visium/"
setwd(filepath_output)
truth <- read.csv("truth.csv", row.names=1)
gsea_nn <- read.csv("nn_scores.csv", row.names=1)
gsea_snn <- read.csv("snn_scores.csv", row.names=1)
gsea_inter <- read.csv("inter_snn_scores.csv", row.names=1)

gsea_nn_spatial <- read.csv("nn_spatial_scores.csv", row.names=1)
gsea_snn_spatial <- read.csv("snn_spatial_scores.csv", row.names=1)
gsea_spatial <- read.csv("spatial_scores.csv", row.names=1)
gsea_union<- read.csv("union_scores.csv", row.names=1)
gsea_alpha <- read.csv("alpha_scores.csv", row.names=1)
gsea_raw <- read.csv("raw.csv", row.names=1)
colnames(gsea_raw) <- "gsea"
truth <- as.vector(truth$x)

#making log-regression-models
log_model <- glm(truth ~ gsea, data = cbind(truth,gsea_raw), family = "binomial")
r_raw <- pR2(log_model)
write.csv(r_raw, "raw_R2.csv")

l <- c("nn", "spatial", "union", "snn", "inter")
for (method in l) {
  temp <- read.csv(paste("~/BA/images Visium/", method, "_scores.csv", sep=""), row.names = 1)
  temp <- cbind(truth, temp)
  alphas <- c("alpha2", "alpha4", "alpha6", "alpha8")
  colnames(temp) <-c("truth", "alpha2", "alpha4", "alpha6", "alpha8")
  write.csv(temp, paste("~/BA/images Visium/", method,"_scores.csv", sep=""))
  model_temp <- list()
  for(i in alphas){
    formula <- paste("truth~",i, sep="")
    model<- glm(formula, data = temp, family = "binomial")
    model_temp[[paste(i, collapse = "_")]] <- pR2(model)
  }
  write.csv(model_temp, paste(method, "_R2.csv", sep=""))
}
l <- c("alpha", "nn_spatial", "snn_spatial")
alphas <- c("alpha2alpha2","alpha2alpha4", "alpha2alpha6", "alpha2alpha8", "alpha4alpha2", "alpha4alpha4", "alpha4alpha6", "alpha4alpha8", "alpha6alpha2", "alpha6alpha4", "alpha6alpha6", "alpha6alpha8", "alpha8alpha2", "alpha8alpha4", "alpha8alpha6", "alpha8alpha8")

for (method in l) {
  temp <- read.csv(paste("~/BA/images Visium/", method, "_scores.csv", sep=""), row.names = 1)
  temp <- cbind(truth, temp)
  colnames(temp) <-c("truth", alphas)
  write.csv(temp, paste("~/BA/images Visium/", method,"_scores.csv", sep=""))
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

pdf("plot_nn.pdf")
plot_quant(gsea_nn[[4]], coordinates, truth)
dev.off()
pdf("plot_snn.pdf")
plot_quant(gsea_snn[[4]], coordinates, truth)
dev.off()
pdf("plot_spatial.pdf")
plot_quant(gsea_spatial[[4]], coordinates, truth)
dev.off()

pdf("plot_union.pdf")
plot_quant(gsea_union[[4]], coordinates, truth)
dev.off()
pdf("plot_inter.pdf")
plot_quant(gsea_inter[[4]], coordinates, truth)
dev.off()
pdf("plot_alpha.pdf")
plot_quant(alpha[,16], coordinates, truth)
dev.off()
pdf("plot_nn_spatial.pdf")
plot_quant(nn_spatial[,16], coordinates, truth)
dev.off()

pdf("plot_snn_spatial.pdf")
plot_quant(snn_spatial[,16], coordinates, truth)
dev.off()
