library(tidyr)
library(reshape2)
library(Seurat)
nn <- read.csv("~/BA/images Visium/nn_scores.csv", row.names=1)[,4]
snn <- read.csv("~/BA/images Visium/snn_scores.csv", row.names=1)[,4]
nn_spatial <- read.csv("~/BA/images Visium/nn_spatial_scores.csv", row.names=1)[,16]
snn_spatial <- read.csv("~/BA/images Visium/snn_spatial_scores.csv", row.names=1)[,16]
union <- read.csv("~/BA/images Visium/union_scores.csv", row.names=1)[,4]
inter <- read.csv("~/BA/images Visium/inter_snn_scores.csv", row.names=1)[,4]
spatial <- read.csv("~/BA/images Visium/spatial_scores.csv", row.names=1)[,4]
alpha <- read.csv("~/BA/images Visium/alpha_scores.csv", row.names=1)[,16]
raw <- read.csv("~/BA/images Visium/raw.csv", row.names=1)[,1]

F1_raw <- F1_quant(raw, truth,raw)
F1_nn <- F1_quant(nn, truth,raw)
F1_snn <- F1_quant(snn, truth,raw)
F1_nn_spatial <- F1_quant(nn_spatial, truth,raw)
F1_snn_spatial <- F1_quant(snn_spatial, truth,raw)
F1_spatial <- F1_quant(spatial, truth,raw)
F1_union <- F1_quant(union, truth,raw)
F1_inter <- F1_quant(inter, truth,raw)
F1_alpha <- F1_quant(alpha, truth,raw)

pred_raw <- ifelse(raw >=F1_raw$threshold[3], 1, 0)
pred_nn <- ifelse(nn >=F1_nn$threshold[5], 1, 0)
pred_snn <- ifelse(snn >=F1_snn$threshold[4], 1, 0)
pred_nn_spatial <- ifelse(nn_spatial >=F1_nn_spatial$threshold[5], 1, 0)
pred_snn_spatial <- ifelse(snn_spatial >=F1_snn_spatial$threshold[4], 1, 0)
pred_union <- ifelse(union >=F1_union$threshold[5], 1, 0)
pred_inter <- ifelse(inter >=F1_inter$threshold[3], 1, 0)
pred_spatial <- ifelse(spatial >=F1_spatial$threshold[4], 1, 0)
pred_alpha <- ifelse(alpha >=F1_alpha$threshold[5], 1, 0)

accuracy <- NULL
preds <- cbind(pred_raw, pred_nn, pred_snn, pred_nn_spatial, pred_snn_spatial, pred_union, pred_inter, pred_spatial,pred_alpha)
for (i in 1:9) {
  conf_mat <- table(truth, preds[,i])
  acc <- sum(diag(conf_mat))/sum(conf_mat)
  accuracy <- append(accuracy, acc)
}
accuracy <- as.matrix(accuracy, nrow=1)
rownames(accuracy) <- c("pred_raw", "pred_nn", "pred_snn", "pred_nn_spatial", "pred_snn_spatial", "pred_union", "pred_inter", "pred_spatial","pred_alpha")
colnames(accuracy) <- "accuracy"
acc2 <- read.csv("~/BA/images Visium/accuracy.csv", row.names=1)





