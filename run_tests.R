library(tidyr)
library(reshape2)
library(Seurat)

#load package (if not installed)
setwd("~/spatialNetSmooth/")
devtools::load_all()
#Run SpatialNetSmooth on Visium-Data
#Set input and output directories
filepath_data= "~/score_smoothing/data/"
filepath_output = "~/BA/images Visium/"
setwd(filepath_data)
#get coordinates
coord_tab = read.csv('./spatial/tissue_positions_list.csv', header=FALSE)

#get count matrix
filename = 'V1_Breast_Cancer_Block_A_Section_2_filtered_feature_bc_matrix.h5'

#Load object
se = Load10X_Spatial(
  './',
  filename = filename
)

#calculate GSEA
seu <- gseaCalc(se)


coordinates <- GetTissueCoordinates(seu, scale="lowres")

#get true labels
library(imager)

file_annot = './spatial/tissue_lowres_image_annotated.png'
im_annot = load.image(file_annot)

der = B(im_annot) %>% liner(thr = 0.58, fill=0) %>%
  bucketfill(x=200, y=200,color="darkorange",sigma=.2) %>%
  bucketfill(x=200, y=450, color="darkorange",sigma=.2) %>%
  bucketfill(x=350, y=330, color="darkorange",sigma=.2) %>%
  bucketfill(x=420, y=300, color="darkorange",sigma=.2) %>%
  bucketfill(x=450, y=350, color="darkorange",sigma=.2) %>%
  bucketfill(x=450, y=400, color="darkorange",sigma=.2) %>%
  bucketfill(x=480, y=420, color="darkorange",sigma=.2) %>% plot()



labeled_gray <- grayscale(der)

plot(labeled_gray)

labeled <- as.matrix(labeled_gray)

tumor <- NULL
coordinates<- GetTissueCoordinates(se, scale="lowres")
for(i in 1:nrow(coordinates)){
  row <- round(coordinates[i,1])
  col <- round(coordinates[i,2])
  istumor <- labeled[row, col]
  if(istumor==0){
    istumor <- 0
  }else if(istumor > 0.75){
    istumor <- 0
  }else{
    istumor <- 1
  }
  tumor <- append(tumor, istumor)
}
truth <- tumor
setwd(filepath_output)
save.csv(truth, "truth.csv")
gsea_raw <- seu@meta.data$gsea_rat_norm
#F1_raw <- F1_quant(gsea_raw, truth, gsea_raw)
pdf("gsea_raw_5.pdf")
plot_quant(gsea_raw, coordinates, truth, 5)
dev.off()
alphas <- c(0.2, 0.4, 0.6, 0.8)
gsea_spatial <- vector(mode = "list", length = 4)
  for (i in 1:4) {
    gsea_spatial[[i]] <- spatial_smooth(seu, a = alphas[i])

  }
spatial <- as.data.frame(gsea_spatial)
colnames(spatial) <- c("alpha=0.2", "alpha=0.4", "alpha=0.6", "alpha=0.8")
write.csv(spatial, "spatial_scores.csv", row.names = T)
spatial <- read.csv("spatial_scores.csv", row.names = 1)

F1_spatial <- vector(mode = "list", length = 4)
for (i in 1:4) {
  F1_spatial[[i]] <- F1_quant(spatial[,i], truth, gsea_raw)
  
}

gsea_nn <- vector(mode = "list", length = 4)
for (i in 1:4) {
  gsea_nn[[i]] <- nn_smooth(seu, a = alphas[i])
  
}
nn <- as.data.frame(gsea_nn)
colnames(nn) <- c("alpha=0.2", "alpha=0.4", "alpha=0.6", "alpha=0.8")
write.csv(nn, "nn_scores.csv", row.names = T)

F1_nn <- vector(mode = "list", length = 4)
for (i in 1:4) {
  F1_nn[[i]] <- F1_quant(nn[,i], truth)
  
}

gsea_snn <- vector(mode = "list", length = 4)
for (i in 1:4) {
  gsea_snn[[i]] <- nn_smooth(seu, a = alphas[i], graph="snn")
  
}
snn <- as.data.frame(gsea_snn)
colnames(snn) <- c("alpha=0.2", "alpha=0.4", "alpha=0.6", "alpha=0.8")
write.csv(snn, "snn_scores.csv", row.names = T)

F1_snn <- vector(mode = "list", length = 4)
for (i in 1:4) {
  F1_snn[[i]] <- F1_quant(snn[,i], truth)
  
}


gsea_nn_spatial <- vector(mode = "list", length = 4)
for (i in 1:4) { # [0.2, 0.2], [0.2,0.4], [0.2, 0.6], ...[0.4, 0.2], [0.4, 0.4]....
  item <- vector(mode = "list", length = 4)
  for (j in 1:4) {
    item[[j]] <- nn_spatial_smooth(seu, a1 = alphas[i],a2 = alphas[j])
  }
  gsea_nn_spatial[[i]] <- item
  
}



nn_spatial<- as.data.frame(gsea_nn_spatial)
colnames(nn_spatial) <- c("a1 = 0.2, a2=0.2","a1 = 0.2, a2=0.4", "a1 = 0.2, a2=0.6", "a1 = 0.2, a2=0.8", "a1 = 0.4, a2=0.2", "a1 = 0.4, a2=0.4", "a1 = 0.4, a2=0.6", "a1 = 0.4, a2=0.8", "a1 = 0.6, a2=0.2", "a1 = 0.6, a2=0.4", "a1 = 0.6, a2=0.6", "a1 = 0.6, a2=0.8", "a1 = 0.8, a2=0.2", "a1 = 0.8, a2=0.4", "a1 = 0.8, a2=0.6", "a1 = 0.8, a2=0.8")
write.csv(nn_spatial, "nn_spatial_scores.csv", row.names = T)

F1_nn_spatial <- vector(mode = "list", length = 16)
for (i in 1:16) {
  F1_nn_spatial[[i]] <- F1_quant(nn_spatial[,i], truth)
  
}

gsea_snn_spatial <- vector(mode = "list", length = 4)
for (i in 1:4) { # [0.2, 0.2], [0.2,0.4], [0.2, 0.6], ...[0.4, 0.2], [0.4, 0.4]....
  item <- vector(mode = "list", length = 4)
  for (j in 1:4) {
    item[[j]] <- nn_spatial_smooth(seu, a1 = alphas[i],a2 = alphas[j], graph="snn")
  }
  gsea_snn_spatial[[i]] <- item
  
}



snn_spatial<- as.data.frame(gsea_snn_spatial)
colnames(snn_spatial) <- c("a1 = 0.2, a2=0.2","a1 = 0.2, a2=0.4", "a1 = 0.2, a2=0.6", "a1 = 0.2, a2=0.8", "a1 = 0.4, a2=0.2", "a1 = 0.4, a2=0.4", "a1 = 0.4, a2=0.6", "a1 = 0.4, a2=0.8", "a1 = 0.6, a2=0.2", "a1 = 0.6, a2=0.4", "a1 = 0.6, a2=0.6", "a1 = 0.6, a2=0.8", "a1 = 0.8, a2=0.2", "a1 = 0.8, a2=0.4", "a1 = 0.8, a2=0.6", "a1 = 0.8, a2=0.8")
write.csv(snn_spatial, "snn_spatial_scores.csv", row.names = T)

F1_snn_spatial <- vector(mode = "list", length = 16)
for (i in 1:16) {
  F1_snn_spatial[[i]] <- F1_quant(snn_spatial[,i], truth)
  
}
gsea_union <- vector(mode = "list", length = 4)
  for (i in 1:4) {
  gsea_union[[i]] <- union_smooth(seu, a = alphas[i])
  
}
union <- as.data.frame(gsea_union)
colnames(union) <- c("a=0.2", "a=0.4", "a=0.6", "a=0.8")
write.csv(union, "union_scores.csv", row.names = T)

F1_union <- vector(mode = "list", length = 4)
for (i in 1:4) {
  F1_union[[i]] <- F1_quant(union[,i], truth)
  
}

gsea_inter <- vector(mode = "list", length = 4)
for (i in 1:4) {
  gsea_inter[[i]] <- inter_smooth(seu, a = alphas[i])
  
}
inter <- as.data.frame(gsea_inter)
colnames(inter) <- c("a=0.2", "a=0.4", "a=0.6", "a=0.8")
write.csv(inter, "inter_scores.csv", row.names = T)

gsea_inter_snn <- vector(mode = "list", length = 4)
for (i in 1:4) {
  gsea_inter_snn[[i]] <- inter_smooth(seu, a = alphas[i], graph="snn")
  
}
inter_snn <- as.data.frame(gsea_inter_snn)
colnames(inter_snn) <- c("a=0.2", "a=0.4", "a=0.6", "a=0.8")
write.csv(inter_snn, "inter_snn_scores.csv", row.names = T)

F1_inter <- vector(mode = "list", length = 4)
for (i in 1:4) {
  F1_inter[[i]] <- F1_quant(inter[,i], truth)
  
}

F1_inter_snn <- vector(mode = "list", length = 4)
for (i in 1:4) {
  F1_inter_snn[[i]] <- F1_quant(inter_snn[,i], truth)
  
}

gsea_alpha <- vector(mode = "list", length = 4)
for (i in 1:4) {
  item <- vector(mode = "list", length = 4)
  for (j in 1:4) {
    item[[j]] <- alpha_nn_spatial_smooth(seu, a = alphas[i], alpha = alphas[j])
  }
  gsea_alpha[[i]] <- item
  
}
alpha <- as.data.frame(gsea_alpha)
colnames(alpha) <- c("a = 0.2, alpha=0.2", "a = 0.2, alpha=0.4", "a = 0.2, alpha=0.6", "a = 0.2, alpha=0.8", "a = 0.4, alpha=0.2", "a = 0.4, alpha=0.4", "a = 0.4, alpha=0.6", "a = 0.4, alpha=0.8", "a = 0.6, alpha=0.2", "a = 0.6, alpha=0.4", "a = 0.6, alpha=0.6", "a = 0.6, alpha=0.8", "a = 0.8, alpha=0.2", "a = 0.8, alpha=0.4", "a = 0.8, alpha=0.6", "a = 0.8, alpha=0.8")
write.csv(alpha, "alpha_scores.csv", row.names = T)

F1_alpha <- vector(mode = "list", length = 16)
for (i in 1:16) {
  F1_alpha[[i]] <- F1_quant(alpha[,i], truth)
  
}
nn_df <- as.data.frame(cbind(F1_nn[[1]]$score, F1_nn[[2]]$score, F1_nn[[3]]$score, F1_nn[[4]]$score))
colnames(nn_df) <- c("a=0.2", "a=0.4", "a=0.6", "a=0.8")
nn_df$Position <- 1:nrow(nn_df)
nn_df$unsmoothed <-F1_raw$score
df_long <- melt(nn_df, id.vars = "Position", variable.name = "Smoothing Parameter", value.name = "Value")  

pdf("F1_nn.pdf")
ggplot(df_long, aes(x = Position, y = Value, color = `Smoothing Parameter`)) +
  geom_point() +
  theme_minimal() +
  labs(title = "F1 scores for nn_smoothing",
       x = "Quantile as Threshold",
       y = "F1 Score")
dev.off()

snn_df <- as.data.frame(cbind(F1_snn[[1]]$score, F1_snn[[2]]$score, F1_snn[[3]]$score, F1_snn[[4]]$score))
colnames(snn_df) <- c("a=0.2", "a=0.4", "a=0.6", "a=0.8")
snn_df$Position <- 1:nrow(snn_df)
snn_df$unsmoothed <-F1_raw$score
df_long <- melt(snn_df, id.vars = "Position", variable.name = "Smoothing Parameter", value.name = "Value")  

pdf("F1_snn.pdf")
ggplot(df_long, aes(x = Position, y = Value, color = `Smoothing Parameter`)) +
  geom_point() +
  theme_minimal() +
  labs(title = "F1 scores for snn_smoothing",
       x = "Quantile as Threshold",
       y = "F1 Score")
dev.off()

spatial_df <- as.data.frame(cbind(F1_spatial[[1]]$score, F1_spatial[[2]]$score, F1_spatial[[3]]$score, F1_spatial[[4]]$score))
colnames(spatial_df) <- c("a=0.2", "a=0.4", "a=0.6", "a=0.8")
spatial_df$Position <- 1:nrow(spatial_df)
spatial_df$unsmoothed <-F1_raw$score

df_long <- melt(spatial_df, id.vars = "Position", variable.name = "Smoothing Parameter", value.name = "Value")  

pdf("F1_spatial.pdf")
ggplot(df_long, aes(x = Position, y = Value, color = `Smoothing Parameter`)) +
  geom_point() +
  theme_minimal() +
  labs(title = "F1 scores for spatial_smoothing",
       x = "Quantile as Threshold",
       y = "F1 Score")
dev.off()

union_df <- as.data.frame(cbind(F1_union[[1]]$score, F1_union[[2]]$score, F1_union[[3]]$score, F1_union[[4]]$score))
colnames(union_df) <- c("a=0.2", "a=0.4", "a=0.6", "a=0.8")
union_df$Position <- 1:nrow(union_df)
union_df$unsmoothed <-F1_raw$score
df_long <- melt(union_df, id.vars = "Position", variable.name = "Smoothing Parameter", value.name = "Value")  

pdf("F1_union.pdf")
ggplot(df_long, aes(x = Position, y = Value, color = `Smoothing Parameter`)) +
  geom_point() +
  theme_minimal() +
  labs(title = "F1 scores for union_smoothing",
       x = "Quantile as Threshold",
       y = "F1 Score")
dev.off()

inter_snn_df <- as.data.frame(cbind(F1_inter_snn[[1]]$score, F1_inter_snn[[2]]$score, F1_inter_snn[[3]]$score, F1_inter_snn[[4]]$score))
colnames(inter_snn_df) <- c("a=0.2", "a=0.4", "a=0.6", "a=0.8")
inter_snn_df$Position <- 1:nrow(inter_snn_df)
inter_snn_df$unsmoothed <-F1_raw$score
df_long <- melt(inter_snn_df, id.vars = "Position", variable.name = "Smoothing Parameter", value.name = "Value")  

pdf("F1_inter_snn.pdf")
ggplot(df_long, aes(x = Position, y = Value, color = `Smoothing Parameter`)) +
  geom_point() +
  theme_minimal() +
  labs(title = "F1 scores for inter_smoothing",
       x = "Quantile as Threshold",
       y = "F1 Score")
dev.off()

inter_df <- as.data.frame(cbind(F1_inter[[1]]$score, F1_inter[[2]]$score, F1_inter[[3]]$score, F1_inter[[4]]$score))
colnames(inter_df) <- c("a=0.2", "a=0.4", "a=0.6", "a=0.8")
inter_df$Position <- 1:nrow(inter_df)
inter_df$unsmoothed <-F1_raw$score
df_long <- melt(inter_df, id.vars = "Position", variable.name = "Smoothing Parameter", value.name = "Value")  

pdf("F1_inter.pdf")
ggplot(df_long, aes(x = Position, y = Value, color = `Smoothing Parameter`)) +
  geom_point() +
  theme_minimal() +
  labs(title = "F1 scores for inter_smoothing",
       x = "Quantile as Threshold",
       y = "F1 Score")
dev.off()

nn_spatial_df <- as.data.frame(cbind(F1_nn_spatial[[1]]$score, F1_nn_spatial[[2]]$score, F1_nn_spatial[[3]]$score, F1_nn_spatial[[4]]$score, F1_nn_spatial[[5]]$score, F1_nn_spatial[[6]]$score, F1_nn_spatial[[7]]$score, F1_nn_spatial[[8]]$score, F1_nn_spatial[[9]]$score, F1_nn_spatial[[10]]$score, F1_nn_spatial[[11]]$score, F1_nn_spatial[[12]]$score, F1_nn_spatial[[13]]$score, F1_nn_spatial[[14]]$score, F1_nn_spatial[[15]]$score, F1_nn_spatial[[16]]$score))
colnames(nn_spatial_df) <- c("a1 = 0.2, a2=0.2","a1 = 0.2, a2=0.4", "a1 = 0.2, a2=0.6", "a1 = 0.2, a2=0.8", "a1 = 0.4, a2=0.2", "a1 = 0.4, a2=0.4", "a1 = 0.4, a2=0.6", "a1 = 0.4, a2=0.8", "a1 = 0.6, a2=0.2", "a1 = 0.6, a2=0.4", "a1 = 0.6, a2=0.6", "a1 = 0.6, a2=0.8", "a1 = 0.8, a2=0.2", "a1 = 0.8, a2=0.4", "a1 = 0.8, a2=0.6", "a1 = 0.8, a2=0.8")
nn_spatial_df$Position <- 1:nrow(nn_spatial_df)
nn_spatial_df$unsmoothed <-F1_raw$score
df_long <- melt(nn_spatial_df, id.vars = "Position", variable.name = "Smoothing Parameters", value.name = "Value")  


pdf("F1_nn_spatial.pdf")
ggplot(df_long, aes(x = Position, y = Value, color = `Smoothing Parameters`)) +
  geom_point() +
  theme_minimal() +
  labs(title = "F1 scores for nn_spatial_smoothing",
       x = "Quantile as Threshold",
       y = "F1 Score")
dev.off()

snn_spatial_df <- as.data.frame(cbind(F1_snn_spatial[[1]]$score, F1_snn_spatial[[2]]$score, F1_snn_spatial[[3]]$score, F1_snn_spatial[[4]]$score, F1_snn_spatial[[5]]$score, F1_snn_spatial[[6]]$score, F1_snn_spatial[[7]]$score, F1_snn_spatial[[8]]$score, F1_snn_spatial[[9]]$score, F1_snn_spatial[[10]]$score, F1_snn_spatial[[11]]$score, F1_snn_spatial[[12]]$score, F1_snn_spatial[[13]]$score, F1_snn_spatial[[14]]$score, F1_snn_spatial[[15]]$score, F1_snn_spatial[[16]]$score))
colnames(snn_spatial_df) <- c("a1 = 0.2, a2=0.2","a1 = 0.2, a2=0.4", "a1 = 0.2, a2=0.6", "a1 = 0.2, a2=0.8", "a1 = 0.4, a2=0.2", "a1 = 0.4, a2=0.4", "a1 = 0.4, a2=0.6", "a1 = 0.4, a2=0.8", "a1 = 0.6, a2=0.2", "a1 = 0.6, a2=0.4", "a1 = 0.6, a2=0.6", "a1 = 0.6, a2=0.8", "a1 = 0.8, a2=0.2", "a1 = 0.8, a2=0.4", "a1 = 0.8, a2=0.6", "a1 = 0.8, a2=0.8")
snn_spatial_df$Position <- 1:nrow(snn_spatial_df)
snn_spatial_df$unsmoothed <-F1_raw$score
df_long <- melt(snn_spatial_df, id.vars = "Position", variable.name = "Smoothing Parameters", value.name = "Value")  


pdf("F1_snn_spatial.pdf")
ggplot(df_long, aes(x = Position, y = Value, color = `Smoothing Parameters`)) +
  geom_point() +
  theme_minimal() +
  labs(title = "F1 scores for snn_spatial_smoothing",
       x = "Quantile as Threshold",
       y = "F1 Score")
dev.off()

alpha_df <- as.data.frame(cbind(F1_alpha[[1]]$score, F1_alpha[[2]]$score, F1_alpha[[3]]$score, F1_alpha[[4]]$score, F1_alpha[[5]]$score, F1_alpha[[6]]$score, F1_alpha[[7]]$score, F1_alpha[[8]]$score, F1_alpha[[9]]$score, F1_alpha[[10]]$score, F1_alpha[[11]]$score, F1_alpha[[12]]$score, F1_alpha[[13]]$score, F1_alpha[[14]]$score, F1_alpha[[15]]$score, F1_alpha[[16]]$score))
colnames(alpha_df) <- c("a = 0.2, alpha=0.2", "a = 0.2, alpha=0.4", "a = 0.2, alpha=0.6", "a = 0.2, alpha=0.8", "a = 0.4, alpha=0.2", "a = 0.4, alpha=0.4", "a = 0.4, alpha=0.6", "a = 0.4, alpha=0.8", "a = 0.6, alpha=0.2", "a = 0.6, alpha=0.4", "a = 0.6, alpha=0.6", "a = 0.6, alpha=0.8", "a = 0.8, alpha=0.2", "a = 0.8, alpha=0.4", "a = 0.8, alpha=0.6", "a = 0.8, alpha=0.8")
alpha_df$Position <- 1:nrow(alpha_df)
alpha_df$unsmoothed <-F1_raw$score
df_long <- melt(alpha_df, id.vars = "Position", variable.name = "Parameters", value.name = "Value")  

pdf("F1_alpha.pdf")
ggplot(df_long, aes(x = Position, y = Value, color = `Parameters`)) +
  geom_point() +
  theme_minimal() +
  labs(title = "F1 scores for alpha_smoothing",
       x = "Quantile as Threshold",
       y = "F1 Score")
dev.off()

best_df <- as.data.frame(cbind(F1_nn[[4]]$score, F1_snn[[4]]$score, F1_inter_snn[[4]]$score, F1_union[[4]]$score, F1_spatial[[4]]$score, F1_nn_spatial[[16]]$score, F1_alpha[[16]]$score, F1_snn_spatial[[16]]$score))
colnames(best_df) <- c("nn", "snn", "intersection", "union", "spatial", "nn_spatial", "alpha", "snn_spatial") 
best_df$Position <- 1:nrow(best_df)
best_df$unsmoothed <-F1_raw$score
df_long <- melt(best_df, id.vars = "Position", variable.name = "Method", value.name = "Value")  

pdf("F1_best.pdf")
ggplot(df_long, aes(x = Position, y = Value, color = `Method`)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Comparison of different smoothing methods for Visium-Dataset",
       x = "Quantile as Threshold",
       y = "F1 Score")
dev.off()


