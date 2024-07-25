library(tidyr)
library(reshape2)
library(VoltRon)
library(dplyr)
library(igraph)
library(GSEABase)
library(pROC)
setwd("~/spatialNetSmooth/")
# F1 nochmal machen 
devtools::load_all()
# count matrix
setwd("~/Documents")
letter <- "H1"
datax <- read.table(paste("./data/ST-cnts/", letter,".tsv.gz", sep=""), header = TRUE, sep = "\t")
entities <- datax$X
datax <- datax[,-1]
datax <- t(as.matrix(datax))
#datax <- LogNormalize(datax, scale.factor = 10000, verbose = TRUE)
# coords 
coords <- read.table(paste("./data/ST-spotfiles/",letter ,"_selection.tsv.gz", sep=""), header = TRUE, sep = "\t")
rownames(coords) <- paste(coords$x, coords$y, sep = "x")
coords2 <- read.table(
  file = paste("./data/meta/",letter, "_labeled_coordinates.tsv", sep=""),
  sep = "\t",
  header = TRUE,
  row.names = 1
)
coords2 <- coords2[order(coords2[,1],decreasing=FALSE),]
coords3 <- coords[order(coords[,3],decreasing=FALSE),]
truth <- array(coords2$label)


fun <- function(x){
  if(x=="invasive cancer"){x <- 1}
  else if(x=="cancer in situ"){x <- 1}
  else{x <- 0}
}
truth <- apply(truth, 1, fun)
truth <- as.matrix(truth, ncol=1)
#filtered_vec <- truth[rownames(truth) %in% entities]
rownames(truth) <- rownames(coords3)
truth <- as.matrix(truth[match(entities, rownames(truth))])
rownames(truth) <- entities

coords <- coords[entities,]
coords <- coords[,c("pixel_x", "pixel_y")]
colnames(coords) <- c("x", "y")
rownames(coords) <- paste("spot", (1:nrow(coords)), sep="_")
rownames(truth) <- rownames(coords)
colnames(datax) <- rownames(coords)
#coords <- read.csv("~/BA/D1_coords.csv", row.names=1)
# images
img <- magick::image_read(paste("./data/ST-imgs/HE/", letter,".jpg", sep=""))
img_info <- magick::image_info(img)


# get ST parameters
scale_param <- img_info$width/6200
params <- list(
  nearestpost.distance = (200*sqrt(2) + 50)*scale_param, # distance to nearest spot
  spot.radius = 50*(scale_param),
  vis.spot.radius = 100*(scale_param))

# make voltron object
temp <- formVoltRon(data = datax, image = img, coords = coords, assay.type = "spot", params = params)
#temp <- flipCoordinates(temp, assay="Custom_spot")
#saveRDS(temp, "~/BA/object_ST.Rds")
#temp <- readRDS("~/BA/object_ST.Rds")
setwd(paste("~/BA/images ST/", letter, sep=""))
seu <- gseaCalc(temp, assay="Custom_spot")
coordinates <- coords
#truth <- tumor
gsea_raw <- as.matrix(Metadata(seu)$gsea_rat_norm)
rownames(gsea_raw) <- rownames(Metadata(seu))
#gsea_raw <- cbind(truth, gsea_raw)
#colnames(gsea_raw) <- c("truth", "gsea")
#write.csv(gsea_raw, paste("~/BA/images ST/", letter, "/gsea_raw.csv", sep=""), row.names=T)
                
F1_raw <- F1_quant(gsea_raw, truth, gsea_raw)
pdf(paste("gsea_raw_",letter, ".pdf", sep=""))
plot_quant(gsea_raw, coordinates, truth, 8, gsea_raw)
dev.off()
alphas <- c(0.2, 0.4, 0.6, 0.8)

gsea_spatial <- vector(mode = "list", length = 4)
for (i in 1:4) {
  gsea_spatial[[i]] <- spatial_smooth(seu, a = alphas[i])
  
}
spatial <- as.data.frame(gsea_spatial)
colnames(spatial) <- c("alpha=0.2", "alpha=0.4", "alpha=0.6", "alpha=0.8")
write.csv(spatial, paste("spatial_scores_", letter, ".csv", sep=""), row.names = T)

roc_spatial <- vector(mode = "list", length = 4)
for (i in 1:4) {
  roc_spatial[[i]] <- roc_quant(truth, spatial[,i])
  
}

gsea_nn <- vector(mode = "list", length = 4)
for (i in 1:4) {
  gsea_nn[[i]] <- nn_smooth(seu, a = alphas[i])
  
}
nn <- as.data.frame(gsea_nn)
colnames(nn) <- c("alpha=0.2", "alpha=0.4", "alpha=0.6", "alpha=0.8")
write.csv(nn, paste("nn_scores_", letter, ".csv", sep=""), row.names = T)

roc_nn <- vector(mode = "list", length = 4)
for (i in 1:4) {
  roc_nn[[i]] <- roc_quant(truth, nn[,i])
  
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
write.csv(nn_spatial, paste("nn_spatial_scores_", letter, ".csv", sep=""), row.names = T)

roc_nn_spatial <- vector(mode = "list", length = 16)
for (i in 1:16) {
  roc_nn_spatial[[i]] <- f1_beta(truth, nn_spatial[,i])
  
}

gsea_union <- vector(mode = "list", length = 4)
for (i in 1:4) {
  gsea_union[[i]] <- union_smooth(seu, a = alphas[i])
  
}
union <- as.data.frame(gsea_union)
colnames(union) <- c("a=0.2", "a=0.4", "a=0.6", "a=0.8")
write.csv(union, paste("union_scores_", letter, ".csv", sep=""), row.names = T)

F1_union <- vector(mode = "list", length = 4)
for (i in 1:4) {
  F1_union[[i]] <- F1_quant(union[,i], truth, gsea_raw)
  
}

#gsea_inter <- vector(mode = "list", length = 4)
#for (i in 1:4) {
#  gsea_inter[[i]] <- inter_smooth(seu, a = alphas[i])
  
#}
#inter <- as.data.frame(gsea_inter)
#colnames(inter) <- c("a=0.2", "a=0.4", "a=0.6", "a=0.8")
#write.csv(inter, "inter_scores.csv", row.names = T)

#F1_inter <- vector(mode = "list", length = 4)
#for (i in 1:4) {
#  F1_inter[[i]] <- F1_quant(inter[,i], truth, gsea_raw)
#  
#}

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
write.csv(alpha, paste("alpha_scores_", letter, ".csv", sep=""), row.names = T)

F1_alpha <- vector(mode = "list", length = 16)
for (i in 1:16) {
  F1_alpha[[i]] <- f1_beta(alpha[,i], truth, gsea_raw)
  
}
nn_df <- as.data.frame(cbind(F1_nn[[1]]$score, F1_nn[[2]]$score, F1_nn[[3]]$score, F1_nn[[4]]$score))
colnames(nn_df) <- c("a=0.2", "a=0.4", "a=0.6", "a=0.8")
nn_df$Position <- 1:nrow(nn_df)
nn_df$unsmoothed <-F1_raw$score
df_long <- melt(nn_df, id.vars = "Position", variable.name = "Smoothing Parameter", value.name = "Value")  

pdf(paste("F1_snn_", letter, ".pdf"))
ggplot(df_long, aes(x = Position, y = Value, color = `Smoothing Parameter`)) +
  geom_point() +
  theme_minimal() +
  labs(title = "F1-beta scores for snn_smoothing",
       x = "Quantile as Threshold",
       y = "F1 Score")
dev.off()

spatial_df <- as.data.frame(cbind(F1_spatial[[1]]$score, F1_spatial[[2]]$score, F1_spatial[[3]]$score, F1_spatial[[4]]$score))
colnames(spatial_df) <- c("a=0.2", "a=0.4", "a=0.6", "a=0.8")
spatial_df$Position <- 1:nrow(spatial_df)
spatial_df$unsmoothed <-F1_raw$score

df_long <- melt(spatial_df, id.vars = "Position", variable.name = "Smoothing Parameter", value.name = "Value")  

pdf(paste("F1_spatial_", letter, ".pdf"))
ggplot(df_long, aes(x = Position, y = Value, color = `Smoothing Parameter`)) +
  geom_point() +
  theme_minimal() +
  labs(title = "F1-beta scores for spatial_smoothing",
       x = "Quantile as Threshold",
       y = "F1 Score")
dev.off()

union_df <- as.data.frame(cbind(F1_union[[1]]$score, F1_union[[2]]$score, F1_union[[3]]$score, F1_union[[4]]$score))
colnames(union_df) <- c("a=0.2", "a=0.4", "a=0.6", "a=0.8")
union_df$Position <- 1:nrow(union_df)
union_df$unsmoothed <-F1_raw$score
df_long <- melt(union_df, id.vars = "Position", variable.name = "Smoothing Parameter", value.name = "Value")  

pdf(paste("F1_union_", letter, ".pdf"))
ggplot(df_long, aes(x = Position, y = Value, color = `Smoothing Parameter`)) +
  geom_point() +
  theme_minimal() +
  labs(title = "F1-beta scores for union_smoothing",
       x = "Quantile as Threshold",
       y = "F1 Score")
dev.off()

#inter_df <- as.data.frame(cbind(F1_inter[[1]]$score, F1_inter[[2]]$score, F1_inter[[3]]$score, F1_inter[[4]]$score))
#colnames(inter_df) <- c("a=0.2", "a=0.4", "a=0.6", "a=0.8")
#inter_df$Position <- 1:nrow(inter_df)
#inter_df$unsmoothed <-F1_raw$score
#df_long <- melt(inter_df, id.vars = "Position", variable.name = "Smoothing Parameter", value.name = "Value")  

#pdf("F1_inter.pdf")
#ggplot(df_long, aes(x = Position, y = Value, color = `Smoothing Parameter`)) +
#  geom_point() +
#  theme_minimal() +
#  labs(title = "F1 scores for inter_smoothing",
#       x = "Quantile as Threshold",
#       y = "F1 Score")
#dev.off()

nn_spatial_df <- as.data.frame(cbind(F1_nn_spatial[[1]]$score, F1_nn_spatial[[2]]$score, F1_nn_spatial[[3]]$score, F1_nn_spatial[[4]]$score, F1_nn_spatial[[5]]$score, F1_nn_spatial[[6]]$score, F1_nn_spatial[[7]]$score, F1_nn_spatial[[8]]$score, F1_nn_spatial[[9]]$score, F1_nn_spatial[[10]]$score, F1_nn_spatial[[11]]$score, F1_nn_spatial[[12]]$score, F1_nn_spatial[[13]]$score, F1_nn_spatial[[14]]$score, F1_nn_spatial[[15]]$score, F1_nn_spatial[[16]]$score))
colnames(nn_spatial_df) <- c("a1 = 0.2, a2=0.2","a1 = 0.2, a2=0.4", "a1 = 0.2, a2=0.6", "a1 = 0.2, a2=0.8", "a1 = 0.4, a2=0.2", "a1 = 0.4, a2=0.4", "a1 = 0.4, a2=0.6", "a1 = 0.4, a2=0.8", "a1 = 0.6, a2=0.2", "a1 = 0.6, a2=0.4", "a1 = 0.6, a2=0.6", "a1 = 0.6, a2=0.8", "a1 = 0.8, a2=0.2", "a1 = 0.8, a2=0.4", "a1 = 0.8, a2=0.6", "a1 = 0.8, a2=0.8")
nn_spatial_df$Position <- 1:nrow(nn_spatial_df)
nn_spatial_df$unsmoothed <-F1_raw$score
df_long <- melt(nn_spatial_df, id.vars = "Position", variable.name = "Smoothing Parameters", value.name = "Value")  


pdf(paste("F1_snn_spatial_", letter, ".pdf"))
ggplot(df_long, aes(x = Position, y = Value, color = `Smoothing Parameters`)) +
  geom_point() +
  theme_minimal() +
  labs(title = "F1-beta scores for snn_spatial_smoothing",
       x = "Quantile as Threshold",
       y = "F1 Score")
dev.off()

alpha_df <- as.data.frame(cbind(F1_alpha[[1]]$score, F1_alpha[[2]]$score, F1_alpha[[3]]$score, F1_alpha[[4]]$score, F1_alpha[[5]]$score, F1_alpha[[6]]$score, F1_alpha[[7]]$score, F1_alpha[[8]]$score, F1_alpha[[9]]$score, F1_alpha[[10]]$score, F1_alpha[[11]]$score, F1_alpha[[12]]$score, F1_alpha[[13]]$score, F1_alpha[[14]]$score, F1_alpha[[15]]$score, F1_alpha[[16]]$score))
colnames(alpha_df) <- c("a = 0.2, alpha=0.2", "a = 0.2, alpha=0.4", "a = 0.2, alpha=0.6", "a = 0.2, alpha=0.8", "a = 0.4, alpha=0.2", "a = 0.4, alpha=0.4", "a = 0.4, alpha=0.6", "a = 0.4, alpha=0.8", "a = 0.6, alpha=0.2", "a = 0.6, alpha=0.4", "a = 0.6, alpha=0.6", "a = 0.6, alpha=0.8", "a = 0.8, alpha=0.2", "a = 0.8, alpha=0.4", "a = 0.8, alpha=0.6", "a = 0.8, alpha=0.8")
alpha_df$Position <- 1:nrow(alpha_df)
alpha_df$unsmoothed <-F1_raw$score
df_long <- melt(alpha_df, id.vars = "Position", variable.name = "Parameters", value.name = "Value")  

pdf(paste("F1_alpha_", letter, ".pdf"))
ggplot(df_long, aes(x = Position, y = Value, color = `Parameters`)) +
  geom_point() +
  theme_minimal() +
  labs(title = "F1-beta scores for alpha_smoothing",
       x = "Quantile as Threshold",
       y = "F1 Score")
dev.off()

best_df <- as.data.frame(cbind(F1_nn[[4]]$score, F1_union[[4]]$score, F1_spatial[[4]]$score, F1_nn_spatial[[16]]$score, F1_alpha[[16]]$score))#removed F1_inter[[4]]$score
colnames(best_df) <- c("snn", "union", "spatial", "snn_spatial", "alpha") #"intersection",
best_df$Position <- 1:nrow(best_df)
best_df$unsmoothed <-F1_raw$score
df_long <- melt(best_df, id.vars = "Position", variable.name = "Method", value.name = "Value")  

pdf(paste("F1_best_", letter, ".pdf"))
ggplot(df_long, aes(x = Position, y = Value, color = `Method`)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Comparison of F1-beta scores of different smoothing methods",
       x = "Quantile as Threshold",
       y = "F1 Score")
dev.off()
pdf(paste("plot_snn_", letter, ".pdf", sep=""))
plot_quant(gsea_nn[[4]], coordinates, truth, 8, gsea_raw)
dev.off()
pdf(paste("plot_spatial_", letter, ".pdf", sep=""))
plot_quant(gsea_spatial[[4]], coordinates, truth, 8, gsea_raw)
dev.off()
pdf(paste("plot_union_", letter, ".pdf", sep=""))
plot_quant(gsea_union[[4]], coordinates, truth, 8, gsea_raw)
dev.off()
#pdf("plot_inter.pdf")
#plot_quant(gsea_inter[[4]], coordinates, truth, 3, gsea_raw)
#dev.off()

pdf(paste("plot_alpha_", letter, ".pdf", sep=""))
plot_quant(alpha[,16], coordinates, truth, 8, gsea_raw)
dev.off()
pdf(paste("plot_snn_spatial", letter, ".pdf", sep=""))
plot_quant(nn_spatial[,16], coordinates, truth, 7, gsea_raw)
dev.off()

pdf("plot_nn_spatial_2.pdf")
plot_quant(nn_spatial[,12], coordinates, truth, 5, gsea_raw)
dev.off()

#compare_df <- as.data.frame(cbind(F1_nn_spatial_2[[12]]$score, F1_nn_spatial[[16]]$score))
#colnames(compare_df) <- c("spatial first", "nn_first")
#compare_df$Position <- 1:nrow(compare_df)
#compare_df$unsmoothed <-F1_raw$score
#df_long <- melt(compare_df, id.vars = "Position", variable.name = "Method", value.name = "Value")  

#pdf("F1_compare.pdf")
#ggplot(df_long, aes(x = Position, y = Value, color = `Method`)) +
#  geom_point() +
#  theme_minimal() +
#  labs(title = "Comparing order of spatial and nn",
#       x = "Quantile as Threshold",
#       y = "F1 Score")
#dev.off()
