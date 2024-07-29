library(tidyr)
library(reshape2)
library(VoltRon)
library(dplyr)
library(igraph)
library(GSEABase)
library(pROC)

#load package (if not installed)
setwd("~/spatialNetSmooth/")
devtools::load_all()
#Set input and output directories
filepath_data= "~/Documents/data/"
filepath_output = "~/BA/images ST/"

#select name of dataset
letter <- "C1"

#read count matrix
datax <- read.table(paste(filepath_data, "/ST-cnts/", letter,".tsv.gz", sep=""), header = TRUE, sep = "\t")
entities <- datax$X
datax <- datax[,-1]
datax <- t(as.matrix(datax))
#datax <- LogNormalize(datax, scale.factor = 10000, verbose = TRUE)
# coords 
#read coordinates
coords <- read.table(paste(filepath_data, "/ST-spotfiles/",letter ,"_selection.tsv.gz", sep=""), header = TRUE, sep = "\t")
rownames(coords) <- paste(coords$x, coords$y, sep = "x")
#get true labels
coords2 <- read.table(
  file = paste(filepath_data, "/meta/",letter, "_labeled_coordinates.tsv", sep=""),
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
img <- magick::image_read(paste(filepath_data,"/ST-imgs/HE/", letter,".jpg", sep=""))
img_info <- magick::image_info(img)


# get ST parameters
scale_param <- img_info$width/6200
params <- list(
  nearestpost.distance = (200*sqrt(2) + 50)*scale_param, # distance to nearest spot
  spot.radius = 50*(scale_param),
  vis.spot.radius = 100*(scale_param))

# make voltron object
temp <- formVoltRon(data = datax, image = img, coords = coords, assay.type = "spot", params = params)

setwd(paste(filepath_output, letter, sep=""))
write.csv(truth, "truth.csv", row.names=T)

#caclulate GSEA
seu <- gseaCalc(temp, assay="Custom_spot")
coordinates <- coords
gsea_raw <- as.matrix(Metadata(seu)$gsea_rat_norm)
rownames(gsea_raw) <- rownames(Metadata(seu))
#gsea_raw <- cbind(truth, gsea_raw)
#colnames(gsea_raw) <- c("truth", "gsea")
#write.csv(gsea_raw, paste("~/BA/images ST/", letter, "/gsea_raw.csv", sep=""), row.names=T)
                

alphas <- c(0.2, 0.4, 0.6, 0.8)

#spatial smoothing
gsea_spatial <- vector(mode = "list", length = 4)
for (i in 1:4) {
  gsea_spatial[[i]] <- spatial_smooth(seu, a = alphas[i])
  
}
spatial <- as.data.frame(gsea_spatial)
colnames(spatial) <- c("alpha=0.2", "alpha=0.4", "alpha=0.6", "alpha=0.8")
write.csv(spatial, paste("spatial_scores_", letter, ".csv", sep=""), row.names = T)


#NN smoothing
gsea_nn <- vector(mode = "list", length = 4)
for (i in 1:4) {
  gsea_nn[[i]] <- nn_smooth(seu, a = alphas[i])
  
}
nn <- as.data.frame(gsea_nn)
colnames(nn) <- c("alpha=0.2", "alpha=0.4", "alpha=0.6", "alpha=0.8")
write.csv(nn, paste("nn_scores_", letter, ".csv", sep=""), row.names = T)


#SNN spatial smoothing
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

#union smoothing
gsea_union <- vector(mode = "list", length = 4)
for (i in 1:4) {
  gsea_union[[i]] <- union_smooth(seu, a = alphas[i])
  
}
union <- as.data.frame(gsea_union)
colnames(union) <- c("a=0.2", "a=0.4", "a=0.6", "a=0.8")
write.csv(union, paste("union_scores_", letter, ".csv", sep=""), row.names = T)



#linear combination smoothing
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


#select best parameters for each method by selecting the column in dataframe
#make spatial plots

truth <- as.vector(truth)
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



