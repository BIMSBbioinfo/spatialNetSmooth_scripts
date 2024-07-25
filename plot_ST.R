nn_spat <- read.csv("~/BA/images ST/best_thresholds_nn_spatial.csv", row.names=1)
raw <- read.csv("~/BA/images ST/best_thresholds_raw.csv", row.names=1)

smooth_nn_spat <- nn_spatial_smooth(seu, a1=0.4, a2=0.8)
smooth_nn_spat <- as.vector(smooth_nn_spat)
truth<- as.vector(truth)
pdf("~/BA/images ST/H1/plot_H1_roc.pdf")
plot_quant(smooth_nn_spat, coordinates, truth)
dev.off()

gsea_raw <- as.vector(gsea_raw)
pdf("~/BA/images ST/H1/plot_H1_raw_roc.pdf")
plot_quant(gsea_raw, coordinates, truth)
dev.off()

accuracy <- nn_spat[2,]
accuracy_raw <- raw[2,]
accuracy_comp <- rbind(accuracy_raw, accuracy)
rownames(accuracy_comp) <- c("unsmoothed", "smoothed")
write.csv(accuracy_comp, "~/BA/images ST/accuracy_comp.csv", row.names=T)
