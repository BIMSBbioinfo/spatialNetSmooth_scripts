library(imager)

file_annot = './data/spatial/tissue_lowres_image_annotated.png'
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
