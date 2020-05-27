library(raster)
library(RStoolbox) # this is for PCA

setwd("C:/lab/") 
snt <- brick("snt_r10.tif")

plot(snt)

# B1 blue
# B2 green
# B3 red
# B4 NIR

# R3 G2 B1
plotRGB(snt,3,2,1, stretch="lin") 
plotRGB(snt,4,3,2, stretch="lin") 

pairs(snt)

### PCA analysis
sntpca <- rasterPCA(snt)
sntpca

summary(sntpca$model)
# 70% of information
plot(sntpca$map) 

plotRGB(sntpca$map, 1, 2, 3, stretch="lin")

# set the moving window
window <- matrix(1, nrow = 5, ncol = 5)
window

sd_snt <- focal(sntpca$map$PC1, w=window, fun=sd)

cl <- colorRampPalette(c('dark blue','green','orange','red'))(100) # 
plot(sd_snt, col=cl)

par(mfrow=c(1,2))
plotRGB(snt,4,3,2, stretch="lin", main="original image") 
plot(sd_snt, col=cl, main="diversity")
