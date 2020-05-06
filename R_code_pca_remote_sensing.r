setwd("C:/lab/")

library(raster)
library(RStoolbox)

# use brick for satellite images
p224r63_2011 <- brick("p224r63_2011_masked.grd")


# b1: blue
# b2: green
# b3: red
# b4: NIR
# b5: SWIR
# b6: thermal infrared
# b7: SWIR
# b8: panchromatic

plotRGB(p224r63_2011, r=5, g=4, b=3, stretch="Lin")

# with ggplot
library("ggplot2")
ggRGB(p224r63_2011,5,4,3)

# do the same for 1988 image
p224r63_1988 <- brick("p224r63_1988_masked.grd")
plotRGB(p224r63_1988, r=5, g=4, b=3, stretch="Lin")

# plot the two images one beside the other
par(mfrow=c(1,2))
plotRGB(p224r63_1988, r=5, g=4, b=3, stretch="Lin")
plotRGB(p224r63_2011, r=5, g=4, b=3, stretch="Lin")

names(p224r63_2011)

# check correlation
dev.off()
plot(p224r63_2011$B1_sre, p224r63_2011$B3_sre)

# PCA
# dicrease the resolution
p224r63_2011_res <- aggregate(p224r63_2011, fact=10)
p224r63_2011_pca <- rasterPCA(p224r63_2011_res)

# check the properties of the pca
p224r63_2011_pca

#link the map to the pca
cl <- colorRampPalette(c('dark grey','grey','light grey'))(100)
plot(p224r63_2011_pca$map, col=cl)

summary(p224r63_2011_pca$model)
pairs(p224r63_2011)

# plotRGB, the names of the component are found in the map properties
plotRGB(p224r63_2011_pca$map, r=1,g=2,b=3, stretch="Lin")

# do the procedure for p224r63_1988
p224r63_1988_res <- aggregate(p224r63_1988, fact=10)
p224r63_1988_pca <- rasterPCA(p224r63_1988_res)
plot(p224r63_1988_pca$map)
summary(p224r63_1988_pca$model)
pairs(p224r63_1988)

# plot the difference in all bands: difference in the pca
difpca <- p224r63_2011_pca$map - p224r63_1988_pca$map
plot(difpca)

cldif <- colorRampPalette(c('blue','black','yellow'))(100)
plot(difpca$PC1,col=cldif)
