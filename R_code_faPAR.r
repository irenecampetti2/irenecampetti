# how to look at chemical cycles from satellites

library(raster)
library(rasterVis)
library(rasterdiv)

plot(copNDVI)
copNDVI <- reclassify(copNDVI,cbind(253:255, NA))
levelplot(copNDVI)

setwd("C:/lab/")

# import the data
faPAR10 <- raster("faPAR10.tif")
levelplot(faPAR10)

# save as a .pdf
pdf("copNDVI.pdf")
levelplot(copNDVI)
dev.off()
pdf("faPAR10.pdf")
levelplot(faPAR10)
dev.off()

####
#regression model between faPAR and NDVI
# amount of heavy metals and erosion
erosion <- c(12,14,16,24,26,40,55,67)
hm <- c(30,100,150,200,260,340,460,600)

plot(erosion, hm, col="red", pch=19, xlab="erosion", ylab="heavy metals", cex=2)
model1 <- lm(hm ~ erosion)
summary(model1)
abline(model1)

library(raster)
library(rasterdiv)
libray(sf)
setwd("C:/lab/")
raster("faPAR10.tif")

plot(faPAR10)
plot(copNDVI)

copNDVI <- reclassify(copNDVI,cbind(253:255,NA),right=TRUE)
random.points <- function(raster,n)
{
lin <- rasterToContour(is.na(raster))
pol <- as(st_union(st_polygonize(st_as_sf(lin))), 'Spatial') # st_union to dissolve geometries
pts <- spsample(pol[1,], n, type = 'random')
}

pts <- random.points(faPAR10,1000)
copNDVIp <- extract(copNDVI, pts)
faPAR10p <- extract(faPAR10,pts)

