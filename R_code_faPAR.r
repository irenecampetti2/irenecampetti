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
