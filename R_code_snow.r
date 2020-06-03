# R_code_snow.r


setwd("C:/lab/")

install.packages("ncdf4")
library(ncdf4)
library(raster)

snowmay <- raster("c_gls_SCE_202005260000_NHEMI_VIIRS_V1.0.1.nc")
cl <- colorRampPalette(c('darkblue','blue','light blue'))(100) 

# Exercise: plot snow cover with the cl palette
plot(snowmay,col=cl)

##### import snow data

setwd("C:/lab/snow")

snow2000 <- raster("snow2000r.tif")
snow2005 <- raster("snow2005r.tif")
snow2010 <- raster("snow2010r.tif")
snow2015 <- raster("snow2015r.tif")
snow2020 <- raster("snow2020r.tif")

par(mfrow=c(2,3))
plot(snow2000,col=cl)
plot(snow2005,col=cl)
plot(snow2010,col=cl)
plot(snow2015,col=cl)
plot(snow2020,col=cl)

# how to plot snow data with lapply
# first of all we make a list of the files we are going to import
rlist <- list.files(pattern="snow")
rlist
import <- lapply(rlist,raster)
snow.multitemp <- stack(import)
plot(snow.multitemp,col=cl)

# let's make a prediction on how the snow is expected to change in 2025
source("prediction.r")

#######################

# set working directory
setwd("C:/lab/snow/")

# Exercise: import the snow cover images altogether
library(raster)
rlist <- list.files(pattern="snow")
import <- lapply(rlist,raster)
snow.multitemp <- stack(import)
cl <- colorRampPalette(c('darkblue','blue','light blue'))(100) 
plot(snow.multitemp, col=cl)

# load name.RData
# import the prediction image

prediction <- raster("predicted.2025.norm.tif")
plot(prediction, col=cl)

# export the output (you made the calculation and you want to send the output to a colleague)
writeRaster(prediction, "final.tif")

# final stack
final.stack <- stack(snow.multitemp, prediction)
plot(final.stack,col=cl)

# export the R graph
pdf("my_final_exciting_graph.pdf")
plot(final.stack, col=cl)
dev.off()

# make a .png of the project
png("my_final_exciting_graph.pdf")
plot(final.stack, col=cl)
dev.off()
