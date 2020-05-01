# R code to view biomass over the world map and calculate changes in ecosystem functions and services

# We are going to use Raster
# Raster is the format of pixels - all data based on pixels

install.packages("rasterdiv") 
install.packages("rasterVis") 
library(rasterdiv)
library(rasterVis)

# Attach the dataset copNDVI
data(copNDVI)

# Plot the dataset
plot(copNDVI)

# We can now make additional graphs on the world map. We should remove some values by using cbind. Remove the values from the 253 to 255
# Put this data as NA value
# Reclassify the data into copNDVI
copNDVI<-reclassify(copNDVI,cbind(253:255,NA)

levelplot(copNDVI)

# Lets change the grain of the images: let's aggregate the NDVI
# The function is: aggregate()
# Increase the pixle dimension by factor of 10: fact=10
# Let's give a new name to the image copNDVI10<-
copNDVI10 <- aggregate(copNDVI, fact=10)

# fact=100 blurry image
copNDVI100<-aggregate(copNDVI,fact=100)


##################
library(ggplot2)

myPalette <- colorRampPalette(c('white','green','dark green'))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(1, 8))

ggR(copNDVI, geom_raster = TRUE) +
scale_fill_gradientn(name = "NDVI", colours = myPalette(100))+
labs(x="Longitude",y="Latitude", fill="")+
#   theme(legend.position = "bottom") +
  NULL
# +
# ggtitle("NDVI")
#################

setwd("C:/lab/")
library(raster)

brick("defor1_.jpg")
defor1<-brick("defor1_.jpg")

brick("defor2_.jpg")
defor2<-brick("defor2_.jpg")

# Bands of landset 
# First band B1: NIR
# Second band B2: RED
# Third band B3: GREEN
plotRGB(defor1, r=1, g=2, b=3, stretch="Lin")
plotRGB(defor2, r=1, g=2, b=3, stretch="Lin")
par(mfrow=c(1,2))

dvi1<-defor1$defor1_.1-defor1$defor1_.2
dvi2<-defor2$defor2_.1-defor2$defor2_.2
# Change the color ramp palletes of the plots
# Assign the color ramp pallet to the name: cl
cl<-colorRampPalette(c("darkblue","yellow","red","black"))(100) 
par(mfrow=c(1,2))

plot(dvi1,col=cl)
plot(dvi2,col=cl)

difdvi<-dvi1-dvi2

dev.off()

# Change the color ramp palletes of the plots. Assign the color ramp pallet to the name: cld
cld<-colorRampPalette(c("blue","white","red"))(100) 

plot(difdvi, col=cld)

# Lets make us of a histogram stretch: "hist"
# This function enhances the 'noise'. 
hist(difdvi)
