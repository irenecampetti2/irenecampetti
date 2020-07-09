# R_code_exam.r

# 1. R code first

# Install the package sp: classes and metods for spatial data
# "" or '' is used for objects outside R
install.packages("sp")

# recall the package sp
library(sp)

# load dataset
data(meuse)

# let's see how the meuse dataset is structured:
meuse

# let's look at the first rows of the dataset:
head(meuse)

# let's plot two variables together
# let's see if zinc concentration is related to that of copper
attach(meuse)
plot(zinc,copper)

# change the color to green -- with col="" or ''
plot(zinc,copper,col="green")
# change the symbol in the plot -- with pch
plot(zinc,copper,col="green",pch=19)
# change symbol's size -- with cex
plot(zinc,copper,col="green",pch=19,cex=2)

##################################
##################################

# 2. R code multipanel

# Install GGally package, used for the function ggpairs()
install.packages("GGally")

# Require the packages we are going to use
library(sp)
library(GGally)

data(meuse) # there is a dataset available named meuse
attach(meuse)

# Exercise: see the names of the variables and plot cadmium versus zinc
# There are two ways to se the names of the variables:
head(meuse) # shows first 6 lines of the dataset
# or View(meuse)
# or meuse
plot(cadmium,zinc)

# add a change in symbol, color, size
plot(cadmium,zinc,pch=15,col="red",cex=2)

# Exercise: make make all the possible pairwise plots of the dataset
# Instead of plotting each pair separately -- plot(x,zinc); plot(x,copper) ecc.. we can use the function pairs()
pairs(meuse) # the output is a matrix of plots

# Use the formula, starting with the symbol ~ . Each term will give a separate variable in the pairs plot
pairs(~ cadmium + copper + lead + zinc, data=meuse)

pairs(meuse[,3:6])

# Exercise: prettify the graph
pairs(meuse[,3:6],pch=8,col="blue",cex=0.5)

# Use ggpairs() function in order to obtain better graphs
ggpairs(meuse[,3:6])

##################################
##################################

# 3. R code spatial
# R code for spatial view of points

library(sp)

data(meuse)

head(meuse)

# set spatial coordinates
coordinates(meuse) = ~x+y

plot(meuse)

# plot zinc concentration, spplot function uses different layers?
spplot(meuse, "zinc")

# Exercise: plot the spatial amount of copper
# main -- allows to add a title
spplot(meuse, "copper", main="Copper concentration")

# create a bubble plot of spatial data, then add title
bubble(meuse, "zinc")
bubble(meuse, "zinc", main="Zinc concentration")

# Exercise: bubble the copper in red
bubble(meuse, "copper", main="Copper concentration", col="red")

### Importing new data

# put the covid_agg.csv file into the lab folder 

# setting the working directory: lab
# for windows
setwd("C:/lab/")

# import the dataset and read it in a table format
# it is possible to name the dataset by --> name <-
covid <- read.table("covid_agg.csv", head=T)

head(covid)

attach(covid)
plot(country, cases)

# plot (covid$country, covid$cases) -- with different orientation of the tick mark labels
plot(country, cases, las=0) # parallel labels
plot(country, cases, las=1) # horizontal labels
plot(country, cases, las=2) # perpendicular labels
plot(country, cases, las=3) # vertical labels

plot(country, cases, las=3, cex.axis=0.5) # cex.axis -- change the size of text for x axis, 50% smaller
plot(country, cases, las=3, cex.axis=0.7)

# install ggplot2 package
install.packages("ggplot2")

# set the working directory
setwd("C:/lab/")

# load dataset
load("spatial.RData")
# list objects in the dataset
ls()

library(ggplot2)  # also require(ggplot2) can be used

# load mpg data
data(mpg)
# show first rows of the dataset
head(mpg)

# to make graphs with ggplot we need three key components: dataset, aesthetich mapping (the variables we want to put in the graph), geometry
ggplot(mpg, aes(x=displ, y=hwy)) + geom_point()
ggplot(mpg, aes(x=displ, y=hwy)) + geom_line()
ggplot(mpg, aes(x=displ, y=hwy)) + geom_polygon()

# in our case the dataset is covid, we want to relate country x cases, and we are using point as geometry
head(covid)
ggplot(covid, aes(x=lon, y=lat, size=cases)) + geom_point()

##################################
##################################

# 4. R code for point pattern analysis: density map

setwd("C:/lab/")

# Install spatstat package and require it
install.packages("spatstat")
library(spatstat)

# Import data, head=T means that there is a column header
covid <- read.table("covid_agg.csv",head=T)

attach(covid)
head(covid)

# Point pattern dataset represented in two dimensions
# The vectors represent longitude and latitude
covids <- ppp(lon, lat, c(-180, 180), c(-90, 90)) # if not using attach() --> ppp(covid$lon, covid$lat, c(-180,180), c(-90,90))

# Estimate density
d <- density(covids)
plot(d)

# Show points at specific coordinates
points(covids)

setwd("C:/lab/")

# Load previous workspace
load(".RData")
ls()

# covids: point pattern
# d: density map

library(spatstat)
plot(d)
points(covids)

# Install and require rgdal package
install.packages("rgdal")
library(rgdal)

# Let's input vector lines (x0y0, x1y1, x2,y2...)
# Import coastlines
coastlines <- readOGR("ne_10m_coastline.shp")

# Add coastlines to the density map
plot(d)
points(covids)
plot(coastlines, add=T) # add=T adds to current plot

# Change the color of the graph
cl <- colorRampPalette(c("yellow", "orange", "red")) (100)
plot(d, col=cl)
points(covids)
plot(coastlines, add=T)

# Exercise: make a colorRampPalette
clr <- colorRampPalette(c("light green", "yellow","orange","violet")) (100)
plot(d, col=clr, main="Densities of covid-19") # main adds a title of the graph
points(covids)
plot(coastlines, add=T)

# New example: abrupt change of colors!
cll <- colorRampPalette(c("light green", "yellow","orange","violet")) (5) # lower number of intermediate colors
plot(d, col=cll, main="Densities of covid-19")
points(covids)
plot(coastlines, add=T)

# Use a higher number of intermediate colors
clr <- colorRampPalette(c("light green", "yellow","orange","violet")) (1000) 
plot(d, col=clr, main="Densities of covid-19")
points(covids)
plot(coastlines, add=T)

# Export as .pdf
pdf("covid_density.pdf")
cl <- colorRampPalette(c("light green", "yellow","orange","violet")) (100)
plot(d, col=cl, main="Densities of covid-19")
points(covids)
plot(coastlines, add=T)
dev.off() # closes the current plot

# Export as .png
png("covid_density.png")
clr <- colorRampPalette(c("light green", "yellow","orange","violet")) (100)
plot(d, col=clr, main="Densities of covid-19")
points(covids)
plot(coastlines, add=T)
dev.off()

##################################
##################################

# 5. R code for multivariate analysis

# Install vegan package for vegetation analysis, then require it
install.packages("vegan")
library(vegan)

# Set working directory for Windows
setwd("C:/lab/")

# Read the dataset and name it biomes
biomes <- read.table("biomes.csv", header=T, sep=",") # header=T tells R that the first row contains the names of the variables of the table
# sep is the separator character

head(biomes) # View(biomes), biomes are also good

# DEtrended CORrespondence ANAlysis decorana
multivar <- decorana(biomes) 
multivar
plot(multivar)

# Biomes types
biomes_types <- read.table("biomes_types.csv", header=T, sep=",")
head(biomes_types)
attach(biomes_types)

ordiellipse(multivar, type, col=1:4, kind = "ehull", lwd=3) # Also col="blue", "red", "green", "black")
# draws lines or polygons for dispersion

ordispider(multivar, type, col=1:4, label= TRUE) # draws a spider diagram in which each point is connected to the group centroid

##################################
##################################

# 6. R code for remote sensing data analysis

# install the packages raster and RStoolbox
install.packages("raster")
install.packages("RStoolbox")

setwd("C:/lab/")

library(raster)

# Create a multi-layer raster object using brick() from a multi-layer file -- RGB colors
p224r63_2011 <- brick("p224r63_2011_masked.grd")

plot(p224r63_2011)

# Bands of Landset
# B1: blue
# B2: green
# B3: red
# B4: NIR

# Create a different color palette and name it cl
cl <- colorRampPalette(c('black','grey','light grey'))(100) # (100) is for the saturation of colors 
plot(p224r63_2011, col=cl)

# Multiframe of different plots
# par() is used to put several graphs in a single plot, mfrow defines the matrix
par(mfrow=c(2,2))

# B1: blue
clb <- colorRampPalette(c('dark blue','blue','light blue'))(100)
plot(p224r63_2011$B1_sre, col=clb)

# Exercise: do the same with the green band B2
clg <- colorRampPalette(c('dark green','green','light green'))(100)
plot(p224r63_2011$B2_sre, col=clg)

# B3: red
clr <- colorRampPalette(c('dark red','red','pink'))(100)   # light red does not exist
plot(p224r63_2011$B3_sre, col=clr)

# B4: NIR
cln <- colorRampPalette(c('red','orange','yellow'))(100)
plot(p224r63_2011$B4_sre, col=cln)


# plotRGB() plots three layers representing the different bandwidths of the electromagnetic spectrum (r,g,b) 
plotRGB(p224r63_2011, r=3, g=2, b=1, stretch="Lin") # stretch increases the contrast of the image, linear stretch
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin")

# Exercise: put NIR on top of the G component
plotRGB(p224r63_2011, r=3, g=4, b=2, stretch="Lin")

#######

setwd("C:/lab/")

# Load the saved workspace
load("rs.RData")
ls()

library(raster)

# Use 1988 data to make a multitemporal comparison
# masked in the image is because there is no data where there was water
# brick() imports the whole package of different bands
p224r63_1988 <- brick("p224r63_1988_masked.grd")
plot(p224r63_1988)

# Exercise: plot in visible RGB 321 both images
p224r63_2011 <- brick("p224r63_2011_masked.grd")
par(mfrow=c(2,1))
plotRGB(p224r63_1988, r=3, g=2, b=1, stretch="Lin")
plotRGB(p224r63_2011, r=3, g=2, b=1, stretch="Lin")

# Exercise: plot in false colour 432 both images
# 4 is for NIR
par(mfrow=c(2,1))
plotRGB(p224r63_1988, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin")

# Enhance the noise
par(mfrow=c(2,1))
plotRGB(p224r63_1988, r=4, g=3, b=2, stretch="Hist") # histogram stretch
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Hist")
# higher noise is observed in 1988, due to humidity


# B1: blue
# B2: green
# B3: red: B3_sre (spectrum reflectance)
# B4: NIR: B4_sre

# Vegetation Index for 2011
dvi2011 <- p224r63_2011$B4_sre - p224r63_2011$B3_sre # $ symbol links elements in R
cl <- colorRampPalette(c('darkorchid3','light blue','lightpink4'))(100)
plot(dvi2011, col=cl)

# Exercise: dvi for 1988
dvi1988 <- p224r63_1988$B4_sre - p224r63_1988$B3_sre
clb <- colorRampPalette(c('darkorchid3','light blue','lightpink4'))(100)
plot(dvi1988, col=clb)

# Difference from one year to another
diff <- dvi2011 - dvi1988
plot(diff)

# Changing grain (dimension of the pixels)
# aggregate() image + factor, which is how much we want to increase the pixels dimension
p224r63_2011res <- aggregate(p224r63_2011, fact=10) # from a pixel of 30 m to a pixel of 300 m
p224r63_2011res100 <- aggregate(p224r63_2011, fact=100)

# Plot by RGB them altogether, the first image and the ones with factor of 10 and factor of 100
par(mfrow=c(3,1))
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011res, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011res100, r=4, g=3, b=2, stretch="Lin")

# By entering the name of the images you get the informations
p224r63_2011
p224r63_2011res
p224r63_2011res100
 
##################################
##################################

# 7. R code for ecosystem functions

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

levelplot(copNDVI) # draws contour plots

# Lets change the grain of the images: let's aggregate the NDVI
# Increase the pixel dimension by factor of 10
# Let's give a new name to the image: copNDVI10
copNDVI10 <- aggregate(copNDVI, fact=10)

# fact=100 blurry image
copNDVI100<-aggregate(copNDVI,fact=100)

##################

setwd("C:/lab/")
library(raster)

# Import RGB images and name them defor1 and defor2
defor1<-brick("defor1_.jpg")
defor2<-brick("defor2_.jpg")

# Bands of landset 
# B1: NIR
# B2: RED
# B3: GREEN
plotRGB(defor1, r=1, g=2, b=3, stretch="Lin")
plotRGB(defor2, r=1, g=2, b=3, stretch="Lin")
                    
par(mfrow=c(1,2))
plotRGB(defor1, r=1, g=2, b=3, stretch="Lin")
plotRGB(defor2, r=1, g=2, b=3, stretch="Lin")

dvi1<-defor1$defor1_.1-defor1$defor1_.2
dvi2<-defor2$defor2_.1-defor2$defor2_.2
                    
# Change the color ramp palletes of the plots
# Assign the color ramp pallet to the name cl
cl<-colorRampPalette(c("darkblue","yellow","red","black"))(100) 
par(mfrow=c(1,2))
plot(dvi1,col=cl)
plot(dvi2,col=cl)

difdvi<- dvi1-dvi2

dev.off()

# Change the color ramp palletes of the plots. Assign the color ramp pallet to the name cld
cld<- colorRampPalette(c("blue","white","red"))(100) 
plot(difdvi, col=cld)

# Use a histogram stretch: "hist"
# This function enhances the noise. 
hist(difdvi)

##################################
##################################
                    
# 8. R code pca remote sensing
                    
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
# Plot RGB
plotRGB(p224r63_2011, r=5, g=4, b=3, stretch="Lin")

# with ggplot, ggRGB() creates ggplot2 plot with RGB
library("ggplot2")
ggRGB(p224r63_2011,5,4,3)

# do the same for 1988 image
p224r63_1988 <- brick("p224r63_1988_masked.grd")
plotRGB(p224r63_1988, r=5, g=4, b=3, stretch="Lin")

# plot the two images one beside the other
par(mfrow=c(1,2))
plotRGB(p224r63_1988, r=5, g=4, b=3, stretch="Lin")
plotRGB(p224r63_2011, r=5, g=4, b=3, stretch="Lin")

names(p224r63_2011) # get the names of the objects

# check correlation
dev.off()
plot(p224r63_2011$B1_sre, p224r63_2011$B3_sre)

# PCA principal component analysis
# decrease the resolution
p224r63_2011_res <- aggregate(p224r63_2011, fact=10)
p224r63_2011_pca <- rasterPCA(p224r63_2011_res)

# check the properties of the pca
p224r63_2011_pca

# link the map to the pca
# make a color ramp palette named cl
cl <- colorRampPalette(c('dark grey','grey','light grey'))(100)
plot(p224r63_2011_pca$map, col=cl)

# summary of the various models                    
summary(p224r63_2011_pca$model)
                    
# make a matrix of scateerplots                    
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
plot(difpca$PC1,col=cldif) # link difpca and PC1

##################################
##################################
          
# 9. R code radiance
                    
# bit example
library(raster)

# create data with 2 columns and 2 rows
toy <- raster(ncol=2, nrow=2, xmn=1, xmx=2, ymn=1, ymx=2)
values(toy) <- c(1.13,1.44,1.55,3.4)

plot(toy)
text(toy, digits=2)

toy2bits <- stretch(toy,minv=0,maxv=3)

# avoid decimal degrees, make use of integer
storage.mode(toy2bits[]) = "integer"
plot(toy2bits)
text(toy2bits, digits=2)

# increase the amount of bits
toy4bits <- stretch(toy,minv=0,maxv=15)
storage.mode(toy4bits[]) = "integer"
plot(toy4bits)
text(toy4bits, digits=2)
 
toy8bits <- stretch(toy,minv=0,maxv=255)
storage.mode(toy8bits[]) = "integer"
plot(toy8bits)
text(toy8bits, digits=2)

# plot altogether
par(mfrow=c(1,4))
plot(toy)
text(toy, digits=2)
plot(toy2bits)
text(toy2bits, digits=2)
plot(toy4bits)
text(toy4bits, digits=2)
plot(toy8bits)
text(toy8bits, digits=2)

#####################

setwd("C:/lab/")
load("faPAR.RData")
library(raster)
library(rasterdiv)

# see data
ls()
faPAR10

#let's see how much space is needed for the 8-bit set
writeRaster(copNDVI, "copNDVI.tif")

#faPAR: levelplot this set
library(rasterVis)
levelplot(faPAR10)
          
##################################
##################################
                    
# 10. R code faPAR
                    
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
# regression model between faPAR and NDVI
# amount of erosion and heavy metals
erosion <- c(12,14,16,24,26,40,55,67)
hm <- c(30,100,150,200,260,340,460,600)

plot(erosion, hm, col="red", pch=19, xlab="erosion", ylab="heavy metals", cex=2)
model1 <- lm(hm ~ erosion) # fitting linear model, used for regression
summary(model1)
abline(model1) # adds a straight line through the current plot

library(raster)
library(rasterdiv)
libray(sf)
setwd("C:/lab/")
                    
faPAR10<-raster("faPAR10.tif")
plot(faPAR10)
plot(copNDVI)

copNDVI <- reclassify(copNDVI,cbind(253:255,NA),right=TRUE)
                    
random.points <- function(raster,n)

pts <- random.points(faPAR10,1000)
copNDVIp <- extract(copNDVI, pts)
faPAR10p <- extract(faPAR10,pts)

##################################
##################################
                    
# 11. R code EBVs
# essential biodiversity variables
                    
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
                    
# make use of NIR
plotRGB(snt,4,3,2, stretch="lin") 

pairs(snt) # scatterplots for each combination of variables 

### PCA analysis
sntpca <- rasterPCA(snt)
sntpca

summary(sntpca$model)
# 70% of information
# plot pca with map
plot(sntpca$map) 

plotRGB(sntpca$map, 1, 2, 3, stretch="lin")

# set the moving window
# create a matrix first
window <- matrix(1, nrow = 5, ncol = 5)
window
# calculate values for neighborhood of focal cells
# focal(linking the main elements, matrix, function) 
sd_snt <- focal(sntpca$map$PC1, w=window, fun=sd)
cl <- colorRampPalette(c('dark blue','green','orange','red'))(100) # 
plot(sd_snt, col=cl)

par(mfrow=c(1,2))
plotRGB(snt,4,3,2, stretch="lin", main="original image") 
plot(sd_snt, col=cl, main="diversity")

##########################

# Cladonia stellaris example
                    
setwd("C:/lab/")

# import the image, RGB layers
library(raster)
clad <- brick("cladonia_stellaris_calaita.JPG")
plotRGB(clad,1,2,3,stretch="lin")

# make a 3x3 matrix
window <- matrix(1, nrow=3, ncol=3)
window

# we want to calculate values for the neighborhood of focal cells --> focal()
                    
pairs(clad)

# PCA
library(RStoolbox)
cladpca <- rasterPCA(clad)
cladpca

summary(cladpca$model)
# 98%

plotRGB(cladpca$map,1,2,3,stretch="lin")

# now we can calculate the values
sd_clad <- focal(cladpca$map$PC1, w=window, fun=sd)

# PC1 aggregate
PC1_agg <- aggregate(cladpca$map$PC1, fact=10)
sd_clad_agg <- focal(PC1_agg, w=window, fun=sd)

# multiframe graph
# plot the calculation
par(mfrow=c(1,2))
cl <- colorRampPalette(c('yellow','violet','black'))(100) 
plot(sd_clad, col=cl)
plot(sd_clad_agg, col=cl)
                    
##################################
##################################
                    
# 12. R code snow
                   
setwd("C:/lab/")

install.packages("ncdf4")
library(ncdf4)
library(raster)

snowmay <- raster("c_gls_SCE_202005260000_NHEMI_VIIRS_V1.0.1.nc")
cl <- colorRampPalette(c('darkblue','blue','light blue'))(100) 

# Exercise: plot snow cover with the cl palette
plot(snowmay,col=cl)

# import snow data
# set the working directory: folder where all the images are
setwd("C:/lab/snow")
                    
# perform raster() and plot() for each image separately --> time consuming way
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

# how to plot snow data with lapply --> fast way
setwd("C:/lab/snow/snow/") # working directory: folder containing only the images
# first of all we make a list of the files we are going to import
rlist <- list.files(pattern="snow") # pattern: picks all the files with "snow" in their names
rlist
# import all the images in the list 
import <- lapply(rlist,raster)
snow.multitemp <- stack(import)
plot(snow.multitemp,col=cl)

# let's make a prediction on how the snow is expected to change in 2025
source("prediction.r")

#######################

# set working directory
setwd("C:/lab/snow/snow")

# Exercise: import the snow cover images altogether
library(raster)
rlist <- list.files(pattern="snow")
import <- lapply(rlist,raster)
snow.multitemp <- stack(import)
cl <- colorRampPalette(c('darkblue','blue','light blue'))(100) 
plot(snow.multitemp, col=cl)

# import the prediction image for 2025 and plot it with color ramp palette
prediction <- raster("predicted.2025.norm.tif")
plot(prediction, col=cl)

# export the output (you made the calculation and you want to send the output to a colleague)
writeRaster(prediction, "final.tif")

# final stack
final.stack <- stack(snow.multitemp, prediction)
plot(final.stack,col=cl)

# export the R graph in .pdf
pdf("my_final_exciting_graph.pdf")
plot(final.stack, col=cl)
dev.off()

# make a .png of the project
png("my_final_exciting_graph.pdf")
plot(final.stack, col=cl)
dev.off()

# plot datas of 2010 and 2020 
plot(snow.multitemp$snow2010r, snow.multitemp$snow2020r)
abline(0,1, col="red")
                    
##################################
##################################
                    
# 13. R code NO2
# air pollution in time
                    
# set the working direcotry
setwd("C:/lab/no2/")
                    
library(raster)

# create a list of the images that we want to import
rlist <- list.files(pattern="EN")

# import images in the list, perform stack and add the color ramp palette
import <- lapply(rlist,raster)
EN <- stack(import)
cl <- colorRampPalette(c('red','orange','yellow'))(100) #
plot(EN, col=cl)

# situation in january and march
par(mfrow=c(1,2))
plot(EN$EN_0001, col=cl)
plot(EN$EN_0013, col=cl)

# RGB space
plotRGB(EN, r=1, g=7, b=13, stretch="lin")

# difference map between the two situations
dif <- EN$EN_0013 - EN$EN_0001
cld <- colorRampPalette(c('blue','white','red'))(100) # 
plot(dif, col=cld)

# quantitative estimate box plot
boxplot(EN)
# to remove the outline
boxplot(EN,outline=F)
# make the box plot horizontal
boxplot(EN,outline=F,horizontal=T,axes=T)

# plot the data of the first image with the data of the last image
plot(EN$EN_0001,EN$EN_0013)
#abline a=0, b=1
abline(0,1,col="red")

##################################
##################################
                    
# 14. R code zoom and crop functions
                    
# set my working directory to the lab folder
setwd("C:/lab/")
                    
# we are going to see two possibilities on how to make use of a specified part of a downloaded data
# we can use a particular part of the data instead of global perspective 
# recall the libraries needed
library(raster)
library(ncdf4)
                    
# import data and name it snow
snow <- raster("c_gls_SCE_202005260000_NHEMI_VIIRS_V1.0.1.nc")
                  
# make a color ramp palette ranging from dark blue to light blue
cl <- colorRampPalette(c('darkblue','blue','light blue'))(100)
                    
# we have to specify which extent we want to zoom in, in the data/image 
# in our case we zoom on Italy --> the first two numbers are ranging in longitude, the last two numbers are ranging in latitude
# assign it the name ext
ext <- c(0, 20, 35, 50)
# use function zoom(), specifyng the image and the extension
zoom(snow, ext=ext)
                    
# use crop() naming it snowitaly
snowitaly <- crop(snow, ext)
plot(snowitaly, col=cl)
                    
# other example with zoom() using drawExtent by clicking on the map and drawing a rectangular that will be our zoom
zoom(snow, ext=drawExtent())
                    
##################################
##################################
                    
# 15. R code interpolation with students' data
                    
setwd("C:/lab/")
# import data, it is a table
# assign it the name inp
library(spatstat)
inp <- read.table("dati_plot55_LAST3.csv", sep=";", head=T) # separation is the symbol ;
head(inp)
attach(inp) # by attaching the dataset, we can work directly in it
# we can make directly use of the columns X and Y
plot(X,Y)
                    
# in order to know the minimum and maximum values make use of summary()
summary(inp)
inppp <- ppp(x=X, y=Y, c(716000,718000),c(4859000,4861000))
# see the name of the columns
names(inp)
                    
# we want to make estimations of the tree coverage --> Canopy.cov
marks(inppp) <- Canopy.cov
# smooth() for interpolation. Estimates data where it has not yet been measured
# assign the name canopy
canopy<-Smooth(inppp)
plot(canopy)
points(inppp, col="green")
# southern part is more dense
marks(inppp) <- cop.lich.mean
lichs <- Smooth(inppp)
plot(lichs)
points(inppp)
       
# plot the maps together to compare them
par(mfrow=c(1,2))
plot(canopy)
points(inppp)
plot(lichs)
points(inppp)              
 
# now plot the three maps together
par(mfrow=c(1,3))
plot(canopy)
points(inppp)
plot(lichs)
points(inppp)
plot(Canopy.cov, cop.lich.mean, col="red", pch=19, cex=2)

dev.off()
                    
#############
                    
# Dati Psammofile
         
# we are going to read the dataset, giving it the name inp.psam
# include header and separator ;
inp.psam<-read.table("dati_psammofile.csv",sep=";", head=T)
attach(inp.psam)
head(inp.psam)
summary(inp.psam)
                    
# by reading the coordinates we can see that we are in the north east of Italy
# X-value: E 
# Y-value: N
plot(E,N)  
inp.psam.ppp <- ppp(x=E,y=N,c(356450,372240),c(5059800,5064150))
                    
# we are going to make use of the variable organic carbon C_org
marks(inp.psam.ppp) <- C_org
                    
# we estimate values for the points we don't have measured values for
C<-Smooth(inp.psam.ppp)
plot(C)
points(inp.psam.ppp)
                    
##################################
##################################
                    
# 16. R code for species distribution modelling
                    
# install the needed packages and recall the libraries
install.packages("sdm")

library(sdm)
library(raster) # for predictors
library(rgdal) # for species

# species data are in the folder "external" by downloading sdm
# import the file
file <- system.file("external/species.shp", package="sdm")

# use the graphical part of the file
species <- shapefile(file)
plot(species)

# species occurrence presence/absence
# 1 for presence, 0 for absence
species$Occurrence
plot(species[species$Occurrence == 1,],col='blue',pch=16)
points(species[species$Occurrence == 0,],col='red',pch=16)

# import environmental variables/predictors
path <- system.file("external", package="sdm")

# make a list of the files
lst <- list.files(path=path,pattern='asc$',full.names = T) 
lst # shows the files listed and their path

# put the predictors altogether with stack function
preds <- stack(lst)
plot(preds)

# change colors
cl <- colorRampPalette(c('blue','orange','red','yellow')) (100)
plot(preds, col=cl)

# show the distribution of Brachipodium rupestris according to each environmental variable
# let's make use of elevation
plot(preds$elevation, col=cl)
points(species[species$Occurrence == 1,], pch=16)

# make use of temperature
plot(preds$temperature, col=cl)
points(species[species$Occurrence == 1,], pch=16)

# make use of precipitation 
plot(preds$precipitation, col=cl)
points(species[species$Occurrence == 1,], pch=16)

# make use of vegetation
plot(preds$vegetation, col=cl)
points(species[species$Occurrence == 1,], pch=16)

# model the sets species and predictors
d <- sdmData(train=species, predictors=preds)
d # see the data

m1 <- sdm(Occurrence ~ elevation + precipitation + temperature + vegetation, data=d, methods = "glm")

p1<- predict(m1, newdata=preds)

# plot species and preds
plot(p1, col=cl)
points(species[species$Occurrence == 1,], pch=16)

s1 <- stack(preds,p1)
plot(s1, col=cl)

                    
             
                    
                    
                    
                    
                    
                    

                    
