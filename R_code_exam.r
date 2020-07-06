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

# Bands of Landsat
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

# 7. R code 

