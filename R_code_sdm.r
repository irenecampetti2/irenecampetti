# R code for species distribution modelling

# install the needed packages and recall the packages
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



