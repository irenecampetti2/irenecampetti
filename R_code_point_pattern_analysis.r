# point pattern analysis: density map

install.packages("spatstat")
library(spatstat)

attach(covid)
head(covid)
covids <- ppp(lon, lat, c(-180, 180), c(-90, 90))

# without attaching the covid set: covids <- ppp(covid$lon, covid$lat, c(-180, 180), c(-90, 90))
d <- density(covids)
plot(d)
points(covids)

setwd("C:/lab/")
load(".RData")
ls()

# covids: point pattern
# d: density map

library(spatstat)
plot(d)
points(covids)

install.packages("rgdal")
library(rgdal)

#let's input vector lines (x0y0, x1y1, x2,y2...)
coastlines <- readOGR("ne_10m_coastline.shp")
plot(d)
points(covids)
plot(coastlines, add=T)

cl <- colorRampPalette(c("yellow", "orange", "red")) (100)
plot(d, col=cl)
points(covids)
plot(coastlines, add=T)

#Exercise: make a colorRampPalette
cl <- colorRampPalette(c("blue", "purple", "pink", "white")) (100)
plot(d, col=cl, main="Densities of covid-19")
points(covids)
plot(coastlines, add=T)

# new example: abrupt change of colors!
cll <- colorRampPalette(c("blue", "purple", "pink", "white")) (5)
plot(d, col=cll, main="Densities of covid-19")
points(covids)
plot(coastlines, add=T)

# Export

pdf("covid_density.pdf")
cl <- colorRampPalette(c("blue", "purple", "pink", "white")) (100)
plot(d, col=cl, main="Densities of covid-19")
points(covids)
plot(coastlines, add=T)
dev.off()
