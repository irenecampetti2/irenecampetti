# R code for spatial view of points

library(sp)

data(meuse)

head(meuse)

# coordinates
coordinates(meuse) = ~x+y

plot(meuse)

spplot(meuse, "zinc")

# Exercise: plot the spatial amount of copper

spplot(meuse, "copper", main="Copper concentration")

bubble(meuse, "zinc")
bubble(meuse, "zinc", main="Zinc concentration")

# Exercise: bubble the copper in red

bubble(meuse, "copper", main="Copper concentration", col="red")

### Importing new data

# put the covid_agg.csv file into the folder lab

# setting the working directory: lab
# windows
setwd("C:/lab/")

covid <- read.table("covid_agg.csv", head=T)

head(covid)

attach(covid)
plot(country, cases)

plot(country, cases, las=0) # parallel labels

plot(country, cases, las=1) # horizontal labels
plot(country, cases, las=2) # perpendicular labels
plot(country, cases, las=3) # vertical labels
plot(country, cases, las=3, cex.axis=0.5)

# ggplot2 package

install.packages("ggplot2")
library(ggplot2)   # also require(ggplot2)

