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

# 2. R code spatial
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
covid <- read.table("covid_agg.csv", head=T)

head(covid)

attach(covid)
plot(country, cases)

# plot (covid$country, covid$cases) -- with different orientation of the tick mark labels
plot(country, cases, las=0) # parallel labels
plot(country, cases, las=1) # horizontal labels
plot(country, cases, las=2) # perpendicular labels
plot(country, cases, las=3) # vertical labels

plot(country, cases, las=3, cex.axis=0.5) # cex.axis -- change the size of labels
plot(country, cases, las=3, cex.axis=0.7)

# install ggplot2 package
install.packages("ggplot2")

# set the working directory
setwd("C:/lab/")

# load dataset
load("spatial.RData")
# list objects in the dataset
ls()

library(ggplot2)  # also require(ggplot2)

# load mpg data
data(mpg)
head(mpg)

#key components: data, aes, geometry
ggplot(mpg, aes(x=displ, y=hwy)) + geom_point()
ggplot(mpg, aes(x=displ, y=hwy)) + geom_line()
ggplot(mpg, aes(x=displ, y=hwy)) + geom_polygon()

head(covid)
ggplot(covid, aes(x=lon, y=lat, size=cases)) + geom_point()
