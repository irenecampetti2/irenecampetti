### multipanel in R the second lecture of monitoring ecosystems

install.packages("sp")
install.packages("GGally")

library(sp) # require(sp) will also work
library(GGally)

data(meuse) # there is a dataset available

attach(meuse)
head(meuse)
plot(cadmium,zinc)
plot(cadmium,zinc,pch=15,col="red",cex=2)

# Ex. make plot all the possible variables

pairs(meuse)
pairs(~ cadmium + copper + lead + zinc, data=meuse)

pairs(meuse[,3:6])

# Ex. prettify the graph

pairs(meuse[,3:6],pch=8,col="blue",cex=0.5)

ggpairs(meuse[,3:6])
