# R code for multivariate analysis

library(vegan)
setwd("C:/lab/")

biomes <- read.table("biomes.csv", header=T, sep=",")
head(biomes) # view(biomes), biomes are also good

# DEtrended CORrespondence ANAlysis
multivar <- decorana(biomes) 
plot(multivar)
multivar

plot(multivar)
biomes_types <- read.table("biomes_types.csv", header=T, sep=",")
head(biomes_types)

attach(biomes_types)
ordiellipse(multivar, type, col=1:4, kind = "ehull", lwd=3)

# Also col="blue", "red", "green", "black")

ordispider(multivar, type, col=1:4, label= TRUE)
