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

