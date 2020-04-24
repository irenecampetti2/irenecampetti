# R code for remote sensing data analysis

install.packages("raster")
install.packages("RStoolbox")

setwd("C:/lab/")
library(raster)
p224r63_2011 <- brick("p224r63_2011_masked.grd")

plot(p224r63_2011)

# B1: blue
# B2: green
# B3: red
# B4: NIR

cl <- colorRampPalette(c('black','grey','light grey'))(100) # 
plot(p224r63_2011, col=cl)

# multiframe of different plots
par(mfrow=c(2,2))
# B1: blue
clb <- colorRampPalette(c('dark blue','blue','light blue'))(100)
plot(p224r63_2011$B1_sre, col=clb)

# Exercise: do the same with the green band

clg <- colorRampPalette(c('dark green','green','light green'))(100)
plot(p224r63_2011$B2_sre, col=clg)

# B3: red
clr <- colorRampPalette(c('dark red','red','pink'))(100)   # light red does not exist
plot(p224r63_2011$B3_sre, col=clr)

# B4: NIR
cln <- colorRampPalette(c('red','orange','yellow'))(100)
plot(p224r63_2011$B4_sre, col=cln)


par(mfrow=c(4,1))
# B1: blue
clb <- colorRampPalette(c('dark blue','blue','light blue'))(100)
plot(p224r63_2011$B1_sre, col=clb)

# Exercise: do the same with the green band

clg <- colorRampPalette(c('dark green','green','light green'))(100)
plot(p224r63_2011$B2_sre, col=clg)

# B3: red
clr <- colorRampPalette(c('dark red','red','pink'))(100)   # light red does not exist
plot(p224r63_2011$B3_sre, col=clr)

# B4: NIR
cln <- colorRampPalette(c('red','orange','yellow'))(100)
plot(p224r63_2011$B4_sre, col=cln)


# plotRGB
plotRGB(p224r63_2011, r=3, g=2, b=1, stretch="Lin")
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin")

# Exercise: put NIR on top of the G component
plotRGB(p224r63_2011, r=3, g=4, b=2, stretch="Lin")

#######

setwd("C:/lab/")
load("rs.RData")
ls()

library(raster)
p224r63_1988 <- brick("p224r63_1988_masked.grd")

plot(p224r63_1988)

# Exercise: plot in visible RGB 321 both images
p224r63_2011 <- brick("p224r63_2011_masked.grd")
par(mfrow=c(2,1))
plotRGB(p224r63_1988, r=3, g=2, b=1, stretch="Lin")
plotRGB(p224r63_2011, r=3, g=2, b=1, stretch="Lin")

# Exercise: plot in false colour 432 both images
par(mfrow=c(2,1))
plotRGB(p224r63_1988, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin")

# Enhance the noise
par(mfrow=c(2,1))
plotRGB(p224r63_1988, r=4, g=3, b=2, stretch="Hist")
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Hist")


# B1: blue
# B2: green
# B3: red
# B4: NIR
dvi2011 <- p224r63_2011$B4_sre - p224r63_2011$B3_sre
cl <- colorRampPalette(c('darkorchid3','light blue','lightpink4'))(100)
plot(dvi2011, col=cl)

# Exercise: dvi for 1988
dvi1988 <- p224r63_1988$B4_sre - p224r63_1988$B3_sre
clb <- colorRampPalette(c('darkorchid3','light blue','lightpink4'))(100)
plot(dvi1988, col=clb)

# Difference from one year to another
diff <- dvi2011 - dvi1988
plot(diff)

# Changing grain
p224r63_2011res <- aggregate(p224r63_2011, fact=10)
p224r63_2011res100 <- aggregate(p224r63_2011, fact=100)

par(mfrow=c(3,1))
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011res, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011res100, r=4, g=3, b=2, stretch="Lin")

p224r63_2011
p224r63_2011res
p224r63_2011res100
  

