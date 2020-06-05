# R_code_no2.r
# set the working direcotry
setwd("C:/lab/no2/")
library(raster)

# create a list of the images that we want to import
rlist <- list.files(pattern="EN")

# import images in the list, perform stack and add the color palette
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
plot(EN$EN_0013,EN$EN_0001)
#abline a=0, b=1
abline(0,1,col="red")




 
