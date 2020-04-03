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
