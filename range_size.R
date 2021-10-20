library(dplyr)
library(ggplot2)
library(alphahull)

d = readxl::read_xlsx("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/localities/zastictus from WAM.xlsx")
d1 = as.data.frame(d[, c("LONGDEC", "LATDEC")])
d2 = d1 %>% mutate_all(as.numeric) %>% unique()

plot(oz)
plot(d2)

ggplot() +
  geom_sf(data = oz, fill = "gray99", color = "gray30", size = 0.15) +
  geom_point(data = d2, 
             aes(x = LONGDEC, y = LATDEC)) 

pt <- SpatialPoints(d2, proj4string=CRS('+proj=longlat +datum=WGS84'))
pt <- spTransform(pt, CRS("+init=epsg:3395"))
b <- rgeos::gBuffer(rgeos::gConvexHull(pt), width=5)
b1 <- spTransform(b, CRS('+proj=longlat +datum=WGS84'))
b2 = spTransform(b, CRS("+init=epsg:3395"))

plot(d2)
plot(b, add = T)
rgeos::gArea(b1)
rgeos::gArea(b2) / 1000 / 1000

library(rgdal)
r = readOGR("~/Dropbox (Personal)/inornatus_gr/Meiri_etal_ranges/GARD1.1_dissolved_ranges/modeled_reptiles.shp")
rnames = r[[1]]
rnames = gsub(" ", "_", rnames)
ct = grep("Ctenotus", rnames)
sizes = rep(NA, length(ct))

for (i in 1:length(ct)) {
  b = r[ct[i], ]
  b2 = spTransform(b, CRS("+init=epsg:3395"))
  sizes[i] = rgeos::gArea(b2) / 1000 / 1000
}

