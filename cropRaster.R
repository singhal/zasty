library(raster)

d = read.csv("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/Ctenotus_zastictus.project3.csv",
             stringsAsFactors = F)
d1 = d[which(d$mtDNA_LINEAGE %in% c("pallasotus", "zastictus")), ]
range(d1$LAT, na.rm = T)
range(d1$LON, na.rm = T)

bbox = c(112, 122, -28, -20)
## confirm that crop box is right
# x = raster("~/Desktop/soils/AWC/AWC_000_005_05_N_P_AU_WAT_D_20140801.tif")
# plot(x)
# polygon(x = c(112, 112, 122, 122), y = c(-27, -20, -20, -27))
# points(d1[,c("LON", "LAT")])

rasters = list.files("~/Desktop/soils/", pattern = "*tif", full.names = T, recursive = T)
for (i in 1:length(rasters)) {
  x = raster(rasters[i])
  cropr = raster::crop(x, bbox)
  newfile = gsub(".tif", "_crop.tif", rasters[i])
  writeRaster(cropr, newfile)
  cat(i, "\n")
}