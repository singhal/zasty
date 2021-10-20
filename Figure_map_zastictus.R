library(readr)
library(rgdal)
library(dplyr)
library(scales)
library(raster)

d = read.csv("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/Ctenotus_zastictus.project4.csv",
             stringsAsFactors = F, na.strings = c("NA", ""))
# catalog number
d1 = read.table("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/localities/0318040-200613084148143.csv",
                sep = "\t", stringsAsFactors = F, header = T)
# hard do do automatically, but one of the tissues is in d3 
# but the record here doesn't have lat longs
d1a = d1[grep("ABTC59808", d1$catalogNumber, invert = T), ]

# catalog number
d2 = read.table("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/localities/records-2021-07-07/records-2021-07-07.csv",
                sep = ",", stringsAsFactors = F, header = T)
# most of d2 is already in d1
# note that d2 includes inaturalist individual
d2a = d2[!d2$catalogNumber %in% d1$catalogNumber, ]

d3 = d[which(d$nDNA_LINEAGE == "zastictus"),]
d3 = d3[which(d3$SAMPLE_ID != "AMR_123129_Ct_zast"),]

# all of these WAM points are in ALA data
d4 = readxl::read_xlsx("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/localities/zastictus from WAM.xlsx")

dd = data.table::rbindlist(list(d1a[, c("decimalLatitude", "decimalLongitude")],
                                d2a[, c("decimalLatitude", "decimalLongitude")],
                                d3[, c("LAT", "LON")]))

wrld = raster::stack("~/Dropbox (Personal)/brazil_gene_flow/data/HYP_HR_SR_W_DR/HYP_HR_SR_W_DR.tif")
oz = crop(wrld, extent(112, 155, -40, -10))

worldmap <- rgdal::readOGR("~/Dropbox (Personal)/brazil_gene_flow/data/ne_10m_land/ne_10m_land.shp")
oze = crop(worldmap, extent(112, 155, -40, -10))
oze2 <- rgeos::gSimplify(oze, 0.05)
oz2 = mask(oz, oze2)

pdf("~/Dropbox (Personal)/publications/Ctenotus_zasticus/figures/Figure1_map.pdf")
plot(NA, xlim = c(112, 155), ylim = c(-40, -10),
     axes = F, xlab = "", ylab = "")
plotRGB(oz2, add =T, maxpixels = ncell(oz))
points(dd[ ,c(2, 1)], pch = 21, bg = alpha("white", 0.5))
dev.off()


