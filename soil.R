library(raster)
library(sp)
library(missMDA)
library(FactoMineR)
library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

rf = list.files("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/soils/",
                pattern = "crop.tif",
               full.names = T, recursive = T)

d = read.csv("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/Ctenotus_zastictus.project3.csv",
             stringsAsFactors = F)
d1 = d[which(d$mtDNA_LINEAGE %in% c("pallasotus", "zastictus")), ]
d1 = d1[d1$SAMPLE_ID != 'AMR_123129_Ct_zast', ]
d2 = d1[complete.cases(d1$LAT), ]
pzpts = SpatialPoints(d2[,c("LON", "LAT")])

bbox = c(112, 122, -28, -20)
wrld = rgdal::readOGR("~/Dropbox (Personal)/brazil_gene_flow/data/ne_10m_land/ne_10m_land.shp")
oz = crop(wrld, bbox)
pts = spsample(oz, n=1000, "random")

s = data.frame(matrix(NA, nrow = length(pts) + nrow(d2), 
                      ncol = length(rf)))
for (i in 1:length(rf)) {
  r = raster(rf[i])
  s[,i] = c(raster::extract(r, pts), raster::extract(r, pzpts))
  cat(i, "\n")
}

s$name = c(paste("random", 1:1000, sep = "."), d2$SAMPLE_ID)
missind = apply(s, 1, function(x){sum(is.na(x))})
# inds that are totally missing
s1 = s[which(missind != 163), ]
missr = apply(s1, 2, function(x){sum(is.na(x))})
s2 = s1[,which(missr != 129)]

nb <- estim_ncpPCA(s2[ ,1:139], method.cv = "Kfold", verbose = TRUE,
                   scale = TRUE)
res.comp <- imputePCA(s2[ ,1:139], ncp = nb$ncp)
x = PCA(res.comp, scale.unit = TRUE)
pccoords = data.frame(x$ind$coord)
pccoords$name = s2$name
pccoords$type = "random"
pccoords[grep("Ct", pccoords$name), "type"] = "pall"
pccoords[grep("zast", pccoords$name), "type"] = "zast"

a = ggplot(pccoords %>% filter(type == "random"), 
       aes(Dim.1, Dim.2)) +
  geom_point(shape = 16, col = "gray80") +
  geom_point(data = pccoords %>% filter(type != "random"),
             aes(Dim.1, Dim.2, fill = type), shape = 21, size = 2) +
  xlab("PC1 axis (32.1%)") + ylab("PC2 axis (25.8%)") +
  scale_fill_manual(values = c("#377eb8", "#e41a1c"))

rf2 = list.files("~/Dropbox (Personal)/macroevolution/eco_IBD_oz/data/AUS_5arc/", pattern = "tif",
                full.names = T)
e = data.frame(matrix(NA, nrow = length(pts) + nrow(d2), 
                      ncol = length(rf2)))
for (i in 1:length(rf2)) {
  r = raster(rf2[i])
  e[,i] = c(raster::extract(r, pts), raster::extract(r, pzpts))
  cat(i, "\n")
}
e$name = c(paste("random", 1:1000, sep = "."), d2$SAMPLE_ID)
missind = apply(e, 1, function(x){sum(is.na(x))})
# inds that are totally missing
e1 = e[which(missind < 20), ]
missr = apply(e1, 2, function(x){sum(is.na(x))})

epca <- estim_ncpPCA(e1[ ,1:27], method.cv = "Kfold", verbose = TRUE,
                   scale = TRUE)
epca2 <- imputePCA(e1[ ,1:27], ncp = epca$ncp)
epca3 = PCA(epca2, scale.unit = TRUE)

eccoords = data.frame(epca3$ind$coord)
eccoords$name = e1$name
eccoords$type = "random"
eccoords[grep("Ct", eccoords$name), "type"] = "pall"
eccoords[grep("zast", eccoords$name), "type"] = "zast"

b = ggplot(eccoords %>% filter(type == "random"), 
       aes(Dim.1, Dim.2)) +
  geom_point(shape = 16, col = "gray80") +
  geom_point(data = eccoords %>% filter(type != "random"),
             aes(Dim.1, Dim.2, fill = type), shape = 21, size = 2) +
  xlab("PC1 axis (48.1%)") + ylab("PC2 axis (20.8%)") +
  scale_fill_manual(values = c("#377eb8", "#e41a1c"))

prow <- plot_grid(
  a + theme(legend.position="none"),
  b + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B"),
  hjust = -1,
  nrow = 1
)


pdf("~/Desktop/soil_env.pdf", height = 3, width = 6)
plot_grid(prow)
dev.off()
