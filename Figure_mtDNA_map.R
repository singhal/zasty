library(dplyr)
library(sf)
library(rnaturalearthdata)
library(rnaturalearth)
library(ggplot2)
library(rgeos)
library(RColorBrewer)

d = read.csv("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/Ctenotus_zastictus.project3.csv",
             stringsAsFactors = F, na.strings = c("NA", ""))

world = ne_countries(scale = "medium", returnclass = "sf")
oz = world[world$admin == "Australia", ]

lins = unique(d$nDNA_LINEAGE)
core = c("duricola", "pallasotus", "piankai", "serventyi", "rhabdotus", "zastictus")
core2 = c("duricola", "pallasotus", "piankai", "zastictus")

d = d[d$SAMPLE_ID != 'AMR_123129_Ct_zast', ]

pdf("~/Desktop/mtDNA_maps.pdf", height = 6, width = 4)
ggplot() +
  geom_sf(data = oz, fill = "gray99", color = "gray30", size = 0.15) +
  geom_point(data = d %>% filter(mtDNA_LINEAGE %in% core), 
             aes(x = LON, y = LAT, fill = mtDNA_LINEAGE), 
             shape = 21, alpha = 0.7,
             size = 2, color = "black", stroke = 0.2) +
  coord_sf(xlim = c(110.5, 156), ylim = c(-46, -8.5)) +
  scale_fill_manual(values = brewer.pal(6, "Paired")) +
  facet_wrap(~mtDNA_LINEAGE, ncol = 2) + 
  theme_void() +
  guides(fill = "none", size = "none", alpha = "none", shape = "none")
dev.off()
