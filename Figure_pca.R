library(adegenet)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

cols = brewer.pal(6, "Paired")
sps = c("duricola", "pallasotus", "piankai", "rhabdotus", 
        "serventyi", "zastictus")

setwd("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/ipyrad_output/zasticus_ingroup_SNPs//")
sampfile = "~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/Ctenotus_zastictus.project4.csv"
indfile =  "zasticus.miss0.7.MAC0.thinned.fam"
snpfile =  "zasticus.miss0.7.MAC0.thinned.snp"

d = read.csv(sampfile, stringsAsFactors = F, na.string = c(""))
inds = read.table(indfile, stringsAsFactors = F, header = F, sep = " ")
inds = data.frame(inds = inds$V2, 
                  sp = d[match(inds$V2, d$SAMPLE_ID), "nDNA_LINEAGE"], 
                  stringsAsFactors = F)
rownames(inds) = inds$inds

names(cols) = sps
inds$col = cols[inds$sp]

##############
# PCA results
##############

t = read.snp(snpfile)
# requires no missing data
# pcasnp <- dudi.pca(t, center=TRUE, scale=FALSE)
pcasnp <- dudi.pco(dist(t), scannf=FALSE, nf=3)

eig.perc <- pcasnp$eig/sum(pcasnp$eig)
pcal = data.frame(pcasnp$li)
pcal$sp = inds$sp

a = ggplot(pcal, aes(A1, A2, fill = sp)) + 
  ggforce::geom_mark_ellipse(aes(x = A1, y = A2, fill = sp), 
                             expand = unit(2, "mm"), size = 0.1, alpha = 0.3) +
  geom_point(shape = 21, size = 2.5) +
  theme_cowplot() +
  theme(legend.position = "none") +
  scale_fill_manual(values = cols) +
  xlab("PC 1 (21.2%)") +
  ylab("PC 2 (15%)")
save_plot("~/Desktop/PCA.pdf", a, base_height = 4, base_width = 4)

# repeat it with other datasets
# results are consistent - zast is often nested
# or barely outside of pallostatus

# repeat it without the highly divergent lineages
divsp = c("rhabdotus", "serventyi", "piankai", "duricola")
dropind = inds[inds$sp %in% divsp, "inds"]
t2 <- t[!indNames(t) %in% dropind]
pcasnp2 <- dudi.pco(dist(t2), scannf=FALSE, nf=3)

eig.perc <- pcasnp2$eig/sum(pcasnp2$eig)
pcal = data.frame(pcasnp2$li)
pcal$sp = inds[match(t2@ind.names, inds$inds), "sp"]
eig.perc

b = ggplot(pcal, aes(A1, A2, fill = sp)) + 
  ggforce::geom_mark_ellipse(aes(x = A1, y = A2, fill = sp), 
                             expand = unit(2, "mm"), size = 0.1, alpha = 0.3) +
  geom_point(shape = 21, size = 2.5) +
  theme_cowplot() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c(cols[2], cols[6]))+ 
  xlab("PC 1 (29.0%)") +
  ylab("PC 2 (20.6%)")
save_plot("~/Desktop/PCA2.pdf", b, base_height = 4, base_width = 4)

