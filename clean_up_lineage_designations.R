d = read.csv("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/Ctenotus_zastictus.project2.csv",
             stringsAsFactors = F)
dd = read.csv("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/mtDNA/Ctenotus_zastictus_mtDNA.6July21.csv",
              stringsAsFactors = F)
dd1 = dd[! (dd$SAMPLE_ID %in% d$SAMPLE_ID), ]

d2 = dplyr::bind_rows(d, dd1)
d2$SEQ = NULL
d2$RESEQ = NULL
d2$SAMPLE_ID.1 = NULL
write.csv(d2, "~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/Ctenotus_zastictus.project3.csv",
          row.names = F)

library(ape)
setwd("~/Dropbox (Personal)/publications/Ctenotus_zasticus/")
d = read.csv("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/Ctenotus_zastictus.project3.csv",
             stringsAsFactors = F)
t = read.tree("data/phylogenies/ctzast1.dropinds.fa.treefile")
t$tip.label = paste(d[match(t$tip.label, d$SAMPLE_ID), "nDNA_LINEAGE"], t$tip.label)

pdf("~/Desktop/nDNA_phylogeny.pdf", height = 20, width = 10)
par(mar = c(0, 0, 0, 0))
plot(t, cex = 0.5)
dev.off()

d1 = d[complete.cases(d$nDNA_LINEAGE), ]
d1[d1$nDNA_LINEAGE != d1$mtDNA_LINEAGE, ]

d = read.csv("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/Ctenotus_zastictus.project3.csv",
             stringsAsFactors = F)
t = read.tree("data/mtDNA/IQtree.cytb_partition.treefile")
t$tip.label = paste(d[match(t$tip.label, d$SAMPLE_ID), "mtDNA_LINEAGE"], t$tip.label)
pdf("~/Desktop/phylogeny.pdf", height = 20, width = 10)
par(mar = c(0, 0, 0, 0))
plot(t, cex = 0.5)
dev.off()