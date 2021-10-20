t = read.tree("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/phylogenies/atlas_morpho_final.treefile")
outs = c("CUMV_14616_Le_dese", "CUMV_14613_Le_dese", "CUMV_14614_Le_dese")
t = root(t, outs, resolve.root = T)
t = drop.tip(t, outs)

setwd("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/")
x = read.csv("Ctenotus_zastictus.project4.csv", stringsAsFactors = F)

keep = c("zastictus", "duricola", "pallasotus",
         "piankai", "serventyi", "rhabdotus", "hanloni",
         "dux", "quattuordecimlineatus2", "atlas", "ariadnae",
         "maryani")
inds = t$tip.label
names(inds) = x[match(inds, x$SAMPLE_ID), "nDNA_LINEAGE"]
keepsp = rep(NA, length(keep))
for (i in 1:length(keep)) {
  curinds = inds[which(names(inds) == keep[i])]
  keepsp[i] = sample(curinds, 1)
}
t1 = keep.tip(t, keepsp)
t1$tip.label = names(inds)[ match(t1$tip.label, inds) ]

d = read.csv("divergence/atlas_group.pi.22July21.csv",
             stringsAsFactors = F)
pi = d[match(t1$tip.label, d$sp1), "pi"]
pi[8] = 0.007
names(pi) = t1$tip.label
pdf("~/Desktop/pi.pdf", height = 3, width = 5)
phytools::dotTree(t1, pi)
dev.off()
