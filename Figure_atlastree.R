library(RColorBrewer)
library(ggtree)

d = read.csv("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/Ctenotus_zastictus.project4.csv",
             stringsAsFactors = F)
cols = rep(brewer.pal(12, "Paired"), 3)

t = treeio::read.iqtree("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/phylogenies/atlas_morpho_final.treefile")

# root and drop
outs = c("CUMV_14616_Le_dese", "CUMV_14613_Le_dese", "CUMV_14614_Le_dese")
drop = c("SAMR_42153_Ct_pant", "CUMV_1457414505_Ct_pant", "WAMR_129335_Ct_pant",
         "WAMR_110668_Ct_duri")
nodeout = phytools::findMRCA(t@phylo, outs)
t1 = treeio::root(t, node = nodeout)
t2 = treeio::drop.tip(t1, outs)
t2 = treeio::drop.tip(t2, drop)

inds = t2@phylo$tip.label
isps = d[match(inds, d$SAMPLE_ID), "nDNA_LINEAGE"]
# sps = c("duricola", "pallasotus", "piankai", "rhabdotus", 
#        "serventyi", "zastictus", "hanloni", "ariadnae",
#        "dux", "quattuordecimlineatus2", "atlas", "grandis")

nodenum = seq(Ntip(t2@phylo) + 1, Ntip(t2@phylo) + Nnode(t2@phylo))
t2@data$upper = FALSE

for (i in nodenum) {
  curlins = isps[ phangorn::Descendants(t2@phylo, i, type = "tips")[[1]]  ]
  curtips = inds[ phangorn::Descendants(t2@phylo, i, type = "tips")[[1]]  ]
  
  if (length(unique(curlins)) > 1 ) {
    t2@data[which(t2@data$node == i), "upper"] = TRUE
  } else if (length(unique(curlins)) == 1) {
    sp = unique(curlins)
    sptips = inds[which(isps == sp)]
    if (setequal(sptips, curtips)) {
      t2@data[which(t2@data$node == i), "upper"] = TRUE
    }
  }
}
t2@data$support = NA
t2@data[ which(t2@data$upper == TRUE & t2@data$SH_aLRT >= 95), "support"] = TRUE
t2@data[ which(t2@data$upper == TRUE & t2@data$SH_aLRT < 95), "support"] = FALSE

a = ggtree(t2) + geom_tiplab(size = 2)
sps = unique(isps)
for (i in 1:length(sps)) {
  tips = inds[ which(isps == sps[i]) ]
  mrca = phytools::findMRCA(t2@phylo, tips)
  a = a + geom_balance(node=mrca, fill= cols[i], color= NA,
                       alpha=0.6, extend= 0.007) +
    geom_cladelabel(mrca, paste("C.", sps[i]), fontface = "italic", offset=0, 
                    barsize = NA, angle=0, offset.text=0.007, 
                    align = T, size = 1)
}
pdf("~/Dropbox (Personal)/publications/Ctenotus_zasticus/figures/FigureSI_atlas.pdf", height = 9, width = 6.5)
a + xlim(0, 0.03) + 
  geom_nodepoint(aes(fill = support, subset = support), shape = 21, size = 2) +
  scale_fill_manual(values = c("white")) +
  theme(legend.position = "none") +
  geom_treescale(x = 0, y = 80, width = 0.005)
dev.off()
