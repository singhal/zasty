library(ape)

d = read.csv("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/Ctenotus_zastictus.project3.csv",
             stringsAsFactors = F)
cols = brewer.pal(6, "Paired")

t = read.tree("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/phylogenies/duricola.phy.treefile")
t$node.label = unlist(lapply(strsplit(t$node.label, "/"), function(x) {x[1]}))

# root and drop
outs = c("CUMV_14366_Ct_dux", "SAMR_36098_Ct_dux", "CUMV_14373_Ct_dux")
t1 = root(t, outs, resolve.root = T)
t2 = drop.tip(t1, outs)
t2$node.label = as.numeric(t2$node.label)

inds = t2$tip.label
isps = d[match(inds, d$SAMPLE_ID), "nDNA_LINEAGE"]
mtsps = d[match(inds, d$SAMPLE_ID), "mtDNA_LINEAGE"]
sps = c("duricola", "pallasotus", "piankai", "rhabdotus", 
        "serventyi", "zastictus")

edgecols = rep("black", length(t$edge.length))
mrcas = rep(NA, length(sps))
for (i in 1:length(sps)) {
  curtips = inds[ which(isps == sps[i]) ]
  mrca = phytools::findMRCA(t2, tips = curtips)
  mrcas[i] = mrca
  desc = phangorn::Descendants(t2, mrca, type = "all")
  edgecols[ which(t2$edge[,2] %in% desc) ] = cols[i]
}

pdf("~/Desktop/phylo.pdf", width = 4, height = 5)
par(mar = c(0, 0, 0, 10), xpd = T)
plot.phylo(t2, edge.color = edgecols, show.tip.label = F)
tiplabels(gsub("_Ct_\\S+", "", t2$tip.label), 1:Ntip(t2), frame = "none", adj = c(0, 0.5), cex = 0.5)
points(rep(0.029, Ntip(t2)), 1:Ntip(t2), pch = 21, bg = cols[ match(isps, sps) ])
points(rep(0.03, Ntip(t2)), 1:Ntip(t2), pch = 21, bg = cols[ match(mtsps, sps) ])
for (i in 1:length(sps)) {
  curtips = inds[ which(isps == sps[i]) ]
  ylim = range(which(t2$tip.label %in% curtips))
  lines(x = c(0.022, 0.022), y = c(ylim[1] + 0.05, ylim[2] - 0.05))
  text(x = 0.0222, y = mean(ylim), 
       paste0("C. ", sps[i]), font = 3, adj = c(0, 0.5), cex = 0.7)
}
for (i in (Ntip(t2) + 1):(Ntip(t2) + Nnode(t2))) {
  curlins = isps[ phangorn::Descendants(t2, i, type = "tips")[[1]]  ]
  curtips = inds[ phangorn::Descendants(t2, i, type = "tips")[[1]]  ]
  
  plotnode = FALSE
  if (length(unique(curlins)) > 1 ) {
    plotnode = TRUE
  } else if (length(unique(curlins)) == 1) {
    sp = unique(curlins)
    sptips = inds[which(isps == sp)]
    if (setequal(sptips, curtips)) {
      plotnode = TRUE
    }
  }
  if (plotnode & t2$node.label[i - Ntip(t2)] > 95) {
   nodelabels("", i, frame = "none", pch = 21, bg = "white", cex = 0.5) 
  }
}
dev.off()