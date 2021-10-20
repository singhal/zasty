require(ggplot2)
require(ggtree)
require(treeio)
require(RColorBrewer)

x <- treeio::read.iqtree("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/atlas_mtDNA/final.rooted.treefile")
d = read.csv("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/Ctenotus_zastictus.project4.csv",
             stringsAsFactors = F)
cols = rep(brewer.pal(12, "Paired"), 3)

# drop bad tips
x1 = treeio::drop.tip(x, c("SAMR_57650_Ct_atla", "SAMR_55210_Ct_atla", "WAMR_102159_Ct_duri", "WAMR_139423_Ct_duri"))
lins = d[match(x1@phylo$tip.label, d$SAMPLE_ID), "mtDNA_LINEAGE"]
names(lins) = x1@phylo$tip.label 

nodenum = seq(Ntip(x1@phylo) + 1, Ntip(x1@phylo) + Nnode(x1@phylo))
x1@data$upper = FALSE
tips = x1@phylo$tip.label
for (i in nodenum) {
  curtips = tips[ phangorn::Descendants(x1@phylo, i, type = "tips")[[1]] ]
  curlins = lins[curtips]
  
  if (length(unique(curlins)) > 1 ) {
    x1@data[which(x1@data$node == i), "upper"] = TRUE
  } else if (length(unique(curlins)) == 1) {
    sp = unique(curlins)
    sptips = names(lins)[which(lins == sp)]
    if (setequal(sptips, curtips)) {
      x1@data[which(x1@data$node == i), "upper"] = TRUE
    }
  }
}


# Mark everyone as "non_selected
x1@data$posterior_val <- FALSE

# In my case, "prob_percent" is the posterior probability.
# Just tag nodes with more than 90 percent:
x1@data$posterior_val[x1@data$SH_aLRT > 95 & x1@data$upper == TRUE] <- TRUE

# label clades
lins = d[match(x1@phylo$tip.label, d$SAMPLE_ID), "mtDNA_LINEAGE"]
names(lins) = x1@phylo$tip.label

a = ggtree(x1)

alllins = unique(lins)
for (i in 1:length(alllins)) {
  tips = names(lins)[which(lins == alllins[i])]
  # if (length(tips) > 1) {
    mrca = phytools::findMRCA(x1@phylo, tips)
    a = a + geom_balance(node=mrca, fill= cols[i], color= NA,
                         alpha=0.6, extend= 0.1) +
      geom_cladelabel(mrca, alllins[i], fontface = "italic", offset=0, 
                      barsize = NA, angle=0, offset.text=0.04, 
                      align = T, size = 0.7)
  # } else {
  #   mrca = which(x1@phylo$tip.label == tips[1]) 
  #   
  #   a + geom_tiplab(aes(subset=(node==mrca), label= alllins[i], 
  #                       fill = cols[i], geom = "label",  # labels not text
  #                       label.padding = unit(0.15, "lines"), # amount of padding around the labels
  #                       label.size = 0))
  # }
  

}

pdf("~/Desktop/mtDNA_tree.pdf", height = 9, width = 4)
a +
  # add in node labels
  geom_nodepoint(aes(subset = posterior_val),  fill = "white", shape = 21, size = 1.5) + 
  xlim(0, 0.8)
dev.off()
