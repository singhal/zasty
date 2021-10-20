require(ggplot2)
require(ggtree)
require(treeio)
require(RColorBrewer)

x <- read.beast("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/mtDNA/Ctenotus_zastictus.node_mutation.mtDNA.tre")
d = read.csv("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/Ctenotus_zastictus.project3.csv",
             stringsAsFactors = F)
cols = brewer.pal(6, "Paired")

# keep only focal species
lins = d[match(x@phylo$tip.label, d$SAMPLE_ID), "mtDNA_LINEAGE"]
names(lins) = x@phylo$tip.label
core = c("duricola", "pallasotus", "piankai", "serventyi", "rhabdotus", "zastictus")
drop = names(lins)[which(!lins %in% core)]
x1 = drop.tip(x, drop)

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
x1@data$posterior_val[x1@data$posterior > 0.95 & x1@data$upper == TRUE] <- TRUE
x1@data$height_0.95_HPD[x1@data$upper == FALSE] <- NA

# label clades
lins = d[match(x1@phylo$tip.label, d$SAMPLE_ID), "mtDNA_LINEAGE"]
names(lins) = x1@phylo$tip.label

a = ggtree(x1)

for (i in 1:length(core)) {
  tips = names(lins)[which(lins == core[i])]
  mrca = phytools::findMRCA(x1@phylo, tips)
  a = a + geom_balance(node=mrca, fill= cols[i], color= NA,
                       alpha=0.6, extend= 0.1) +
    geom_cladelabel(mrca, core[i], fontface = "italic", offset=0, 
                    barsize = NA, angle=0, offset.text=0.5, 
                    align = T, size = 1)
}

pdf("~/Desktop/mtDNA_tree.pdf", height = 6, width = 4)
# add in bars # add in posterior
a = a + geom_range("height_0.95_HPD", color='gray20', size= 1.5, alpha=.3) +
  # add in node labels
  geom_nodepoint(aes(subset = posterior_val)) +
  # add in scale bar
  theme_tree2()
a1 = revts(a)
a1 + xlim(-10, 2.5)
dev.off()