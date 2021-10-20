library(ape)
library(ggplot2)
library(cowplot)
library(dplyr)
theme_set(theme_cowplot())
setwd("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/")

x = read.csv("Ctenotus_zastictus.project4.csv", stringsAsFactors = F)
x = x[complete.cases(x$mtDNA_LINEAGE), ]
x1 = x[which(x$mtDNA_LINEAGE == x$nDNA_LINEAGE | !complete.cases(x$nDNA_LINEAGE)), ]
mt = read.dna("atlas_mtDNA/Ctenotus_zastictus_mtDNA.cytb.2Sept21.fasta.aln",
              format = "fasta")
mt1 = mt[which(attr(mt, "dimnames")[[1]] %in% x1$SAMPLE_ID), ]
mtdist = dist.dna(mt1, model = "T92", pairwise.deletion = T, as.matrix = T)
x2 = x1[x1$SAMPLE_ID %in% attr(mt1, "dimnames")[[1]], ]

keep = c("zastictus", "duricola", "pallasotus",
         "piankai", "serventyi", "rhabdotus", "hanloni",
         "dux", "quattuordecimlineatus2", "atlas", "ariadnae",
         "maryani")

spcombo = combn(keep, 2)
mtdists = data.frame(sp1 = spcombo[1, ],
                     sp2 = spcombo[2, ],
                     mtdiv = NA,
                     stringsAsFactors = F)
for (i in 1:nrow(mtdists)) {
    sp1 = mtdists[i, "sp1"]
    sp2 = mtdists[i, "sp2"]
    inds1 = x2[x2$mtDNA_LINEAGE == sp1, "SAMPLE_ID"]
    inds2 = x2[x2$mtDNA_LINEAGE == sp2, "SAMPLE_ID"]
    tmpdist = mtdist[inds1,inds2]
    mtdists[i, "mtdiv"] = mean(tmpdist)
}
mtdist2 = mtdists %>% group_by(sp1) %>% 
  slice_min(n=1, order_by = mtdiv) %>%
  ungroup() %>% 
  mutate(rank = dense_rank(-desc(mtdiv)))

d = read.csv("divergence/atlas_group.dxy.8Sept21.csv",
             stringsAsFactors = F)
d = d %>% filter(sp1 %in% keep)
d1 = d %>% group_by(sp1) %>% 
  slice_min(n=1, order_by = dxy) %>%
  ungroup() %>% 
  mutate(rank = dense_rank(-desc(dxy)))
d1$zast = "no"
d1[grep("zast", d1$sp1), "zast"] = "yes"
a = ggplot(d1, aes(rank, dxy)) +
  geom_jitter(aes(fill = zast), shape = 21, size = 2) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylab(expression(d[xy])) +
  scale_fill_manual(values = c("#377eb8", "#e41a1c"))
  

f = read.csv("divergence/atlas_group.fst.8Sept21.csv",
             stringsAsFactors = F)
f = f %>% filter(sp1 %in% keep)
f1 = f %>% group_by(sp1) %>% 
  slice_min(n=1, order_by = fst) %>%
  ungroup() %>% 
  mutate(rank = dense_rank(-desc(fst)))
f1$zast = "no"
f1[grep("zast", f1$sp1), "zast"] = "yes"
b = ggplot(f1, aes(rank, fst)) +
  geom_jitter(aes(fill = zast), shape = 21, size = 2) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylab(expression(F[ST])) +
  scale_fill_manual(values = c("#377eb8", "#e41a1c"))

ab = plot_grid(a, b, labels = c("A", "B"))
save_plot("~/Desktop/fst_dxy.pdf", ab,
          ncol = 2, base_height = 2.75, base_width = 4)

f2 = inner_join(f1, d1, by = c("sp1" = "sp1"))
f3 = left_join(f2, mtdist2, by = c("sp1" = "sp1"))
a = ggplot(f3, aes(dxy, fst)) +
  geom_point(aes(fill = zast.x), shape = 21, size = 2) +
  theme_classic() +
  theme(legend.position = "none") +
  ylab(expression(F[ST])) +
  xlab(expression("nuclear " ~ d[xy])) +
  scale_fill_manual(values = c("#377eb8", "#e41a1c"))
b = ggplot(f3, aes(dxy, mtdiv)) +
  geom_point(aes(fill = zast.x), shape = 21, size = 2) +
  theme_classic() +
  theme(legend.position = "none") +
  xlab(expression(nuclear ~ d[xy])) +
  ylab(expression(mtDNA ~ d[xy])) +
  scale_fill_manual(values = c("#377eb8", "#e41a1c"))
ggplot(f3, aes(fst, mtdiv)) +
  geom_point(aes(fill = zast.x), shape = 21, size = 2) +
  theme_classic() +
  theme(legend.position = "none") +
  xlab(expression(F[ST])) +
  ylab(expression(mtDNA ~ d[xy])) +
  scale_fill_manual(values = c("#377eb8", "#e41a1c"))
ab = plot_grid(a, b, labels = c("A", "B"))
save_plot("../figures/Figure6_dxy_fst.png", ab,
          ncol = 2, base_height = 2.75, base_width = 3.5)
