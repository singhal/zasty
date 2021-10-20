cols = brewer.pal(6, "Paired")
sporder = c("zastictus", "pallasotus", "piankai", "duricola")
d = read.csv("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/Ctenotus_zastictus.project3.csv",
             stringsAsFactors = F)

adir = "~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/admixture/nodiv/"

pdf("~/Dropbox (Personal)/publications/Ctenotus_zasticus/figures/Figure4_admix.pdf",
    height = 6, width = 4)
layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8), ncol = 2))
for (n in 2:4) {
  file = paste(adir, "duricola_nodiv.miss0.7.MAC0.thinned.", n, ".Q", sep = "")
  k = read.table(file, stringsAsFactors = F)
  inds = read.table(paste(adir, "duricola_nodiv.miss0.7.MAC0.thinned.fam", sep = ""))
  inds = data.frame(inds = inds$V2, 
                    sp = d[match(inds$V2, d$SAMPLE_ID), "nDNA_LINEAGE"], 
                    stringsAsFactors = F)
  rownames(k) = paste(inds$sp, inds$inds)
  ksp = inds$sp
  k = k[order(match(ksp, sporder)), ]
  
  par(mar=c(0, 0, 0, 0), xpd=FALSE)
  plot(NULL, ylim=c(0, 1.5), xlim=c(1, nrow(k)), 
       axes = FALSE, xlab="", ylab="")
  for (i in 1:nrow(k)) {
    for (j in 1:ncol(k)) {
      if (j == 1) {
        start = 0
      } else {
        start = sum(k[i, 1:(j-1)])
      }
      end = start + k[i, j]
      polygon(x = c(i, i  + 1, i + 1, i), 
              y = c(start, start, end, end),
              col=cols[j], border=NA)
    }
  }
  
  k$label = gsub(" \\S+$", "", rownames(k))
  labels = unique(k$label)
  for (i in 1:length(labels)) {
    xmin = min(which(k$label == labels[i])) + 0.2
    xmax = max(which(k$label == labels[i])) + 1 - 0.2
    lines(y = c(1.02, 1.02), x = c(xmin, xmax))
    xloc = xmin + 0.3 * (xmax - xmin)
    graphics::text(labels[i], y = 1.07, x = xloc,
                   srt = 45, adj = c(0, 0.5))
  }
}


bdir = "~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/admixture_sim/"

for (n in 2:4) {
  file = paste(bdir, "admix_sim.miss0.7.MAC2.thinned.subsample.", n, ".Q", sep = "")
  k = read.table(file, stringsAsFactors = F)
  inds = read.table(paste(bdir, "admix_sim.miss0.7.MAC2.thinned.subsample.fam", sep = ""))
  inds = data.frame(inds = inds$V2, 
                    sp = d[match(inds$V2, d$SAMPLE_ID), "nDNA_LINEAGE"], 
                    stringsAsFactors = F)
  rownames(k) = paste(inds$sp, inds$inds)
  ksp = inds$sp
  k = k[order(match(ksp, sporder)), ]
  
  par(mar=c(0, 0, 0, 0), xpd=FALSE)
  plot(NULL, ylim=c(0, 1.5), xlim=c(1, nrow(k) + 1), 
       axes = FALSE, xlab="", ylab="")
  for (i in 1:nrow(k)) {
    for (j in 1:ncol(k)) {
      if (j == 1) {
        start = 0
      } else {
        start = sum(k[i, 1:(j-1)])
      }
      end = start + k[i, j]
      polygon(x = c(i, i  + 1, i + 1, i), 
              y = c(start, start, end, end),
              col=cols[j], border=NA)
    }
  }
  
  k$label = gsub(" \\S+$", "", rownames(k))
  labels = unique(k$label)
  for (i in 1:length(labels)) {
    xmin = min(which(k$label == labels[i])) + 0.2
    xmax = max(which(k$label == labels[i])) + 1 - 0.2
    lines(y = c(1.02, 1.02), x = c(xmin, xmax))
    xloc = xmin + 0.3 * (xmax - xmin)
    graphics::text(labels[i], y = 1.07, x = xloc,
                   srt = 45, adj = c(0, 0.5))
  }
}
dev.off()