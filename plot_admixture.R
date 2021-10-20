# https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
cols = c('#e6194b', '#3cb44b', '#ffe119', 
         '#4363d8', '#f58231', '#911eb4', 
         '#46f0f0', '#f032e6', '#bcf60c', 
         '#fabebe', '#008080', '#e6beff', 
         '#9a6324', '#fffac8', '#800000', 
         '#aaffc3', '#808000', '#ffd8b1', 
         '#000075', '#808080', '#ffffff', 
         '#000000')

setwd("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/admixture/nodiv/")
sampfile = "../Ctenotus_zastictus.project2.csv"

d = read.csv(sampfile, stringsAsFactors = F)

pdf("~/Desktop/test.pdf", width = 10, height = 4)
for (n in 2:10) {
  file = paste("duricola_nodiv.miss0.7.MAC0.thinned.", n, ".Q", sep = "")
  k = read.table(file, stringsAsFactors = F)
  inds = read.table(paste("duricola_nodiv.miss0.7.MAC0.thinned.fam", sep = ""))
  inds = data.frame(inds = inds$V2, 
                    sp = d[match(inds$V2, d$SAMPLE_ID), "LINEAGE"], 
                    stringsAsFactors = F)
  rownames(k) = paste(inds$sp, inds$inds)
  k = k[order(rownames(k)), ]
  
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
    graphics::text(labels[i], y = 1.05, x = mean(c(xmin, xmax)),
                   srt = 90, adj = c(0, 0.5))
  }
}
dev.off()