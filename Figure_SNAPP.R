library(ape)

t = read.nexus("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/bfd_snapp/snapp/snap.tre")

cols = brewer.pal(6, "Paired")
sps = c("duricola", "pallasotus", "piankai", "rhabdotus", 
        "serventyi", "zastictus")

boxlabel<-function(x,y,text,cex=1,bg="transparent",offset=0){
  w<-strwidth(text)*cex*1.1
  h<-strheight(text)*cex*1.4
  os<-offset*strwidth("W")*cex
  rect(x+os,y-0.9*h,x+w+os,y+0.9*h,col=bg,border=0)
  text(x,y,text,pos=4,offset=offset,font=3)
}

pdf("~/Desktop/snapp.pdf", height = 2, width = 3)
par(mar = c(0, 0, 0, 5), xpd = T)
ape::plot.phylo(t, show.tip.label = F)
pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
N<-Ntip(tree)
for (i in 1:length(t$tip.label)) {
  col = cols[which(t$tip.label == sps[i])]
  cursp = paste("C. ", t$tip.label[i], sep = "")
  boxlabel(pp$xx[i],pp$yy[i],cursp,bg = alpha(col, 0.5))
}
for (i in (Ntip(t) + 1):(Ntip(t) + Nnode(t))) {
  nodelabels("", i, pch = 21, bg = "white", frame = "none")
}
dev.off()