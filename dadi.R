process <- function(x) {
  d = read.table(x, stringsAsFactors = F, sep = "\t", header = T)
  num = gsub(".*Number_", "", x)
  num = gsub("\\..*", "", num)
  names(d)[7] = "params"
  d$num = num
  d$param_num = length(strsplit(d$params[1], ",")[[1]])
  return(d)
}

f = list.files("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/dadi//", pattern = "optimized", full.names = T)
f = f[grep("vic", f, invert = T)]
f = f[grep("founder", f, invert = T)]
f = f[grep("anc", f, invert = T)]
x = lapply(f, process)
d = do.call("rbind", x)
d2 = d %>% group_by(Model, num) %>% slice_max(log.likelihood)
View(d2)
unique(d2$Model)
d3 = d2%>% group_by(Model) %>% slice_max(log.likelihood)
