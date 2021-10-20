require(alphahull)
# require(rgeos)
require(rgdal)
require(maps)
require(mapdata)
require(maptools)
require(rJava)
require(dismo)
require(raster)
require(sp)
data(wrld_simpl)

# source some of pascal's helper scripts
pdir = '/Users/singhal/Dropbox (Personal)/scripts/eco_IBD/pascal_scripts/'
source(paste(pdir, 'ah2sp.R', sep=""))
source(paste(pdir, 'ENM2range.R', sep=""))

# location of raster files
cur_rast_dir = '/Users/singhal/Dropbox (Personal)/macroevolution/eco_IBD_oz/data/AUS_5arc/'

# min_pts
min_pts = 3
# threshold for ENM
# anything in the 99% of the distribution 
# for suitability will be treated as part of the range
threshold <- 0.05

# out dir
out_dir = '/Users/singhal/Dropbox (Personal)/publications/Ctenotus_zasticus/data/niches/'

# ind data 
d = read.csv("~/Dropbox (Personal)/publications/Ctenotus_zasticus/data/Ctenotus_zastictus.project3.csv",
             stringsAsFactors = F)
sps = c("duricola", "pallasotus", "piankai", "rhabdotus", "serventyi", "zastictus")
d1 = d[which(d$mtDNA_LINEAGE %in% sps), ]
d2 = d1[complete.cases(d1$LAT), ]

get_cur_env_data <- function(rast_dir) {
	#load rasters
	setwd(rast_dir)
	env <- stack(list.files(pattern='.tif$'))
	names(env) <- gsub('aus5min_', '', names(env))
	
	# only select bioclim variables
	# don't have historical data for other variables
	env <- env[[paste('bio',1:19,sep='')]]

	#mask so that all rasters have same NA values
	for (i in 2:19) {
		env[[paste('bio',i,sep='')]] <- mask(env[[paste('bio',i,sep='')]], env[['bio1']])
	}
	
	return(env)
}



run_maxent <- function(cur_env, occ) {
	occ = SpatialPointsDataFrame(coords=occ[,c('LON','LAT')],
	                             data=occ, proj4string=CRS('+proj=longlat +datum=WGS84'))
	
	# don't want to run maxent if all the coordinates
	# are in the same grid - no information there
	if (nrow(occ) < 10) {
		# number unique cells
		num_cells = length(unique(cellFromXY(cur_env[[1]], occ)))
		if (num_cells == 1) {
			# will halt things in the next step
			occ = occ[1,]
		}
	}
	
	# only run MaxEnt if you have a certain number of points
	if (nrow(occ) >= min_pts) {
		# figure out what occurrence points are associated with missing data
		x <- extract(cur_env[[c('bio1')]], occ)
		occ <- occ[which(complete.cases(x)),]
		
		# actually run maxent!
		xm <- maxent(x=cur_env, p=occ, args=c("randomseed=true"))

		# thins variables so no overfitting	
		perm <- xm@results
		# gets permutation significance results
		perm <- perm[grep('permutation', rownames(perm)),]
		names(perm) <- gsub('.permutation.importance', '', names(perm))
		# selects only those that are significant above 5
		topVar <- names(which(perm > 5))
		# cannot build a model with a single variable
		# so if there is only one top variable, then add the second
		if (length(topVar) == 1) { 
			topVar <- c(topVar, names(sort(perm, decreasing=TRUE))[2])
		}

		# build the model again with just the topvar
		xm <- maxent(x=cur_env[[topVar]], p=occ, args=c("randomseed=true"))
		# make current prediction
		cur_pred <- predict(xm, cur_env[[topVar]], progress='text')
		
		#least training presence threshold
		# presenceProbs <- extract(cur_pred, occ)
		# thresh <- quantile(presenceProbs,threshold)
		# converts 0 - 1 habitat suitability to binary presence / absence
		# px2 = cur_pred
		# px2[which(values(cur_pred) < thresh)] <- NA
		# px2[which(values(cur_pred) >= thresh)] <- 1
	
		# master <- binaryToPolygon(cur_pred, projection = proj4string(occ), buffer=10000, occ, dropNoOcc=TRUE)
		results = list(cur_pred, xm)
		
		return(results)
	} else {
		# returns NA if too few points
		# or if all points in the same grid cell
		return(list(NA, NA))
	}
}

cur_env = get_cur_env_data(cur_rast_dir)

sps2 = c("pallasotus")
res = vector("list", length(sps2))
for (i in 1:length(sps2)) {
	occ = d2[d2$mtDNA_LINEAGE == sps2[i], ]
	res[[i]] = run_maxent(cur_env, occ)
}
names(res) = sps2

pall = d2[d2$mtDNA_LINEAGE == "pallasotus", ]
pall = SpatialPointsDataFrame(coords=pall[,c('LON','LAT')],
                              data=pall, proj4string=CRS('+proj=longlat +datum=WGS84'))

zast = d2[d2$mtDNA_LINEAGE == "zastictus", ]
zast = zast[zast$SAMPLE_ID != "AMR_123129_Ct_zast", ]
zast = SpatialPointsDataFrame(coords=zast[,c('LON','LAT')],
                             data=zast, proj4string=CRS('+proj=longlat +datum=WGS84'))
csps = c("pallasotus")
pdf("~/Desktop/ENM.pdf", height = 2.5, width = 3)
par(mar = c(0, 0, 0, 3))
for (i in 1:length(csps)) {
	# make map
	plot(res[[csps[i]]][[1]], axes = F, box = F,
	     legend.args = list(text = 'suitability', side = 4, 
	                        font = 2, line = 2.5, cex = 0.8),
	     legend.shrink=0.75,
	     axis.args=list(at= c(0, 0.25, 0.5, 0.75, 0.98),
	                    labels=c(0, 0.25, 0.5, 0.75, 1), 
	                    cex.axis=0.6, tck = -0.5))
  occ = d2[d2$mtDNA_LINEAGE == csps[i], ]
  occ = SpatialPointsDataFrame(coords=occ[,c('LON','LAT')],
                                data=occ, proj4string=CRS('+proj=longlat +datum=WGS84'))
  points(occ, pch = 3, cex = 0.4)
  points(zast, pch = 21, col = "red", cex = 0.8)
  cat(extract(res[[csps[i]]][[1]], zast), "\n")
}
dev.off()

mean(extract(res[[csps[i]]][[1]], pall), na.rm = T)
