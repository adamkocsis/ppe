# Script written to illustrate the use of GPMs in Paleobiology
# Ádám T. Kocsis, 2026-04-15, Erlangen

library(rgplates) # tectonic reconstructions, requires GPLATES desktop app
library(chronosphere) # data acquisition
library(divDyn) # geologic timescale and data binning
library(icosa) # icosahedral gridding
library(rampage) # color ramp
library(vegan) # distance metric
library(viridisLite) # color ramp
library(predicts) # MaxEnt wrapper
library(smoothr) # resampling shapes

# working directory
setwd("/mnt/sky/Dropbox/WorkSpace/2025-06-10_PPE/")

################################################################################
# Data preparation
################################################################################

dir.create("data/chronosphere", showWarnings=FALSE)

# PALEOMAP model
PALEOMAP <- fetch("paleomap", "model", datadir="data/chronosphere", ver="v19o_r1c")

# PALEOMAP DEM
dems <- fetch("paleomap", "dem", res=0.1, datadir="data/chronosphere", ver="v24221")

# PALEOMAP Paleocoastlines
pc<- fetch("paleomap", "paleocoastlines", datadir="data/chronosphere", ver="7")

# BRIDGE SST reconstruction HadCM3- interpolated to work with Paleocoastlines
ssinter<- fetch("SOM-kocsis-provinciality", datadir="data/chronosphere", ver="v1.0")

# PBDB data
pbdb<- fetch("pbdb", datadir="data/chronosphere", ver="20260412")


# the hexagonl grid to be used
hex <- icosa::hexagrid(deg=8, sf=TRUE)

################################################################################
# PBDB Data preparation (crude)
# taxomomic subset
dat <- pbdb[pbdb$class=="Trilobita", ]

# get the devonian bunch
# Stratigraphy
data(stages)
data(keys)

# B. the stg entries (lookup)
stgMin <- divDyn::categorize(dat[,"early_interval"],keys$stgInt)
stgMax <- divDyn::categorize(dat[,"late_interval"],keys$stgInt)

stgMin <- as.numeric(stgMin)
stgMax <- as.numeric(stgMax)

# empty container
dat$stg <- rep(NA, nrow(dat))
# select entries, where
stgCondition <- c(
# the early and late interval fields indicate the same stg
	which(stgMax==stgMin),
# or the late_intervarl field is empty
	which(stgMax==-1))
dat$stg[stgCondition] <- stgMin[stgCondition]

# Lower devonian subset (Emsian stage)
ld <- c(31)
mid <- (max(stages$bottom[ld])+min(stages$top[ld]))/2

# the occurrencs
occs <- dat[which(dat$stg%in%ld), ]

# closest product date 410 Ma
prod <- divDyn::matchtime(as.numeric(names(ssinter)),mid)

# stage mid coordinates
midCoords <- rgplates::reconstruct(occs[, c("lng", "lat")], age=mid, model=PALEOMAP)
colnames(midCoords) <- c("mid_plng", "mid_plat")

# product  coordinates
prodCoords <- rgplates::reconstruct(occs[, c("lng", "lat")], age=prod, model=PALEOMAP)
colnames(prodCoords) <- c("prod_plng", "prod_plat")

# join them
occs <- cbind(occs,midCoords, prodCoords)

# Get the reconsturctions (mid age and product)
midPlates<- reconstruct("static_polygons",  age=mid, model=PALEOMAP)
prodPlates<- reconstruct("static_polygons",  age=prod, model=PALEOMAP)

################################################################################
# a. Visualization of coordinates (plates + coastlines)
################################################################################

# projection used on non long-lat maps (robinson)
proj <- "ESRI:54030"

# colors
bgAll <- "#0044DDbb"
colAll <- "#00000099"
colFocal <- "#FF4400bb"
bgFocal <- "#00FF00BB"
colFocal <- "#ffffff99"

#original
colorAll <- "#0044DDbb"
colorFocal <- "#FF4400bb"

# occurrences for mid
occsMid <- occs[!is.na(occs$mid_plng), ]

# the edge of the map!
me <- rgplates::mapedge()

# Visualize
dir.create("export/gpmuse/", showWarnings=FALSE)
png("export/gpmuse/basics.png", width=2000, height=1000, pointsize=24, bg=NA)
	par(mai=rep(0.1, 4))
	plot(me, col="white")
	plot(midPlates$geometry, col="#BBBBBBBB", border=NA, add=TRUE)
	# the points (all trilobites)
	plot(me, add=TRUE, border="black")
	segments(x0=-180, y0=seq(-90, 90, 15), x1=180, y1=seq(-90, 90, 15), lty=2, lwd=2, col="gray80")
	segments(x0=-180, y0=c(-45, 0, 45), x1=180, y1=c(-45, 0, 45), lty=2, lwd=2, col="gray30")
	points(unique(occsMid[, c("mid_plng", "mid_plat")]), col=colAll,bg=bgAll, pch=21, cex=1.5)
	plot(me, border="black", lwd=4, add=TRUE, col=NA)
dev.off()


# the Projected versions
	# reprojected edge of the map
	meProj <- mapedge(crs=proj)
	# reprojected static polygons
	midPlatesProj <- st_transform(smoothr::densify(midPlates), crs=proj)

	# occurrence reprojected coordinates: add to occsMid
	occsMidProj <- st_as_sf(occsMid[, c("collection_no", "mid_plng", "mid_plat")], coords=c("mid_plng", "mid_plat"), crs="WGS84")
	occsMidProj <- st_transform(occsMidProj, crs=proj)
	midProjCoords <- st_coordinates(occsMidProj)
	colnames(midProjCoords) <- c("mid_plng_proj", "mid_plat_proj")
	occsMid<- cbind(occsMid, midProjCoords)

	# reprojet the graticules
	grat <- data.frame()
	gratLats <- seq(-90,90, 15)
	for(i in 1:length(gratLats)){
		ps <- seq(-180, 180, 1)
		oneLat<- cbind(long=ps, lat=rep(gratLats[i], length(ps)))
		grat<- rbind(grat, oneLat)
	}
	gratProj <- st_as_sf(grat, coords=c("long", "lat"), crs="WGS84")
	gratProj <- st_transform(gratProj, crs=proj)

	# Visualize
	png("export/gpmuse/basicProj.png", width=2000, height=1000, pointsize=24, bg=NA)
		par(mai=rep(0.1, 4))
		plot(meProj, col="white")
		plot(midPlatesProj$geometry, col="#BBBBBBBB", border=NA, add=TRUE)
		plot(hex, col=NA, border="gray60", add=TRUE, crs=proj)
		for(i in seq(-90, 90, 15)) lines(st_coordinates(gratProj)[grat[,2]==i,],  lty=2, lwd=2, col="gray75")
		for(i in c(-45,45)) lines(st_coordinates(gratProj)[grat[,2]==i,], lty=2, lwd=2, col="gray30")
		lines(st_coordinates(gratProj)[grat[,2]==0,], lwd=4, col="gray30")
		points(unique(occsMid[, c("mid_plng_proj", "mid_plat_proj")]), col=colAll,bg=bgAll, pch=21, cex=1.5)
		plot(meProj, border="black", lwd=4, add=TRUE, col=NA)
	dev.off()

################################################################################
# A. Range visuaziation
################################################################################

# focal genus
focalGen <- "Kettneraspis"

# photo from - need to redraw
# https://www.fossilera.com/fossils/1-65-kettneraspis-trilobite-with-long-occipital-spectacular-prep?srsltid=AfmBOoo8ZQhAKxKupcFtj45E6DzpLtRHnTSu-31r2-t8CQYFos8NyIf3

# the color of the occupied cell
colCell <- "#DD330066"

# 1. icosa grid occupancy

# add the located cell
occsMid$cell <- icosa::locate(hex, occsMid[, c("mid_plng", "mid_plat")])

# the occurrences of the focal taxon
focalMid <- occsMid[which(occsMid$genus==focalGen), ]

png("export/gpmuse/onerange.png", width=2000, height=1000, pointsize=24, bg=NA)
	par(mai=c(0.1,0.1, 1, 0.1))
	plot(me, col="white")
	plot(midPlates$geometry, col="#BBBBBBBB", border=NA, add=TRUE)
	plot(hex, col=NA, border="gray60", add=TRUE)
	plot(hex, unique(focalMid$cell), col=colCell, border=substr(colCell, 1,7), add=TRUE, lwd=3)
	points(unique(occsMid[, c("mid_plng", "mid_plat")]), col=colAll,bg=bgAll, pch=21, cex=1.5)
	points(unique(focalMid[, c("mid_plng", "mid_plat")]), col=colFocal,bg=bgFocal, pch=23, cex=2)
	mtext(side=3, line=0.5, text=paste0(focalGen, ", cell occupancy = ", length(unique(focalMid$cell))), cex=2)
	plot(me, border="black", lwd=4, add=TRUE, col=NA)
dev.off()


# the reprojected version
png("export/gpmuse/onerangeProj.png", width=2000, height=1000, pointsize=24, bg=NA)
	par(mai=c(0.1,0.1, 1, 0.1))
	plot(meProj, col="white")
	plot(midPlatesProj$geometry, col="#BBBBBBBB", border=NA, add=TRUE)
	plot(hex, col=NA, border="gray60", add=TRUE, crs=proj)
	plot(hex, unique(focalMid$cell), col=colCell, border=substr(colCell, 1,7), add=TRUE, lwd=3, crs=proj)
	points(unique(occsMid[, c("mid_plng_proj", "mid_plat_proj")]), col=colAll,bg=bgAll, pch=21, cex=1.5)
	points(unique(focalMid[, c("mid_plng_proj", "mid_plat_proj")]), col=colFocal,bg=bgFocal, pch=23, cex=2)
#	mtext(side=3, line=0.5, text=paste0(focalGen, ", cell occupancy = ", length(unique(focalMid$cell))), cex=2)
	plot(meProj, border="black", lwd=4, add=TRUE, col=NA)
dev.off()

# 2. histogram of ranges
# calculate histogram of cell occupancies
occDist <- tapply(X=occsMid$cell, INDEX=occsMid$genus, FUN=function(x) length(unique(x)) )

# plot
png("export/gpmuse/rangeHist.png", width=1000, height=1000, pointsize=24)
	hist(occDist, main="Distribution of cell occupancies of genera", xlab="Cell Occupancy", col=colorAll, border="white")
	box()
dev.off()


################################################################################
# B. Latitudinal patterns
################################################################################

## # Trilobite occurrences from the Emsian stage
# collection table
collLat <- unique(occsMid[, c("collection_no", "mid_plat")])

# Collection-level alpha richness
alpha <- tapply(X=occsMid$genus, INDEX=occsMid$collection_no, FUN=function(x) length(levels(factor(x))))

collLat$alpha <- alpha[as.character(collLat$collection_no)]
## plot(y=collLat$mid_plat,x=collLat$alpha, col="gray", pch=16, ylim=c(-90, 90), ylab="Latitude", xlab="Collection alpha richness (Unstandardized)", )

## library(mgcv)
## collMod <- gam(alpha ~ s(mid_plat), data = collLat)
## newLat <- seq(-90, 90,0.1)
## gamMod <- predict(collMod, newdata=data.frame(mid_plat=newLat))
## lines(x=gamMod, y=newLat, col="darkred")
## # latOccDensity <- hist(occsMid$mid_plat, breaks=seq(-90, 90, 5))
## latCollDensity <- hist(collLat$mid_plat, breaks=seq(-90, 90, 5))


# Total richness in grid cell as a function of latitude
gridAlpha <- tapply(X=occsMid$genus, INDEX=occsMid$cell, FUN=function(x) length(levels(factor(x))))
gridCenters <- as.data.frame(centers(hex))
gridCenters$alpha <- gridAlpha[rownames(gridCenters)]
gridCenters <- gridCenters[!is.na(gridCenters$alpha), ]

# visualize
# plot(y=gridCenters$lat,x=gridCenters$alpha, col="gray", pch=16, ylim=c(-90, 90), ylab="Latitude", xlab="Richness in grid cells (Unstandardized)", )

## library(mgcv)
## gridMod <- gam(alpha ~ s(lat), data = gridCenters)
## newLat <- seq(-90, 90,0.1)
## gamMod <- predict(gridMod, newdata=data.frame(lat=newLat))
## lines(x=gamMod, y=newLat, col="darkred")
# latOccDensity <- hist(occsMid$mid_plat, breaks=seq(-90, 90, 5))

# create latitudinal bins
bin <- 15
boundaries <- seq(-90, 90, bin)
#abline(h=boundaries,lty=2, col="gray60")

# the midpoints
mids <- (boundaries[2:length(boundaries)-1]+ boundaries[2:length(boundaries)])/2

# what to visualize
latMeans <- rep(NA, length(boundaries)-1) # mean

# which latitudinal bins do the cells belong? (order is based on boundaries!)
whichBin <- cut(gridCenters$lat, boundaries, labels=FALSE)

# calculate the mean of the bins
binMeans <- tapply(INDEX=whichBin, X=gridCenters$alpha, mean, na.rm=TRUE)

# the latitudinal averages
latMeans[as.numeric(names(binMeans))] <- binMeans
#lines(x=latMeans, y=mids, type="o")

# the number of grid cells occupied
# gridHist<- hist(gridCenters$lat, breaks=boundaries,plot=FALSE)
# lines(y=gridHist$mids, x=gridHist$counts*3, type="o", col="red")

# occurrence density
latOccDensity <- hist(occsMid$mid_plat, breaks=boundaries, plot=FALSE)
# lines(y=latOccDensity$mids, x=latOccDensity$counts/10, type="o", col="blue")

# collection density
latCollDensity <- hist(collLat$mid_plat, breaks=boundaries, plot=FALSE)
# lines(y=latCollDensity$mids, x=latCollDensity$counts/2, type="o", col="red")


png("export/gpmuse/latitude.png", width=800, height=2000, pointsize=24)
	par(mai=c(2.1, 1, 2.1, 1))
	plot(
		y=gridCenters$lat,x=gridCenters$alpha, col="gray", pch=16, ylim=c(-90, 90),
		ylab="", xlab="Richness in grid cells (Unstandardised)", axes=FALSE, yaxs="i", cex.lab=1.5)
	axis(2, at=seq(-90, 90, 30), label=seq(-90, 90, 30))
	axis(1)
	axis(3, at=seq(0, 160, 20), label=seq(0, 160, 20)*10, col="blue", col.ticks="blue", col.axis="blue")
	lines(y=latOccDensity$mids[latOccDensity$counts>0], x=latOccDensity$counts[latOccDensity$counts>0]/10, type="o", col="blue", pch=16, cex=1.5)
	lines(x=latMeans, y=mids, type="o", pch=16, cex=1.5)
	abline(h=boundaries,lty=2, col="gray60")
	box()
	mtext(text="Occurrence records in latitudinal bins", col="blue", side=3, line=3, cex=1.5)
dev.off()

################################################################################
# C. Climate model use
################################################################################
# grab sst reconstruction closest to the stage mid
sst <- ssinter[as.character(prod)]

# There is an interpolation artifact at the southermost coast of Gondwana, leading
# to erroneous -10 degrees SST. These are corrected here manually
values(sst)[which(values(sst)< -2)] <- -2

# repojrect occurrence data for this target age
occsProd <- occs[!is.na(occs$prod_plng), ]
occsProdProj <- sf::st_as_sf(occsProd[, c("collection_no", "prod_plng", "prod_plat")], coords=c("prod_plng", "prod_plat"), crs="WGS84")
occsProdProj <- sf::st_transform(occsProdProj, crs=proj)
prodProjCoords <- sf::st_coordinates(occsProdProj)
colnames(prodProjCoords) <- c("prod_plng_proj", "prod_plat_proj")
occsProd<- cbind(occsProd, prodProjCoords)

# the paleocoastlines
coast <- pc[as.character(prod),"coast"]
margin <- pc[as.character(prod),"margin" ]

# the occurrences
occsProd <- occsProd[which(occsProd$stg==31), ]
focalProd <- occsProd[which(occsProd$genus==focalGen), ]

png("export/gpmuse/sst.png", width=2000, height=1000, pointsize=24, bg=NA)
	plot(sst, col=gradinv(256), box=TRUE, axes=FALSE)
	plot(margin, col="#BBBBBB88", add=TRUE, border=NA)
	plot(coast, col="gray90", add=TRUE, border="gray90", lwd=3)
	points(unique(occsProd[, c("prod_plng", "prod_plat")]),  bg=bgAll,col=colAll, pch=21, cex=1.5, lwd=2)
	points(unique(focalProd[, c("prod_plng", "prod_plat")]), bg=bgFocal, col=colFocal, pch=23, cex=2, lwd=3)
	plot(me, border="black", lwd=4, add=TRUE, col=NA)
dev.off()


# reprojected version of sst
crs(sst) <- "WGS84"
sstProj <- terra::project(sst, proj)

# reproject paleocoastlines
coastProj <- sf::st_transform(smoothr::densify(coast), crs=proj)
marginProj <- sf::st_transform(smoothr::densify(margin), crs=proj)

png("export/gpmuse/sstProj.png", width=2000, height=1000, pointsize=24, bg=NA)
	plot(sstProj, col=gradinv(256), box=FALSE, axes=FALSE, ylim=c(st_bbox(meProj)[c('ymin', 'ymax')]))
	plot(marginProj, col="#BBBBBB88", add=TRUE, border=NA)
	plot(coastProj, col="gray90", add=TRUE, border="gray90", lwd=3)
	points(unique(occsProd[, c("prod_plng_proj", "prod_plat_proj")]),  bg=bgAll,col=colAll, pch=21, cex=1.5, lwd=2)
	points(unique(focalProd[, c("prod_plng_proj", "prod_plat_proj")]), bg=bgFocal, col=colFocal, pch=23, cex=2, lwd=3)
	plot(meProj, border="black", lwd=4, add=TRUE, col=NA)
dev.off()

# histogram
focalSST <- terra::extract(sst, focalProd[, c("prod_plng", "prod_plat")])
focalSST <- focalSST[,"sst_inter_1deg_400.nc"]
#hist(focalSST, breaks=seq(-10, 40, 2.5), xlim=c(-5, 35))

allSST <- terra::extract(sst, occsProd[, c("prod_plng", "prod_plat")])
allSST <- allSST[,"sst_inter_1deg_400.nc"]
#hist(allSST[,"sst_inter_1deg_400.nc"], breaks=seq(-10, 40, 2.5), xlim=c(-5, 35))

# calculate the densities of these
allDensity <- density(allSST[!is.na(allSST)])
focalDensity <- density(focalSST[!is.na(focalSST)])

#need to make Kettnerapsis italic - later!
png("export/gpmuse/temphist.png", width=1000, height=1000, pointsize=24, bg=NA)
	plot(NULL, NULL, xlim=c(-2, 40), ylim=c(0,max(c(allDensity$y, focalDensity$y))*1.05), yaxs="i", xaxs="i", axes=FALSE, ylab="", xlab="Mean Annual SST (°C)")
	usr <- par()$usr
	rect(ybottom=usr[3], ytop=usr[4], xleft=usr[1], xright=usr[2], col="white", border="black", lwd=3)
	polygon(allDensity$x, allDensity$y, col=colorAll)
	polygon(c(min(focalDensity$x)-1, focalDensity$x), c(0,focalDensity$y), col=colorFocal)
	axis(1)
	legend("topleft", bty="n", legend=c("All trilobites", focalGen), fill=c(colorAll, colorFocal), cex=1.5, inset=c(0.1,0.1))
dev.off()


################################################################################
# Bioregionalization (mid)
################################################################################

# a bunch of random colors
load("data/allHex.RData")

# plotting the script is here
source("code/methods/plots.R")

# get rid of those where either the genus or cell is missing
occsMidUse <- occsMid[!is.na(occsMid$genus) & !is.na(occsMid$cell),]

# contingency matrix
cont <- table(occsMidUse$cell, occsMidUse$genus)

# incidence
cont[cont>1] <- 1

# Method 1. Compositional similarity
distmat <- vegan::vegdist(cont, method="jaccard")

# clustering
cluster <- hclust(distmat, "ward.D2")

# plot this
plot(cluster)
h <- 1.1
abline(h=h, col="red")

# cutting the dendrogram-> membership vector
mem <- cutree(cluster, h=h)

# long-lat version
png("export/gpmuse/bioregionalization.png", width=2000, height=1000, pointsize=24, bg=NA)
	par(mai=rep(0.1, 4))
	plot(me, col="white")
	biogeoplot(mem=mem, bg=midPlates$geometry, gri=hex, ylim=c(-90, 90), add=TRUE, alpha="66")
	plot(me, border="black", lwd=4, col=NA, add=TRUE)
dev.off()

# the reprojected version
png("export/gpmuse/bioregionalizationProj.png", width=2000, height=1000, pointsize=24, bg=NA)
	par(mai=rep(0.1, 4))
	plot(meProj, col="white")
	biogeoplot(mem=mem, bg=midPlates$geometry, gri=hex, ylim=c(-90, 90), add=TRUE, alpha="66", crs=proj)
	plot(meProj, border="black", lwd=4, col=NA, add=TRUE)
dev.off()

################################################################################
# DEM (prod)
################################################################################
# select the appropriate DEm
emsianDEM<- dems[as.character(prod)]

# plot longlat
png("export/gpmuse/dem.png", width=2000, height=1000, pointsize=24, bg=NA)
	par(mai=rep(0.1, 4))
	plot(emsianDEM, col=paleomap$col, breaks=paleomap$breaks, legend=FALSE, axes=FALSE)
	plot(me, col=NA, border="black", lwd=4,add=TRUE)
dev.off()

# plot reprojected
png("export/gpmuse/demProj.png", width=2000, height=1000, pointsize=24, bg=NA)
	par(mai=rep(0.1, 4))
	plot(project(emsianDEM, proj), col=paleomap$col, breaks=paleomap$breaks, legend=FALSE, axes=FALSE)
	plot(meProj, col=NA, border="black", lwd=4,add=TRUE)
dev.off()

################################################################################
# Maxent?
################################################################################
# get the monthly data - coarser resolution!
allClim <- list.files("data/400_teXPc")

# select ocean surface layers
allClim <- allClim[grep("pfc", allClim)]

# stack of rasters
allMonths <- terra::rast(paste0("data/400_teXPc/", allClim) )

# get only the temperature variable
allTemp <- allMonths[[which(varnames(allMonths)=="temp_mm_uo")]]

# check how it relates to the resampled version above
## resst <- resample(mean(allTemp), sst)
## plot(resst-sst)

# rotate by 180 deg
allTemp <- terra::rotate(allTemp)

# calculate: annual mean, annual minimum, annual maximum
minTemp <- min(allTemp)
maxTemp <- max(allTemp)
meanTemp <- mean(allTemp)
layers <- c(meanTemp, minTemp, maxTemp)

# the occurrence record coordinate data
dat <- focalProd[, c("prod_plng", "prod_plat")]

# if you want to show points
#points(unique(focalProd[, c("prod_plng", "prod_plat")]), bg=bgFocal, col=colFocal, pch=23, cex=2, lwd=3)

# do a quick, dirty maxent with the 3 predictors
# Latches on the rJava-run Maxent jar
xf <- predicts::MaxEnt(layers, dat) # loads of misssing values...

# prediction suitability
pred <- predicts::predict(xf, layers)

# longlat plot
# plot(pred)

# set CRS for reprojection
crs(pred) <- "WGS84"

# mask out the non-shelf area
masked <- terra::mask(pred,sf::st_as_sf(margin))
## plot(masked)


# grab calculated threshold (based on training set!)
res <- read.csv(file.path(xf@path, "maxentResults.csv"))
tr <- res[1,"Maximum.training.sensitivity.plus.specificity.Cloglog.threshold"]

# do thresholding
one <- masked>tr
terra::values(one)[!terra::values(one)] <- NA

# create threshold-basd color ramp
colDat <- data.frame(color=rev(viridisLite::viridis(3)),z=c(max(values(masked), na.rm=TRUE), tr,0 ))
suitCols <- rampage::expand(colDat, n=256)

# plot longlat
png("export/gpmuse/maxent_suitability.png", width=2000, height=1000, pointsize=24, bg=NA)
	plot(masked, axes=FALSE, xlim=c(-180, 180), ylim=c(-90, 90), col=suitCols$col, breaks=suitCols$breaks, type="continuous")
	plot(margin, add=TRUE, col="gray", border=NA)
	plot(coast, col="black", add=TRUE)
	plot(masked, axes=FALSE, add=TRUE, legend=FALSE, col=suitCols$col, breaks=suitCols$breaks)
	plot(me, lwd=4, col=NA, border="black", add=TRUE)
	points(unique(focalProd[, c("prod_plng", "prod_plat")]), bg=bgFocal, col=colFocal, pch=23, cex=2, lwd=3)
dev.off()

# plot reprojected
png("export/gpmuse/maxent_suitabilityProj.png", width=2000, height=1000, pointsize=24, bg=NA)
	plot(project(masked, proj), axes=FALSE, xlim=st_bbox(meProj)[c("xmin", "xmax")], ylim=st_bbox(meProj)[c("ymin", "ymax")], col=suitCols$col, breaks=suitCols$breaks, type="continuous")
	plot(meProj, col="white", add=TRUE)
	plot(marginProj, add=TRUE, col="gray", border=NA)
	plot(coastProj, col="black", add=TRUE)
	plot(project(masked, proj), axes=FALSE, add=TRUE, legend=FALSE, col=suitCols$col, breaks=suitCols$breaks)
	plot(meProj, lwd=4, col=NA, border="black", add=TRUE)
	points(unique(focalProd[, c("prod_plng_proj", "prod_plat_proj")]), bg=bgFocal, col=colFocal, pch=23, cex=2, lwd=3)
dev.off()
