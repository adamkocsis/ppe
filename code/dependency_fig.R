# Script written to illustrate the dependency of paleogeographic products
# Ádám T. Kocsis, 2026-06-17, Erlangen
 
library(rgplates)
library(chronosphere)
library(rampage)
library(smoothr) # resampling shapes


setwd("/mnt/sky/Dropbox/WorkSpace/2025-06-10_PPE/")


# The cambrian
age <- 500
model <- "PALEOMAP"
proj <- "ESRI:54009"
meProj <- rgplates::mapedge(crs=proj)

dir.create("export", showWarnings=FALSE)
dir.create("export/dependency", showWarnings=FALSE)


################################################################################
# A. Global plate model
mod <- fetch("paleomap", "model", ver="v3-GPlates")
poly <- reconstruct("static_polygons", age=age, model=mod)
polyer <- smoothr::densify(poly, n=10)
coast <- reconstruct("coastlines", age=age, model=model)
coastProj <- st_transform(coast, crs=proj)
polyProj <- st_transform(polyer, crs=proj)



png("export/dependency/gpm.png", width=3000, height=1500, pointsize=24, bg=NA)
	par(mai=rep(0.1,4))
	plot(meProj, col="#1A6BB0", border="gray90")
	sphereshade(left="#0f3f67", right="#1A6BB0", crs=proj)
	plot(polyProj$geometry, border=NA, col="gray90", add=TRUE)
	plot(coastProj$geometry, border="black", col="gray90", add=TRUE, lwd=2)
	plot(meProj, col=NA, lwd=2, border="gray90", add=TRUE)
dev.off()

png("export/dependency/alter/gpm2.png", width=3000, height=1500, pointsize=24, bg=NA)
	par(mai=rep(0.1,4))
	plot(meProj, col="#1A6BB0", border="gray90")
	plot(coastProj, border="black", col="gray90", add=TRUE)
	sphereshade(left="#000000", right="#000000", crs=proj, left.alpha=0.25, right.alpha=0)
	plot(meProj, col=NA, lwd=2, border="gray90", add=TRUE)
dev.off()
################################################################################
# B. DEM
# for coloring the rasters
dems <- fetch("paleomap", "dem", res=0.1, data="data/chronosphere")
dem <- dems[as.character(age)]

# for better visualization - upscale to ensure good geometry
demRe<- resample(dem, rast(res=0.05))
# and then project
demProj <- project(demRe, proj, mask=TRUE, use_gdal=FALSE)

# a smoother topographic map
topocol <- expand(topos$ptolemy, n=512)

png("export/dependency/dem.png", width=3000, height=1500, pointsize=24, bg=NA)
	par(mai=rep(0.1,4))
	plot(demProj, col=topocol$col, breaks=topocol$breaks, legend=FALSE, smooth=TRUE, axes=FALSE)
	sphereshade(left="#000000", right="#000000", crs=proj, left.alpha=0.25, right.alpha=0)
	plot(meProj, col=NA, border="white", lwd=3, add=TRUE)
dev.off()



################################################################################
# C. Weathering and degassing latitude bands
at <- c(30,-30)
dev <- 8 

weather <- rast(res=0.1)
cells <- xyFromCell(weather, 1:ncell(weather))
oneDensity <- dnorm(x=cells[,2], mean=at[1], sd=dev)
twoDensity <- dnorm(x=cells[,2], mean=at[2], sd=dev)
densities <- oneDensity + twoDensity 

# the vlaues
values(weather) <- densities

weatherProj <- project(weather, proj)


base <- c("#ffffff", "#008e00")
pal <- colorRampPalette(base)
cols <- pal(256)

# version 2 - transparency
## cols<- rgb(red=0, green=142, blue=0, alpha=0:255, maxColorValue=255)

png("export/dependency/weathering.png", width=3000, height=1500, pointsize=24, bg=NA)
	plot(weatherProj, col=cols,  axes=FALSE, legend=FALSE)
	plot(coastProj, border="black", col="#99999966", add=TRUE)
	sphereshade(left="#000000", right="#000000", crs=proj, left.alpha=0.25, right.alpha=0)
	plot(meProj, add=TRUE)
dev.off()

################################################################################
# D. Climate - mean annual temperature
gmst<- chronosphere::fetch("paleomap", "gmst", datadir="data/chronosphere")


landProj <- demProj
values(landProj) <-NA
values(landProj)[values(demProj)>0] <- 1

one <- gmst["500"]
crs(one) <- "WGS84"
oneRe <- resample(one, rast(res=0.5))
gmstProj<- project(oneRe, proj)
png("export/dependency/climate.png", width=3000, height=1500, pointsize=24, bg=NA)
	plot(gmstProj, col=gradinv(256), axes=FALSE, legend=FALSE)
#	plot(coastProj, border="black", col="#99999988", add=TRUE)
	plot(landProj, col="#99999988", axes=FALSE, legend=FALSE, add=TRUE)
	sphereshade(left="#000000", right="#000000", crs=proj, left.alpha=0.15, right.alpha=0)
	plot(meProj, add=TRUE, lwd=3)
dev.off()



################################################################################
# E. Circulation
# FROM 
# https://www.paleo.bristol.ac.uk/ummodel/scripts/papers/Valdes_et_al_2021.html
# https://www.paleo.bristol.ac.uk/ummodel/users/Valdes_et_al_2021/new2/
library(ncdf4)
nc <- nc_open("data/500_texPwo/teXPwo.pfclann.nc")
names(nc$var)

u <- rotate(rast("data/500_texPwo/teXPwo.pfclann.nc", subds="ucurrTot_mm_dpth"))
v <- rotate(rast("data/500_texPwo/teXPwo.pfclann.nc", subds="vcurrTot_mm_dpth"))

# the strength of the currents
magnitude <- sqrt(u^2 + v^2)


reMag <- resample(magnitude, rast(res=1))
plot(reMag)

magnitudeProj <- project(reMag, proj)
plot(magnitudeProj)
plot(meProj, add=TRUE)


# project and resmaple
uProj <- resample(u,rast(res=1))
vProj <- resample(v,rast(res=1))


xy <- xyFromCell(uProj, 1:ncell(uProj))
df <- cbind(xy, u=values(uProj), v=values(vProj))
colnames(df) <- c("lon", "lat", "u", "v")

# The GG stands for Garbage, Garbage...
# Ekkora kókányolást még életemben nem csináltam plottal...

# Step 1. Draw long-lat raster with streamlines, export to png
png("export/dependency/alter/current_base.png", width=5000, height=2500, pointsize=60, bg=NA)

par(mar=rep(0,4))
(g <- ggplot(df, aes(lon, lat)) + #theme_bw()+
	theme(panel.background = element_rect(fill='transparent'), panel.border = element_blank(), panel.grid.major = element_blank(),plot.background = element_rect(fill = "gray"),
		panel.grid.minor = element_blank(), axis.line = element_blank(), axis.text.x=element_blank(),
		axis.text.y=element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
	geom_spatraster(data = magnitude , show.legend=FALSE) +
	 scale_fill_grass_c(palette = "oranges")+
    geom_streamline(aes(dx = dlon(u, lat), dy = dlat(v)), L = 40 , skip=12, arrow.angle=20, linewidth = 2)  )
dev.off()


# Step 2. Read in the long-lat png as RGB channels
library(png)
block <- readPNG("export/dependency/alter/current_base.png")
rgb <- c(rast(block[,,1]), rast(block[,,2]), rast(block[,,3]))
# check!
plotRGB(rgb[[1:3]]*255)

# Step 3. Crop the margins off the plot...

# the valid plo without the crap margins
theplot <- ext(c(xmin=238, xmax=4762, ymin=151, ymax=2349))
cro <- crop(rgb, theplot)
ext(cro) <- ext(rast())
crs(cro) <- "WGS84"

croProj <- project(cro, proj)

# Step 4: profit!
png("export/dependency/alter/current.png", width=3000, height=1500, pointsize=24, bg=NA)
	plotRGB(croProj*255)
#	plot(landProj, axes=FALSE,legend=FALSE, col="gray", add=TRUE)
	sphereshade(left="#000000", right="#000000", crs=proj, left.alpha=0.15, right.alpha=0)
	plot(mapedge(crs=proj), add=TRUE, lwd=3)
dev.off()


# ...or:
# Simplest solution: crop with imagemagick, project with gplates...
# convert -crop 4524x2198+238+151 current_base.png current_longlat.png
# the projection this way is somewhat nicer


################################################################################
# F. Sea level

# create a color palette
seacols <- mako(10)[4:10]
seapal <- colorRampPalette(seacols)(30)
coldf <-data.frame(
	color=c(seapal[c(1,1)], seapal, "#ffffffff", "#ffffffff"),
	z=c(-10000, -250, seq(-200, +200, length.out=length(seapal)), 250, 10000)
)
seaColors <- expand(coldf, n=1000)



sl <- dem
## values(sl)[values(dem) > 200] <- NA
## values(sl)[values(dem) < -200] <- NA
plot(sl, col=mako(6))

sl <- resample(sl, rast(res=0.1))
slproj <- project(sl, proj)
png("export/dependency/sealevel.png", width=3000, height=1500, pointsize=24, bg=NA)
	plot(slproj,  col=seaColors$col,breaks=seaColors$breaks, axes=FALSE, legend=FALSE)
	sphereshade(left="#000000", right="#000000", crs=proj, left.alpha=0.1, right.alpha=0)
	plot(mapedge(crs=proj), add=TRUE, lwd=3)
dev.off()



writeRaster(slproj, file="export/dem500moll.tif")

################################################################################
pbdb <- fetch("pbdb", ser="occs4", ver="20260412", datadir="data/chronosphere")
# approximate!

library(divDyn)

# Based on the mas this is the Miaolingian approximately
miao <- pbdb[pbdb$max_ma>=497 & pbdb$min_ma>497 & pbdb$max_ma<=506.5 & pbdb$min_ma<506.5,  ]
data(keys)
miao$lith<-divDyn::categorize(miao$lithology1,keys$lith)

cambColl <- unique(miao[, c("collection_no", "lng","lat", "lith")])
pcoords <- reconstruct(cambColl[, c("lng", "lat")], age=age, model=model)
colnames(pcoords) <- c("plong", "plat")
cambColl <- cbind(cambColl, pcoords )





png("export/dependency/paleoecology.png", width=3000, height=1500, pointsize=24, bg=NA)
	par(mai=rep(0.1,4))
	PlotOccs(x=cambColl, map=polyer, crs="ESRI:4326", plng="plong", plat="plat", 
		points.cex=4, col="white", symbol=NULL,
		symbol.col=NULL, map.bgdamp=FALSE, coloredsub=TRUE)
dev.off()

################################################################################
# Bioregionalization
 
# a bunch of random colors
load("data/allHex.RData")

# plotting the script is here
source("code/methods/plots.R")

miaoCoords <- merge(miao, cambColl[, c("collection_no", "plong", "plat")], by="collection_no")

library(icosa)
# very crude example of regionalization
hex <- hexagrid(deg=5, sf=TRUE)
miaoCoords$cell<- locate(hex,miaoCoords[, c("plong", "plat")] )

# contingency
cont <- table(miaoCoords$cell, miaoCoords$genus)

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

png("export/dependency/bioregionalization.png", width=3000, height=1500, pointsize=24, bg=NA)
	par(mai=rep(0.1,4))
	plot(meProj, col="white")
	biogeoplot(mem=mem, bg=poly$geometry, gri=hex, ylim=c(-90, 90), add=TRUE,
		alpha="66", crs=proj, gri.border="#ffffff33", labels=FALSE, gri.lwd=5)
	sphereshade(left="#000000", right="#000000", crs=proj, left.alpha=0.15, right.alpha=0)
dev.off()

################################################################################
# Evaporites
boucot <- read.csv("data/AJB Paleozoic v9.csv")

# thej
lower <- boucot[boucot$Period=="Cambrian" & boucot$LMU=="Lower", ]
lower <- cbind(lower,reconstruct(lower[, c("LONG", "LAT")], age=age, model=mod))

lower <- lower[!is.na(lower$paleolong) & !is.na(lower$paleolat), ]
lowerSF <- st_as_sf(lower, coords=c("paleolong", "paleolat"), crs="WGS84")
lowerProj <- st_transform(lowerSF, proj)


# Preciptation - evaporation
evap <- rotate(rast("data/500_texPwo/teXPw_precipevap_ann_fsy.nc"))
crs(evap) <- "WGS84"
reEvap <- resample(evap, rast(res=0.5))
reEvapProj <- project(reEvap, proj)

# make a land mask for this
landMask <- magnitude
values(landMask) <- NA
values(landMask)[is.na(values(magnitude))] <- TRUE
plot(landMask)
# get rid of the polar artifact
corr <- mask(landMask, st_sf(mapedge(ymin=70)))
values(landMask)[!is.na(values(corr))] <- NA

reLandMask <- resample(landMask, rast(res=0.5))
reLandMaskProj <- project(reLandMask, proj)
plot(reLandMaskProj)

toPlot <- c("Evaporites", "Gypsum", "Halite", "Anhydrites")

png("export/dependency/evaporites.png", width=3000, height=1500, pointsize=24, bg=NA)
	par(mai=rep(0.1,4))
	plot(reEvapProj, axes=FALSE, legend=FALSE, col=inferno(256))
	plot(polyProj$geometry, border="#ffffff44",lwd=3,  add=TRUE)
	plot(reLandMaskProj, axes=FALSE, add=TRUE, col="#ffffff66")
	for(i in 1:length(toPlot)){	
		plotThis<- lowerProj[lowerProj$Lithology==toPlot[i], ]
		plot(plotThis$geometry, pch=(21:24)[i], col="black", bg=turbo(4)[i], add=TRUE, lwd=2, cex=5)
	}
	plot(meProj, col=NA, lwd=2, border="gray90", add=TRUE)
	sphereshade(left="#000000", right="#000000", crs=proj, left.alpha=0.15, right.alpha=0)
dev.off()

################################################################################
# Miaolingian clastics
library(icosa)
library(rampage)
library(terra)



# estimate the density based on these
carb <- cambColl[cambColl$lith=="carbonate", c("plong", "plat")]
colnames(carb) <- c("long", "lat")
sil <- cambColl[cambColl$lith=="siliciclastic", c("plong", "plat")]
colnames(sil) <- c("long", "lat")

# a grid for the density estimation
gr <- hexagrid(deg=15)
CellCount <- function(x) table(x$cell)
carbO <- grapply(x=carb, out=rast(),y=gr,  iter=500, FUN=CellCount, miss=0.05)
silO <- grapply(x=sil, out=rast(),y=gr,  iter=500, FUN=CellCount, miss=0.05)


# to be used as a mask in case
polyProper <- poly[st_geometry_type(poly$geometry)=="MULTIPOLYGON",]


# create color ramps
redDF <-data.frame(
	color=colorRampPalette(c("#ffffff", gradinv(7)[6]))(7),
	z=c(-8,-2, -0.5, 0, +0.5,+2, 8)
)
reder<- expand(redDF, 256)

blueDF <-data.frame(
	color=colorRampPalette(c("#ffffff", gradinv(7)[2]))(7),
	z=c(-8,-2, -0.5, 0, +0.5,+2, 8)
)
bluer<- expand(blueDF, 256)


# mask the two rasters with 
polyCarbo<- mask(carbO, polyProper)
plot(log(polyCarbo), col=reder$col, breaks=reder$breaks, legend=FALSE)

polySilo<- mask(silO, polyProper)
plot(log(polySilo), col=bluer$col, breaks=bluer$breaks, legend=FALSE)


# And plot them!

png("export/dependency/carbonates.png", width=3000, height=1500, pointsize=24, bg=NA)
PlotLithology(
	x=cambColl[cambColl$lith=="carbonate", ],
	ras=polyCarbo, log=TRUE, proj=proj,
	plng="plong", plat="plat",
	points.cex=5, col=gradinv(7)[7], coloramp=reder )
dev.off()




png("export/dependency/siliciclasitics.png", width=3000, height=1500, pointsize=24, bg=NA)
PlotLithology(
	x=cambColl[cambColl$lith=="siliciclastic", ],
	ras=polySilo, log=TRUE, proj=proj,
	plng="plong", plat="plat",
	points.cex=5, col=gradinv(7)[1], coloramp=bluer , pch=23)
dev.off()





