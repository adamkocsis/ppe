# Script written to illustrate the dependency of paleogeographic products
# Ádám T. Kocsis, 2026-06-17, Erlangen
 
library(rgplates) # relies on both the local GPlates app, as well as the GWS (accessed locally)
library(chronosphere)
library(rampage)
library(smoothr) # resampling shapes
library(terra) # rasters
library(rampage) # color ramps
library(viridisLite) # directly manipulate netcdf
library(divDyn) # geologic timescale and data binning
library(icosa) # icosahedral gridding
library(vegan) # dissimilarity metric

# for the streamlines only
library(ncdf4) # directly manipulate netcdf
library(ggplot2) # grid driver
library(metR) # streamlines
library(tidyterra) # background
library(png) # reading it back

# set working directory
setwd("/mnt/sky/Dropbox/WorkSpace/2025-06-10_PPE/")

# source some plotting candies
source("code/methods/plots.R")

# The Cambrian mid age
age <- 500

# This code relies on GWS 1.0.0.
# If you want to run this locally, you can run the following shell command:
# docker run -d --rm -p 18000:80 gplates/gws:v1.0.0

# And then in R run
rgplates::setgws("http://localhost:18000/")

# after this you can turn 
model <- "PALEOMAP"

# generally used map projection (Mollweide)
proj <- "ESRI:54009"
meProj <- rgplates::mapedge(crs=proj)

# create export directories
dir.create("export", showWarnings=FALSE)
dir.create("export/dependency", showWarnings=FALSE)

################################################################################
# A. Global plate model

# ensure directory is there
dir.create("data/chronosphere", showWarnings=FALSE)

# To be used with the GPlates app. This version is the same as in the GWS.
mod <- chronosphere::fetch("paleomap", "model", datadir="data/chronosphere",  ver="v3-GPlates" )

# static polygons (not in GWS)
poly <- rgplates::reconstruct("static_polygons", age=age, model=mod)

# densify to make it look nicer in projections
polyer <- smoothr::densify(poly, n=10)

# get the modern coastlines (from GWS)
coast <- rgplates::reconstruct("coastlines", age=age, model=model)

# project both to target
coastProj <- sf::st_transform(coast, crs=proj)
polyProj <- sf::st_transform(polyer, crs=proj)


# Two different versions plotted
png("export/dependency/gpm.png", width=3000, height=1500, pointsize=24, bg=NA)
	par(mai=rep(0.1,4))
	plot(meProj, col="#1A6BB0", border="gray90")
	sphereshade(left="#0f3f67", right="#1A6BB0", crs=proj)
	plot(polyProj$geometry, border=NA, col="gray90", add=TRUE)
	plot(coastProj$geometry, border="black", col="gray90", add=TRUE, lwd=2)
	plot(meProj, col=NA, lwd=2, border="gray90", add=TRUE)
dev.off()

################################################################################
# B. Paleomap DEM

# for coloring the rasters
dems <- chronosphere::fetch("paleomap", "dem", res=0.1, data="data/chronosphere")
dem <- dems[as.character(age)]

# for better visualization - upscale to ensure good geometry
demRe<- terra::resample(dem, terra::rast(res=0.05))

# and then project
demProj <- terra::project(demRe, proj, mask=TRUE, use_gdal=FALSE)

# a smoother topographic map
data(topos, package="rampage")
topocol <- rampage::expand(topos$ptolemy, n=512)

png("export/dependency/dem.png", width=3000, height=1500, pointsize=24, bg=NA)
	par(mai=rep(0.1,4))
	plot(demProj, col=topocol$col, breaks=topocol$breaks, legend=FALSE, smooth=TRUE, axes=FALSE)
	sphereshade(left="#000000", right="#000000", crs=proj, left.alpha=0.25, right.alpha=0)
	plot(meProj, col=NA, border="white", lwd=3, add=TRUE)
dev.off()



################################################################################
# C. Weathering and degassing latitude bands - with fabricated data

# define laitudes
at <- c(30,-30)

# standard deviaion (zone width)
dev <- 8 

# define 
weather <- terra::rast(res=0.1)
cells <- terra::xyFromCell(weather, 1:terra::ncell(weather))
oneDensity <- dnorm(x=cells[,2], mean=at[1], sd=dev)
twoDensity <- dnorm(x=cells[,2], mean=at[2], sd=dev)
densities <- oneDensity + twoDensity 

# the vlaues
terra::values(weather) <- densities

# project to target 
weatherProj <- terra::project(weather, proj)

# make a color ramp
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

# The GMST of C Scotese
gmst<- chronosphere::fetch("paleomap", "gmst", ver="scotese02a_v21321", datadir="data/chronosphere")

# Get the values on land
landProj <- demProj
terra::values(landProj) <-NA
terra::values(landProj)[terra::values(demProj)>0] <- 1

# the Camrian one
one <- gmst[as.character(age)]

# enforce missing CRS
terra::crs(one) <- "WGS84"

# resample and project
oneRe <- terra::resample(one, terra::rast(res=0.5))
gmstProj<- terra::project(oneRe, proj)

# render 
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

# open to identify layers for velocities 
nc <- ncdf4::nc_open("data/500_texPwo/teXPwo.pfclann.nc")
names(nc$var)

# load layers and rotate
u <- terra::rotate(terra::rast("data/500_texPwo/teXPwo.pfclann.nc", subds="ucurrTot_mm_dpth"))
v <- terra::rotate(terra::rast("data/500_texPwo/teXPwo.pfclann.nc", subds="vcurrTot_mm_dpth"))

# the strength of the currents
magnitude <- sqrt(u^2 + v^2)

# resample magnitude
reMag <- resample(magnitude, terra::rast(res=1))

# and project
magnitudeProj <- terra::project(reMag, proj)

# preliminary plotting
## plot(magnitudeProj)
## plot(meProj, add=TRUE)

# project and resmaple velocities separately
uProj <- terra::resample(u,terra::rast(res=1))
vProj <- terra::resample(v,terra::rast(res=1))


# change to a data.frame so it can work with the monster below
xy <- terra::xyFromCell(uProj, 1:terra::ncell(uProj))
df <- cbind(xy, u=values(uProj), v=values(vProj))
colnames(df) <- c("lon", "lat", "u", "v")

# Ekkora kókányolást még életemben nem csináltam plottal... mekkora egy gagyi sz*r ez?
# Step 1. Draw long-lat raster with streamlines, export to png
dir.create("export/dependency/alter", showWarnings=FALSE)
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
block <- png::readPNG("export/dependency/alter/current_base.png")
rgb <- c(terra::rast(block[,,1]), terra::rast(block[,,2]), terra::rast(block[,,3]))
# check!
terra::plotRGB(rgb[[1:3]]*255)

# Step 3. Crop the margins off the plot... wow
# the valid plo without the crap margins
theplot <- terra::ext(c(xmin=238, xmax=4762, ymin=151, ymax=2349))
cro <- terra::crop(rgb, theplot)
ext(cro) <- terra::ext(terra::rast())
crs(cro) <- "WGS84"

croProj <- terra::project(cro, proj)

# Step 4: profit!
png("export/dependency/current_gplates.png", width=3000, height=1500, pointsize=24, bg=NA)
	terra::plotRGB(croProj*255)
#	plot(landProj, axes=FALSE,legend=FALSE, col="gray", add=TRUE)
	sphereshade(left="#000000", right="#000000", crs=proj, left.alpha=0.15, right.alpha=0)
	plot(rgplates::mapedge(crs=proj), add=TRUE, lwd=3)
dev.off()


# ...or:
# Simplest solution: crop with imagemagick, project with gplates...
# convert -crop 4524x2198+238+151 current_base.png current_longlat.png
# the projection this way is somewhat nicer... this was used in the paper.


################################################################################
# F. Sea level
# create a color palette
seacols <- viridisLite::mako(10)[4:10]
seapal <- colorRampPalette(seacols)(30)
coldf <-data.frame(
	color=c(seapal[c(1,1)], seapal, "#ffffffff", "#ffffffff"),
	z=c(-10000, -250, seq(-200, +200, length.out=length(seapal)), 250, 10000)
)
seaColors <- rampage::expand(coldf, n=1000)



# copy over the dem raster
sl <- dem

# make it nicer
sl <- terra::resample(sl, terra::rast(res=0.1))
slproj <- project(sl, proj)
png("export/dependency/sealevel.png", width=3000, height=1500, pointsize=24, bg=NA)
	plot(slproj,  col=seaColors$col,breaks=seaColors$breaks, axes=FALSE, legend=FALSE)
	sphereshade(left="#000000", right="#000000", crs=proj, left.alpha=0.1, right.alpha=0)
	plot(rgplates::mapedge(crs=proj), add=TRUE, lwd=3)
dev.off()



## writeRaster(slproj, file="export/dem500moll.tif")

################################################################################
# get the PBDB to illustrate Paleoecology,  etc.
pbdb <- chronosphere::fetch("pbdb", ser="occs4", ver="20260412", datadir="data/chronosphere")
# approximate!

# Based on the this is the Miaolingian approximately (corresponding to 500Ma)
miao <- pbdb[pbdb$max_ma>=497 & pbdb$min_ma>497 & pbdb$max_ma<=506.5 & pbdb$min_ma<506.5,  ]

# get lithology info 
data(keys, package="divDyn")
miao$lith<-divDyn::categorize(miao$lithology1,keys$lith)

# reconstruct paleocoords
cambColl <- unique(miao[, c("collection_no", "lng","lat", "lith")])
pcoords <- rgplates::reconstruct(cambColl[, c("lng", "lat")], age=age, model=model)
colnames(pcoords) <- c("plong", "plat")
cambColl <- cbind(cambColl, pcoords )

# illustrate paleoecology
png("export/dependency/paleoecology.png", width=3000, height=1500, pointsize=24, bg=NA)
	par(mai=rep(0.1,4))
	PlotOccs(x=cambColl, map=polyer, crs=proj, plng="plong", plat="plat", 
		points.cex=4, col="white", symbol=NULL,
		symbol.col=NULL, map.bgdamp=FALSE, coloredsub=TRUE)
dev.off()

################################################################################
# Bioregionalization
 
# a bunch of random colors
load("data/allHex.RData")

# join collections together with occurrences
miaoCoords <- merge(miao, cambColl[, c("collection_no", "plong", "plat")], by="collection_no")

# very crude example of regionalization
hex <- icosa::hexagrid(deg=5, sf=TRUE)
miaoCoords$cell<- icosa::locate(hex,miaoCoords[, c("plong", "plat")] )

# contingency
cont <- table(miaoCoords$cell, miaoCoords$genus)

# incidence
cont[cont>1] <- 1

# Method 1. Compositional similarity
distmat <- vegan::vegdist(cont, method="jaccard")

# clustering
cluster <- hclust(distmat, "ward.D2")

# plot this
h <- 1.1
## plot(cluster)
## abline(h=h, col="red")

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

# Select the appropriate records 
lower <- boucot[boucot$Period=="Cambrian" & boucot$LMU=="Lower", ]

# and reconstruct paleogeography
lower <- cbind(lower,rgplates::reconstruct(lower[, c("LONG", "LAT")], age=age, model=mod))

lower <- lower[!is.na(lower$paleolong) & !is.na(lower$paleolat), ]
lowerSF <- sf::st_as_sf(lower, coords=c("paleolong", "paleolat"), crs="WGS84")
lowerProj <- sf::st_transform(lowerSF, proj)

# Preciptation - evaporation
# same source as the circulation
evap <- terra::rotate(terra::rast("data/500_texPwo/teXPw_precipevap_ann_fsy.nc"))
terra::crs(evap) <- "WGS84"

# resample and project
reEvap <- terra::resample(evap, terra::rast(res=0.5))
reEvapProj <- terra::project(reEvap, proj)

# make a land mask for this
landMask <- magnitude
terra::values(landMask) <- NA
terra::values(landMask)[is.na(terra::values(magnitude))] <- TRUE

# example check
plot(landMask)

# get rid of the polar artifact
corr <- terra::mask(landMask, sf::st_sf(mapedge(ymin=70)))
terra::values(landMask)[!is.na(terra::values(corr))] <- NA

# resample and project...
reLandMask <- terra::resample(landMask, terra::rast(res=0.5))
reLandMaskProj <- terra::project(reLandMask, proj)
plot(reLandMaskProj)

# wha facies to plot?
toPlot <- c("Evaporites", "Gypsum", "Halite", "Anhydrites")

png("export/dependency/evaporites.png", width=3000, height=1500, pointsize=24, bg=NA)
	par(mai=rep(0.1,4))
	plot(reEvapProj, axes=FALSE, legend=FALSE, col=viridisLite::inferno(256))
	plot(polyProj$geometry, border="#ffffff44",lwd=3,  add=TRUE)
	plot(reLandMaskProj, axes=FALSE, add=TRUE, col="#ffffff66")
	for(i in 1:length(toPlot)){	
		plotThis<- lowerProj[lowerProj$Lithology==toPlot[i], ]
		plot(plotThis$geometry, pch=(21:24)[i], col="black", bg=viridisLite::turbo(4)[i], add=TRUE, lwd=2, cex=5)
	}
	plot(meProj, col=NA, lwd=2, border="gray90", add=TRUE)
	sphereshade(left="#000000", right="#000000", crs=proj, left.alpha=0.15, right.alpha=0)
dev.off()

################################################################################
# Miaolingian clastics

# estimate the density based on these
# carbonates
carb <- cambColl[cambColl$lith=="carbonate", c("plong", "plat")]
colnames(carb) <- c("long", "lat")
# siliciclastics
sil <- cambColl[cambColl$lith=="siliciclastic", c("plong", "plat")]
colnames(sil) <- c("long", "lat")

# a grid for the smoothing/density estimation 
gr <- icosa::hexagrid(deg=15)

# function to apply in every cell
CellCount <- function(x) table(x$cell)
carbO <- icosa::grapply(x=carb, out=terra::rast(),y=gr,  iter=500, FUN=CellCount, miss=0.05)
silO <- icosa::grapply(x=sil, out=terra::rast(),y=gr,  iter=500, FUN=CellCount, miss=0.05)

# to be used as a mask in case for land
polyProper <- poly[sf::st_geometry_type(poly$geometry)=="MULTIPOLYGON",]

# two different color ramps
# create color ramps
redDF <-data.frame(
	color=colorRampPalette(c("#ffffff", rampage::gradinv(7)[6]))(7),
	z=c(-8,-2, -0.5, 0, +0.5,+2, 8)
)
reder<- rampage::expand(redDF, 256)

blueDF <-data.frame(
	color=colorRampPalette(c("#ffffff", rampage::gradinv(7)[2]))(7),
	z=c(-8,-2, -0.5, 0, +0.5,+2, 8)
)
bluer<- rampage::expand(blueDF, 256)

# mask the two rasters with 
polyCarbo<- terra::mask(carbO, polyProper)
plot(log(polyCarbo), col=reder$col, breaks=reder$breaks, legend=FALSE)

polySilo<- terra::mask(silO, polyProper)
plot(log(polySilo), col=bluer$col, breaks=bluer$breaks, legend=FALSE)

# And plot them!
png("export/dependency/carbonates.png", width=3000, height=1500, pointsize=24, bg=NA)
PlotLithology(
	x=cambColl[cambColl$lith=="carbonate", ],
	ras=polyCarbo, log=TRUE, proj=proj,
	plng="plong", plat="plat",
	points.cex=5, col=rampage::gradinv(7)[7], coloramp=reder )
dev.off()


png("export/dependency/siliciclasitics.png", width=3000, height=1500, pointsize=24, bg=NA)
PlotLithology(
	x=cambColl[cambColl$lith=="siliciclastic", ],
	ras=polySilo, log=TRUE, proj=proj,
	plng="plong", plat="plat",
	points.cex=5, col=rampage::gradinv(7)[1], coloramp=bluer , pch=23)
dev.off()





