library(rgplates) # requires at least 0.6.1!
library(chronosphere)

setwd("/mnt/sky/Dropbox/WorkSpace/2025-06-10_PPE/")

# get the TorsvikCocks model (rgplates::platemodel class object)
source("code/methods/TorsvikCocks2025.R")
source("code/methods/TorsvikCocks2017.R") #outputs TC2017

# get the model data from the chronosphere
MERDITH2021 <- chronosphere::fetch("GPlates", "MERDITH2021", datadir="data/chronosphere")
CAO2024 <- chronosphere::fetch("GPlates", "CAO2024", datadir="data/chronosphere")
PALEOMAP <- chronosphere::fetch("paleomap", "model", datadir="data/chronosphere")
MULLER2022 <- chronosphere::fetch("GPlates", "MULLER2022", datadir="data/chronosphere") # uses 1000_0_rotfile_MantleOpt.rot!

# The chosen projection
proj <- "ESRI:54030"

# the edge of the map
meProj <- rgplates::mapedge(crs=proj)

# use this for a local instance of the GWS
# setgws("http://localhost:18000/", check=FALSE)

# GWS - based
platesTC2017mantle <- rgplates::reconstruct("static_polygons", age=500, model=TC2017, anchor=0)
platesTC2017mag <- rgplates::reconstruct("static_polygons", age=500, model=TC2017, anchor=1)
platesTC2025 <- rgplates::reconstruct("coastlines", age=500, model=TC2025)
microTC2025 <- rgplates::reconstruct("microcontinents", age=500, model=TC2025)
#platesMERDITH2021 <- rgplates::reconstruct("coastlines", age=500, model=MERDITH2021)
platesMERDITH2021 <- rgplates::reconstruct("coastlines", age=500, model="MERDITH2021")
platesMULLER2022 <- rgplates::reconstruct("coastlines", age=500, model=MULLER2022)

# Gives the same result as the web???
## platesMULLER2022gws <- rgplates::reconstruct("coastlines", age=500, model="MULLER2022")
## plot(platesMULLER2022gws$geometry, col="#FF000088", border=NA)
## plot(platesMULLER2022$geometry, col="#0000FF88", add=TRUE, border=NA)
# project
platesTC2017mantleProj <- st_transform(platesTC2017mantle, crs=proj)
platesTC2017magProj <- st_transform(platesTC2017mag, crs=proj)
platesTC2025Proj <- st_transform(platesTC2025, crs=proj)
microTC2025Proj <- st_transform(microTC2025, crs=proj)
platesMERDITH2021Proj <- st_transform(platesMERDITH2021, crs=proj) # requires a bit of interpolation!
platesMULLER2022Proj <- st_transform(platesMULLER2022, crs=proj) # requires a bit of interpolation!

# trond model versions
plot(meProj)
plot(platesTC2017mantleProj$geometry, col="#44444466", add=TRUE, border=NA)
plot(platesTC2017magProj$geometry, col="#99999966", add=TRUE, border=NA)
plot(platesTC2025Proj$geometry, col="#FF000066", add=TRUE, border=NA)
plot(microTC2025Proj$geometry, col="#FF000066", add=TRUE, border=NA)


sphereshade <- function(n=180,res=5, left="#ffffff", right="#aaaaaa", crs="EPSG:4326"){
	breaks <- seq(-180, 180, length.out=n+1)

	colFun <- colorRampPalette(c(left, right))
	cols <- colFun(n)
	for(i in 1:n){
		oneSlice <- mapedge(xmin=breaks[i],xmax=breaks[i+1],crs=crs)
		plot(oneSlice, border=NA, col=cols[i], add=TRUE)
	}

}


################################################################################
colMerdith2021 <-"#880000"
colMuller2022 <-"#000088"
colTCmag <- "#00FF00"
colTCmantle <- "#008844"
alpha <- "66"


# comparison
dir.create("export", showWarnings=FALSE)
dir.create("export/compare/", showWarnings=FALSE)

# 0. Two trond models
png("export/compare/TC.png", width=3000, height=1500, pointsize=24, bg="transparent")
par(mai=rep(0.1,4))
plot(meProj)
sphereshade(crs=proj, right="#ffffff", left="#eeeeee")
plot(platesTC2017mantleProj$geometry, col=paste0(colTCmantle, alpha), add=TRUE, border=NA)
plot(platesTC2017magProj$geometry, col=paste0(colTCmag, alpha), add=TRUE, border=NA)
legend("topright", fill=paste0(c(colTCmag, colTCmantle),alpha), legend=c("TorsvikCocks2017 (PMRF)", "TorsvikCocks2017 (MRF)"), bty="n", cex=2, inset=c(0.18,0.1))
plot(meProj, add=TRUE, lwd=7)
dev.off()


# 1. only palemagnetic reference
png("export/compare/mag.png", width=3000, height=1500, pointsize=24, bg="transparent")
par(mai=rep(0.1,4))
plot(meProj)
sphereshade(crs=proj, right="#ffffff", left="#eeeeee")
plot(platesTC2017magProj$geometry, col=paste0(colTCmag,alpha), add=TRUE, border=NA)
plot(platesMERDITH2021Proj$geometry, col=paste0(colMerdith2021,alpha), add=TRUE, border=NA)
plot(meProj, add=TRUE, border="darkred", lwd=10)
legend("topright", fill=paste0(c(colTCmag, colMerdith2021),alpha), legend=c("TorsvikCocks2017 (PMRF)", "MERDITH2021"), bty="n", cex=2, inset=c(0.18,0.1))
dev.off()


# only mantle reference
png("export/compare/mantle.png", width=3000, height=1500, pointsize=24, bg="transparent")
par(mai=rep(0.1,4))
plot(meProj)
sphereshade(crs=proj, right="#ffffff", left="#eeeeee")
plot(platesTC2017mantleProj$geometry, col=paste0(colTCmantle,alpha), add=TRUE, border=NA)
plot(platesMULLER2022Proj$geometry, col=paste0(colMuller2022,alpha), add=TRUE, border=NA)
legend("topright", fill=paste0(c(colTCmantle, colMuller2022),alpha), legend=c("TorsvikCocks2017 (MRF)", "MULLER2022"), bty="n", cex=2, inset=c(0.18,0.1))
plot(meProj, add=TRUE, lwd=7)
dev.off()
