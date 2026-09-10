# Script written to illustrate the effects of temporal mismatch
# Ádám T. Kocsis, 2026-04-15, Erlangen
library(rgplates)
library(chronosphere)
library(divDyn)
library(icosa)

# set working directory
setwd("/mnt/sky/Dropbox/WorkSpace/2025-06-10_PPE/")

# This example uses only the GWS-based implemetation
# use this for a local instance of the GWS
# rgplates::setgws("http://localhost:18000/")

# ge the coastlines
tcCoast<- rgplates::reconstruct("coastlines", age=400, anchor=1 , model="TorsvikCocks2017")
pmCoast<- rgplates::reconstruct("coastlines",  age=400, model="PALEOMAP")

# The Paleobiology Database
pbdb<- chronosphere::fetch("pbdb", ser="occs4", ver="20260412", datadir="data/chronosphere")

# select trilobites
dat <- pbdb[which(pbdb$class=="Trilobita"), ]

# get the devonian bunch
# Stratigraphy
data(keys, package="divDyn")

# B. the stg entries (lookup) - as in kocsis2019divDyn
stgMin<-divDyn::categorize(dat[,"early_interval"],keys$stgInt)
stgMax<-divDyn::categorize(dat[,"late_interval"],keys$stgInt)

stgMin<-as.numeric(stgMin)
stgMax<-as.numeric(stgMax)

# empty container
dat$stg <- rep(NA, nrow(dat))
# select entries, where
stgCondition <- c(
# the early and late interval fields indicate the same stg
	which(stgMax==stgMin),
# or the late_intervarl field is empty
	which(stgMax==-1))
dat$stg[stgCondition] <- stgMin[stgCondition]

# the lower Devonian
lowerDevonian <- dat[which(dat$stg%in%c(29,30,31)), ]

# reconstruct with PALEOMAP
pmCoords <- rgplates::reconstruct(lowerDevonian[, c("lng", "lat")], age=400, model="PALEOMAP")
colnames(pmCoords) <- c("pm_lng", "pm_lat")

#  reconstructTorsvik and Cocks 2017
tcCoords <- rgplates::reconstruct(lowerDevonian[, c("lng", "lat")], age=400, model="TorsvikCocks2017", anchor=1)
colnames(tcCoords) <- c("tc_lng", "tc_lat")
lowerDevonian  <- cbind(lowerDevonian,pmCoords, tcCoords)

# Assign to grid cells
hex <- icosa::hexagrid(deg=10)
hex <- icosa::newsf(hex, res=300) # increase resolution of sf object

# locate coordinaes
lowerDevonian$pm_cell<- icosa::locate(hex, lowerDevonian[, c("pm_lng", "pm_lat")])
lowerDevonian$tc_cell<- icosa::locate(hex, lowerDevonian[, c("tc_lng", "tc_lat")])


# Find Australia - based on coordinaes grabbed from the maps
# Paleomap
pmAustrliaPoints <- sf::st_as_sf(as.data.frame(cbind(long=160, lat=-14.5)), coords=c("long", "lat"), crs="WGS84")
pmCoast2 <- pmCoast
pmCoast2$index <- 1:nrow(pmCoast2)
pmJoin <- sf::st_join(pmAustrliaPoints, sf::st_make_valid(pmCoast2), join=sf::st_intersects)

# TorsvikCocks
tcAustrliaPoints <- sf::st_as_sf(as.data.frame(cbind(long=132, lat=-18.5)), coords=c("long", "lat"), crs="WGS84")
tcCoast2 <- tcCoast
tcCoast2$index <- 1:nrow(tcCoast2)
sf::sf_use_s2(FALSE)
tcJoin <- sf::st_join(tcAustrliaPoints, sf::st_make_valid(tcCoast2), join=sf::st_intersects)
sf::sf_use_s2(TRUE)

# set the projection
proj <- "ESRI:54009"

################################################################################
# Actual example: Trimerus with use mollweide!

#' @param x occurrence data frame
#' @param gen Genus name character string.
#' @param model Model abbreviation (tc, pm)
#' @param coast Sf for the modern coastline reconstruction
#' @param proj EPSG id of projection.
#' @param gr Icosa grid object.
#' @param highlight.index Which polygon is to be highlighted?
#' @param highlight.col with what color filling?
#' @param highlight.border with what color border?
#' @param dir where should be the plots put?
RangePlot <- function(x, gen="", model, coast, proj="ESRI:54009", gr=hex, highlight.index=NULL,
	highlight.col="#ffffd6ff", highlight.border="#c0b100ff", dir="export/mismatch"){

	# create directory
	dir.create("export", showWarnings=FALSE)
	dir.create(dir, showWarnings=FALSE)

	# get the data
	genDat <-x[grepl(gen, x$genus), ]

	# plot a reconstruction
	base <- sf::st_transform(coast, proj)

	coordColumns <- paste0(model, "_", c("lng", "lat"))
	# omit missing values
	genDat <- genDat[!is.na(genDat[, coordColumns[1]]), ]
	# cooredinates
	occs <- sf::st_as_sf(genDat[,
		c("collection_no", coordColumns)], coords=coordColumns, crs="WGS84")	
	occs <- sf::st_transform(occs, proj)

	# maximum great circle distances!
	png(paste0(dir, "/",gen,"_",model, ".png" ),
		height=800, width=1500, bg=NA)
		par(mai=rep(0.2, 4))
		plot(rgplates::mapedge(crs=proj), col="white", reset=FALSE,border="gray30", lwd=3)
		plot(base, col="gray", border=NA, add=TRUE)
		plot(gr, border="gray90", col=NA, add=TRUE, crs=proj)
		if(!is.null(highlight.index))
			plot(base[highlight.index,], col=highlight.col, border=highlight.border, add=TRUE, lwd=2)

		plot(hex, unique(genDat[, paste0(model, "_cell")]), border="gray90",
			col="#0088AA77", add=TRUE, crs=proj)
		plot(occs, col="red", pch=3, lwd=3, add=TRUE)
		plot(rgplates::mapedge(crs=proj), col=NA, add=TRUE,border="gray30", lwd=3)

	dev.off()

}

# make the plots
RangePlot(x=lowerDevonian, proj=proj, gen="Trimerus", model="pm",
	coast=pmCoast, highlight.index=pmJoin$index, highlight.border="#c0b100ff")
RangePlot(x=lowerDevonian, proj=proj, gen="Trimerus", model="tc",
	coast=tcCoast, highlight.index=tcJoin$index, highlight.border="#c0b100ff")
