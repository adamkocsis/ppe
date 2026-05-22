library(rgplates)
library(chronosphere)

setwd("/mnt/sky/Dropbox/WorkSpace/2025-06-10_PPE/")

# Get the models!
# get the TorsvikCocks model (rgplates::platemodel class object)
source("code/methods/TorsvikCocks2017.R") #outputs TC2017
PALEOMAP <- fetch("paleomap", "model", datadir="data/chronosphere")

# compare australia
tcCoast<- reconstruct("coastlines", age=400, anchor=1 , model="TorsvikCocks2017")
pmCoast<- reconstruct("coastlines",  age=400, model="PALEOMAP")


# for the modern basemap
pmModern<- reconstruct("coastlines",  age=0, model="PALEOMAP")

png("export/modern.png", width=4000, height=2000)
plot(pmModern$geometry, col="gray", border=NA)
dev.off()


dir.create("export/ranges", showWarnings=FALSE)
dir.create("export/ranges/tc400", showWarnings=FALSE)

png("export/ranges/tc400base.png", width=1000, height=1000)
plot(tcCoast$geometry, col="gray", border=NA, xlim=c(115, 155), ylim=c(-40, 10))
abline(v=seq(90, 150, 15), lty=2, col="gray80")
abline(h=seq(-45, 15, 15), lty=2, col="gray80")
dev.off()

dir.create("export/ranges/pm400", showWarnings=FALSE)

png("export/ranges/pm400base.png", width=1000, height=1000)
plot(pmCoast$geometry, col="gray", border=NA, xlim=c(125, 175), ylim=c(-40, 10))
abline(v=seq(90, 150, 15), lty=2, col="gray80")
abline(h=seq(-45, 15, 15), lty=2, col="gray80")
dev.off()

pbdb<- fetch("pbdb", datadir="data/chronosphere")
dat <- pbdb[pbdb$class=="Trilobita", ]

# get the devonian bunch
# Stratigraphy
library(divDyn)
data(keys)

# B. the stg entries (lookup)
stgMin<-categorize(dat[,"early_interval"],keys$stgInt)
stgMax<-categorize(dat[,"late_interval"],keys$stgInt)

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

# categorize entries as they are in the lookup table
tenMin <- categorize(dat[ ,"early_interval"], keys$tenInt) 
tenMax <- categorize(dat[ ,"late_interval"], keys$tenInt)

# Convert to simple numeric values

tenMin <- as.numeric(tenMin) 
tenMax <- as.numeric(tenMax)

# Define an empty vector

dat$ten <- rep(NA, nrow(dat))

# Assign 10 Myr year bin numbers to each row in myOrdBrachDataforAffinities

# Select entries where
tenCondition <- c(
	# the early and late interval fields indicate the same bin 
	which(tenMax==tenMin),
	# or the late_interval field is empty
	which(tenMax==-1))

# in these entries, use the bin indicated by the early_interval
dat$ten[tenCondition] <- tenMin[tenCondition]



emsian <- dat[which(dat$stg==31), ]
emsian <- dat[which(dat$ten==13), ]
emsian <- dat[which(dat$stg%in%c(29,30,31)), ]

# reconstruct the paleomap
pmCoords <- reconstruct(emsian[, c("lng", "lat")], age=400, model="PALEOMAP")
colnames(pmCoords) <- c("pm_lng", "pm_lat")
tcCoords <- reconstruct(emsian[, c("lng", "lat")], age=400, model="TorsvikCocks2017", anchor=1)
colnames(tcCoords) <- c("tc_lng", "tc_lat")
emsian  <- cbind(emsian,pmCoords, tcCoords)

allGen <- unique(emsian$genus)

dir.create("export/ranges/pm400_stg29_30_31", showWarnings=FALSE)
for(i in 1:length(allGen)){
	gen <- allGen[i]
	genDat <-emsian[which(emsian$genus==gen), ]

	png(paste0("export/ranges/pm400_stg29_30_31/", gsub(" ", "", gen), ".png"), width=2000, height=1000, pointsize=30)
	plot(pmCoast, col="gray", border=NA, main=gen, reset=FALSE)
	points(genDat$pm_lng, genDat$pm_lat, col="red", pch=3, lwd=3)
	dev.off()


}

## Ceratocephala
## Cryphaspis
## Gravicalymene
## Leonaspis
## Maurotarion
## Phacops
## Proetus
## Radiaspis
## Scutellum

# 29-31
# Crotealocephalus
# Gravicalymene
# Proetus (sl)


library(icosa)
hex <- hexagrid(deg=10)
hex <- newsf(hex, res=300)

emsian$pm_cell<- locate(hex, emsian[, c("pm_lng", "pm_lat")])
emsian$tc_cell<- locate(hex, emsian[, c("tc_lng", "tc_lat")])


## gen <- "Ceratocephala"
## gen <- "Phacops"
## gen <- "Leonaspis"



dir.create("export/ranges/lowerdev/", showWarnings=FALSE)


tasmanGen <- c(
	"Acastella",
	"Belenopyge",
	"Calymene",
	"Ceratocephala",
	"Cheirurus",
	"Chotecops",
	"Coniproetus",
	"Cornuproetus",
	"Crotalocephalus",
	"Cyphaspis	",
	"Dalmanites",
	"Dicranurus",
	"Gravicalymene",
	"Harpes",
	"Harpidella",
	"Kainops",
	"Koneprusia",
	"Leonaspis",
	"Nephranomma",
	"Odontochile",
	"Paciphacops",
	"Phacops",
	"Proetus",
	"Prokops",
	"Radiaspis",
	"Scutellum",
	"Trimerus",
	"Zlichovaspis"
)

# except otarion and otarionella
for(i in 1:length(tasmanGen)){
	gen<-tasmanGen[i]
	genDat <-emsian[grepl(gen, emsian$genus), ]

	png(paste0("export/ranges/lowerdev/", gsub(" ", "", gen), "_pm.png"), width=2000, height=1000, pointsize=30)
	par(mai=rep(0.1,4))
	plot(pmCoast, col="gray", border=NA, reset=FALSE)
	plot(hex, border="gray90", col=NA, add=TRUE)
	plot(hex, unique(genDat$pm_cell), border="gray90", col="#0088AA77", add=TRUE)
	points(genDat$pm_lng, genDat$pm_lat, col="red", pch=3, lwd=3)
	dev.off()

	png(paste0("export/ranges/lowerdev/", gsub(" ", "", gen), "_tc.png"), width=2000, height=1000, pointsize=30)
	par(mai=rep(0.1,4))
	plot(tcCoast, col="gray", border=NA, reset=FALSE)
	plot(hex, border="gray90", col=NA, add=TRUE)
	plot(hex, unique(genDat$tc_cell), border="gray90", col="#0088AA77", add=TRUE)
	points(genDat$tc_lng, genDat$tc_lat, col="red", pch=3, lwd=3)
	dev.off()
}


################################################################################
# Actual example: Trimerus - use mollweide!

## gen <- "Trimerus"
## x <- emsian
## model <- "pm"
## coast <- pmCoast


NicePlot <- function(x, gen="", model, coast, proj="ESRI:54009", gr=hex, highlight.index=NULL, highlight.col="#ffffd6ff", highlight.border="#c0b100ff"){
	dir.create("export/elements/", showWarnings=FALSE)
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

	png(paste0("export/elements/",gen,"_",model, ".png" ),
		height=800, width=1500, bg=NA)
	par(mai=rep(0.2, 4))
	plot(mapedge(crs=proj), col="white", reset=FALSE,border="gray30", lwd=3)
	plot(base, col="gray", border=NA, add=TRUE)
	plot(hex, border="gray90", col=NA, add=TRUE, crs=proj)
	if(!is.null(highlight.index))
		plot(base[highlight.index,], col=highlight.col, border=highlight.border, add=TRUE, lwd=2)

	plot(hex, unique(genDat[, paste0(model, "_cell")]), border="gray90", col="#0088AA77", add=TRUE, crs=proj)
	plot(occs, col="red", pch=3, lwd=3, add=TRUE)
	plot(mapedge(crs=proj), col=NA, add=TRUE,border="gray30", lwd=3)
	
	dev.off()


}

# Find Australia
# Paleomap
pmAustrliaPoints <- st_as_sf(as.data.frame(cbind(long=160, lat=-14.5)), coords=c("long", "lat"), crs="WGS84")
pmCoast2 <- pmCoast
pmCoast2$index <- 1:nrow(pmCoast2)
pmJoin <- sf::st_join(pmAustrliaPoints, st_make_valid(pmCoast2), join=st_intersects)
# TorsvikCocks
tcAustrliaPoints <- st_as_sf(as.data.frame(cbind(long=132, lat=-18.5)), coords=c("long", "lat"), crs="WGS84")
tcCoast2 <- tcCoast
tcCoast2$index <- 1:nrow(tcCoast2)
sf_use_s2(FALSE)
tcJoin <- sf::st_join(tcAustrliaPoints, st_make_valid(tcCoast2), join=st_intersects)
sf_use_s2(TRUE)

# the projections
proj <- "ESRI:54030"
proj <- "ESRI:54009"

NicePlot(x=emsian, proj=proj, gen="Trimerus", model="pm", coast=pmCoast, highlight.index=pmJoin$index, highlight.border="#c0b100ff")
NicePlot(x=emsian, proj=proj, gen="Trimerus", model="tc", coast=tcCoast, highlight.index=tcJoin$index, highlight.border="#c0b100ff")
