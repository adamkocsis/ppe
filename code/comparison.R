# Script written to illustrate the use of GPMs in Paleobiology
# Ádám T. Kocsis, 2026-04-15, Erlangen
library(rgplates) 
library(chronosphere)
library(smoothr)

# working direcory
setwd("/mnt/sky/Dropbox/WorkSpace/2025-06-10_PPE/")

# source some plotting candies
source("code/methods/plots.R")

# The chosen projection
proj <- "ESRI:54030"

# the edge of the map
meProj <- rgplates::mapedge(crs=proj)

# This uses only the GWS-based implemetation
# use this for a local instance of the GWS
setgws("http://localhost:18000/")

# GWS - based
TC2017mantle <- rgplates::reconstruct("coastlines", age=500, model="TorsvikCocks2017", anchor=0)
TC2017mag <- rgplates::reconstruct("coastlines", age=500, model="TorsvikCocks2017", anchor=1)
MERDITH2021 <- rgplates::reconstruct("coastlines", age=500, model="MERDITH2021")
MULLER2022 <- rgplates::reconstruct("coastlines", age=500, model="MULLER2022")

# increase point count for smoother projectsion
TC2017mantle <- smoothr::densify(TC2017mantle, n=10)
TC2017mag <- smoothr::densify(TC2017mag, n=10)
MERDITH2021 <- smoothr::densify(MERDITH2021, n=10)
MULLER2022 <- smoothr::densify(MULLER2022, n=10)

# project all
TC2017mantleProj <- sf::st_transform(TC2017mantle, crs=proj)
TC2017magProj <- sf::st_transform(TC2017mag, crs=proj)
MERDITH2021Proj <- sf::st_transform(MERDITH2021, crs=proj) 
MULLER2022Proj <- sf::st_transform(MULLER2022, crs=proj) 


# set the colors
colMerdith2021 <-"#880000"
colMuller2022 <-"#000088"
colTCmag <- "#00FF00"
colTCmantle <- "#008844"
alpha <- "66"


# comparison
dir.create("export", showWarnings=FALSE)
dir.create("export/compare/", showWarnings=FALSE)

# 1. Two TorsvikCocks models
png("export/compare/TC.png", width=3000, height=1500, pointsize=24, bg="transparent")
par(mai=rep(0.1,4), family="mono")
plot(meProj)
sphereshade(crs=proj, right="#ffffff", left="#eeeeee")
plot(TC2017mantleProj$geometry, col=paste0(colTCmantle, alpha), add=TRUE, border=NA)
plot(TC2017magProj$geometry, col=paste0(colTCmag, alpha), add=TRUE, border=NA)
legend("topright", fill=paste0(c(colTCmag, colTCmantle),alpha), legend=c("TorsvikCocks2017 (PMRF)", "TorsvikCocks2017 (MRF)"), bty="n", cex=2, inset=c(0.22,0.05))
plot(meProj, add=TRUE, lwd=7)
dev.off()


# 2. only palemagnetic reference
png("export/compare/mag.png", width=3000, height=1500, pointsize=24, bg="transparent")
par(mai=rep(0.1,4), family="mono")
plot(meProj)
sphereshade(crs=proj, right="#ffffff", left="#eeeeee")
plot(TC2017magProj$geometry, col=paste0(colTCmag,alpha), add=TRUE, border=NA)
plot(MERDITH2021Proj$geometry, col=paste0(colMerdith2021,alpha), add=TRUE, border=NA)
plot(meProj, add=TRUE, border="darkred", lwd=10)
legend("topright", fill=paste0(c(colTCmag, colMerdith2021),alpha), legend=c("TorsvikCocks2017 (PMRF)", "MERDITH2021"), bty="n", cex=2, inset=c(0.22,0.05))
dev.off()


# 3. only mantle reference
png("export/compare/mantle.png", width=3000, height=1500, pointsize=24, bg="transparent")
par(mai=rep(0.1,4), family="mono")
plot(meProj)
sphereshade(crs=proj, right="#ffffff", left="#eeeeee")
plot(TC2017mantleProj$geometry, col=paste0(colTCmantle,alpha), add=TRUE, border=NA)
plot(MULLER2022Proj$geometry, col=paste0(colMuller2022,alpha), add=TRUE, border=NA)
legend("topright", fill=paste0(c(colTCmantle, colMuller2022),alpha), legend=c("TorsvikCocks2017 (MRF)", "MULLER2022"), bty="n", cex=2, inset=c(0.22,0.05))
plot(meProj, add=TRUE, lwd=7)
dev.off()
