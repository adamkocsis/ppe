library(rgplates)
library(chronosphere)

setwd("/mnt/sky/Dropbox/WorkSpace/2025-06-10_PPE/")

# get the TorsvikCocks model (rgplates::platemodel class object)
source("code/methods/TorsvikCocks2025.R")


# get the model data from the chronosphere
MERDITH2021 <- chronosphere::fetch("GPlates", "MERDITH2021", datadir="data/chronosphere")
CAO2024 <- chronosphere::fetch("GPlates", "CAO2024", datadir="data/chronosphere")
PALEOMAP <- chronosphere::fetch("paleomap", "model", datadir="data/chronosphere")

# The chosen projection
proj <- "ESRI:54030"


me <- rgplates::mapedge()
meProj <- st_transform(me, proj)

# GWS - based
platesTC2017_mantle_gws <- rgplates::reconstruct("static_polygons", age=400, model="TorsvikCocks2017", anchor=0)
platesTC2017_mag_gws <- rgplates::reconstruct("static_polygons", age=400, model="TorsvikCocks2017", anchor=1)

# project
platesTC2017_mantle_gwsProj <- st_transform(platesTC2017_mantle_gws, crs=proj)
platesTC2017_mag_gwsProj <- st_transform(platesTC2017_mag_gws, crs=proj)

# Two trond models (to be upgraded to 2025
plot(meProj)
plot(platesTC2017_mantle_gwsProj$geometry, col="#44444466", add=TRUE, border=NA)
plot(platesTC2017_mag_gwsProj$geometry, col="#99999966", add=TRUE, border=NA)

# only mantle reference

# only mantle reference
