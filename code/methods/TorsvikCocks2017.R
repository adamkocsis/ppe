
# where they are
dir <- "data/TorsvikCocks2017"

# features
features<- c("static_polygons"=file.path(dir,"Torsvik_Cocks_2016_Terranes.gpml"))

# The model
TC2017 <- platemodel(rotation=file.path(dir,"Torsvik_Cocks_HybridRotationFile.rot"), features=features)
