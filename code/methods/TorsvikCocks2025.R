
feat <- data.frame(
	feature_collection=file.path("data/TorsvikCocks2025", c("CEED6_MICROCONTINENTS.shp","PHAB2024_LAND.shp", "PHAB_Exposed_Land_2023.shp" )),
	from=1000,
	to=0

)
rownames(feat) <- c("microcontinents", "coastlines", "land")

# the platemodel
TC2025 <- platemodel(
	rotation=c("data/TorsvikCocks2025/PHAB2023_ROTATION_ENGINE.rot"),
	features=feat
)

## coastsMantle <- reconstruct("coastlines", age=400, model=pt, anchor=0)
## microsMantle <- reconstruct("microcontinents", age=400, model=pt, anchor=0)
## coastsMag <- reconstruct("coastlines", age=400, model=pt, anchor=1)
## microsMag <- reconstruct("microcontinents", age=400, model=pt, anchor=1)



#plot(coastsMantle$geometry, col="#FF000044", border=NA)
#plot(coastsMag$geometry, col="#0000FF44", border=NA, add=TRUE)
