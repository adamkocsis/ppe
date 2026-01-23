#' Plotting biogeographic membership
#' Adam T. Kocsis (Erlangen, 2020-06-17)
#' CC-BY 4.0
#' @param mem Membership vector.
#' @param cols color vector.
#' @param bg Spatial object, background.
#' @param alpha alpha values of region colors
#' @param labels should the labels be plotted
#' @param gri icosa grid used for plotting
biogeoplot <- function(mem, cols=allHex, bg=land$geometry, alpha="99", labels=TRUE, gri=gr, crs="EPSG:4326",denser=TRUE,...){
	# empty vector for the colors
	member <- rep(NA, nrow(gri@faces))

	# every entry corresponds to a face, ordered in the grid
	names(member) <- rownames(gri@faces)

	# color every entry
	reorder <- cols[mem] # implies: names(reorder) <- names(mem)

	# assign colors to appropriate face
	member[names(mem)] <- paste0(reorder, alpha)

	# plot empty background
	if(crs!="EPSG:4326"){
		if(denser){
			bg <- smoothr::densify(bg)
		}
		bg <- st_transform(bg, crs)
	}
	plot(bg, col="#BBBBBBBB",border=NA, ...)

	# plot colors
	plot(gri, col=member, add=TRUE, border="gray60", crs=crs)

	if(labels){
		# centroids reordered to match the memberhsip vector
		cent <- centers(gri)[names(mem),]
		if(crs!="EPSG:4326"){
			centSF <- st_as_sf(as.data.frame(cent), coords=c("long", "lat"), crs="WGS84")
			centSF <- st_transform(centSF, crs=proj)
			cent[, c("long", "lat")] <- st_coordinates(centSF)
		}
		# plot the membership
		text(x=cent[,1], y=cent[,2], label=mem, cex=0.6)
	}
}
