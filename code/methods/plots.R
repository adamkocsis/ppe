#' Plotting biogeographic membership
#' Adam T. Kocsis (Erlangen, 2020-06-17)
#' CC-BY 4.0
#' @param mem Membership vector.
#' @param cols color vector.
#' @param bg Spatial object, background.
#' @param alpha alpha values of region colors
#' @param labels should the labels be plotted
#' @param gri icosa grid used for plotting
biogeoplot <- function(mem, cols=allHex, bg=land$geometry, alpha="99", labels=TRUE, gri=gr,gri.border="gray60", gri.lwd=1,crs="EPSG:4326",denser=TRUE,...){
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
	plot(gri, col=member, add=TRUE, border=gri.border, crs=crs, lwd=gri.lwd)

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

#' Plot a paleogeographic map with PBDB occurrences
#'
#' A function by Adam to make a single nice map
#' 
#' @param x Occurrence data frame
#' @param map A single sf object for the background
#' @param crs CRS string for the entire projection.
#' @param plng Paleolongitude column to be used from x
#' @param plat Paleolatitude column to be used from x
#' @param main Top row title
#' @param graticules.lwd Line width of graticules
#' @param points.cex cex of points
#' @param pbdb.col Color of the PBDB
#' @param gbdb.col Color of the GBDB
#' @param map.bg Ocean backgroudn color 
#' @param symbol Symbol of the geological period 
#' @param smybol.col Color of the geological period symbol
#' @param map.bgdamp Should the background be damped?
#' @param coloredsub Should the legend be a colored text - TRUE is very sensitive to resizing of the plot!
PlotOccs <- function(x, map, crs="ESRI:4326", plng="plng510", plat="plat510", 
	points.cex=1, col="red", map.bg="#1A6BB0", symbol=NULL,
	symbol.col=NULL, map.bgdamp=FALSE, coloredsub=TRUE){

	# get rid of missing data
	x <- x[!is.na(x[, plng]) | !is.na(x[, plat]),]

	# grab the collections
	collections <- unique(x[c("collection_no", plng, plat)])

	# project the map
	mapProj <- sf::st_transform(map, crs)

	# transform coordinates
	sfcolls <- sf::st_as_sf(collections, coords=c(plng,plat), crs="WGS84")
	slcProj <- sf::st_transform(sfcolls, crs)

	# make a new plot
	me <- rgplates::mapedge(crs=crs)
	plot(me, reset=FALSE, col=NA)
	sphereshade(left="#0f3f67", right="#1A6BB0", crs=crs)
	if(map.bgdamp) {
		plot(me, add=TRUE, col="#ffffff88")
	}
	plot(mapProj$geometry, add=TRUE, col="gray")
	plot(slcProj$geometry, col="black", bg=col, add=TRUE, pch=21, cex=points.cex, lwd=3)
	plot(me, reset=FALSE, col=NA, add=TRUE)
	
}

#' Plot a Carbonates 
#'
#' A function by Adam to make a single nice map
#' 
#' @param x Occurrence data frame
#' @param map A single sf object for the background
#' @param crs CRS string for the entire projection.
#' @param plng Paleolongitude column to be used from x
#' @param plat Paleolatitude column to be used from x
#' @param main Top row title
#' @param graticules.lwd Line width of graticules
#' @param points.cex cex of points
#' @param pbdb.col Color of the PBDB
#' @param gbdb.col Color of the GBDB
#' @param map.bg Ocean backgroudn color 
#' @param symbol Symbol of the geological period 
#' @param smybol.col Color of the geological period symbol
#' @param map.bgdamp Should the background be damped?
#' @param coloredsub Should the legend be a colored text - TRUE is very sensitive to resizing of the plot!
PlotLithology <- function(x, ras, log=TRUE, proj="ESRI:4326", plng="plng510", plat="plat510",
	points.cex=1, col , coloramp, pch=21){

	# get rid of missing data
	x <- x[!is.na(x[, plng]) | !is.na(x[, plat]),]

	# grab the collections
	collectionsIn <- unique(x[c("collection_no", plng, plat)])

	# project the map
#	ras <- resample(ras, rast(res=0.05))
	mapProj <- project(ras, y=proj)

	# transform coordinates
	sfcolls <- st_as_sf(collectionsIn, coords=c(plng,plat), crs="WGS84")
	slcProj <- st_transform(sfcolls, proj)

	# make a new plot
	me <- rgplates::mapedge(crs=proj)
	if(log){
		plot(log(mapProj), col=coloramp$col, breaks=coloramp$breaks, legend=FALSE, axes=FALSE)
		plot(me, col="white", add=TRUE)
		plot(log(mapProj), col=coloramp$col, breaks=coloramp$breaks, legend=FALSE, axes=FALSE, add=TRUE)
	}else{
		plot(mapProj, col=coloramp$col, breaks=coloramp$breaks, legend=FALSE, axes=FALSE)
	}
	plot(slcProj$geometry, col="black", bg=col, add=TRUE, pch=pch, cex=points.cex, lwd=3)
	sphereshade(left="#000000", right="#000000", crs=proj, left.alpha=0.15, right.alpha=0)
	plot(me, reset=FALSE, col=NA, add=TRUE)
	
}

#' Plot shading over a geographic map
#'
#' An eye-candy utility function for 3d-effect of global maps
#'
#' The function plots a gradient on a geographic projection, either as a color or as an alpha gradient, the function relies on \link{\code{mapedge}} for the calculation of bands.
#' @param n The number of color bands.
#' @param left RGB (no alpha!) color for the left side of the gradient, input to \link[grDevices]{\code{colorRampPalette}}.
#' @param right RGB (no alpha!) color for the right side of the gradient, input to \link[grDevices]{\code{colorRampPalette}}.
#' @param crs A coordinate reference system string to be passed to \link{\code{mapedge}}.
#' @param left.alpha The alpha value (0-1) used for the left side of the plot.
#' @param right.alpha The alpha value (0-1) used for the right side of the plot.
#' @param west The \code{xmin} argument of \link{\code{mapedge}}.
#' @param east The \code{xmax} argument of \link{\code{mapedge}}.
#' @export
#' @examples
#' # basic use
#' me <- mapedge()
#' plot(me)
#' sphereshade()
#'
#' # In Robinson projection
#' me <- mapedge(crs="ESRI:54030")
#' plot(me)
#' sphereshade(left="#0f3f67", right="#1A6BB0", crs="ESRI:54030")
sphereshade <- function(n=180, left="#ffffff", right="#aaaaaa", crs="EPSG:4326", left.alpha=1, right.alpha=1, west=-180, east=180){
	if(length(n)!=1 | !is.numeric(n)) stop("The number of bands ('n') has to be a single positive integer.")
	if(length(left.alpha)!=1 | length(right.alpha)!=1) stop("You can only provide a single alpha value.")
	# create breakpoints
	breaks <- seq(west, east, length.out=n+1)

	# create a color ramp
	colFun <- grDevices::colorRampPalette(c(left, right))
	cols <- colFun(n)

	# transparency is needed
	alphas <- seq(left.alpha, right.alpha, length.out=n)
	alphaVals <- as.character(as.hexmode(round(alphas*255)))
	alphaVals[nchar(alphaVals)==1] <- paste0("0", alphaVals[nchar(alphaVals)==1])


	# plot them sequentiallly
	for(i in 1:n){
		oneSlice <- mapedge(xmin=breaks[i],xmax=breaks[i+1],crs=crs)
		plot(oneSlice, border=NA, col=paste0(cols[i], as.character(alphaVals[i])), add=TRUE)
	}

}
