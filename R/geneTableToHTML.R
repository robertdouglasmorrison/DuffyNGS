# geneTableToHTML.R

# turn a set of sampleIDs and a data.frame into an HTML table with clickable gene pictures

`geneTableToHTMLandPlots` <- function( geneDF, sampleIDset, colorset, htmlFile="geneTable.html",
				results.path=".", html.path=".", N=100, makePlots=TRUE,
				genesToPlot=NULL, tailWidth=1500, minYmax=5, title="", label="", useLog=FALSE,
				geneNameColumn="GENE_ID", PLOT.FUN=NULL, ...){

	# make sure we can open a drawing device
	if ( ! (capabilities()[ "png"] )) stop( " Unable to open PNG device...")

	# get all the needed WIG data
	if ( makePlots && is.null(PLOT.FUN)) {
		cat( "\n\nPreparing to make Gene Plots...")
		WIG.defaults()
		WIGlist <- loadMultipleWIGs( sampleIDset, results.path=results.path)
	} else {
		WIGlist <- NULL
	}

	# make a subfolder for the plots
	pngFolder <- linkPngFolder <- paste( "pngPlots", sep=".")
	pngFolder <- file.path( html.path, pngFolder)
	if ( ! file.exists( pngFolder)) dir.create( pngFolder, recursive=T, showWarnings=F)

	if ( ! is.null( geneDF)) {

		# limit the HTML table and plots to the top N genes
		N <- min( N, nrow(geneDF))
		geneDF <- geneDF[ 1:N, ]
		genes <- geneDF$GENE_ID
	
		# OK make the HTML file
		htmlFile <- file.path( html.path, htmlFile)
		table2html( geneDF, fileout=htmlFile, title=title, linkPaths=linkPngFolder)
	}

	# now make those individual plots
	if ( makePlots) {
		if ( ! is.null( genesToPlot)) genes <- genesToPlot

		makeAllMultiWIGgenePlots( WIGlist, colorset, genes, path=pngFolder, tailWidth=tailWidth, 
				label=label, geneNameColumn=geneNameColumn, useLog=useLog, minYmax=minYmax,
				PLOT.FUN=PLOT.FUN, ...)

		# now save some global data for adding extra gene plots later...
		GPLOT_WIGlist <<- WIGlist
		GPLOT_colorset <<- colorset
		GPLOT_path <<- pngFolder
		GPLOT_tail <<- tailWidth
		GPLOT_label <<- label
		GPLOT_gname <<- geneNameColumn
	}

	return()
}


`extraGenesToHTMLandPlots` <- function( genesToPlot=NULL) {

	# turn a set of geneIDs into extra gene plots to supplement clickable HTML tables
	if ( is.null( genesToPlot)) return()
	cat( "\nMaking 'extra' Gene Plots... ")

	cat( "     Skipping.. Too slow for now...")
	return()

	# we can skip over plots that are already done, to save a bit of time
	nIn <- length( genesToPlot)
	doneGenes <- sub( ".png", "", dir( GPLOT_path), fixed=T)
	genesToPlot <- sort( setdiff( genesToPlot, doneGenes))
	nNow <- length( genesToPlot)
	cat( "     skipping ", (nIn-nNow), "already made...")
	if ( length( genesToPlot) < 1) return()

	# there is a chance that these are not true GENE_IDs, but instead other types of names...
	if ( GPLOT_gname != "GENE_ID") {
		
		gmap <- getCurrentGeneMap()
		where <- match( genesToPlot, gmap[[ GPLOT_gname]], nomatch=0)
		useGeneName <- rep( NA, times=length(genesToPlot))
		useGeneName[ where > 0] <- gmap$GENE_ID[ where]
		genesToPlot <- useGeneName[ ! is.na( useGeneName)]
	}

	# now make those individual plots
	makeAllMultiWIGgenePlots( GPLOT_WIGlist, GPLOT_colorset, genesToPlot, path=GPLOT_path, 
			tailWidth=GPLOT_tail, label=GPLOT_label, geneNameColumn=GPLOT_gname)
	
	return()
}



`geneTableToHTMLandPlotsCleanup` <- function() {

	# remove storage that made stuff easier
	try( rm( GPLOT_WIGlist, GPLOT_colorset, GPLOT_path, GPLOT_tail, GPLOT_label, GPLOT_gname, 
			envir=.GlobalEnv))
	return()
}
