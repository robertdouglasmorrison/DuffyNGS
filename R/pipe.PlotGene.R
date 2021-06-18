# pipe.PlotGene.R

`pipe.PlotGene` <- function( sampleIDs, genes, annotationFile="Annotation.txt",
				optionsFile="Options.txt", results.path=NULL, targetID=NULL,
				colorColumn="Color", PLOT.FUN=NULL, plotFormat=c("","png","pdf"), 
				plot.path="Gene_Plots", keepShortGeneName=NULL, pileup.col=c(4,2,1), 
				png.width=1000, png.height=700, pdf.width=10, pdf.height=7, ...) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	if( is.null( targetID)) targetID <- getOptionValue( optT, "targetID", notfound="HsPf", verbose=F)
	setCurrentTarget( targetID)
	speciesSet <- getCurrentTargetSpecies()
	curSpecies <- getCurrentSpecies()
	on.exit( setCurrentSpecies( curSpecies))

	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	plotFormat <- match.arg( plotFormat)
	if (plotFormat != "") {
		if ( ! file.exists( plot.path)) dir.create( plot.path, recursive=T, showWarnings=F)
	}

	# turn the list of genes into the list of species to visit
	NG <- length( genes)
	gspecies <- rep.int( NA, NG)
	gptr <- rep.int( 0, NG)
	for ( speciesID in speciesSet) {
		setCurrentSpecies( speciesID)
		gmap <- getCurrentGeneMap()
		where <- match( genes, gmap$GENE_ID, nomatch=0)
		gspecies[ where > 0] <- speciesID
		gptr[ where > 0] <- where[ where > 0]

		where <- match( genes, gmap$NAME, nomatch=0)
		gspecies[ where > 0] <- speciesID
		gptr[ where > 0] <- where[ where > 0]
	}

	# at this point, we should know most all of them
	notfound <- which( is.na( gspecies))
	if ( length(notfound) > 0) {
		cat( "\nSome Gene names not found: ", genes[ notfound])
	}

	# make sure we know all these sampleIDs
	annT <- readAnnotationTable( annotationFile)
	doSamples <- intersect( sampleIDs, annT$SampleID)
	if ( length( doSamples) < 1) {
		cat( "\nNo Samples found that match annotation: ", annT$SampleID)
		return()
	}
	if ( length( doSamples) < length( sampleIDs)) {
		cat( "\nSome Samples not found: ", setdiff( sampleIDs, doSamples))
	}
	doMultiWIG <- ( length( doSamples) > 1)
	wigColors <- annT[[ colorColumn]][ match( doSamples, annT$SampleID)]

	# visit those species...
	doSpecies <- sort( unique( gspecies[ ! is.na(gspecies)]))
	for (speciesID in doSpecies) {

		# get the WIG data we need
		setCurrentSpecies(speciesID)
		gmap <- getCurrentGeneMap()
		prefix <- getCurrentSpeciesFilePrefix()
		if ( length( doSpecies) > 1) cat( "\n", speciesID, ": \n", sep="")

		if (doMultiWIG) {
			WIGlist <- loadMultipleWIGs( doSamples, speciesID, results.path)
		} else {
			wigfile <- paste( doSamples, prefix, "WIG.rda", sep=".")
			wigfile <- file.path( results.path, "wig", wigfile)
			load( wigfile)
		}

		# now visit each gene in this species
		# use chromosomal order for speed
		mygenes <- which( gspecies == speciesID)
		mygenes <- mygenes[ order( gptr[mygenes])]

		for ( ig in mygenes) {
			g <- gmap$GENE_ID[ gptr[ ig]]
			gname <- g
			if ( ! is.null( keepShortGeneName)) gname <- shortGeneName( g, keep=as.numeric(keepShortGeneName))

			if (doMultiWIG) {
				plotMultiWIGgene( WIGlist, colors=wigColors, gene=g, PLOT.FUN=PLOT.FUN, ...)
			} else {
				plotWIGgene( wiggles, gene=g, col=pileup.col, ...)
			}
			dev.flush()

			if ( plotFormat == "png") {
				plotfile <- paste( gname, "png", sep=".")
				plotfile <- file.cleanSpecialCharactersFromFileName( plotfile)
				plotfile <- file.path( plot.path, plotfile)
				dev.print( png, plotfile, width=png.width, height=png.height, "bg"="white")
			}
			if ( plotFormat == "pdf") {
				plotfile <- paste( gname, "pdf", sep=".")
				plotfile <- file.cleanSpecialCharactersFromFileName( plotfile)
				plotfile <- file.path( plot.path, plotfile)
				dev.print( pdf, plotfile, width=pdf.width, height=pdf.height, "bg"="white")
			}
			if (length(mygenes) > 1) cat( "  ",gname)
		}
	}
	if ( NG > 1) cat( "\n")
		
	return( )
}


`pipe.PlotRegion` <- function( sampleIDs, seqIDs, starts=NULL, stops=NULL, annotationFile="Annotation.txt",
				optionsFile="Options.txt", results.path=NULL, targetID=NULL,
				colorColumn="Color", plotFormat=c("","png","pdf"), 
				plot.path="Gene_Plots", keepShortGeneName=NULL, pileup.col=c(4,2,1), 
				png.width=1000, png.height=700, pdf.width=10, pdf.height=7, ...) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	if( is.null( targetID)) targetID <- getOptionValue( optT, "targetID", notfound="HsPf", verbose=F)
	setCurrentTarget( targetID)
	speciesSet <- getCurrentTargetSpecies()
	curSpecies <- getCurrentSpecies()
	on.exit( setCurrentSpecies( curSpecies))

	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	plotFormat <- match.arg( plotFormat)
	if (plotFormat != "") {
		if ( ! file.exists( plot.path)) dir.create( plot.path, recursive=T, showWarnings=F)
	}

	# turn the list of chromosomes into the list of species to visit
	NC <- length( seqIDs)
	cspecies <- rep.int( NA, NC)
	cptr <- rep.int( 0, NC)
	for ( speciesID in speciesSet) {
		setCurrentSpecies( speciesID)
		gmap <- getCurrentGeneMap()
		where <- match( seqIDs, gmap$SEQ_ID, nomatch=0)
		cspecies[ where > 0] <- speciesID
		cptr[ where > 0] <- where[ where > 0]
	}

	# at this point, we should know most all of them
	notfound <- which( is.na( cspecies))
	if ( length(notfound) > 0) {
		cat( "\nSome SeqID names not found: ", seqIDs[ notfound])
	}

	# make sure we know all these sampleIDs
	annT <- readAnnotationTable( annotationFile)
	doSamples <- intersect( sampleIDs, annT$SampleID)
	if ( length( doSamples) < 1) {
		cat( "\nNo Samples found that match annotation: ", annT$SampleID)
		return()
	}
	if ( length( doSamples) < length( sampleIDs)) {
		cat( "\nSome Samples not found: ", setdiff( sampleIDs, doSamples))
	}
	doMultiWIG <- ( length( doSamples) > 1)
	wigColors <- annT[[ colorColumn]][ match( doSamples, annT$SampleID)]

	# visit those species...
	doSpecies <- sort( unique( cspecies[ ! is.na(cspecies)]))
	for (speciesID in doSpecies) {

		# get the WIG data we need
		setCurrentSpecies(speciesID)
		gmap <- getCurrentGeneMap()
		prefix <- getCurrentSpeciesFilePrefix()
		if ( length( doSpecies) > 1) cat( "\n", speciesID, ": \n", sep="")

		if (doMultiWIG) {
			WIGlist <- loadMultipleWIGs( doSamples, speciesID, results.path)
		} else {
			wigfile <- paste( doSamples, prefix, "WIG.rda", sep=".")
			wigfile <- file.path( results.path, "wig", wigfile)
			load( wigfile)
		}

		# now visit each chromosome in this species
		# use chromosomal order for speed
		mychrs <- which( cspecies == speciesID)

		for ( ic in mychrs) {
			s <- gmap$SEQ_ID[ cptr[ ic]]
			mystart <- starts[1]
			mystop <- stops[1]
			wh <- match( s, seqIDs)
			if ( !is.null(starts) && length(starts) >= wh) mystart=starts[wh]
			if ( !is.null(stops) && length(stops) >= wh) mystop=stops[wh]

			if (doMultiWIG) {
				plotMultiWIGregion( WIGlist, colors=wigColors, seqid=s, position=mystart, end=mystop, ...)
			} else {
				plotWIGregion( wiggles, seqid=s, col=pileup.col, position=mystart, end=mystop, ...)
			}
			dev.flush()

			if ( plotFormat == "png") {
				plotfile <- paste( s, "png", sep=".")
				plotfile <- file.cleanSpecialCharactersFromFileName( plotfile)
				plotfile <- file.path( plot.path, plotfile)
				dev.print( png, plotfile, width=png.width, height=png.height, "bg"="white")
			}
			if ( plotFormat == "pdf") {
				plotfile <- paste( s, "pdf", sep=".")
				plotfile <- file.cleanSpecialCharactersFromFileName( plotfile)
				plotfile <- file.path( plot.path, plotfile)
				dev.print( png, plotfile, width=pdf.width, height=pdf.height, "bg"="white")
			}
			if (length(mychrs) > 1) cat( "  ", s)
		}
	}
	if ( NC > 1) cat( "\n")
		
	return( )
}

