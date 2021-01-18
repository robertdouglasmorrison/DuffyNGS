# pipe.CellTypeProportions.R -- do various ways of calling immune cell subset proportions from gene transcription

`pipe.CellTypeProportions` <- function( sampleID, annotationFile="Annotation.txt", optionsFile="Options.txt", 
				speciesID=getCurrentSpecies(), results.path=NULL, recalculate=FALSE, 
				verbose=TRUE) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	annT <- readAnnotationTable( annotationFile)

	# set up for this species
	if ( speciesID != getCurrentSpecies()) {
		setCurrentSpecies(speciesID)
	}
	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	celltype.path <- file.path( results.path, "CellTypeProportions", sampleID)
	if ( ! file.exists( celltype.path)) dir.create( celltype.path, recursive=TRUE)

	# we will use the transcriptome as our input data
	transcriptFile <- file.path( results.path, "transcript", paste( sampleID, prefix, "Transcript.txt", sep="."))
	if ( ! file.exists( transcriptFile)) {
		cat( "\nError: required transcriptome result file not found: ", transcriptFile)
		return( NULL)
	}
	
	# set up if needed
	verifyCellTypeSetup()
	cellTypeColors <- getCellTypeColors()
	N_CellTypes <- length( cellTypeColors)

	# name for the file of results
	celltypeDetailsFile <- file.path( celltype.path, paste( sampleID, prefix, "CellTypeProportions.csv", sep="."))

	if (recalculate || ! file.exists(celltypeDetailsFile)) {
		cat( "\nProcessing Sample:  ", sampleID, "\n")

		# we have 2 different functions, using 2 different fit tools.  Do them all and average.
		# and the deconvolution tool seems to see log2 transsformed data very differently
		myColor <- annT$Color[ match( sampleID, annT$SampleID)]
		if ( is.null(myColor) || is.na(myColor) || myColor == "") myColor <- "orchid1"

		pcts1 <- pcts2 <- pcts3 <- pcts4 <- pcts5 <- pcts6 <- rep.int( NA, N_CellTypes)
		names(pcts1) <- names(pcts2) <- names(pcts3) <- names(pcts4) <- names(pcts5) <- names(pcts6) <- names(cellTypeColors)
		
		# 1)  Steepest Descent of the 28-dimensional immune profile
		cat( "\n1. Cell Type Profile: Fit by Steepest Descent:\n")
		ans1 <- fitCellTypeProfileFromFile( f=transcriptFile, sid=sampleID, col=myColor, max.iterations=200, 
						makePlot='final', plot.path=celltype.path, algorithm='steep')
		if ( ! is.null(ans1)) pcts1 <- ans1$CellProportions

		# 2)  Nonlinear Least Squares of the 28-dimensional immune profile
		cat( "\n2. Cell Type Profile: Fit by Nonlinear Least Squares (NLS):\n")
		ans2 <- fitCellTypeProfileFromFile( f=transcriptFile, sid=sampleID, col=myColor, 
						makePlot='final', plot.path=celltype.path, algorithm='nls')
		if ( ! is.null(ans2)) pcts2 <- ans2$CellProportions

		# 3)  Transcriptome Deconvolution, by NLS using the 'port' algorithm...
		cat( "\n3. Cell Type Deconvolution:  Fit RPKM by Nonlinear Least Squares (NLS):\n")
		ans3 <- fileSet.TranscriptDeconvolution( files=transcriptFile, fids=sampleID, algorithm="port",
						useLog=FALSE, plot=TRUE, plot.path=celltype.path, plot.col=cellTypeColors)
		if ( ! is.null(ans3)) pcts3 <- ans3$BestFit
		
		# 4)  Transcriptome Deconvolution, by NLS using the 'GenSA' simulated annealing algorithm...
		cat( "\n4. Cell Type Deconvolution:  Fit RPKM by Simulated Annealing (GenSA):\n")
		ans4 <- fileSet.TranscriptDeconvolution( files=transcriptFile, fids=sampleID, algorithm="GenSA",
						useLog=FALSE, plot=TRUE, plot.path=celltype.path, plot.col=cellTypeColors)
		if ( ! is.null(ans4)) pcts4 <- ans4$BestFit
		
		# 5)  Transcriptome Deconvolution, by NLS using the 'port' algorithm...
		cat( "\n5. Cell Type Deconvolution:  Fit Log2(RPKM) by Nonlinear Least Squares (NLS):\n")
		cat( "\n           Not numerically stable..  Skip this for now..")
		#ans5 <- fileSet.TranscriptDeconvolution( files=transcriptFile, fids=sampleID, algorithm="port",
		#				useLog=TRUE, plot=TRUE, plot.path=celltype.path, plot.col=cellTypeColors)
		#if ( ! is.null(ans5)) pcts5 <- ans5$BestFit
		
		# 6)  Transcriptome Deconvolution, by NLS using the 'GenSA' simulated annealing algorithm...
		cat( "\n6. Cell Type Deconvolution:  Fit Log2(RPKM) by Simulated Annealing (GenSA):\n")
		ans6 <- fileSet.TranscriptDeconvolution( files=transcriptFile, fids=sampleID, algorithm="GenSA",
						useLog=TRUE, plot=TRUE, plot.path=celltype.path, plot.col=cellTypeColors)
		if ( ! is.null(ans6)) pcts6 <- ans6$BestFit

		# done with all 4 methods, now merge and average
		cellM <- matrix( NA, nrow=N_CellTypes, ncol=6)
		colnames(cellM) <- c( "Fit.SteepDescent", "Fit.NLS", "Deconv.NLS", "Deconv.GenSA", "Deconv.Log2.NLS", "Deconv.Log2.GenSA")
		cellM[ ,1] <- pcts1
		cellM[ ,2] <- pcts2
		cellM[ ,3] <- pcts3
		cellM[ ,4] <- pcts4
		cellM[ ,5] <- pcts5
		cellM[ ,6] <- pcts6
		
		cellMean <- apply( cellM, 1, mean, na.rm=T)
		cellMean <- round( cellMean * 100 / sum(cellMean), digits=3)
		
		# leave out the NLS of Log2 data for now...
		cellAns <- data.frame( "CellType"=names(pcts1), "Final.Proportions"=cellMean, round(cellM[,-5],digits=3), stringsAsFactors=F)
		
		write.table( cellAns, celltypeDetailsFile, sep=",", quote=T, row.names=F)
	} else {
		cat( "\nUsing Previously Processed Sample:  ", sampleID, "\n")
		cellAns <- read.csv( celltypeDetailsFile, as.is=T)
	}

	out <- data.frame( "SampleID"=sampleID, cellAns, stringsAsFactors=F)
	return( invisible( out))
}


`pipe.CellTypeCompare` <- function( sampleIDset, groups, levels=sort(unique(as.character(groups))), 
				annotationFile="Annotation.txt", optionsFile="Options.txt", 
				speciesID=getCurrentSpecies(), results.path=NULL, cellProportionsColumn="Final.Proportions", 
				minPerGroup=3, test=t.test, plot=TRUE, plot.mode=c("bars", "auto", "pie", "lines"),
				label="", wt.fold=1, wt.pvalue=2, min.percent=0.1, verbose=TRUE, ...) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	annT <- readAnnotationTable( annotationFile)

	# set up for this species
	if ( speciesID != getCurrentSpecies()) {
		setCurrentSpecies(speciesID)
	}
	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	celltype.path <- file.path( results.path, "CellTypeProportions")
	if ( ! file.exists( celltype.path)) {
		cat( "\nExpected folder of Cell Type Proportion calls not found: ", celltype.path)
		return(NULL)
	}

	# set up if needed
	verifyCellTypeSetup()
	cellTypeColors <- getCellTypeColors()
	N_CellTypes <- length( cellTypeColors)

	# gather the cell type proportions data from every sample
	NS <- length( sampleIDset)
	if ( length(groups) != NS) stop( "Length of 'groups' must equal number of samples")
	ctpM <- matrix( NA, nrow=N_CellTypes, ncol=NS)
	colnames(ctpM) <- sampleIDset
	rownames(ctpM) <- names( cellTypeColors)

	for ( i in 1:NS) {
		sid <- sampleIDset[i]

		# name for  file of results
		my.celltype.path <- file.path( celltype.path, sid)
		celltypeDetailsFile <- file.path( my.celltype.path, paste( sid, prefix, "CellTypeProportions.csv", sep="."))
		if ( ! file.exists(celltypeDetailsFile)) {
			cat( "\nCell Type Proportions results not found for sample: ", sid, "|", celltypeDetailsFile)
			next
		}
		cellAns <- read.csv( celltypeDetailsFile, as.is=T)
		myCellPcts <- cellAns[[ cellProportionsColumn]]
		if ( is.null( myCellPcts)) {
			cat( "\nError: No column of cell type proportions file has name '", cellProportionColumn, "'", sep="")
			return(NULL)
		}

		# stash these values
		ctpM[ , i] <- myCellPcts
	}

	# now with this matrix of cell type percentages, call the comparison tool
	plot.mode <- match.arg( plot.mode)
	ans <- compareTranscriptProportions( ctpM, groups=groups, levels=levels, col=cellTypeColors, 
					minPerGroup=minPerGroup, test=test, plot=plot, plot.mode=plot.mode, label=label, 
					wt.fold=wt.fold, wt.pvalue=wt.pvalue, min.percent=min.percent, ...)

	# if we did plot, decide what to call it...
	levelString <- paste( levels, collapse=".v.")
	plotFile <- paste( "CellTypeCompare_", levelString, "_", toupper(plot.mode), ".Plot.png", sep="")
	dev.print( png, file.path( celltype.path, plotFile), width=900)
	plotFile <- paste( "CellTypeCompare_", levelString, "_", toupper(plot.mode), ".Plot.pdf", sep="")
	dev.print( pdf, file.path( celltype.path, plotFile), width=12)

	out <- vector( mode="list")
	out[[1]] <- ctpM
	names(out)[1] <- "Proportions.Matrix"

	if ( length(levels) < 3) {
		out[[2]] <- ans
		names(out)[2] <- "Comparison.Results"
	} else {
		for (j in 1:length(ans)) {
			out[[j+1]] <- ans[[j]]
			names(out)[j+1] <- names(ans)[j]
		}
	}

	return( invisible( out))
}

