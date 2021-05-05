# pipe.CellTypeProportions.R -- do various ways of calling immune cell subset proportions from gene transcription

`pipe.CellTypeProportions` <- function( sampleID, annotationFile="Annotation.txt", optionsFile="Options.txt", 
				speciesID=getCurrentSpecies(), results.path=NULL, geneUniverse=NULL,
				recalculate=c("none", "all", "missing", "profile", "deconvolution"), 
				makePlots=c( "final", "none", "all"),
				verbose=TRUE) {

	if ( length(sampleID) > 1) {
		cat( "\nWarning: expected a single sample ID..")
		sampleID <- sampleID[1]
	}

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
	cellAns <- NULL
	pcts1 <- pcts2 <- pcts3 <- pcts4 <- pcts5 <- pcts6 <- pcts7 <- pcts8 <- rep.int( NA, N_CellTypes)
	names(pcts1) <- names(pcts2) <- names(pcts3) <- names(pcts4) <- names(pcts5) <- names(pcts6) <- names(pcts7) <-  names(pcts8) <- names(cellTypeColors)
	TestNames <- c( "Fit.SteepDescent", "Fit.NLS", "Fit.GenSA", "Deconv.NLS", "Deconv.GenSA", "Deconv.Log2.NLS", "Deconv.Log2.GenSA", "Deconv.SteepDescent")
	
	recalculate <- match.arg( recalculate)
	makePlots <- match.arg( makePlots)
	deconvPlot <- (makePlots %in% c( "final", "all"))

	if ( file.exists(celltypeDetailsFile)) {
		cellAns <- read.csv( celltypeDetailsFile, as.is=T)
		if ( TestNames[1] %in% colnames(cellAns)) pcts1 <- as.numeric( cellAns[[ TestNames[1]]])
		if ( TestNames[2] %in% colnames(cellAns)) pcts2 <- as.numeric( cellAns[[ TestNames[2]]])
		if ( TestNames[3] %in% colnames(cellAns)) pcts3 <- as.numeric( cellAns[[ TestNames[3]]])
		if ( TestNames[4] %in% colnames(cellAns)) pcts4 <- as.numeric( cellAns[[ TestNames[4]]])
		if ( TestNames[5] %in% colnames(cellAns)) pcts5 <- as.numeric( cellAns[[ TestNames[5]]])
		if ( TestNames[6] %in% colnames(cellAns)) pcts6 <- as.numeric( cellAns[[ TestNames[6]]])
		if ( TestNames[7] %in% colnames(cellAns)) pcts7 <- as.numeric( cellAns[[ TestNames[7]]])
		if ( TestNames[8] %in% colnames(cellAns)) pcts8 <- as.numeric( cellAns[[ TestNames[8]]])
	}
	
	if (recalculate != "none" || is.null(cellAns)) {
		cat( "\nProcessing Sample:  ", sampleID, "\n")

		# we have 2 different functions, using 2-4 different fit tools & settings.  Do them all and average.
		# Note that the deconvolution tool seems to see log2 transformed data very differently
		myColor <- annT$Color[ match( sampleID, annT$SampleID)]
		if ( is.null(myColor) || is.na(myColor) || myColor == "") myColor <- "orchid1"

		# 1)  Steepest Descent of the 28-dimensional immune profile
		if ( is.null(cellAns) || (recalculate %in% c("all","profile")) || (recalculate == "missing" && all(is.na(pcts1)))) {
			cat( "\n1. Cell Type Profile: Fit by Steepest Descent:\n")
			ans1 <- fitCellTypeProfileFromFile( f=transcriptFile, sid=sampleID, col=myColor, max.iterations=200, 
						makePlots=makePlots, plot.path=celltype.path, algorithm='steep',
						geneUniverse=geneUniverse)
			if ( ! is.null(ans1)) pcts1 <- ans1$CellProportions
		}

		# 2)  Nonlinear Least Squares of the 28-dimensional immune profile
		if ( is.null(cellAns) || (recalculate %in% c("all","profile")) || (recalculate == "missing" && all(is.na(pcts2)))) {
			cat( "\n2. Cell Type Profile: Fit by Nonlinear Least Squares (NLS):\n")
			ans2 <- fitCellTypeProfileFromFile( f=transcriptFile, sid=sampleID, col=myColor, 
						makePlots=makePlots, plot.path=celltype.path, algorithm='nls',
						geneUniverse=geneUniverse)
			if ( ! is.null(ans2)) pcts2 <- ans2$CellProportions
		}

		# 3)  GenSA of the 28-dimensional immune profile
		if ( is.null(cellAns) || (recalculate %in% c("all","profile")) || (recalculate == "missing" && all(is.na(pcts3)))) {
			cat( "\n3. Cell Type Profile: Fit by Simulated Annealing (GenSA):\n")
			ans3 <- fitCellTypeProfileFromFile( f=transcriptFile, sid=sampleID, col=myColor, 
						makePlots=makePlots, plot.path=celltype.path, algorithm='GenSA',
						geneUniverse=geneUniverse)
			if ( ! is.null(ans3)) pcts3 <- ans3$CellProportions
		}

		# 4)  Transcriptome Deconvolution, by NLS using the 'port' algorithm...
		if ( is.null(cellAns) || (recalculate %in% c("all","deconvolution")) || (recalculate == "missing" && all(is.na(pcts4)))) {
			cat( "\n4. Cell Type Deconvolution:  Fit RPKM by Nonlinear Least Squares (NLS):\n")
			ans4 <- fileSet.TranscriptDeconvolution( files=transcriptFile, fids=sampleID, algorithm="port",
						useLog=FALSE, plot=deconvPlot, plot.path=celltype.path, plot.col=cellTypeColors,
						geneUniverse=geneUniverse, verbose=F)
			if ( ! is.null(ans4)) pcts4 <- ans4$BestFit
		}
		
		# 5)  Transcriptome Deconvolution, by NLS using the 'GenSA' simulated annealing algorithm...
		if ( is.null(cellAns) || (recalculate %in% c("all","deconvolution")) || (recalculate == "missing" && all(is.na(pcts5)))) {
			cat( "\n5. Cell Type Deconvolution:  Fit RPKM by Simulated Annealing (GenSA):\n")
			ans5 <- fileSet.TranscriptDeconvolution( files=transcriptFile, fids=sampleID, algorithm="GenSA",
						useLog=FALSE, plot=deconvPlot, plot.path=celltype.path, plot.col=cellTypeColors,
						geneUniverse=geneUniverse, verbose=F)
			if ( ! is.null(ans5)) pcts5 <- ans5$BestFit
		}
		
		# 6)  Transcriptome Deconvolution, by NLS using the 'port' algorithm...
		if ( is.null(cellAns) || (recalculate %in% c("all","deconvolution")) || (recalculate == "missing" && all(is.na(pcts6)))) {
			cat( "\n6. Cell Type Deconvolution:  Fit Log2(RPKM) by Nonlinear Least Squares (NLS):\n")
			cat( "           Not numerically stable..  Skip this method for now..\n")
			#ans6 <- fileSet.TranscriptDeconvolution( files=transcriptFile, fids=sampleID, algorithm="port",
			#				useLog=TRUE, plot=deconvPlot, plot.path=celltype.path, plot.col=cellTypeColors,,
			#				geneUniverse=geneUniverseverbose=F)
			#if ( ! is.null(ans6)) pcts6 <- ans6$BestFit
		}
		
		# 7)  Transcriptome Deconvolution, by NLS using the 'GenSA' simulated annealing algorithm...
		if ( is.null(cellAns) || (recalculate %in% c("all","deconvolution")) || (recalculate == "missing" && all(is.na(pcts7)))) {
			cat( "\n7. Cell Type Deconvolution:  Fit Log2(RPKM) by Simulated Annealing (GenSA):\n")
			ans7 <- fileSet.TranscriptDeconvolution( files=transcriptFile, fids=sampleID, algorithm="GenSA",
						useLog=TRUE, plot=deconvPlot, plot.path=celltype.path, plot.col=cellTypeColors,
						geneUniverse=geneUniverse, verbose=F)
			if ( ! is.null(ans7)) pcts7 <- ans7$BestFit
		}
		
		# 8)  Transcriptome Deconvolution, by Steepest Descent...
		if ( is.null(cellAns) || (recalculate %in% c("all","deconvolution")) || (recalculate == "missing" && all(is.na(pcts8)))) {
			cat( "\n8. Cell Type Deconvolution:  Fit RPKM by Steepest Descent:\n")
			ans8 <- fileSet.TranscriptDeconvolution( files=transcriptFile, fids=sampleID, algorithm="steep",
						useLog=FALSE, plot=deconvPlot, plot.path=celltype.path, plot.col=cellTypeColors,
						geneUniverse=geneUniverse, verbose=F)
			if ( ! is.null(ans8)) pcts8 <- ans8$BestFit
		}
	}
	
	# done with all methods, now merge and average
	cellM <- matrix( NA, nrow=N_CellTypes, ncol=8)
	colnames(cellM) <- TestNames
	cellM[ ,1] <- pcts1
	cellM[ ,2] <- pcts2
	cellM[ ,3] <- pcts3
	cellM[ ,4] <- pcts4
	cellM[ ,5] <- pcts5
	cellM[ ,6] <- pcts6
	cellM[ ,7] <- pcts7
	cellM[ ,8] <- pcts8
		
	cellMean <- apply( cellM, 1, mean, na.rm=T)
	cellMean <- round( cellMean * 100 / sum(cellMean), digits=3)
		
	# leave out the NLS of Log2 data for now...
	cellAns <- data.frame( "CellType"=names(cellTypeColors), "Final.Proportions"=cellMean, round(cellM[,-6],digits=3), stringsAsFactors=F)
	write.table( cellAns, celltypeDetailsFile, sep=",", quote=T, row.names=F)

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

	# package up the final details
	out <- vector( mode="list")
	out[[1]] <- ctpM
	names(out)[1] <- "Proportions.Matrix"
	outFile <- file.path( celltype.path, "All.CellType.Proportion.Details.csv")
	write.csv( ctpM, outFile)

	if ( length(levels) < 3) {
		out[[2]] <- ans
		names(out)[2] <- "Comparison.Results"
		outFile <- file.path( celltype.path, paste( "CellTypeCompare_", levelString, "_Details.csv", sep=""))
		write.csv( ans, outFile, row.names=F)
	} else {
		for (j in 1:length(ans)) {
			out[[j+1]] <- ans[[j]]
			names(out)[j+1] <- names(ans)[j]
			outFile <- file.path( celltype.path, paste( "CellTypeCompare_", names(ans)[j], "_Details.csv", sep=""))
			write.csv( ans[[j]], outFile, row.names=F)
		}
	}

	return( invisible( out))
}

