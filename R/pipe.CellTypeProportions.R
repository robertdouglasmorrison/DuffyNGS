# pipe.CellTypeProportions.R -- do various ways of calling immune cell subset proportions from gene transcription

`pipe.CellTypeProportions` <- function( sampleID, annotationFile="Annotation.txt", optionsFile="Options.txt", 
				speciesID=getCurrentSpecies(), results.path=NULL, geneUniverse=NULL, genePercentage=NULL,
				recalculate=c("none", "all", "missing", "profile", "deconvolution", "nls", "GenSA","steep"), 
				mode=c("weighted.mean","average","best"), makePlots=c( "final", "none", "all"),
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
	# grab the gene expression values for later RMS deviation calcuations
	unitsColumn <- getExpressionUnitsColumn( optionsFile, verbose=verbose)
	units <- sub( "_[MU]$", "", unitsColumn)
	sampleMatrix <- expressionFileSetToMatrix( transcriptFile, sampleID)
	
	# set up the current defined cell type details
	CellTypeSetup( reload=TRUE)
	
	# modify the shape of the cell type data to reflect what we are fitting against
	reshapeCellTypeMatrix( f=transcriptFile, intensityColumn=unitsColumn, verbose=verbose)
	
	# get the details of the cell type data
	cellTypeColors <- getCellTypeColors()
	cellTypeNames <- names(cellTypeColors)
	N_CellTypes <- length( cellTypeColors)
	reference <- getCellTypeReference()
	cellTypeMatrix <- getCellTypeMatrix()

	# let the caller reduce the gene universe by percentage instead of by naming genes
	if ( ! is.null( genePercentage)) {
		if ( ! is.null( geneUniverse)) {
			cat( "\n  Warning: only one of 'geneUniverse' or 'genePercentage' can be non-NULL")
			return( NULL)
		}
		genePercentage <- as.numeric( genePercentage)
		if ( genePercentage > 100) {
			cat( "\n  Warning: 'genePercentage' must be in the range (0 to 1) or (2 to 100)")
			return( NULL)
		} else if ( genePercentage > 1.0) {
			genePercentage <- genePercentage / 100
		}
		# get the target matrix, and decide which genes we will keep, least DE get discarded
		de <- apply( cellTypeMatrix, 1, function(x) diff( range( x, na.rm=T)) / max( mean(x,na.rm=T,1)))
		ord <- order( de, decreasing=T)
		# decide how many genes we keep
		nGenesKeep <- round( nrow(cellTypeMatrix) * genePercentage)
		# that gives us our universe of genes
		keep <- ord[ 1:nGenesKeep]
		geneUniverse <- rownames(cellTypeMatrix)[ keep]
		cat( "\nPercentage of genes to use in fitting =", genePercentage*100, "\nThus reducing to a universe of", nGenesKeep, "most DE genes.\n")
	}

	# name for the file of results
	celltypeDetailsFile <- file.path( celltype.path, paste( sampleID, prefix, reference, "Proportions.csv", sep="."))
	celltypePlotFile <- file.path( celltype.path, paste( sampleID, prefix, reference, "Proportions.Final.pdf", sep="."))
	cellAns <- NULL
	pcts1 <- pcts2 <- pcts3 <- pcts4 <- pcts5 <- pcts6 <- rep.int( NA, N_CellTypes)
	names(pcts1) <- names(pcts2) <- names(pcts3) <- names(pcts4) <- names(pcts5) <- names(pcts6) <- cellTypeNames
	TestNames <- c( "Fit.SteepDescent", "Fit.NLS", "Fit.GenSA", "Deconv.NLS", "Deconv.GenSA", "Deconv.SteepDescent")
	
	recalculate <- match.arg( recalculate)
	makePlots <- match.arg( makePlots)
	deconvPlot <- (makePlots %in% c( "final", "all"))

	if ( file.exists(celltypeDetailsFile)) {
		cellAns <- read.csv( celltypeDetailsFile, as.is=T)
		# make sure existing data is in the currently expected order
		cellWhere <- match( cellTypeNames, cellAns$CellType)
		cellAns <- cellAns[ cellWhere, ]
		if ( TestNames[1] %in% colnames(cellAns)) pcts1 <- as.numeric( cellAns[[ TestNames[1]]])
		if ( TestNames[2] %in% colnames(cellAns)) pcts2 <- as.numeric( cellAns[[ TestNames[2]]])
		if ( TestNames[3] %in% colnames(cellAns)) pcts3 <- as.numeric( cellAns[[ TestNames[3]]])
		if ( TestNames[4] %in% colnames(cellAns)) pcts4 <- as.numeric( cellAns[[ TestNames[4]]])
		if ( TestNames[5] %in% colnames(cellAns)) pcts5 <- as.numeric( cellAns[[ TestNames[5]]])
		if ( TestNames[6] %in% colnames(cellAns)) pcts6 <- as.numeric( cellAns[[ TestNames[6]]])
	}
	
	if (recalculate != "none" || is.null(cellAns)) {
		cat( "\nProcessing Sample:  ", sampleID, "\n")

		# we have 2 different functions, using 2-4 different fit tools & settings.  Do them all and average.
		# Note that the deconvolution tool seems to see log2 transformed data very differently
		myColor <- annT$Color[ match( sampleID, annT$SampleID)]
		if ( is.null(myColor) || is.na(myColor) || myColor == "") myColor <- "orchid1"

		# 1)  Steepest Descent of the 28-dimensional immune profile
		if ( is.null(cellAns) || (recalculate %in% c("all","profile","steep")) || (recalculate == "missing" && all(is.na(pcts1)))) {
			cat( "\n1. Cell Type Profile: Fit Dimension Proportions by Steepest Descent:\n")
			ans1 <- fitCellTypeProfileFromFile( f=transcriptFile, sid=sampleID, col=myColor, max.iterations=200, 
						makePlots=makePlots, plot.path=celltype.path, algorithm='steep',
						geneUniverse=geneUniverse, verbose=verbose)
			if ( ! is.null(ans1)) {
				pcts1 <- ans1$CellProportions
				cellWhere <- match( cellTypeNames, names(pcts1))
				pcts1 <- pcts1[ cellWhere]
			}
		}

		# 2)  Nonlinear Least Squares of the 28-dimensional immune profile
		if ( is.null(cellAns) || (recalculate %in% c("all","profile","nls")) || (recalculate == "missing" && all(is.na(pcts2)))) {
			cat( "\n2. Cell Type Profile: Fit Dimension Proportions by Nonlinear Least Squares (NLS):\n")
			ans2 <- fitCellTypeProfileFromFile( f=transcriptFile, sid=sampleID, col=myColor, 
						makePlots=makePlots, plot.path=celltype.path, algorithm='nls',
						geneUniverse=geneUniverse, verbose=verbose)
			if ( ! is.null(ans2)) {
				pcts2 <- ans2$CellProportions
				cellWhere <- match( cellTypeNames, names(pcts2))
				pcts2 <- pcts2[ cellWhere]
			}
		}

		# 3)  GenSA of the 28-dimensional immune profile
		if ( is.null(cellAns) || (recalculate %in% c("all","profile","GenSA")) || (recalculate == "missing" && all(is.na(pcts3)))) {
			cat( "\n3. Cell Type Profile: Fit Dimension Proportions by Simulated Annealing (GenSA):\n")
			ans3 <- fitCellTypeProfileFromFile( f=transcriptFile, sid=sampleID, col=myColor, 
						makePlots=makePlots, plot.path=celltype.path, algorithm='GenSA',
						geneUniverse=geneUniverse, verbose=verbose)
			if ( ! is.null(ans3)) {
				pcts3 <- ans3$CellProportions
				cellWhere <- match( cellTypeNames, names(pcts3))
				pcts3 <- pcts3[ cellWhere]
			}
		}

		# 4)  Transcriptome Deconvolution, by NLS using the 'port' algorithm...
		if ( is.null(cellAns) || (recalculate %in% c("all","deconvolution","nls")) || (recalculate == "missing" && all(is.na(pcts4)))) {
			cat( "\n4. Cell Type Deconvolution:  Fit Gene", units, "by Nonlinear Least Squares (NLS):\n")
			ans4 <- fileSet.TranscriptDeconvolution( files=transcriptFile, fids=sampleID, algorithm="port",
						useLog=FALSE, plot=deconvPlot, plot.path=celltype.path, plot.col=cellTypeColors,
						geneUniverse=geneUniverse, verbose=verbose)
			if ( ! is.null(ans4)) {
				# deconvolution tool returns weights.  Trun them to percentages
				pcts4 <- ans4$BestFit * 100
				cellWhere <- match( cellTypeNames, rownames(pcts4))
				pcts4 <- pcts4[ cellWhere, 1]
			}
		}
		
		# 5)  Transcriptome Deconvolution, by NLS using the 'GenSA' simulated annealing algorithm...
		if ( is.null(cellAns) || (recalculate %in% c("all","deconvolution","GenSA")) || (recalculate == "missing" && all(is.na(pcts5)))) {
			cat( "\n5. Cell Type Deconvolution:  Fit Gene", units, "by Simulated Annealing (GenSA):\n")
			ans5 <- fileSet.TranscriptDeconvolution( files=transcriptFile, fids=sampleID, algorithm="GenSA",
						useLog=FALSE, plot=deconvPlot, plot.path=celltype.path, plot.col=cellTypeColors,
						geneUniverse=geneUniverse, verbose=verbose)
			if ( ! is.null(ans5)) {
				# deconvolution tool returns weights.  Trun them to percentages
				pcts5 <- ans5$BestFit * 100
				cellWhere <- match( cellTypeNames, rownames(pcts5))
				pcts5 <- pcts5[ cellWhere, 1]
			}
		}
		
		# 6)  Transcriptome Deconvolution, by Steepest Descent...
		if ( is.null(cellAns) || (recalculate %in% c("all","deconvolution","steep")) || (recalculate == "missing" && all(is.na(pcts6)))) {
			cat( "\n6. Cell Type Deconvolution:  Fit Gene", units, "by Steepest Descent:\n")
			ans6 <- fileSet.TranscriptDeconvolution( files=transcriptFile, fids=sampleID, algorithm="steep",
						useLog=FALSE, plot=deconvPlot, plot.path=celltype.path, plot.col=cellTypeColors,
						geneUniverse=geneUniverse, verbose=verbose)
			if ( ! is.null(ans6)) {
				# deconvolution tool returns weights.  Trun them to percentages
				pcts6 <- ans6$BestFit * 100
				cellWhere <- match( cellTypeNames, rownames(pcts6))
				pcts6 <- pcts6[ cellWhere, 1]
			}
		}
	}
	
	# done with all methods, now merge
	cellM <- matrix( NA, nrow=N_CellTypes, ncol=6)
	colnames(cellM) <- TestNames
	cellM[ ,1] <- pcts1
	cellM[ ,2] <- pcts2
	cellM[ ,3] <- pcts3
	cellM[ ,4] <- pcts4
	cellM[ ,5] <- pcts5
	cellM[ ,6] <- pcts6

	# catch any that are all zero
	isAllZero <- which( apply( cellM, 2, function(x) all( x == 0, na.rm=T)))
	if ( length(isAllZero)) cellM[ , isAllZero] <- NA

	# let's independently eavaluate the RMS fit for every method
	rmsdAns <- rmsDeviation( obsMatrix=sampleMatrix, calcMatrix=cellTypeMatrix, calcPcts=cellM, geneUniverse=geneUniverse)
	if (verbose) {
		cat( "\nGene expression RMS Deviation by method:\n")
		print( rmsdAns)
	}

	# final answer is either the average or the one method with lowest RMS
	# final answer will now be 'proportion', not 'percentage'.  Not summing to 100% anymore.
	mode <- match.arg( mode)
	if ( mode == "weighted.mean") {
		wts <- 1 / ((rmsdAns^2) / min(rmsdAns^2))
		names(wts) <- names(rmsdAns)
		if (verbose) {
			cat( "\nFit Weights by method:\n")
			print( wts);
		}
		cellMean <- apply( cellM, 1, weighted.mean, w=wts, na.rm=T)
		cellMean <- round( cellMean, digits=5)
		if (verbose) cat( "\n  Final answer is weighted mean of all methods, using relative RMSD as weights")
	} else if ( mode == "average") {
		wts <- rep.int( 1, 6)
		names(wts) <- names(rmsdAns)
		if (verbose) {
			cat( "\nFit Weights by method:\n")
			print( wts);
		}
		cellMean <- apply( cellM, 1, mean, na.rm=T)
		cellMean <- round( cellMean, digits=5)
		if (verbose) cat( "\n  Final answer is average of all methods")
	} else {
		best <- which.min(rmsdAns)
		wts <- rep.int( 0, 6)
		wts[ best] <- 1
		names(wts) <- names(rmsdAns)
		if (verbose) {
			cat( "\nFit Weights by method:\n")
			print( wts);
		}
		cellMean <- round( cellM[ , best], digits=5)
		if (verbose) cat( "\n  Final answer is method with best (lowest) RMSD: ", names(rmsdAns)[best])
	}
		
	cellAns <- data.frame( "CellType"=cellTypeNames, "Final.Proportions"=cellMean, round(cellM,digits=3), stringsAsFactors=F)
	write.table( cellAns, celltypeDetailsFile, sep=",", quote=T, row.names=F)

	# make a plot image of the final average
	if ( makePlots != "none") {
		tmpM <- cbind( cellM, cellMean)
		rownames(tmpM) <- cellTypeNames
		colnames(tmpM) <- sub( "\\.", "\n", c( colnames(cellM), "Final.Proportions"))
		plotTranscriptProportions( tmpM, col=cellTypeColors, label=paste( sampleID, ":  Final average of all methods", sep=""))
		printPlot( celltypePlotFile, width=12, height=9)
	}

	out <- data.frame( "SampleID"=sampleID, cellAns, stringsAsFactors=F)
	return( invisible( out))
}


`rmsDeviation` <- function( obsMatrix, calcMatrix, calcPcts, geneUniverse=NULL, DEBUG=FALSE) {

	# measure the Root Mean Square (RMS) deviation of each set of fitted cell type percentages
	# first get the set of genes to use
	obsGenes <- shortGeneName( rownames( obsMatrix), keep=1)
	calcGenes <- rownames( calcMatrix)
	useGenes <- sort( intersect( obsGenes, calcGenes))
	if ( ! is.null( geneUniverse)) {
		useGenes <- intersect( useGenes, geneUniverse)
	}
	NG <- length( useGenes)
	whereObs <- match( useGenes, obsGenes)
	whereCalc <- match( useGenes, calcGenes)

	# force the reference to scale match the observed, to be exactly what the deconvolution tools saw
	NcellTypes <- nrow(calcPcts)
	Nmethods <- ncol( calcPcts)
	obsValue <- obsMatrix[ whereObs, 1]
	obsSum <- sum( obsValue, na.rm=T)
	calcValues <- calcMatrix[ whereCalc, ]
	calcSums <- apply( calcValues, 2, sum, na.rm=T)
	calcScales <- obsSum / calcSums
	for ( j in 1:NcellTypes) calcValues[ , j] <- calcValues[ ,j] * calcScales[j]

	if (DEBUG) {
		cat( "\n  Debug:  Total Observed Expression: ", obsSum)
		cat( "\n  Debug:  Total Calculated Expression: \n")
		print( apply( calcValues, 2, sum, na.rm=T))
	}

	# now calc the deviation for each proportions answer
	rmsOut <- rep.int( NA, Nmethods)
	names(rmsOut) <- colnames(calcPcts)
	for ( j in 1:Nmethods) {
		if ( all( is.na( calcPcts[ ,j]))) next
		# build a calculated transcriptome
		thisCalcValue <- rep.int( 0, NG)
		for ( i in 1:NcellTypes) {
			# turn the given percentages back to weights
			thisCellValue <- calcValues[ , i] * calcPcts[ i, j] / 100
			thisCellValue[ is.na( thisCellValue)] <- 0
			thisCalcValue <- thisCalcValue + thisCellValue
		}
		vDiff <- obsValue - thisCalcValue
		rms <- sqrt( mean( vDiff * vDiff))
		rmsOut[j] <- round( rms, digits=3)
	}
	return( rmsOut)
}


`pipe.CellTypeCompare` <- function( sampleIDset, groups=sampleIDset, levels=sort(unique(as.character(groups))), 
				annotationFile="Annotation.txt", optionsFile="Options.txt", 
				speciesID=getCurrentSpecies(), results.path=NULL, cellProportionsColumn="Final.Proportions", 
				prefix=NULL, 
				minPerGroup=3, test=t.test, plot=TRUE, plot.mode=c("bars", "auto", "pie", "lines"),
				label="", wt.fold=1, wt.pvalue=2, min.percent=0.1, 
				significance.scaling=TRUE, verbose=TRUE, ...) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	annT <- readAnnotationTable( annotationFile)

	# set up for this species
	if ( speciesID != getCurrentSpecies()) {
		setCurrentSpecies(speciesID)
	}
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	celltype.path <- file.path( results.path, "CellTypeProportions")
	if ( ! file.exists( celltype.path)) {
		cat( "\nExpected folder of Cell Type Proportion calls not found: ", celltype.path)
		return(NULL)
	}

	# set up if needed
	CellTypeSetup()
	cellTypeColors <- getCellTypeColors()
	N_CellTypes <- length( cellTypeColors)
	reference <- getCellTypeReference()

	# gather the cell type proportions data from every sample
	NS <- length( sampleIDset)
	if ( length(groups) != NS) stop( "Length of 'groups' must equal number of samples")
	# allow passsing in of explicit prefix to allow mixed organism copmares..
	if ( is.null( prefix)) {
		prefix <- getCurrentSpeciesFilePrefix()
	}
	prefix <- rep( prefix, length.out=NS)

	ctpM <- matrix( NA, nrow=N_CellTypes, ncol=NS)
	colnames(ctpM) <- sampleIDset
	rownames(ctpM) <- names( cellTypeColors)

	for ( i in 1:NS) {
		sid <- sampleIDset[i]

		# name for  file of results
		my.celltype.path <- file.path( celltype.path, sid)
		celltypeDetailsFile <- file.path( my.celltype.path, paste( sid, prefix[i], reference, "Proportions.csv", sep="."))
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

		# stash these values, allowing some to be missing
		myCellNames <- cellAns$CellType
		where <- match( rownames(ctpM), myCellNames)
		ctpM[ , i] <- myCellPcts[ where]
	}

	# now with this matrix of cell type percentages, call the comparison tool
	plot.mode <- match.arg( plot.mode)
	ans <- compareTranscriptProportions( ctpM, groups=groups, levels=levels, col=cellTypeColors, 
					minPerGroup=minPerGroup, test=test, plot=plot, plot.mode=plot.mode, label=label, 
					wt.fold=wt.fold, wt.pvalue=wt.pvalue, min.percent=min.percent, 
					significance.scaling=significance.scaling, ...)

	# if we did plot, decide what to call it...
	dev.type <- getPlotDeviceType( optT)
	levelString <- paste( levels, collapse=".v.")
	if ( length(levels) > 3) levelString <- paste( paste( levels[1:3], collapse=".v."), ".v.etc", sep="")
	plotFile <- paste( "Compare", prefix[1], reference, "_", levelString, "_", toupper(plot.mode), ".Plot.", dev.type, sep="")
	printPlot( file.path( celltype.path, plotFile))

	# package up the final details
	out <- vector( mode="list")
	out[[1]] <- ctpM
	names(out)[1] <- "Proportions.Matrix"
	outFile <- file.path( celltype.path, paste( "All", prefix[1], reference, "Proportion.Details.csv", sep="."))
	write.csv( ctpM, outFile)

	if ( length(levels) < 3) {
		out[[2]] <- ans
		names(out)[2] <- "Comparison.Results"
		outFile <- file.path( celltype.path, paste( prefix[1], ".", reference, ".Compare_", levelString, "_Details.csv", sep=""))
		write.csv( ans, outFile, row.names=F)
	} else {
		for (j in 1:length(ans)) {
			out[[j+1]] <- ans[[j]]
			names(out)[j+1] <- names(ans)[j]
			outFile <- file.path( celltype.path, paste( prefix[1], ".", reference, ".Compare_", names(ans)[j], "_Details.csv", sep=""))
			write.csv( ans[[j]], outFile, row.names=F)
		}
	}

	return( invisible( out))
}

