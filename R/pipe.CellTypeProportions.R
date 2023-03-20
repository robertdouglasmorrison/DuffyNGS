# pipe.CellTypeProportions.R -- do various ways of calling immune cell subset proportions from gene transcription

`pipe.CellTypeProportions` <- function( sampleID, annotationFile="Annotation.txt", optionsFile="Options.txt", 
				speciesID=getCurrentSpecies(), results.path=NULL, geneUniverse=NULL, genePercentage=NULL,
				recalculate=c("none", "all", "missing", "profile", "deconvolution", "nls", "GenSA","steep"), 
				mode=c("weighted.mean","average","best"), makePlots=c( "final", "none", "all"),
				verbose=TRUE, ...) {

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
	celltypeStatsFile <- file.path( celltype.path, paste( sampleID, prefix, reference, "Fit.Statistics.csv", sep="."))
	celltypePlotFile <- file.path( celltype.path, paste( sampleID, prefix, reference, "Proportions.Final.pdf", sep="."))
	pcts1 <- pcts2 <- pcts3 <- pcts4 <- pcts5 <- pcts6 <- rep.int( NA, N_CellTypes)
	pcts7 <- pcts8 <- pcts9 <- pcts10 <- pcts11 <- pcts12 <- rep.int( NA, N_CellTypes)
	names(pcts1) <- names(pcts2) <- names(pcts3) <- names(pcts4) <- names(pcts5) <- names(pcts6) <- cellTypeNames
	names(pcts7) <- names(pcts8) <- names(pcts9) <- names(pcts10) <- names(pcts11) <- names(pcts12) <- cellTypeNames
	TestNames <- c( "Fit.Steep", "Fit.NLS", "Fit.GenSA", "Deconv.NLS", "Deconv.GenSA", "Deconv.Steep",
			"Fit.Steep.QN", "Fit.NLS.QN", "Fit.GenSA.QN", "Deconv.NLS.QN", "Deconv.GenSA.QN", "Deconv.Steep.QN")
	
	recalculate <- match.arg( recalculate)
	makePlots <- match.arg( makePlots)
	deconvPlot <- (makePlots %in% c( "final", "all"))

	cellAns <- NULL
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
		if ( TestNames[7] %in% colnames(cellAns)) pcts7 <- as.numeric( cellAns[[ TestNames[7]]])
		if ( TestNames[8] %in% colnames(cellAns)) pcts8 <- as.numeric( cellAns[[ TestNames[8]]])
		if ( TestNames[9] %in% colnames(cellAns)) pcts9 <- as.numeric( cellAns[[ TestNames[9]]])
		if ( TestNames[10] %in% colnames(cellAns)) pcts10 <- as.numeric( cellAns[[ TestNames[10]]])
		if ( TestNames[11] %in% colnames(cellAns)) pcts11 <- as.numeric( cellAns[[ TestNames[11]]])
		if ( TestNames[12] %in% colnames(cellAns)) pcts12 <- as.numeric( cellAns[[ TestNames[12]]])
	}
	
	if (recalculate != "none" || is.null(cellAns)) {
		cat( "\nProcessing Sample:  ", sampleID, "\n")

		# we have 2 different methods, using 3 different fit tools, and 2 different reference normalization modes.  Do them all and average.
		myColor <- annT$Color[ match( sampleID, annT$SampleID)]
		if ( is.null(myColor) || is.na(myColor) || myColor == "") myColor <- "orchid1"

		# 1)  Steepest Descent of the 27-dimensional immune profile
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

		# 2)  Nonlinear Least Squares of the 27-dimensional immune profile
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

		# 3)  GenSA of the 27-dimensional immune profile
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
	# save a copy of the original cell type matrix
	cellTypeMatrix.orig <- cellTypeMatrix

	# now modify the shape of the cell type data to reflect the actual transcriptome we are fitting against
	reshapeCellTypeMatrix( f=transcriptFile, intensityColumn=unitsColumn, verbose=verbose)
	cellTypeMatrix <- getCellTypeMatrix()
	cellTypeMatrix.qn <- cellTypeMatrix
	
	# now call all 6 tools again, using this normalized reference
	if (recalculate != "none" || is.null(cellAns)) {

		# 7)  Steepest Descent of the 27-dimensional immune profile
		if ( is.null(cellAns) || (recalculate %in% c("all","profile","steep")) || (recalculate == "missing" && all(is.na(pcts7)))) {
			cat( "\n7. Cell Type Profile: Fit QN Dimension Proportions by Steepest Descent:\n")
			ans7 <- fitCellTypeProfileFromFile( f=transcriptFile, sid=sampleID, col=myColor, max.iterations=200, 
						makePlots=makePlots, plot.path=celltype.path, algorithm='steep',
						geneUniverse=geneUniverse, verbose=verbose)
			if ( ! is.null(ans7)) {
				pcts7 <- ans7$CellProportions
				cellWhere <- match( cellTypeNames, names(pcts7))
				pcts7 <- pcts7[ cellWhere]
			}
		}

		# 8)  Nonlinear Least Squares of the 27-dimensional immune profile
		if ( is.null(cellAns) || (recalculate %in% c("all","profile","nls")) || (recalculate == "missing" && all(is.na(pcts8)))) {
			cat( "\n8. Cell Type Profile: Fit QN Dimension Proportions by Nonlinear Least Squares (NLS):\n")
			ans8 <- fitCellTypeProfileFromFile( f=transcriptFile, sid=sampleID, col=myColor, 
						makePlots=makePlots, plot.path=celltype.path, algorithm='nls',
						geneUniverse=geneUniverse, verbose=verbose)
			if ( ! is.null(ans8)) {
				pcts8 <- ans8$CellProportions
				cellWhere <- match( cellTypeNames, names(pcts8))
				pcts8 <- pcts8[ cellWhere]
			}
		}

		# 9)  GenSA of the 27-dimensional immune profile
		if ( is.null(cellAns) || (recalculate %in% c("all","profile","GenSA")) || (recalculate == "missing" && all(is.na(pcts9)))) {
			cat( "\n9. Cell Type Profile: Fit QN Dimension Proportions by Simulated Annealing (GenSA):\n")
			ans9 <- fitCellTypeProfileFromFile( f=transcriptFile, sid=sampleID, col=myColor, 
						makePlots=makePlots, plot.path=celltype.path, algorithm='GenSA',
						geneUniverse=geneUniverse, verbose=verbose)
			if ( ! is.null(ans9)) {
				pcts9 <- ans9$CellProportions
				cellWhere <- match( cellTypeNames, names(pcts9))
				pcts9 <- pcts9[ cellWhere]
			}
		}

		# 10)  Transcriptome Deconvolution, by NLS using the 'port' algorithm...
		if ( is.null(cellAns) || (recalculate %in% c("all","deconvolution","nls")) || (recalculate == "missing" && all(is.na(pcts10)))) {
			cat( "\n10. Cell Type Deconvolution:  Fit QN Gene", units, "by Nonlinear Least Squares (NLS):\n")
			ans10 <- fileSet.TranscriptDeconvolution( files=transcriptFile, fids=sampleID, algorithm="port",
						useLog=FALSE, plot=deconvPlot, plot.path=celltype.path, plot.col=cellTypeColors,
						geneUniverse=geneUniverse, verbose=verbose)
			if ( ! is.null(ans10)) {
				# deconvolution tool returns weights.  Trun them to percentages
				pcts10 <- ans10$BestFit * 100
				cellWhere <- match( cellTypeNames, rownames(pcts10))
				pcts10 <- pcts10[ cellWhere, 1]
			}
		}
		
		# 11)  Transcriptome Deconvolution, by NLS using the 'GenSA' simulated annealing algorithm...
		if ( is.null(cellAns) || (recalculate %in% c("all","deconvolution","GenSA")) || (recalculate == "missing" && all(is.na(pcts11)))) {
			cat( "\n11. Cell Type Deconvolution:  Fit QN Gene", units, "by Simulated Annealing (GenSA):\n")
			ans11 <- fileSet.TranscriptDeconvolution( files=transcriptFile, fids=sampleID, algorithm="GenSA",
						useLog=FALSE, plot=deconvPlot, plot.path=celltype.path, plot.col=cellTypeColors,
						geneUniverse=geneUniverse, verbose=verbose)
			if ( ! is.null(ans11)) {
				# deconvolution tool returns weights.  Trun them to percentages
				pcts11 <- ans11$BestFit * 100
				cellWhere <- match( cellTypeNames, rownames(pcts11))
				pcts11 <- pcts11[ cellWhere, 1]
			}
		}
		
		# 12)  Transcriptome Deconvolution, by Steepest Descent...
		if ( is.null(cellAns) || (recalculate %in% c("all","deconvolution","steep")) || (recalculate == "missing" && all(is.na(pcts12)))) {
			cat( "\n12. Cell Type Deconvolution:  Fit QN Gene", units, "by Steepest Descent:\n")
			ans12 <- fileSet.TranscriptDeconvolution( files=transcriptFile, fids=sampleID, algorithm="steep",
						useLog=FALSE, plot=deconvPlot, plot.path=celltype.path, plot.col=cellTypeColors,
						geneUniverse=geneUniverse, verbose=verbose)
			if ( ! is.null(ans12)) {
				# deconvolution tool returns weights.  Trun them to percentages
				pcts12 <- ans12$BestFit * 100
				cellWhere <- match( cellTypeNames, rownames(pcts12))
				pcts12 <- pcts12[ cellWhere, 1]
			}
		}
	}
	# restore the non-normalized cell type details after all QN steps are done
	CellTypeSetup( reload=TRUE)
	
	# done with all methods, now merge
	cellM <- matrix( NA, nrow=N_CellTypes, ncol=12)
	colnames(cellM) <- TestNames
	cellM[ ,1] <- pcts1
	cellM[ ,2] <- pcts2
	cellM[ ,3] <- pcts3
	cellM[ ,4] <- pcts4
	cellM[ ,5] <- pcts5
	cellM[ ,6] <- pcts6
	cellM[ ,7] <- pcts7
	cellM[ ,8] <- pcts8
	cellM[ ,9] <- pcts9
	cellM[ ,10] <- pcts10
	cellM[ ,11] <- pcts11
	cellM[ ,12] <- pcts12

	# catch any that are all zero
	isAllZero <- which( apply( cellM, 2, function(x) all( x == 0, na.rm=T)))
	if ( length(isAllZero)) cellM[ , isAllZero] <- NA

	# let's independently eavaluate the RMS and COD fit for every method
	evalAns1 <- rms.cod.evaluation( obsMatrix=sampleMatrix, calcMatrix=cellTypeMatrix.orig, calcPcts=cellM[ ,1:6], geneUniverse=geneUniverse)
	evalAns2 <- rms.cod.evaluation( obsMatrix=sampleMatrix, calcMatrix=cellTypeMatrix.qn, calcPcts=cellM[ ,7:12], geneUniverse=geneUniverse)
	rmsdAns <- c( evalAns1$rmsd, evalAns2$rmsd)
	codAns <- c( evalAns1$cod, evalAns2$cod)
	if (verbose) {
		cat( "\nGene expression RMS Deviation by method:\n")
		print( rmsdAns)
		cat( "Gene expression COD (Coefficient of Determination) by method:\n")
		print( codAns)
	}

	# final answer is either the average or the one method with lowest RMS
	# final answer will now be actual 'proportion', not 'percentage'.  Not summing to 100% anymore.
	mode <- match.arg( mode)
	if ( mode == "weighted.mean") {
		# since COD is always in 0..1, it is a better metric than RMS for doing the weights
		cod <- codAns
		cod[ cod <= 0] <- 0
		cod[ is.na(cod)] <- 0
		# Note:  there is a small chance that ALL CoD can be negative.  If so, revert to using normalized RMS deviation instead
		if ( all( cod == 0)) {
			cat( "\n  Warning:  all CoD are negative, reverting to using normalized RMSD instead.")
			nrms <- min(rmsdAns,na.rm=T) / rmsdAns
			nrms[ nrms <= 0] <- 0
			nrms[ is.na(nrms)] <- 0
			cod <- nrms
		}
		codRange <- range(cod[cod>0])
		if ( diff(codRange) > 0.05) {
			# if large enough spread, treat worst as zero and prorate all others back to original scale
			codTmp <- (cod - codRange[1])/diff(codRange)
			codTmp[ codTmp < 0] <- 0
			wts <- codTmp ^ 2
		} else {
			# if minimal spread, none get zeroed
			wts <- (cod / codRange[2]) ^ 2
		}
		wts <- round(wts, digits=4)
		names(wts) <- names(codAns)
		if (verbose) {
			cat( "\nFit Weights by 'Weighted Mean':\n")
			print( wts);
		}
		cellMean <- apply( cellM, 1, weighted.mean, w=wts, na.rm=T)
		cellMean <- round( cellMean, digits=5)
		if (verbose) cat( "\nFinal answer is weighted mean of all methods, using relative CoD squared as weights")
	} else if ( mode == "average") {
		wts <- rep.int( 1, 12)
		names(wts) <- names(codAns)
		if (verbose) {
			cat( "\nFit Weights by 'Average':\n")
			print( wts);
		}
		cellMean <- apply( cellM, 1, mean, na.rm=T)
		cellMean <- round( cellMean, digits=5)
		if (verbose) cat( "\nFinal answer is average of all methods")
	} else {
		best <- which.max(codAns)[1]
		wts <- rep.int( 0, 12)
		wts[ best] <- 1
		names(wts) <- names(codAns)
		if (verbose) {
			cat( "\nFit Weights by 'Best':\n")
			print( wts);
		}
		cellMean <- round( cellM[ , best], digits=5)
		if (verbose) cat( "\nFinal answer is method with best (highest) CoD: ", names(codAns)[best])
	}
		
	cellAns <- data.frame( "CellType"=cellTypeNames, "Final.Proportions"=cellMean, round(cellM,digits=3), stringsAsFactors=F)
	write.table( cellAns, celltypeDetailsFile, sep=",", quote=T, row.names=F)
	statsM <- matrix( c( rmsdAns, codAns, wts), nrow=12, ncol=3)
	rownames(statsM) <- names(rmsdAns)
	statsAns <- data.frame( "Statistic"=c("RMSD", "CoD", "Weight"), t(statsM), stringsAsFactors=F)
	write.table( statsAns, celltypeStatsFile, sep=",", quote=T, row.names=F)

	# make a plot image of the final average
	if ( makePlots != "none") {
		tmpM <- cbind( cellM, cellMean)
		rownames(tmpM) <- cellTypeNames
		colnames(tmpM) <- sub( "\\.", "\n", c( colnames(cellM), "Final.Proportions"))
		plotAns <- plotTranscriptProportions( tmpM, col=cellTypeColors, label=paste( sampleID, ":  Final average of all methods", sep=""), ...)
		text( c(-0.05,plotAns$bar.centers), rep.int(101.5,14), c("CoD=", as.character(round(codAns,dig=3)),""), cex=0.65, col=1, xpd=NA)
		text( c(-0.05,plotAns$bar.centers), rep.int(-1.5,14), c("Wt=", as.character(round(wts,dig=3)),""), cex=0.65, col=1, xpd=NA)
		dev.flush()
		printPlot( celltypePlotFile, width=12, height=9)
	}

	out <- data.frame( "SampleID"=sampleID, cellAns, stringsAsFactors=F)
	return( invisible( out))
}


`rms.cod.evaluation` <- function( obsMatrix, calcMatrix, calcPcts, geneUniverse=NULL, DEBUG=FALSE) {

	# measure the Root Mean Square (RMS) deviation and the CoD of each set of fitted cell type percentages
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
	rmsOut <- codOut <- rep.int( NA, Nmethods)
	names(rmsOut) <- names(codOut) <- colnames(calcPcts)
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
		resSq <- (obsValue - thisCalcValue) ^ 2
		rms <- sqrt( mean( resSq))
		rmsOut[j] <- round( rms, digits=3)
		SSres <- sum( resSq)
		obsSq <- (obsValue - mean(obsValue)) ^ 2
		SStot <- sum( obsSq)
		cod <- 1 - (SSres/SStot)
		codOut[j] <- round( cod, digits=3)
	}
	return( list( "rmsd"=rmsOut, "cod"=codOut))
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
	plotFile <- paste( "Compare.", prefix[1], ".", reference, "_", levelString, "_", toupper(plot.mode), ".Plot.", dev.type, sep="")
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

