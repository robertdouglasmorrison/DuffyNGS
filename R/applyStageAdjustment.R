# applyStageAdjustent.R


# take the set of all transcripts and perform a stage adjustment to put them on similar stages
`pipe.StageAdjustment` <- function(  sampleIDset=NULL, newResultsPath=NULL, annotationFile="Annotation.txt", 
				optionsFile="Options.txt", rate=2, 
				tolerance=1, max.iterations=300, label="", legend.cex=1,
				scaleFile="./stageAdjustments.txt", verbose=TRUE) {

	optionTable <- readOptionsTable( optionsFile)
	oldResultsPath <- getOptionValue( optionTable, "results.path", notfound=".")
	cat( "\nReading transcripts from folder:   ", oldResultsPath)

	if ( is.null( newResultsPath)) {
		newResultsPath <- file.path( dirname(oldResultsPath), 
					paste( "adjusted", basename(oldResultsPath), sep="_"))
	}
	cat( "\nWriting adusted results to folder: ", newResultsPath,"\n")
	if ( ! file.exists( newResultsPath)) dir.create( newResultsPath, recursive=TRUE, showWarnings=FALSE)

	annT <- readAnnotationTable( annotationFile)
	myColors <- annT$Color
	mySamples <- annT$SampleID
	if ( ! is.null( sampleIDset)) {
		keep <- which( mySamples %in% sampleIDset)
		mySamples <- mySamples[keep]
		myColors <- myColors[keep]
	}

	mySpecies <- getCurrentSpecies()
	myPrefix <- getCurrentSpeciesFilePrefix()

	myPFfiles <- paste( mySamples, myPrefix, "Transcript.txt", sep=".")
	myPFfiles <- file.path( oldResultsPath, "transcript", myPFfiles)
	for ( thisfile in myPFfiles) {
		if ( ! file.exists( thisfile)) stop( paste( "Transcript file not found: ", thisfile))
	}

	# set up for plotting and adjusting the transcripts
	fileSet <- myPFfiles
	fidSet <- mySamples
	fcolors <- myColors

	verifyLifeCycleSetup()

	savemai <- par( "mai")
	par( "mai"=c( 1.4, 0.9, 0.9, 0.4))

	plotLifeCycleStageFromFileSet( fileSet, fidSet, fcolors=fcolors, intensityColumn="RPKM_M",
		label=paste( "Pre-Adjustment:  ", label), legend.cex=legend.cex)
	dev.print( png, file.path( ".", "stagePlot_Original.png"), width=900, height=600)
	
	ans <- adjustLifeCycleStageFileSet( fileSet, fidSet, fcolors=fcolors, intensityColumn="RPKM_M",
		label=paste( "Post-Adjustment:  ", label), rate=rate, tolerance=tolerance, max.iterations=max.iterations,
		legend.cex=legend.cex)
	dev.print( png, file.path( ".", "stagePlot_Adjusted.png"), width=900, height=600)

	par( "mai"=savemai)

	scaleFactors <- ans$scaleFactors
	scaleDF <- data.frame( "GENE_ID"=rownames( scaleFactors), scaleFactors)
	cat( "\nWriting Stage Adjustment scale factors to file: ", scaleFile)
	write.table( scaleDF, file=scaleFile, sep="\t", quote=FALSE, row.names=TRUE)

	whoDone <- applyStageAdjustment( scaleFile, sampleIDs=fidSet, newResultsPath=newResultsPath, 
			annotationFile=annotationFile, optionsFile=optionsFile, verbose=verbose)
	if ( length(whoDone) < length( fileSet)) stop( "Not all datasets successfully adjusted...")

	cat("\n\nStage Adjustment Done.\n")
	return()
}


# take a gene scale factor matrix from doing a life cycle stage adjustment,
# and apply it to all affected WIG datasets

`applyStageAdjustment` <- function( scaleFile, sampleIDs=colnames(scaleFile), newResultsPath=".", 
				annotationFile="Annotation.txt", optionsFile="Options.txt", verbose=TRUE) {

	if (verbose) cat( "\n\nApplying Stage Adjustments ",
		"\n  Scaling file:               ", scaleFile, 
		"\n  Destination Results folder: ", newResultsPath, "\n")

	optT <- readOptionsTable( optionsFile)
	oldResultsPath <- getOptionValue( optT, "results.path", notfound=".", verbose=verbose)
	annT <- readAnnotationTable( annotationFile)

	f <- scaleFile
	if ( ! file.exists( f)) {
		warning( "\nStage Scaling file not found:  ", f, "\n")
		return( vector())
	}

	# make a parallel folder for results
	wigPath <- file.path( newResultsPath, "wig")
	if ( ! file.exists( wigPath)) dir.create( wigPath, showWarnings=FALSE)
	transPath <- file.path( newResultsPath, "transcript")
	if ( ! file.exists( transPath)) dir.create( transPath, showWarnings=FALSE)

	# get the matrix of gene scale factors for each sample
	# note that the sampleIDs may not match the colnames exactly...
	scaling <- read.delim( f, as.is=TRUE)
	if ( "GENE_ID" %in% colnames(scaling)) {
		mySamples <- setdiff( colnames( scaling), "GENE_ID")
		myGenes <- scaling$GENE_ID
		mySampleIDs <- setdiff( sampleIDs, "GENE_ID")
		colOffset <- 1
	} else {
		mySamples <- colnames( scaling)
		myGenes <- rownames( scaling)
		mySampleIDs <- sampleIDs
		colOffset <- 0
	}

	# we may have been given a subset to adjust...
	if (length(mySampleIDs) < length(mySamples)) {
		# make the SIDs look like they had the colname substitutions
		tmpSIDs <- gsub( "-", ".", mySampleIDs, fixed=T)
		tmpSIDs <- gsub( "+", ".", tmpSIDs, fixed=T)
		tmpSIDs <- gsub( "%", ".", tmpSIDs, fixed=T)
		tmpSIDs <- gsub( "(", ".", tmpSIDs, fixed=T)
		tmpSIDs <- gsub( ")", ".", tmpSIDs, fixed=T)
		tmpSIDs <- sub( "(^[0-9])", "X\\1", tmpSIDs, fixed=F)
		whereScale <- match( tmpSIDs, mySamples, nomatch=0)
		if (any( whereScale == 0)) {
			cat( "\nUnable to match some SamplesIDs to scaling columns:")
			who <- which( whereScale == 0)
			cat( "\nSampleIDs:   ", mySampleIDs[who])
			cat( "\nLooked for:  ", tmpSIDs[who])
			stop( "Fixed in 'applyStageAdjustment()'")
		}
		scaling <- scaling[ , whereScale+colOffset, drop=FALSE]
		mySamples <- mySamples[ whereScale]
		mySampleIDs <- mySampleIDs[ whereScale > 0]
	}

	# order into chromosomal order for faster scaling
	myPrefix <- getCurrentSpeciesFilePrefix()
	gmap <- getCurrentGeneMap()
	gOrder <- base::match( myGenes, gmap$GENE_ID, nomatch=NA)
	ord <- base::order( gOrder)
	scaling <- scaling[ ord, , drop=FALSE]
	myGenes <- myGenes[ ord]


	`stageAdjustOneSample` <- function( iSample, verbose=TRUE) {

		# get that WIG object
		thisColumn <- mySamples[ iSample]
		thisSampleID <- mySampleIDs[ iSample]
		wigFile <- paste( thisSampleID, myPrefix, "WIG.rda", sep=".")
		wigFullFile <- file.path( oldResultsPath, "wig", wigFile)
		if ( ! file.exists( wigFullFile)) {
			cat( "WIG_file not found... ", wigFullFile, "\nSkipping: ", thisSampleID)
			return("")
		}
		load( wigFullFile)

		# get that column of scale factors
		stageColumn <- base::match( thisColumn, colnames( scaling))
		myScaleFactors <- scaling[ , stageColumn]

		# apply the scaling
		if (verbose) cat( "\nScaling:  ", thisSampleID)
		wiggles <- stageAdjustOneWIG( wiggles, myGenes, myScaleFactors,
				newPath=newResultsPath, oldPath=oldResultsPath, verbose=verbose)
		wigOutFile <- file.path( newResultsPath, "wig", wigFile)
		save( wiggles, file=wigOutFile)

		# remake the transcript..
		if (verbose) cat( "\nRemake Transcript:  ", thisSampleID, "\n")
		sampleID <- thisSampleID
		readSense <- getReadSense( sampleID, annotationFile)
		originalID <- originalSamplePairID( sampleID, annotationFile)
		useBothStrands <- ! getAnnotationTrue( annT, originalID, "StrandSpecific", notfound=TRUE)
		keepIntergenics <- getAnnotationTrue( annT, originalID, "KeepIntergenics", notfound=FALSE)
		fileOutTrans <- paste( sampleID, myPrefix, "Transcript.txt", sep=".")
		fileOutTrans <- file.path( newResultsPath, "transcript", fileOutTrans)
		trans <- calcWigTranscriptome( wiggles, useBothStrands=useBothStrands, 
					keepIntergenics=keepIntergenics, fileout=fileOutTrans)
		return( thisSampleID)
	}

	cat( "\nAdjusting all samples..\n")
	adjustedSamples <- multicore.lapply( 1:length(mySamples), FUN=stageAdjustOneSample, verbose=verbose)
	cat("\nDone.\n")
	return( adjustedSamples)
}


`stageAdjustOneWIG` <- function( wiggles, genes, scales, newPath, oldPath, verbose=TRUE) {

	gmap <- getCurrentGeneMap()
	curSeq <- ""
	curWigChunk <- NULL

	# make a new copy, and adjust all the file path names
	newWIG <- wiggles
	oldFiles <- wiggles$Info$FileName
	newWIG$Info$FileName <- sub( oldPath, newPath, oldFiles, fixed=T)
	oldSubFolder <- wiggles$SubWigFolder
	newWIG$SubWigFolder <- sub( oldPath, newPath, oldSubFolder, fixed=T)

	# there is a chance that some Seqs have no genes to adjust, so force
	# move all wiggle chunks to make sure...
	smap <- getCurrentSeqMap()
	for ( seqid in smap$SEQ_ID) {
		curWigChunk <- WIG_getWigglesOneSeq( wiggles, seqid, verbose=F)
		newWIG <- WIG_updateWigglesOneSeq( newWIG, seqid, curWigChunk, errorNotFound=FALSE)
	}
	curWigChunk <- NULL

	# we will change the read count of some genes, but we need to end up with the 
	# unchanged genes still having the same RPKM when we're done.
	oldReadCnt <- newWIG$Info$TotalReads
	oldUniqueReadCnt <- newWIG$Info$UniqueReads
	netDeltaUniqueCnts <- netDeltaCnts <- 0

	# adjust the wiggle tracks for each gene...
	for( i in 1:length(genes)) {
		if ( is.na(scales[i])) next
		if ( scales[i] == 1) next
		if ( scales[i] == 0) next
		gptr <- base::match( genes[i], gmap$GENE_ID, nomatch=0)
		if ( gptr == 0) next
		thisSeq <- gmap$SEQ_ID[gptr]

		# get that seq's wiggle tracks
		if ( thisSeq != curSeq) {
			# save last seq's new data
			if ( ! is.null( curWigChunk)) {
				curWigChunk$Plus <- thisPlus
				curWigChunk$Minus <- thisMinus
				curWigChunk$PlusUnique <- thisPlusUnique
				curWigChunk$MinusUnique <- thisMinusUnique
				newWIG <- WIG_updateWigglesOneSeq( newWIG, curSeq, curWigChunk)
				if (verbose) cat( "  ", curSeq)
				curWigChunk <- NULL
			}
			curWigChunk <- WIG_getWigglesOneSeq( wiggles, thisSeq)
			if ( is.null( curWigChunk)) next
			thisPlus <- curWigChunk$Plus
			thisMinus <- curWigChunk$Minus
			thisPlusUnique <- curWigChunk$PlusUnique
			thisMinusUnique <- curWigChunk$MinusUnique
			curSeq <- thisSeq
		}

		# get the gene edges for this gene
		start <- gmap$POSITION[gptr]
		stop <- gmap$END[gptr]

		# the bsse depth tables are not fixed, so we need to find this region in each
		if ( ! is.null( thisPlus)) {
			midpts <- (thisPlus$START + thisPlus$STOP) /2
			who <- which( midpts >= start & midpts <= stop)
			if( length(who) > 0) {
				thisPlus$DEPTH[who] <- thisPlus$DEPTH[who] * scales[i]
			}
		}
		if ( ! is.null( thisMinus)) {
			midpts <- (thisMinus$START + thisMinus$STOP) /2
			who <- which( midpts >= start & midpts <= stop)
			if( length(who) > 0) {
				thisMinus$DEPTH[who] <- thisMinus$DEPTH[who] * scales[i]
			}
		}
		if ( ! is.null( thisPlusUnique)) {
			midpts <- (thisPlusUnique$START + thisPlusUnique$STOP) /2
			who <- which( midpts >= start & midpts <= stop)
			if( length(who) > 0) {
				thisPlusUnique$DEPTH[who] <- thisPlusUnique$DEPTH[who] * scales[i]
			}
		}
		if ( ! is.null( thisMinusUnique)) {
			midpts <- (thisMinusUnique$START + thisMinusUnique$STOP) /2
			who <- which( midpts >= start & midpts <= stop)
			if( length(who) > 0) {
				thisMinusUnique$DEPTH[who] <- thisMinusUnique$DEPTH[who] * scales[i]
			}
		}
	}

	# store the last set back...
	if ( ! is.null( curWigChunk)) {
		curWigChunk$Plus <- thisPlus
		curWigChunk$Minus <- thisMinus
		curWigChunk$PlusUnique <- thisPlusUnique
		curWigChunk$MinusUnique <- thisMinusUnique
		newWIG <- WIG_updateWigglesOneSeq( newWIG, curSeq, curWigChunk)
		curWigChunk <- NULL
	}

	return( newWIG)
}

