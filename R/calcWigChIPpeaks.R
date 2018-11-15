# calcWigChIPpeaks.R

# turn a Wiggles object into a ChIP peaks list

`calcWigChIPpeaks` <- function( WIG, fileout, doPeakSearch=TRUE, peak.type="auto", cutoff.medians=3, 
				canonical.width=50, use.multicore=TRUE, controlPeaksFile=NULL,
				scaledReadCount=NULL, watermark.gene=NULL, 
				verbose=!interactive(), visualize=interactive() ) {

	setCurrentSpecies( WIG$Info$Species)
	seqMap <- getCurrentSeqMap()
	geneMap <- getCurrentGeneMap()
	sampleID <- WIG$Info$SampleID

	require(ROC)

	cat( "\n\nFinding ChIP Peaks: \t", sampleID)
	cat( "\nSpecies:              \t", getCurrentSpecies())
	cat( "\nAlignment files: \n  ", paste( WIG$Info$FileName, collapse="\n  "), "\n", sep="")

	# visit every chromosome
	allPeaks <- data.frame()

	if (doPeakSearch) {

	    totalReads <- WIG$Info$TotalReads
	    scaleFactor <- 1
	    if ( ! is.null( scaledReadCount)) scaleFactor <- scaledReadCount / totalReads

	    for( iseq in 1:nrow( seqMap)) {

		seqID <- seqMap$SEQ_ID[iseq]
		cat( "\n", seqID) 

		# get the wiggle tracks for this chromosome
		wiggleChunk <- WIG_getWigglesOneSeq( WIG, seqID)
		wigP <- wiggleChunk$Plus
		wigM <- wiggleChunk$Minus
		if ( ! ("Combo" %in% names(wiggleChunk))) {
			cat( "  Combining Plus & Minus strand..")
			wigC <- merge.baseDepthTables( wigP, wigM)

			# update the WIG data...
			cat( "  Resaving WIG data..")
			if ( ! is.null( wiggleChunk$PlusUnique)) wiggleChunk$PlusUnique <- NULL
			if ( ! is.null( wiggleChunk$MinusUnique)) wiggleChunk$MinusUnique <- NULL
			wiggleChunk$Combo <- wigC
			WIG <- WIG_updateWigglesOneSeq( WIG, seqID, wiggleChunk)
			cat( "  Done.\n")
		} else {
			wigC <- wiggleChunk$Combo
		}

		if ( ! is.null( watermark.gene)) {
			cat( "  Subtracting watermark..")
			wigP <- subtractWatermark( wigP, watermark.gene)
			wigM <- subtractWatermark( wigM, watermark.gene)
			wigC <- subtractWatermark( wigC, watermark.gene)
		}

		# make it look like what the peak picker wants
		cat( "  Vectorizing Wiggle Data..")
		vecP <- baseDepthTableToVector(wigP)
		wigP <- list( "x"=as.integer(names(vecP)), "y"=as.numeric( vecP) * scaleFactor, "strand"="Plus")
		rm(vecP)
		vecM <- baseDepthTableToVector(wigM)
		wigM <- list( "x"=as.integer(names(vecM)), "y"=as.numeric( vecM) * scaleFactor, "strand"="Minus")
		rm(vecM)
		vecC <- baseDepthTableToVector(wigC)
		wigC <- list( "x"=as.integer(names(vecC)), "y"=as.numeric( vecC) * scaleFactor, "strand"="Combo")
		rm(vecC)
		gc()


		if (visualize) {
			plotPath <- file.path( dirname( fileout), "ChIP.plots")
			if ( ! file.exists( plotPath)) dir.create( plotPath, recursive=T)
			plot.file <- file.path( plotPath, "rawPeak.png")
		} else {
			plot.file <- NULL
		}

		# prep file names
		rocfileP <- sub( "txt$", paste( seqID, "Plus.ROC.png", sep="."), fileout)
		rocfileM <- sub( "txt$", paste( seqID, "Minus.ROC.png", sep="."), fileout)
		rocfileC <- sub( "txt$", paste( seqID, "Combo.ROC.png", sep="."), fileout)
		pkfileP <- sub( "txt$", paste( seqID, "Plus.txt", sep="."), fileout)
		pkfileM <- sub( "txt$", paste( seqID, "Minus.txt", sep="."), fileout)
		pkfileC <- sub( "txt$", paste( seqID, "Combo.txt", sep="."), fileout)
		scorefile <- sub( "txt$", paste( seqID, "ScoreDetails.txt", sep="."), fileout)

		cat( "  Picking peaks..")
		if ( use.multicore) {
			wlist <- list( wigP, wigM, wigC)
			ans <- multicore.lapply( wlist, FUN=findPeaks, peak.type=peak.type, 
					cutoff.medians=cutoff.medians, canonical.width=canonical.width, 
					plot.file=plot.file, roc.file=fileout, visualize=visualize,
					controlPeaksFile=controlPeaksFile)
			ansP <<- ans[[1]]$peaks
			detailsP <<- ans[[1]]$scoreDetails
			cat( "   N_Plus: ", nrow(ansP), "\n")
			ansM <<- ans[[2]]$peaks
			detailsM <<- ans[[2]]$scoreDetails
			cat( "   N_Minus: ", nrow(ansM), "\n")
			ansC <<- ans[[3]]$peaks
			detailsC <<- ans[[3]]$scoreDetails
			cat( "   N_Combo: ", nrow(ansC), "\n")
		} else {
			ans <- findPeaks( wigP, peak.type=peak.type, cutoff.medians=cutoff.medians, 
						canonical.width=canonical.width, plot.file=plot.file,
						roc.file=rocfileP, visualize=visualize,
						controlPeaksFile=controlPeaksFile)
			ansP <<- ans$peaks
			detailsP <<- ans$scoreDetails
			cat( "   N_Plus: ", nrow(ansP), "\n")
			dev.print( png, rocfileP, width=700, height=720, bg='white')

			ans <- findPeaks( wigM, peak.type=peak.type, cutoff.medians=cutoff.medians, 
						canonical.width=canonical.width, plot.file=plot.file,
						roc.file=rocfileM, visualize=visualize,
						controlPeaksFile=controlPeaksFile)
			ansM <<- ans$peaks
			detailsM <<- ans$scoreDetails
			cat( "   N_Minus: ", nrow(ansM), "\n")
			dev.print( png, rocfileM, width=700, height=720, bg='white')

			ans <- findPeaks( wigC, peak.type=peak.type, cutoff.medians=cutoff.medians, 
						canonical.width=canonical.width, plot.file=plot.file,
						roc.file=rocfileC, visualize=visualize,
						controlPeaksFile=controlPeaksFile)
			ansC <<- ans$peaks
			detailsC <<- ans$scoreDetails
			cat( "   N_Combo: ", nrow(ansC), "\n")
			dev.print( png, rocfileC, width=700, height=720, bg='white')
		}
		write.table( ansP, pkfileP, sep="\t", quote=F, row.names=F, na="")
		write.table( ansM, pkfileM, sep="\t", quote=F, row.names=F, na="")
		write.table( ansC, pkfileC, sep="\t", quote=F, row.names=F, na="")

		scoreDetails <- rbind( detailsP, detailsM, detailsC)
		ord <- order( scoreDetails$loc)
		scoreDetails <- scoreDetails[ord, ]
		scoreDetails$loc <- round( scoreDetails$loc)
		write.table( scoreDetails, scorefile, sep="\t", quote=F, row.names=F)
	    }
	} else {
	    cat( "\nUsing already-picked peaks...")
	}

	# wrap it all together
	out <- summarizeWigChIPpeaks( WIG, fileout, canonical.width=canonical.width,
				controlPeaksFile=controlPeaksFile, scaledReadCount=scaledReadCount,
				watermark.gene=watermark.gene)

	return( out)
}


`summarizeWigChIPpeaks` <- function( WIG, fileout, canonical.width=50, controlPeaksFile=NULL,
					scaledReadCount=NULL, watermark.gene=NULL) {

	setCurrentSpecies( WIG$Info$Species)
	seqMap <- getCurrentSeqMap()
	geneMap <- getCurrentGeneMap()
	sampleID <- WIG$Info$SampleID
	totalReads <- WIG$Info$TotalReads

	# visit every chromosome
	allPeaks <- data.frame()

	for( iseq in 1:nrow( seqMap)) {
		seqID <- seqMap$SEQ_ID[iseq]
		pkfileP <- sub( "txt$", paste( seqID, "Plus.txt", sep="."), fileout)
		pkfileM <- sub( "txt$", paste( seqID, "Minus.txt", sep="."), fileout)
		pkfileC <- sub( "txt$", paste( seqID, "Combo.txt", sep="."), fileout)

		cat( "\n\nMerging 'ChIP' triplets...")
		ansP <- read.delim( pkfileP, as.is=T)
		ansM <- read.delim( pkfileM, as.is=T)
		ansC <- read.delim( pkfileC, as.is=T)
		ans <- mergeChIPpeaks( ansP, ansM, ansC)

		allPeaks <- rbind( allPeaks, ans)

		scorefile <- sub( "txt$", paste( seqID, "ScoreDetails.txt", sep="."), fileout)
		scoreDetails <- read.delim( scorefile, as.is=T)
	}

	# the combo peaks can be messy...  Re-fit using the forward and reverse peaks edges...
	ans <- refitComboChIPpeaks( WIG, allPeaks, canonical.width=canonical.width,
					controlPeaksFile=controlPeaksFile,
					scaledReadCount=scaledReadCount, watermark.gene=watermark.gene)
	allPeaks <- ans$peaks
	refitDetails <- ans$scoreDetails

	# with all 3 peaks ready, make a final score
	ans <- scoreThreePeaks( allPeaks, canonical.width=canonical.width,
					controlPeaksFile=controlPeaksFile)
	allPeaks <- ans$peaks
	tripletDetails <- ans$scoreDetails

	scoreDetails <- rbind( scoreDetails, refitDetails, tripletDetails)
	ord <- order( scoreDetails$loc)
	scoreDetails <- scoreDetails[ ord, ]
	scoreDetails$loc <- round( scoreDetails$loc)
	write.table( scoreDetails, scorefile, sep="\t", quote=F, row.names=F)

	out <- finalizeChIPpeaks( WIG, allPeaks, scaledReadCount=scaledReadCount)
	write.table( out, file=fileout, sep="\t", quote=FALSE, row.names=FALSE)
	cat( "\nWrote ChIP peaks file:  \t", fileout)
	cat( "\nN_Peaks:  ", nrow(out))
	
	return( out)
}


`mergeChIPpeaks` <- function( dfP, dfM, dfC) {

	# given peaks from Plus, Minus and Combined, find the real ones
	dfP$center <- as.integer( round( dfP$center))
	dfM$center <- as.integer( round( dfM$center))
	dfC$center <- as.integer( round( dfC$center))

	ord <- order( dfP$height, decreasing=T)
	dfP <- dfP[ ord, ]
	colnames(dfP) <- gsub( "\\.", "", colnames(dfP))
	NP <- nrow(dfP)
	cat( "\nN_PlusStrand:   ", NP)

	ord <- order( dfM$height, decreasing=T)
	dfM <- dfM[ ord, ]
	colnames(dfM) <- gsub( "\\.", "", colnames(dfM))
	NM <- nrow(dfM)
	cat( "\nN_MinusStrand:  ", NM)

	ord <- order( dfC$height, decreasing=T)
	dfC <- dfC[ ord, ]
	colnames(dfC) <- gsub( "\\.", "", colnames(dfC))
	# only keep the 'real' ones
	dfC <- subset( dfC, height > 0)
	NC <- nrow(dfC)
	cat( "\nN_Combo:        ", NC)

	if ( any( c(NP, NM, NC) == 0)) stop( "Failed to find enough good peaks")

	outC <- outP <- outM <- data.frame()
	
	# make a small 'empty' entry
	empty <- dfP[ 1, ]
	for (j in 1:ncol(empty)) empty[ ,j] <- NA
	empty$type <- ""
	empty$status <- "notDetected"
	empty$pvalue <- 1

	# keep a tally of who is still available
	availP <- 1:NP
	centerP <- dfP$center[ availP]
	highP <- dfP$height[ availP]
	widthP <- dfP$width[ availP]
	availM <- 1:NM
	centerM <- dfM$center[ availM]
	highM <- dfM$height[ availM]
	widthM <- dfM$width[ availM]

	# visit each 'combo' peak in desending order and find evidence for both kids
	for ( ic in 1:NC) {
		thisMean <- dfC$center[ic]
		thisSD <- dfC$width[ic]
		thisHigh <- dfC$height[ic]

		# since the parent peak has the worse chance of being the ideal shape, try
		# to base the child calls on the children themselves

		# the child peak must be at least a small amount away from the parent center
		# and no more than 1.5 widths away
		thisHalfLo <- thisMean - abs(thisSD*2.5)
		thisHalfSmallHi <- thisMean + abs(thisSD*0.5)
		thisHalfHi <- thisMean + abs(thisSD*2.5)
		thisHalfSmallLo <- thisMean - abs(thisSD*0.5)

		# find what kids fall under me, take the biggest (first)
		gotP <- gotM <- FALSE
		myCenP <- myCenM <- myWidP <- myWidM <- NA

		hitP <- which( centerP > thisHalfLo & centerP < thisHalfSmallHi)
		hitM <- which( centerM > thisHalfSmallLo & centerM < thisHalfHi)
		if ( length(hitP) > 0) {
			hitP <- hitP[1]
			myCenP <- centerP[ hitP]
			myWidP <- abs(widthP[ hitP])
			gotP <- TRUE
		}
		if ( length(hitM) > 0) {
			hitM <- hitM[1]
			myCenM <- centerM[ hitM]
			myWidM <- abs(widthM[ hitM])
			gotM <- TRUE
		}

		# verify these kids have the orientation and separation we expect
		# invalid results get an empty entry, otherwise keep and note 'used'
		# Change:  we will do this later via scoring, let them stay picked for now...
		#if ( gotP && gotM) {
		#	deltaCenter <- myCenM - myCenP
		#	avgWidth <- min( c( myWidP, myWidM))
		#	# too close or wrong direction
		#	if ( deltaCenter < (avgWidth/2)) gotM <- gotP <- FALSE
		#	# too far apart
		#	if ( deltaCenter > (avgWidth*6)) gotM <- gotP <- FALSE
		#}

		# if we passed, use those kids...
		if ( ! gotP) {
			outP <- rbind( outP, empty)
		} else {
			outP <- rbind( outP, dfP[ availP[ hitP], ])
			availP <- setdiff( availP, availP[ hitP])
			centerP <- dfP$center[ availP]
			highP <- dfP$height[ availP]
			widthP <- dfP$width[ availP]
		}
		if ( ! gotM) {
			outM <- rbind( outM, empty)
		} else {
			outM <- rbind( outM, dfM[ availM[ hitM], ])
			availM <- setdiff( availM, availM[ hitM])
			centerM <- dfM$center[ availM]
			highM <- dfM$height[ availM]
			widthM <- dfM$width[ availM]
		}

		# we always keep the combo
		outC <- rbind( outC, dfC[ic, ])
	}

	# combine what we want
	colnames(outP) <- paste( "F", colnames(outP), sep="")
	colnames(outM) <- paste( "R", colnames(outM), sep="")
	colnames(outC) <- paste( "C", colnames(outC), sep="")

	out <- cbind( outC, outP, outM, stringsAsFactors=F)
	rownames(out) <- 1:nrow(out)
	return( out)
}


`finalizeChIPpeaks` <- function( WIG, tbl, scaledReadCount=NULL) {

	# we have the final set of peaks...

	#turn volumes into normalized reads
	totalReads <- WIG$Info$TotalReads
	readLen <- WIG$Info$ReadLength
	scaleFac <- (1000000 / totalReads) / readLen

	# we need to 'undo' the effect of scaling all samples to get these
	# values to reflect what was really in the dataset.
	scaleScaleFactor <- 1
	if ( ! is.null( scaledReadCount)) scaleScaleFactor <- totalReads / scaledReadCount
	tbl$Fvpm <- tbl$Fvolume * scaleFac * scaleScaleFactor
	tbl$Rvpm <- tbl$Rvolume * scaleFac * scaleScaleFactor
	tbl$VPM <- tbl$Fvpm + tbl$Rvpm
	Fresid <- tbl$Fresidual * scaleFac * scaleScaleFactor
	Rresid <- tbl$Rresidual * scaleFac * scaleScaleFactor

	# add some measures of the overall result
	tbl$Model_Error <- as.percent( (Fresid + Rresid)/tbl$VPM)

	# select the columns we need to keep
	outColumns <- c( "VPM", "P_Value", "Score",
			"Cheight", "Ccenter", "Ctimes", 
			"Ftype", "Rtype", "Fpvalue", "Rpvalue", "Fscore", "Rscore", 
			"Fcenter", "Rcenter", "Fheight", "Rheight", "Fwidth", "Rwidth", 
			"Model_Error",  
			"Cpvalue", "Cscore", "Cwidth", "Cfloor", "Ctype", "Cstatus", "Cstart", "Cstop", 
			"Fvpm", "Ffloor", "Fstatus", "Fstart", "Fstop", "Ftimes",
			"Rvpm", "Rfloor", "Rstatus", "Rstart", "Rstop", "Rtimes")
	theseColumns <- match( outColumns, colnames(tbl), nomatch=0)
	otherColumns <- setdiff( 1:ncol(tbl), theseColumns)
	out <- tbl[ , c( theseColumns, otherColumns)]
	ord <- order( out$P_Value, -out$Score, -out$Cheight)
	out <- out[ ord, ]

	# assign the proximal gene for each
	cat( "\nAssigning proximal gene names...")
	geneName <- selectProximalGene( out$Ccenter)

	out <- data.frame( "GENE_ID"=geneName, out, stringsAsFactors=F)
	rownames(out) <- 1:nrow(out)

	return(out)
}


`selectProximalGene` <- function( centers) {

	N <- length(centers)
	outName <- rep( "", times=N)

	# get the genes, and their coding direction
	gmap <- subset( getCurrentGeneMap(), REAL_G == TRUE)
	NG <- nrow(gmap)
	gmapNames <- gmap$GENE_ID
	starts <- gmap$POSITION
	stops <- gmap$END
	strand <- gmap$STRAND
	isMinus <- which( strand == "-")
	starts[ isMinus] <- gmap$END[ isMinus]
	stops[ isMinus] <- gmap$POSITION[ isMinus]

	bestGenePtr <- findInterval( centers, gmap$POSITION)

	# visit each center and see exactly how it lays wrt coding genes
	for ( j in 1:N) {
		thisptr <- bestGenePtr[j]

		# special case:  before the first, and at the end
		if (thisptr < 1) {
			outName[j] <- gmapNames[1]
			next
		}
		if (thisptr >= NG) {
			outName[j] <- gmapNames[NG]
			next
		}

		# is it inside?
		if ( centers[j] >= gmap$POSITION[thisptr] & centers[j] <= gmap$END[thisptr]) {
			outName[j] <- gmapNames[thisptr]
			next
		}

		# we are between 2 genes
		signPrev <- strand[ thisptr]
		signNext <- strand[ thisptr + 1]
		if ( all( c( signPrev, signNext) == "+")) {
			outName[j] <- gmapNames[thisptr+1]
			next
		}
		if ( all( c( signPrev, signNext) == "-")) {
			outName[j] <- gmapNames[thisptr]
			next
		}
		if ( signPrev == "+" && signNext == "-") {
			dPrev <- centers[j] - stops[thisptr]
			dNext <- stops[ thisptr+1] - centers[j]
			useptr <- if ( dPrev < dNext) thisptr else (thisptr+1)
			outName[j] <- gmapNames[useptr]
			next
		}
		if ( signPrev == "-" && signNext == "+") {
			dPrev <- centers[j] - starts[thisptr]
			dNext <- starts[ thisptr+1] - centers[j]
			useptr <- if ( dPrev < dNext) thisptr else (thisptr+1)
			outName[j] <- gmapNames[useptr]
			next
		}
		if ( any( c( signPrev, signNext) == "")) {
			dPrev <- centers[j] - gmap$END[thisptr]
			dNext <- gmap$POSITION[ thisptr+1] - centers[j]
			useptr <- if ( dPrev < dNext) thisptr else (thisptr+1)
			outName[j] <- gmapNames[useptr]
			next
		}

		# there times that one may be blank, like non-coding RNAs
		# that's all choices
		cat( "\nInvalid choice! ", j, centers[j], thisptr, gmapNames[thisptr], 
				signPrev, signNext)
	}

	return( outName)
}


`makeControlPeaksFile` <- function( sampleSet, path, fileout="controlPeaks.csv") {

	if ( ! file.exists( path)) {
		cat( "\nPath to existing Control Peak results files not found: ", path)
		return(NULL)
	}

	prefix <- getCurrentSpeciesFilePrefix()
	smap <- getCurrentSeqMap()

	#sampName <- Fscore <- Rscore <- Cscore <- Sscore <- vector()
	#Fcenter <- Rcenter <- Ccenter <- Scenter <- vector()
	out <- data.frame()

	for ( s in sampleSet) {
		smlDF <- data.frame()
		for (seqID in smap$SEQ_ID) {
			smlF <- smlR <- smlC <- smlM <- data.frame()

			f <- paste( s, prefix, "ChIPpeaks", seqID, "Plus.txt", sep=".")
			f <- file.path( path, s, f)
			if ( ! file.exists( f)) {
				cat( "\nControl Peaks file not found: ", f)
			} else {
				tbl <- read.delim( f, as.is=T)
				tmpF <- tbl$score
				tmpFc <- tbl$center
				smlF <- data.frame( "Sample"=s, "Class"="Plus", "Center"=round(tmpFc), "Score"=tmpF,
						stringsAsFactors=F)
			}

			f <- paste( s, prefix, "ChIPpeaks", seqID, "Minus.txt", sep=".")
			f <- file.path( path, s, f)
			if ( ! file.exists( f)) {
				cat( "\nControl Peaks file not found: ", f)
			} else {
				tbl <- read.delim( f, as.is=T)
				tmpR <- tbl$score
				tmpRc <- tbl$center
				smlR <- data.frame( "Sample"=s, "Class"="Minus", "Center"=round(tmpRc), "Score"=tmpR,
						stringsAsFactors=F)
			}

			f <- paste( s, prefix, "ChIPpeaks", seqID, "Combo.txt", sep=".")
			f <- file.path( path, s, f)
			if ( ! file.exists( f)) {
				cat( "\nControl Peaks file not found: ", f)
			} else {
				tbl <- read.delim( f, as.is=T)
				tmpC <- tbl$score
				tmpCc <- tbl$center
				smlC <- data.frame( "Sample"=s, "Class"="Combo", "Center"=round(tmpCc), "Score"=tmpC,
						stringsAsFactors=F)
			}

			f <- paste( s, prefix, "ChIPpeaks", "txt", sep=".")
			f <- file.path( path, s, f)
			if ( ! file.exists( f)) {
				cat( "\nControl Peaks file not found: ", f)
			} else {
				tbl <- read.delim( f, as.is=T)
				tmpM <- tbl$Score
				tmpMc <- tbl$Ccenter
				smlM <- data.frame( "Sample"=s, "Class"="Final", "Center"=round(tmpMc), "Score"=tmpM,
						stringsAsFactors=F)
			}

			# take the 3 strands of peaks and the final meta-score and keep them all
			smlDF <- rbind( smlDF, smlF, smlR, smlC, smlM)
			cat( "\nUsing as control:  ", s)
		}
		ord <- order( smlDF$Center)
		smlDF <- smlDF[ ord, ]
		rownames( smlDF) <- 1:nrow(smlDF)

		out <- rbind( out, smlDF)
	}

	write.table( out, fileout, sep=",", quote=T, row.names=F, na="")
	cat( "\nWrote Control Peaks file: ", fileout, "\nN_Peaks: ", nrow(out))
	return( fileout)
}



`refitComboChIPpeaks` <- function( WIG, allPeaks, canonical.width=50, controlPeaksFile=NULL,
				scaledReadCount=NULL, watermark.gene=NULL) {

	setCurrentSpecies( WIG$Info$Species)
	seqMap <- getCurrentSeqMap()
	sampleID <- WIG$Info$SampleID
	totalReads <- WIG$Info$TotalReads
	scaleFactor <- 1
	if ( ! is.null( scaledReadCount)) scaleFactor <- scaledReadCount / totalReads

	cat( "\n\nRe-Fitting Combo Peaks...")

	# visit every chromosome
	allOut <- allPeaks

	for( iseq in 1:nrow( seqMap)) {

		seqID <- seqMap$SEQ_ID[iseq]
		# get the wiggle tracks for this chromosome
		wiggleChunk <- WIG_getWigglesOneSeq( WIG, seqID)
		wigC <- wiggleChunk$Combo

		if ( ! is.null( watermark.gene)) {
			cat( "  Subtracting watermark..")
			wigC <- subtractWatermark( wigC, watermark.gene)
		}
		# make it look like what the peak picker wants
		cat( "  Vectorizing Combo Peak Data..")
		vecC <- baseDepthTableToVector(wigC)
		allX <- as.integer(names(vecC))
		allY <- as.numeric(vecC) * scaleFactor
		rm( vecC)
		medianY <- median( allY[ allY > 0])

		for ( i in 1:nrow( allPeaks)) {

			# use the F and R edges to clip the C region
			mystart <- as.integer( min( allPeaks$Fstart[i], allPeaks$Rstart[i]))
			mystop <- as.integer( max( allPeaks$Fstop[i], allPeaks$Rstop[i]))
			if ( any( is.na( c( mystart, mystop)))) next
			mypts <- mystart : mystop
			who <- which( allX %in% mypts)
			if ( length( who) >= length( mypts)/2) {
				myXs <- allX[who]
				myYs <- allY[who]
			} else {
				# the X,Y data is too sparse...  fill it in...
				where <- findInterval( mypts, allX, all.inside=T)
				myYs <- allY[ where]
				myXs <- mypts
			}
			nPts <- diff( range( myXs)) + 1
			# extend the range with some fake floor
			newXlows <- myXs - nPts - 10
			newXhighs <- myXs + nPts + 10
			myNewXs <- c( newXlows, myXs, newXhighs)
			myfloorIn <- allPeaks$Ffloor[i] + allPeaks$Rfloor[i]
			myNewYs <- c( rep(myfloorIn, length(newXlows)), myYs, rep(myfloorIn,length(newXhighs)))
			
			#cat( "\n", i, allPeaks$Ccenter[i], mystart, mystop, nPts)

			ans <- findOneBestPeak( myNewXs, myNewYs, peak.type="gaussian", fit.floor=TRUE, visualize=FALSE)
			if ( is.null( ans)) next

			allOut$Ccenter[i] <- myCenter <- ans$center
			allOut$Cwidth[i] <- ans$width
			allOut$Cheight[i] <- myHeight <- ans$height
			# ignore the floor, as there wasn't any given to the fitter
			myFloor <- max( ans$floor, min( allPeaks$Ffloor[i], allPeaks$Rfloor[i]))
			allOut$Cfloor[i] <- myFloor
			allOut$Cvolume[i] <- ans$volume
			allOut$Cresidual[i] <- ans$residual
			allOut$Cdrift[i] <- myCenter - median(mypts)
			myStart <- min( ans$x)
			myStop <- max( ans$x)
			biggerTail <- max( (myCenter - myStart), (myStop - myCenter))
			allOut$Cstart[i] <- myCenter - biggerTail
			allOut$Cstop[i] <- myCenter + biggerTail
			allOut$Ctimes[i] <- (myFloor + myHeight) / medianY
			allOut$Ctype[i] <- ans$type
			if ( i %% 200 == 0) cat( "\n", i, allPeaks$Ctype[i], allPeaks$Cwidth[i], 
						allOut$Ctype[i], allOut$Cwidth[i])
		}

		# with new fits, we now need to re-score and re-Pvalue too..
		cat( "\nRe-scoreing...")
		tmp <- data.frame( "center"=allOut$Ccenter, "height"=allOut$Cheight, "width"=allOut$Cwidth, 
				"floor"=allOut$Cfloor, "drift"=allOut$Cdrift)

		scoreAns <- scorePeaks( tmp, canonical.width * 2, cutoff=medianY, strand="refitCombo")
		allOut$Cscore <- scoreAns$score

		if ( ! is.null( controlPeaksFile)) {
			cat( "\nRe-do P-values...")
			ctrlTbl <- read.csv( controlPeaksFile, as.is=T)
			ctrlTbl <- subset( ctrlTbl, Class == "Combo")
			ctrlScore <- ctrlTbl$Score
			ctrlScore <- ctrlScore[ ! is.na(ctrlScore)]
			Nctrl <- length(ctrlScore)
			useCut <- quantile( ctrlScore, 0.95, na.rm=T)

			# pvalue is based on percent of randoms that were better
			pval <- sapply( 1:nrow(allOut), function(ip) {
					n <- sum( ctrlScore >= allOut$Cscore[ip])
					return( n / Nctrl)
				})
			allOut$Cpvalue <- pval
			allOut$Cstatus[ allOut$Cscore >= useCut] <- "good"
		}
	}

	return( list( "peaks"=allOut, "scoreDetails"=scoreAns))
}


`scoreThreePeaks` <- function( allPeaks, canonical.width=50, controlPeaksFile=NULL) {

	# we have the 3 separate peaks with their scores and P-values

	# now try to incorporate any measurements of how they should be relative to each other
	N <- nrow(allPeaks)
	outScore <- rep( 0, times=N)
	ht <- wd <- dr <- rel <- rep( NA, times=N)

	prevScore <- prevPval <- rep( NA, times=N)

	for ( i in 1:N) {
		avgScore <- mean( c( allPeaks$Fscore[i], allPeaks$Rscore[i], allPeaks$Cscore[i]))
		prevScore[i] <- avgScore
		avgPval <- p.combine( c( allPeaks$Fpvalue[i], allPeaks$Rpvalue[i], allPeaks$Cpvalue[i]))
		prevPval[i] <- avgPval
		if ( is.na(avgScore)) next

		# the kids should be separate...  ideal is about 2 widths... but be a bit laxer...
		sep <- allPeaks$Rcenter[i] - allPeaks$Fcenter[i]
		if ( sep < 1) sep <- 1
		wTerm <- 1.0
		lowerWid <- canonical.width*0.6
		upperWid <- (canonical.width*3)
		if ( sep < lowerWid) wTerm <- (sep/lowerWid)
		if ( sep > upperWid) wTerm <- (upperWid / sep)
		wd[i] <- wTerm

		# the kid peaks should be about the same height
		hF <- allPeaks$Fheight[i]
		hR <- allPeaks$Rheight[i]
		hC <- allPeaks$Cheight[i]
		if ( any( is.na( c(hF,hR,hC)))) next
		if ( any( c(hF,hR,hC) < 1)) next
		hMean <- mean( c( hF, hR))
		relTerm <- sqrt( min( c( hF, hR) / hMean))
		if ( hC < hMean) relTerm <- relTerm * sqrt(hC/hMean)
		rel[i] <- relTerm

		# that's all for now...
		outScore[i] <- avgScore * wTerm * relTerm
	}

	allOut <- allPeaks
	allOut$Score <- outScore

	# bundle up the score facts to look like the other scoring tool...
	out2 <- data.frame( "loc"=allPeaks$Ccenter, "strand"=rep("3-tuple",N), "ht"=ht, "wd"=wd,
				"dr"=dr, "rel"=rel, "score"=outScore, stringsAsFactors=F)

	if ( ! is.null( controlPeaksFile)) {
		cat( "\nRe-do P-values...")
		ctrlTbl <- read.csv( controlPeaksFile, as.is=T)
		ctrlTbl <- subset( ctrlTbl, Class == "Final")
		ctrlScore <- ctrlTbl$Score
		ctrlScore <- ctrlScore[ ! is.na(ctrlScore)]
		Nctrl <- length(ctrlScore)
		useCut <- quantile( ctrlScore, 0.95, na.rm=T)

		# pvalue is based on percent of randoms that were better
		pval <- sapply( 1:nrow(allOut), function(ip) {
				n <- sum( ctrlScore >= allOut$Score[ip])
				return( n / Nctrl)
			})
		allOut$P_Value <- pval
	} else {

		# we don't have the random peaks to use...
		ord <- order( prevScore)
		prevScore <- prevScore[ ord]
		prevPval <- prevPval[ ord]
		if ( length( drops <- which( is.na( prevScore))) > 0) {
			prevScore <- prevScore[ -drops]
			prevPval <- prevPval[ -drops]
		}
		where <- findInterval( allOut$Score, prevScore, all.inside=TRUE)
		myP <- prevPval[ where]
		#myP <- mapply( allPeaks$Fpvalue, allPeaks$Rpvalue, allPeaks$Cpvalue, FUN=max)
		#myP[ is.na(myP)] <- 1.0
		allOut$P_Value <- myP
	}

	return( list( "peaks"=allOut, "scoreDetails"=out2))
}


`subtractWatermark` <- function( wigStrand, geneID) {

	gmap <- getCurrentGeneMap()
	gptr <- match( geneID, gmap$GENE_ID, nomatch=0)
	if ( gptr == 0) {
		cat( "\nNo watermark subtracted...  Gene not found: ", geneID)
		return( wigStrand)
	}

	gBeg <- gmap$POSITION[ gptr]
	gEnd <- gmap$END[ gptr]

	x <- wigStrand$START
	y <- as.numeric( wigStrand$DEPTH)
	newY <- y
	overallMedian <- median( y, na.rm=T)

	waterBeg <- findInterval( gBeg, x, all.inside=T)
	waterEnd <- findInterval( gEnd, x, all.inside=T)
	waterY <-y[ waterBeg:waterEnd]
	waterMedian <- median( waterY, na.rm=T) 

	remove <- waterMedian - overallMedian
	waterY <- waterY - remove
	myFloor <- overallMedian
	waterY <- ifelse( waterY < myFloor, myFloor, waterY)
	newY[ waterBeg:waterEnd] <- waterY

	# put back this slightly lowered Depth data
	wigStrand$DEPTH <- newY
	return( wigStrand)
}
