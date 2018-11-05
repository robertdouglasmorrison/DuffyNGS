# calcWigRIPpeaks.R

# turn a Wiggles object into a RIP peaks list

`calcWigRIPpeaks` <- function( WIG, fileout, doPeakSearch=TRUE, peak.type="auto", cutoff.medians=3, 
				min.width=50, max.width=1000, use.multicore=FALSE, controlPeaksFile=NULL,
				scaledReadCount=NULL, verbose=!interactive(), visualize=interactive() ) {

	setCurrentSpecies( WIG$Info$Species)
	seqMap <- getCurrentSeqMap()
	geneMap <- getCurrentGeneMap()
	sampleID <- WIG$Info$SampleID

	require(ROC)

	cat( "\n\nFinding RIP Peaks: \t", sampleID)
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

		# make it look like what the peak picker wants
		cat( "  Vectorizing Wiggle Data..")
		vecP <- baseDepthTableToVector(wigP)
		wigP <- list( "x"=as.integer(names(vecP)), "y"=as.numeric( vecP) * scaleFactor, "strand"="Plus")
		rm(vecP)
		vecM <- baseDepthTableToVector(wigM)
		wigM <- list( "x"=as.integer(names(vecM)), "y"=as.numeric( vecM) * scaleFactor, "strand"="Minus")
		rm(vecM)
		gc()

		if (visualize) {
			plotPath <- file.path( dirname( fileout), "RIP.plots")
			if ( ! file.exists( plotPath)) dir.create( plotPath, recursive=T)
			plot.file <- file.path( plotPath, "rawPeak.png")
		} else {
			plot.file <- NULL
		}

		# prep file names
		rocfileP <- sub( "txt$", paste( seqID, "Plus.ROC.png", sep="."), fileout)
		rocfileM <- sub( "txt$", paste( seqID, "Minus.ROC.png", sep="."), fileout)
		pkfileP <- sub( "txt$", paste( seqID, "Plus.txt", sep="."), fileout)
		pkfileM <- sub( "txt$", paste( seqID, "Minus.txt", sep="."), fileout)
		scorefile <- sub( "txt$", paste( seqID, "ScoreDetails.txt", sep="."), fileout)

		cat( "  Picking peaks..")
		if ( use.multicore) {
			wlist <- list( wigP, wigM)
			ans <- multicore.lapply( wlist, FUN=findRIPpeaks, peak.type=peak.type, 
					cutoff.medians=cutoff.medians, min.width=min.width, max.width=max.width, 
					plot.file=plot.file, roc.file=fileout, visualize=visualize,
					controlPeaksFile=controlPeaksFile)
			ansP <<- ans[[1]]$peaks
			detailsP <<- ans[[1]]$scoreDetails
			cat( "   N_Plus: ", nrow(ansP), "\n")
			ansM <<- ans[[2]]$peaks
			detailsM <<- ans[[2]]$scoreDetails
			cat( "   N_Minus: ", nrow(ansM), "\n")
		} else {
			ans <- findRIPpeaks( wigP, peak.type=peak.type, cutoff.medians=cutoff.medians, 
						min.width=min.width, max.width=max.width, plot.file=plot.file,
						roc.file=rocfileP, visualize=visualize,
						controlPeaksFile=controlPeaksFile)
			ansP <<- ans$peaks
			detailsP <<- ans$scoreDetails
			cat( "   N_Plus: ", nrow(ansP), "\n")
			dev.print( png, rocfileP, width=700, height=720, bg='white')

			ans <- findRIPpeaks( wigM, peak.type=peak.type, cutoff.medians=cutoff.medians, 
						min.width=min.width, max.width=max.width, plot.file=plot.file,
						roc.file=rocfileM, visualize=visualize,
						controlPeaksFile=controlPeaksFile)
			ansM <<- ans$peaks
			detailsM <<- ans$scoreDetails
			cat( "   N_Minus: ", nrow(ansM), "\n")
			dev.print( png, rocfileM, width=700, height=720, bg='white')
		}
		write.table( ansP, pkfileP, sep="\t", quote=F, row.names=F, na="")
		write.table( ansM, pkfileM, sep="\t", quote=F, row.names=F, na="")

		scoreDetails <- rbind( detailsP, detailsM)
		ord <- order( scoreDetails$loc)
		scoreDetails <- scoreDetails[ord, ]
		scoreDetails$loc <- round( scoreDetails$loc)
		write.table( scoreDetails, scorefile, sep="\t", quote=F, row.names=F)
	    }
	} else {
	    cat( "\nUsing already-picked peaks...")
	}

	# wrap it all together
	out <- summarizeWigRIPpeaks( WIG, fileout, min.width=max.width, max.width=max.width, 
				controlPeaksFile=controlPeaksFile, scaledReadCount=scaledReadCount)

	return( out)
}


`summarizeWigRIPpeaks` <- function( WIG, fileout, min.width=50, max.width=1000, controlPeaksFile=NULL,
					scaledReadCount=NULL) {

	setCurrentSpecies( WIG$Info$Species)
	seqMap <- getCurrentSeqMap()
	geneMap <- getCurrentGeneMap()
	sampleID <- WIG$Info$SampleID
	totalReads <- WIG$Info$TotalReads

	# visit every chromosome
	allPeaks <- allDetails <- data.frame()

	for( iseq in 1:nrow( seqMap)) {
		seqID <- seqMap$SEQ_ID[iseq]
		pkfileP <- sub( "txt$", paste( seqID, "Plus.txt", sep="."), fileout)
		pkfileM <- sub( "txt$", paste( seqID, "Minus.txt", sep="."), fileout)

		cat( "\n\nMerging 'RIP' peaks from 'Plus' and 'Minus' strands..")
		ansP <- read.delim( pkfileP, as.is=T)
		ansM <- read.delim( pkfileM, as.is=T)
		ans <- mergeRIPpeaks( ansP, ansM)

		allPeaks <- rbind( allPeaks, ans)

		scorefile <- sub( "txt$", paste( seqID, "ScoreDetails.txt", sep="."), fileout)
		scoreDetails <- read.delim( scorefile, as.is=T)
		allDetails <- rbind( allDetails, scoreDetails)
	}

	ord <- order( allDetails$loc)
	allDetails <- allDetails[ ord, ]
	allDetails$loc <- round( allDetails$loc)
	write.table( allDetails, scorefile, sep="\t", quote=F, row.names=F)

	out <- finalizeRIPpeaks( WIG, allPeaks, scaledReadCount=scaledReadCount)
	write.table( out, file=fileout, sep="\t", quote=FALSE, row.names=FALSE)
	cat( "\nWrote RIP peaks file:  \t", fileout)
	cat( "\nN_Peaks:  ", nrow(out))
	
	return( out)
}


`mergeRIPpeaks` <- function( dfP, dfM) {

	# given RIP peaks from Plus, Minus, join into a final set
	# there is no need to "pair pos & neg" as in ChIP data.  Just keep both and note Strand
	NP <- nrow(dfP)
	cat( "\nN_PlusStrand:   ", NP)
	dfP$center <- as.integer( round( dfP$center))
	outP <- dfP
	outP$strand <- "Plus"

	NM <- nrow(dfM)
	cat( "\nN_MinusStrand:  ", NM)
	dfM$center <- as.integer( round( dfM$center))
	outM <- dfM
	outM$strand <- "Minus"

	out <- rbind( outP, outM)
	ord <- order( out$score, decreasing=T)
	out <- out[ ord, ]
	if ( nrow(out)) rownames(out) <- 1:nrow(out)
	return( out)
}


`finalizeRIPpeaks` <- function( WIG, tbl, scaledReadCount=NULL) {

	# we have the final set of peaks...

	#turn volumes into normalized reads
	totalReads <- WIG$Info$TotalReads
	readLen <- WIG$Info$ReadLength
	scaleFac <- (1000000 / totalReads) / readLen

	# we need to 'undo' the effect of scaling all samples to get these
	# values to reflect what was really in the dataset.
	scaleScaleFactor <- 1
	if ( ! is.null( scaledReadCount)) scaleScaleFactor <- totalReads / scaledReadCount
	tbl$VPM <- tbl$volume * scaleFac * scaleScaleFactor
	myResid <- tbl$residual * scaleFac * scaleScaleFactor

	# add some measures of the overall result
	tbl$Model_Error <- as.percent( myResid / tbl$VPM)

	# select the columns we need to keep
	# fix some column names on the fly
	colnames(tbl)[ colnames(tbl) == "p.value"] <- "P_Value"
	colnames(tbl)[ colnames(tbl) == "score"] <- "Score"
	colnames(tbl)[ colnames(tbl) == "height"] <- "Height"
	colnames(tbl)[ colnames(tbl) == "center"] <- "Center"
	colnames(tbl)[ colnames(tbl) == "type"] <- "Type"
	colnames(tbl)[ colnames(tbl) == "strand"] <- "Strand"
	colnames(tbl)[ colnames(tbl) == "start"] <- "Start"
	colnames(tbl)[ colnames(tbl) == "stop"] <- "Stop"
	colnames(tbl)[ colnames(tbl) == "width"] <- "Width"
	colnames(tbl)[ colnames(tbl) == "times"] <- "Times"
	colnames(tbl)[ colnames(tbl) == "floor"] <- "Floor"
	outColumns <- c( "VPM", "P_Value", "Score", "Height", "Center", "Type", "Strand",
			"Start", "Stop", "Width", "Times", "Floor", "Model_Error")
	theseColumns <- match( outColumns, colnames(tbl), nomatch=0)
	otherColumns <- setdiff( 1:ncol(tbl), theseColumns)
	out <- tbl[ , c( theseColumns, otherColumns)]
	ord <- order( out$P_Value, -out$Score, -out$VPM, -out$Height)
	out <- out[ ord, ]

	# assign the proximal gene for each
	cat( "\nAssigning proximal gene names...")
	geneName <- selectProximalGene( out$Center)

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
	SAV <<- bestGenePtr

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

			f <- paste( s, prefix, "RIPpeaks", seqID, "Plus.txt", sep=".")
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

			f <- paste( s, prefix, "RIPpeaks", seqID, "Minus.txt", sep=".")
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

			f <- paste( s, prefix, "RIPpeaks", seqID, "Combo.txt", sep=".")
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

			f <- paste( s, prefix, "RIPpeaks", "txt", sep=".")
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

