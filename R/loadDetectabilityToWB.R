# loadDetectabilityToWB.R

# load up the WiggleBin structure that captures short read uniqueness

`loadDetectabilityToWB` <- function( filein, WB=NULL, speciesID=getCurrentSpecies(), 
					type=c( "selfUnique", "otherDetectable"), 
					otherSpeciesID="", readSize=32) {

	cat( "\nReading in file: ", filein)
	conIn <- file( filein, open="r")
	chunkSize <- 1000000

	modeText <- "Building"
	if ( ! is.null( WB)) modeText <- "Adding to"

	type <- match.arg( type)
	if ( type == 'selfUnique') {
		cat( "\n", modeText, "  'selfUnique' WiggleBin for ", speciesID, "...")
	} else {
		cat( "\n", modeText, "  'detectable' WiggleBin for ", speciesID, " by ", otherSpeciesID, "...")
	}

	if ( is.null( WB)) {
		setCurrentSpecies( speciesID)
		wb <- newWB( speciesID)
		wb$Info$FileName <- filein
	} else {
		wb <- WB
		wb$Info$FileName <- base::append( wb$Info$FileName, filein)
	}
	
	uniqueCounts <- wb$BinData
	nReadTotal <- 0
	looping <- TRUE
	while ( looping) {

	hasHeader <- (nReadTotal == 0)
	try(  tmp <- read.delim( file=conIn, header=hasHeader, nrows=chunkSize, 
			colClasses=c( "character", "integer", "integer"), as.is=TRUE))
	if ( nrow(tmp) == 0) break
	if ( nrow(tmp) < chunkSize) looping <- FALSE

	colnames( tmp) <- c( "SEQ_ID", "START", "END")
	seqids <- tmp$SEQ_ID
	froms <- tmp$START
	tos <- tmp$END
	nNow <- nrow( tmp)
	nReadTotal <- nReadTotal + nNow

	curSeq <- ""
	binsetptr <- NULL

	oneBinSize <- WB_getBinSize()

	for( i in 1:nrow(tmp)) {

		# this is one region of unique reads, all overlapping by exactly one base
		iseq <- seqids[i]
		if ( iseq != curSeq) {
			if ( !is.null( binsetptr)) uniqueCounts[[ binsetptr]] <- thisBinSet
			binsetptr <- WB_getBinSetPtrFromSeqID( wb, iseq)
			curSeq <- iseq
			cat( "\n", curSeq, "\tbinPtr=", binsetptr)
			thisBinSet <- uniqueCounts[[ binsetptr]]
			nBins <- nrow(thisBinSet)
			cat( "\tnBins=",nBins)
			allBins <- 1:nBins
		}
		from <- froms[i]
		to <- tos[i]
		# determine how many fixed width reads are here
		nreads <- (to - from) - readSize + 2

		# increment the bin count of all touched bins for each read
		for ( ir in 1:nreads) {
			#firstBin <- WB_base2bin( from + ir - 1)
			firstBin <- ceiling( ( from + ir - 1) / oneBinSize)
			#lastBin <- WB_base2bin( to - nreads + ir)
			lastBin <- ceiling( ( to - nreads + ir) / oneBinSize)
			if ( firstBin < 1 || lastBin > nBins) {
				cat( "\nbad entry: ", base::unlist( tmp[ i, ]), "\t(first,last,Nbins)= ", firstBin, lastBin, nBins)
				next
			}
			x <- firstBin:lastBin
			thisBinSet[ x, WB_UNIQUE_CNT] <- thisBinSet[ x, WB_UNIQUE_CNT] + 1
		}
	}
	uniqueCounts[[ binsetptr]] <- thisBinSet
	cat( " ", nReadTotal)
	} #while

	# show some stats
	cat( "\nDone.")
	close( conIn)
	if ( type == "selfUnique") {
		cat( "\nN_'selfUnique' read regions in file: ", nReadTotal)
	} else {
		cat( "\nN_'otherDetected' read regions in file: ", nReadTotal)
	}

	maxUniqueReadsPerWigBin <- (WB_getBinSize() - readSize + 1) + (( readSize - 1) * 2)

	# now turn the second column into pcts
	for ( i in 1:length( uniqueCounts)) {
		tmp <- uniqueCounts[[i]][ , WB_UNIQUE_CNT] / maxUniqueReadsPerWigBin
		uniqueCounts[[i]][ , WB_UNIQUE_PCT]  <- tmp
	}

	# its ready
	wb$BinData <- uniqueCounts

	# save a binary copy
	if ( type == "selfUnique") {
		fileName <- paste( getCurrentSpeciesFilePrefix(), type, "WB.rda", sep=".")
	} else {
		fileName <- paste( getCurrentSpeciesFilePrefix(), "by", getOtherSpeciesFilePrefix( otherSpeciesID), 
					"detectable", "WB.rda", sep=".")
	}
	assign( "WB", value=wb)
	save( WB, file=fileName)
	cat( "\nWrote binary WiggleBin detectability object: ", fileName,"\n")

	return( WB)
}

