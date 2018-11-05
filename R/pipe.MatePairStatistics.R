# pipe.MatePairStatistics.R

# look at the metrics of aligned paired-end mates

`pipe.MatePairStatistics` <- function( sampleID,  annotationFile="Annotation.txt", optionsFile="Options.txt",
			maxReads=NULL, maxInsertSize=4000, verbose=TRUE) {

	# look at BAM file details
	optT <- readOptionsTable( optionsFile)
	results.path <- getOptionValue( optT, "results.path", notfound=".")
	genomicFile <- file.path( results.path, "align", paste( sampleID, "genomic.bam", sep="."))
	if ( ! file.exists( genomicFile)) {
		genomicFile <- file.path( results.path, "align", paste( sampleID, "genomic.sorted.bam", sep="."))
	}
	if ( ! file.exists( genomicFile)) {
		cat( "\nNo genomic BAM file found for sample:  ", sampleID)
		return(NULL)
	}

	readBufferSize <- as.integer( getOptionValue( optT, "readBufferSize", notfound="1000000"))
	if ( !is.null(maxReads)) readBufferSize <- min( maxReads, readBufferSize)

	conGenome <- bamReader( genomicFile)
	headerGenome <- getHeader( conGenome)
	genomeRefData <- getRefData( conGenome)

	# set up to do a buffer at a time
	pairSizes <- pairSpecies <- unpairSpecies <- vector();
	nReads <- nPairs <- nUnpaired <- 0;
	hasMore <- TRUE

	# grab a buffer
	repeat {
		if ( ! hasMore) break
		if ( !is.null(maxReads) && nReads >= maxReads) break
		cat( "\nReadBAM..")
		chunk <- getNextChunk( conGenome, n=readBufferSize, alignedOnly=TRUE)
		nNow <- size(chunk)
		if ( nNow < 1) break
		if ( nNow < readBufferSize) hasMore <- FALSE
		nReads <- nReads + nNow
		cat( "  N_Aligned: ", prettyNum( as.integer( nReads), big.mark=","))

		# do it as a data frame to interate over the file less times
		chunkDF <- as.data.frame( chunk)
		myRefID <- chunkDF$refid
		myInsertSize <- chunkDF$insertsize
		# the SeqID field is in spliceMap notation  'seqID::geneID::spliceID'
		mySIDs <- refID2seqID( myRefID, refData=genomeRefData)
		mySpeciesIDs <- getSpeciesFromSeqID( mySIDs)

		# now gather up what we want
		isPair <- which( myInsertSize > 0 & myInsertSize <= maxInsertSize)
		isUnpaired <- which( myInsertSize == 0 | myInsertSize > maxInsertSize)
		nPairs <- nPairs + length( isPair)
		nUnpaired <- nUnpaired + length( isUnpaired)

		# keep insert size by species
		pairSpecies <- c( pairSpecies, mySpeciesIDs[isPair])
		pairSizes <- c( pairSizes, myInsertSize[isPair])
		unpairSpecies <- c( unpairSpecies, mySpeciesIDs[isUnpaired])
		rm( chunk, chunkDF)
		gc()
	} # end of each buffer...

	bamClose( conGenome)

	# tally results
	cat( "\nN_Aligns: ", nReads, "\n")
	tbl <- as.percent( tblAll <- c("Paired"=nPairs*2, "Unpaired"=nUnpaired), big.value=nReads)
	print( tbl)
	cat( "\nMate Pairs by Species:\n")
	tbl <- as.percent( tblPair <- table( pairSpecies), big.value=nPairs)
	print( tbl)
	cat( "\nUnpaired Reads by Species:\n")
	tbl <- as.percent( tblUnpair <- table( unpairSpecies), big.value=nUnpaired)
	print( tbl)

	# look at the density curve by species
	ans <- tapply( pairSizes, factor( pairSpecies), density)
	bigY <- 0
	for (densAns in ans) bigY <- max( bigY, max(densAns$y))

	plot( 1,1, type='n', xlim=c(0,maxInsertSize), ylim=c(0,bigY), xlab="Insert Size (bp)", ylab="Density")
	speciesSet <- sort( unique( pairSpecies))
	colorAns <- colorBySpecies( speciesSet)
	for (i in 1:length(ans)) {
		densAns <- ans[[i]]
		lines( densAns, lwd=2, col=colorAns$colors[i])
	}
	legend( "topright", names(colorAns$legend.colors), lwd=3, col=colorAns$legend.colors)

	return( list( "Overview"=tblAll, "Pairs"=tblPair, "Unpaired"=tblUnpair))
}
