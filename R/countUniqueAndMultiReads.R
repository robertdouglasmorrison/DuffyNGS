# countUniqueAndMultiReads.R


`countUniqueAndMultiReads` <- function( filein, readBufferSize=1000000, verbose=TRUE) {

	# this tool ASSUMES raw unsorted BAM data!!

	if ( ! file.exists( filein)) {
		cat( "\nBAM file not found: ", filein)
		return( list( "Unique"=0, "Multi"=0))
	}

	con <- bamReader( filein)
	
	# read in the alignment file one buffer at a time
	hasMore <- TRUE
	nUniq <- nMult <- 0
	repeat {
		if ( ! hasMore) break
		if (verbose) cat( ".")
		chunk <- getNextChunk( con, n=readBufferSize, alignedOnly=TRUE)
		nNow <- size(chunk)
		if ( nNow < 1) break
		if ( nNow < readBufferSize) hasMore <- FALSE

		whoUniq <- which.unique.align(chunk)
		whoMult <- which.multi.read(chunk)

		nUniq <- nUniq + length( whoUniq)
		nMult <- nMult + length( whoMult)
	} # end of each buffer...

	bamClose( con)
	return( list( "Unique"=nUniq, "Multi"=nMult))
}


`countReadsBySpecies` <- function( filein, readBufferSize=1000000, verbose=TRUE) {

	# this tools is FINE with SORTED BAM files..!!

	# count up the reads in a BAM file by their organism
	speciesSet <- getCurrentTargetSpecies()
	NSP <- length( speciesSet)
	nReads <- rep.int(0, NSP)
	names( nReads) <- speciesSet

	if ( ! file.exists( filein)) {
		cat( "\nBAM file not found: ", filein)
		return( nReads)
	}

	con <- bamReader( filein)
	refData <- getRefData( con)
	
	# pre fetch the set of all possible seqIDs, to speed up  the calls
	allSeqIDs <- ptrSeqIDs <- vector()
	curSpecies <- getCurrentSpecies()
	on.exit( setCurrentSpecies( curSpecies))
	for ( i in 1:NSP) {
		setCurrentSpecies( speciesSet[i])
		smap <- getCurrentSeqMap()
		ns <- nrow(smap)
		allSeqIDs <- c( allSeqIDs, smap$SEQ_ID)
		ptrSeqIDs <- c( ptrSeqIDs, rep.int( i, ns))
	}

	# read in the alignment file one buffer at a time
	hasMore <- TRUE
	repeat {
		if ( ! hasMore) break
		if (verbose) cat( ".")
		chunk <- getNextChunk( con, n=readBufferSize, alignedOnly=TRUE)
		nNow <- size(chunk)
		if ( nNow < 1) break
		if ( nNow < readBufferSize) hasMore <- FALSE

		# we want reads, not alignments!
		whoPrim <- which( ! secondaryAlign(chunk))

		# look at those SeqIDs
		refIDs <- refID(chunk)[whoPrim]
		seqIDs <- refID2seqID( refIDs, refData=refData)
		where <- match( seqIDs, allSeqIDs, nomatch=0)
		myPtrs <- ptrSeqIDs[ where]
		for ( i in 1:NSP) {
			nNow <- sum( myPtrs == i)
			nReads[i] <- nReads[i] + nNow
		}

	} # end of each buffer...

	bamClose( con)
	return( as.list( nReads))
}


`countReadsBySequence` <- function( filein, readBufferSize=1000000, verbose=TRUE) {

	# count up the reads in a BAM file by their chromosome hits
	if ( ! file.exists( filein)) {
		cat( "\nBAM file not found: ", filein)
		return( nReads)
	}
	con <- bamReader( filein)
	refData <- getRefData( con)
	
	# get the set of possible sequence IDs from the BAM file
	seqSet <- refData$SN
	NSEQ <- length( seqSet)
	nReads <- rep.int(0, NSEQ)
	names(nReads) <- seqSet
	
	# read in the alignment file one buffer at a time
	hasMore <- TRUE
	repeat {
		if ( ! hasMore) break
		if (verbose) cat( ".")
		chunk <- getNextChunk( con, n=readBufferSize, alignedOnly=TRUE)
		nNow <- size(chunk)
		if ( nNow < 1) break
		if ( nNow < readBufferSize) hasMore <- FALSE
		# look at those SeqIDs
		refIDs <- refID(chunk)
		seqIDs <- refID2seqID( refIDs, refData=refData)
		where <- match( seqIDs, seqSet, nomatch=0)
		for ( i in 1:NSEQ) {
			nNow <- sum( where == i)
			nReads[i] <- nReads[i] + nNow
		}
	} # end of each buffer...

	bamClose( con)
	return( nReads)
}
