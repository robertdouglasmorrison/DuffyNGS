# spliceJunctionLookup.R

# get the SpliceMap rows that go with each splice alignment ID

`spliceJunctionLookup` <-
function( spliceIDs, spliceMapPath=bowtie2Par("IndexPath"), spliceMapPrefix="spliceMap") {

	# make sure the path to compressed spliceMaps is ok
	info <- file.info( spliceMapPath)
	if ( is.na( info$isdir) || ( ! info$isdir)) stop("'spliceMapPath' directory not found.")
	
	# turn the splice alignmentIDs back to seqID, geneID, spliceID
	# the IDs come in as:    seqid::geneid::spliceid
	terms <- base::unlist( strsplit( spliceIDs, split="::", fixed=TRUE))
	nT <- length(terms)
	seqIDset <- terms[ seq( 1, nT, by=3) ]
	geneIDset <- terms[ seq( 2, nT, by=3) ]
	nID <- length( seqIDset)

	# make storage for the results
	gOut <- splOut <- sOut <- sizOut <- posOut <- endOut <- vector( mode="character", length=nID)

	# factor by the seqID, so we visit each compressed spliceMap once
	seqFactor <- factor( seqIDset)
	rowPtrs <- tapply( 1:nID, INDEX=seqFactor, FUN=NULL)

	for ( i in 1:nlevels( seqFactor)) {
		thisSeqID <- levels( seqFactor)[i]
		theseRows <- which( rowPtrs == i)

		# get that file... creates an object called 'smallDF'
		thisfilename <- base::paste( spliceMapPrefix, thisSeqID, "rda", sep=".")
		thisfile <- file.path( spliceMapPath, thisfilename)
		if ( ! file.exists( thisfile)) stop( paste( "Required compressed spliceMap file not found:  ", thisfile))
		load( file=thisfile, envir=environment())

		# make a matching alignmentID to check against
		smlID <- base::paste( smallDF$SEQ_ID, smallDF$GENE_ID, smallDF$SPLICE_ID, sep="::")
		hits <- base::match( spliceIDs[ theseRows], smlID, nomatch=0)
		outRows <- theseRows[ hits > 0]
		fromRows <- hits[ hits > 0]

		# move the data
		gOut[ outRows] <- smallDF$GENE_ID[ fromRows]
		splOut[ outRows] <- smallDF$SPLICE_ID[ fromRows]
		sOut[ outRows] <- smallDF$SEQ_ID[ fromRows]
		sizOut[ outRows] <- smallDF$SIZES[ fromRows]
		posOut[ outRows] <- smallDF$POSITIONS[ fromRows]
		endOut[ outRows] <- smallDF$ENDS[ fromRows]
	}

	# check for missing splices
	whoBlank <- which( gOut == "")
	if ( length( whoBlank) > 0) {
		missing <- base::table( spliceIDs[whoBlank])
		cat( "\nSome Splices missing from splice junction maps:\n")
		print( head( sort( missing, decreasing=TRUE), 100))
	}

	out <- data.frame( gOut, splOut, sOut, sizOut, posOut, endOut, stringsAsFactors=FALSE)
	colnames( out) <- SPLICEMAP_COLUMNS
	rownames( out) <- 1:nrow( out)
	return( out)
}

