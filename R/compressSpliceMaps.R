# compressSpliceMaps.R

# combine onr or more splice maps for decoding splice junction alignments

`compressSpliceMaps` <- function( outPath=NULL, spliceMaps, spliceMapPrefix="SpliceMap") {

	if ( is.null( outPath)) stop( "Missing required 'outPath' argument for writing compressed 'spliceMap' files.")

	# check and/or make directory
	if ( ! file.exists( outPath)) dir.create( outPath, recursive=TRUE)

	# loop over all spliceMaps
	for( f in spliceMaps) {
		cat( "\nCompressing:   ", f, "\n")
		mydf <- read.delim( file=f, as.is=TRUE)

		# make sure it looks like a spliceMap
		if ( ! all( colnames(mydf) == SPLICEMAP_COLUMNS)) {
			warning( paste( "File does not have SpliceMap column names...skipping:  ", f))
			next
		}

		# factor by seqID
		seqFactor <- factor( mydf$SEQ_ID)
		rowPtrs <- tapply( 1:nrow(mydf), INDEX=seqFactor, FUN=NULL)

		for ( i in 1:nlevels(seqFactor)) {
			seqName <- levels( seqFactor)[ i]
			myrows <- which( rowPtrs == i)
			smallDF <- mydf[ myrows, ]
			outfile <- paste( spliceMapPrefix, seqName, "rda", sep=".")
			cat( "\n", outfile, "\tN_splices: ", nrow(smallDF))
			save( smallDF, file=file.path( outPath, outfile))
		}

	}
	cat( "\nDone.\n")
}
