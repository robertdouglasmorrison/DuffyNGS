# syntheticReads.R -- make all possible overlapping reads from a region of a chromosome

`syntheticReads` <- function( seqID, seqDNA, type=c("fasta", "fastq"), readSize=32, 
			firstBase=1, nReads=100000, outfile=NULL) {

	# set up to get all N-mers for this chunk, watching for end of chromosome
	nBase <- length( seqDNA)
	sizeM1 <- readSize - 1
	veryLastStart <- nBase - sizeM1
	if ( firstBase < 1 || firstBase > veryLastStart) {
		stop( "Invalid base range: 'firstBase'")
	}
	lastBase <- firstBase + nReads - 1
	if ( lastBase > veryLastStart) lastBase <- veryLastStart
	fromSet <- firstBase:lastBase
	toSet <- fromSet + sizeM1

	# make the little sequences
	who <- base::mapply( FUN=`:`, fromSet, toSet, SIMPLIFY=FALSE)
	dna <- sapply( who, FUN=function(x) base::paste( seqDNA[x], collapse=""))

	# make the readID names
	nams <- base::paste( seqID, "::", fromSet, sep="")

	# format it...
	type <- match.arg( type)
	if ( type == "fasta") {
		outText <- base::paste( ">", nams, "\n", dna, sep="")
	} else {
		score <- base::paste( rep("I", times=readSize), collapse="")
		outText <- base::paste( "@", nams, "\n", dna, "\n+\n", score, sep="")
	}

	if ( ! is.null( outfile)) {
		con <- file( outfile, open="w")
		writeLines( outText, con=con)
		close( con)	
		return( list( "nReads"=length(fromSet), "reads"=dna, "text"=NULL))
	} else {
		return( list( "nReads"=length(fromSet), "reads"=dna, "text"=outText))
	}
}


`syntheticReadsToFile` <- function( seqID, seqDNA, outfile, type=c("fasta", "fastq"), readSize=32, verbose=T) {

	# set up to make a file of reads from a chunk of DNA
	type <- match.arg( type)
	nBase <- nchar( seqDNA[1])
	if ( nBase < readSize) stop( "expected a character string of genomic DNA")
	sizeM1 <- readSize - 1
	veryLastStart <- nBase - sizeM1

	myDNA <- strsplit( seqDNA, split="")[[1]]

	# make the little sequences in chunks
	conOut <- gzfile( outfile, open="wt")
	chunkSize=1000000
	nReadsThisSeq <- 0

	for( ib in seq( 1, veryLastStart, by=chunkSize)) {

		# build a small .fasta file of N-mers
		if (verbose) cat( "\nmaking reads..")
		ans <- syntheticReads( seqID=seqID, seqDNA=myDNA, type=type, readSize=readSize, 
				firstBase=ib, nReads=chunkSize, outfile=NULL)
		nReads <- ans$nReads
		nReadsThisSeq <- nReadsThisSeq + nReads
		if (verbose) cat( "  writing reads..")
		writeLines( ans$text, conOut)
		if (verbose) cat( "  N = ", nReadsThisSeq)
	}
	close( conOut)
	if (verbose) cat( "\nWrote file:  ", outfile, "\nN_Reads: ", nReadsThisSeq)
}
