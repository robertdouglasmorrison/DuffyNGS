# pipe.GenomicAlign.R

# run a sample's  .fastq file thru the Genomic Alignment pipeline

`pipe.GenomicAlign` <- function( inputFastqFile, sampleID, optionsFile="Options.txt", 
				asMatePairs=FALSE, rawReadCount=NULL, verbose=TRUE) {

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nStarting pipe 'GenomicAlign' on Sample:     ", sampleID, 
			"\n\nInput Fastq file:  ", inputFastqFile,"\n")
	}
	gc()
	
	startT <- proc.time()

	# prep the filenames...
	finalHit <- paste( sampleID, "genomic.bam", sep=".")
	finalNohit <- getNotGenomicFastqFileNames( sampleID, asMatePairs)
	file.delete( c( finalHit, finalNohit))

	# make sure there are reads...
	nReadsIn <- getFileLineCount( inputFastqFile[1], sampleID, verbose=TRUE, mode="quick") / 4
	if ( nReadsIn < 1) stop( paste( "pipe.GenomicAlign:  missing or empty input file:  ", inputFastqFile[1]))

	optT <- readOptionsTable( optionsFile)
	resultsPath <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	alignIndexPath <- getOptionValue( optT, "bowtie2Index.path")
	alignIndex <- getOptionValue( optT, "GenomicIndex")
	alignPolicy <- getOptionValue( optT, "GenomicAlignmentPolicy")
	readBufferSize <- as.numeric( getOptionValue( optT, "readBufferSize", notfound="1000000"))
	maxMultiHits <- as.numeric( getOptionValue( optT, "maxMultiHits", notfound=10))

	maxFastqReads <- as.numeric( getOptionValue( optT, "maxFastqReads", notfound=""))
	if ( is.na( maxFastqReads)) maxFastqReads <- NULL

	finalHit <- file.path( resultsPath, "align", finalHit)
	file.delete( finalHit)

	# if doing pairs, the no hit prefix looks like the unpaired format
	nohitPrefix <- finalNohit
	if ( asMatePairs) nohitPrefix <- getNotGenomicFastqFileNames( sampleID, asMatePairs=FALSE)

	# do the genomic alignment
	ans <- fastqToBAM( inputFastqFile, finalHit, k=maxMultiHits, sampleID=sampleID,
			optionsFile=optionsFile, noHitsFile=nohitPrefix, 
			alignIndex=alignIndex, alignPolicy=alignPolicy, asMatePairs=asMatePairs,
			maxReads=maxFastqReads, verbose=verbose)
	nReadsIn <- ans$RawReads
	nrUnique <- ans$UniqueReads
	nrMulti <- ans$MultiReads
	nrNohit <- ans$NoHitReads
	nrAligned <- nrUnique + nrMulti

	if ( nReadsIn > 0) quickFileLineCountRecord( inputFastqFile, sampleID, lineCount=nReadsIn*4, readCount=nReadsIn)
	if ( nrNohit >= 0) quickFileLineCountRecord( finalNohit, sampleID, lineCount=nrNohit*4, readCount=nrNohit)
	if ( nrAligned >= 0) quickFileLineCountRecord( finalHit, sampleID, lineCount=nrAligned, readCount=nrAligned)

	# we may not have known how many reads there really were until now...
	if (is.null(rawReadCount)) rawReadCount <- nReadsIn

	# if no reads did align at all, redelete the created BAM file
	if ( nrAligned < 1) file.delete( finalHit)


	makeGenomicSummary <- function() {

	out <- vector()
	out <- base::append( out, base::paste( "\nFinished pipe 'GenomicAlign':     \t Results: "))
	out <- base::append( out, "\n")
	out <- base::append( out, base::paste( "\nInput File:            \t", inputFastqFile))
	out <- base::append( out, base::paste( "\nRaw Reads:             \t", prettyNum( nReadsIn, format="d", 
			big.mark=",", width=12)))
	out <- base::append( out, base::paste( "\n\nAlignment Policy:      \t", alignPolicy))
	out <- base::append( out, base::paste( "\n\nNot Genomic File:      \t", finalNohit))
	out <- base::append( out, base::paste( "\nN_NotGenomic Reads:    \t", prettyNum( nrNohit, format="d", 
			big.mark=",", width=12), "\t", as.percent( nrNohit, big.value=nReadsIn),
			"\t", as.percent( nrNohit, big.value=rawReadCount)))
	out <- base::append( out, base::paste( "\n\nHit Genomic File:      \t", finalHit))
	out <- base::append( out, base::paste( "\nN_Unique Reads:        \t", prettyNum( nrUnique, format="d", 
			big.mark=",", width=12), "\t", as.percent( nrUnique, big.value=nReadsIn),
			"\t", as.percent( nrUnique, big.value=rawReadCount)))
	out <- base::append( out, base::paste( "\nN_Multi Reads:         \t", prettyNum( nrMulti, format="d", 
			big.mark=",", width=12), "\t", as.percent( nrMulti, big.value=nReadsIn),
			"\t", as.percent( nrMulti, big.value=rawReadCount)))
	
	grandTotal <- nrUnique + nrMulti
	out <- base::append( out, base::paste( "\n\nTotal Genomic Reads:   \t", prettyNum( grandTotal, 
			format="d", big.mark=",", width=12), "\t", as.percent( grandTotal, 
			big.value=nReadsIn), "\t", as.percent( grandTotal, big.value=rawReadCount)))

	return( out)
	}  # end of 'makeResultText()'


	myTime <- elapsedProcTime( startT, proc.time(), N=nReadsIn)

	outText <- makeGenomicSummary()
	cat( outText)

	allSummary <- c( outText)
	summary.path <- file.path( resultsPath, "summary")
	if ( ! file.exists( summary.path)) dir.create( summary.path, recursive=T)
	summaryFile <- file.path( summary.path, paste( sampleID, "genomic.Summary.txt", sep="."))
	writeLines( allSummary, con=summaryFile, sep="")

	cat( "\n\nTiming Stats: \n")
	print( myTime)
	gc()

	return( list( "RawReads"=nReadsIn, "UniqueReads"=nrUnique, "MultiReads"= nrMulti, 
			"NoHitReads"=nrNohit, "Time"=myTime))
}

