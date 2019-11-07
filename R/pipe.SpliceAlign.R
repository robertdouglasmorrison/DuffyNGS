# pipe.SpliceAlign.R

# run a sample's  .fastq file thru the Splice Alignment pipeline

`pipe.SpliceAlign` <- function( inputFastqFile, sampleID, annotationFile="Annotation.txt", 
				optionsFile="Options.txt", rawReadCount=NULL, 
				asMatePairs=FALSE, verbose=TRUE) {

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nStarting pipe 'SpliceAlign' on Sample:     ", sampleID, 
			"\n\nInput Fastq file:  ", inputFastqFile,"\n")
	}
	
	startT <- proc.time()
	gc()

	# prep the filenames...  for splice-based files
	spliceHit <- paste( sampleID, "splice.bam", sep=".")
	finalHit <- paste( sampleID, "splice.converted.bam", sep=".")
	finalNohit <- getNotSpliceFastqFileNames( sampleID, asMatePairs)
	file.delete( c( spliceHit, finalHit, finalNohit))

	optT <- readOptionsTable( optionsFile)
	resultsPath <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	spliceHit <- file.path( resultsPath, "splicing", spliceHit)
	finalHit <- file.path( resultsPath, "splicing", finalHit)
	file.delete( c( spliceHit, finalHit, finalNohit))

	# make sure there are reads...
	nReadsIn <- getFileLineCount( inputFastqFile[1], sampleID, mode="quick", verbose=TRUE) / 4
	if ( nReadsIn < 1) cat( "\npipe.SpliceAlign:  missing or empty input file:  ", inputFastqFile[1])

	# we need a copy of the genomic bam file for its ref data
	genomeFile <- paste( sampleID, "genomic.bam", sep=".")
	genomeFile <- file.path( resultsPath, "align", genomeFile)

	ans <- NULL
	nrUnique <- nrNohit <- nrMulti <- naBad <- 0

	# get the alignment index info
	alignIndexPath <- getOptionValue( optT, "bowtie2Index.path")
	alignIndex <- getOptionValue( optT, "SpliceIndex")
	alignPolicy <- getOptionValue( optT, "SpliceAlignmentPolicy")
	maxMultiHits <- as.numeric( getOptionValue( optT, "maxMultiHits", notfound=10))
	readBufferSize <- as.numeric( getOptionValue( optT, "readBufferSize", notfound=1000000))

	# get the map filename prefix and expected read strand info
	spliceMapPrefix <- getOptionValue( optT, "SpliceMapPrefix", notfound="spliceMap")

	# get the read sense of this set of reads
	readSense <- getReadSense( sampleID, annotationFile)

	# if doing pairs, the no hit prefix looks like the unpaired format
	nohitPrefix <- finalNohit
	if ( asMatePairs) nohitPrefix <- getNotSpliceFastqFileNames( sampleID, asMatePairs=FALSE)

	if ( nReadsIn > 0) {
		ans <- fastqToBAM( inputFastqFile, spliceHit, k=maxMultiHits, sampleID=sampleID,
				optionsFile=optionsFile, annotationFile=annotationFile, noHitsFile=nohitPrefix, 
				alignIndex=alignIndex, alignPolicy=alignPolicy, 
				asMatePairs=asMatePairs, verbose=verbose)
		nReadsIn <- ans$RawReads
		nrUnique <- ans$UniqueReads
		nrMulti <- ans$MultiReads
		nrAligned <- nrUnique + nrMulti
		nrNohit <- ans$NoHitReads

		if ( nrNohit >= 0) quickFileLineCountRecord( finalNohit, sampleID, lineCount=nrNohit*4, readCount=nrNohit)
		if ( nrAligned >= 0) quickFileLineCountRecord( spliceHit, sampleID, lineCount=nrAligned)
		
		# in case no splices found, re-delete that created BAM file
		if ( nrAligned < 1) file.delete( spliceHit)
	}

	summarize <- function() {

		txt <- vector()
		txt <- base::append( txt, base::paste( "\nFinished pipe 'SpliceAlign':     \t Results: \n"))
		txt <- base::append( txt, base::paste( "\n\nInput File:           \t", inputFastqFile))
		txt <- base::append( txt, base::paste( "\nRaw Reads:            \t", prettyNum( nReadsIn, 
				format="d", big.mark=",")))
		txt <- base::append( txt, base::paste( "\n\nAlignment Policy:     \t", alignPolicy[1]))
		txt <- base::append( txt, base::paste( "\n\nNotSplice File:       \t", finalNohit))
		txt <- base::append( txt, base::paste( "\nN_NotSplice Reads:    \t", prettyNum( nrNohit, 
				format="d", big.mark=","), "\t", 
				as.percent( nrNohit, big.value=nReadsIn)))
		txt <- base::append( txt, base::paste( "\n\nHit Splice File:      \t", finalHit))
		txt <- base::append( txt, base::paste( "\nN_Unique Reads:       \t", prettyNum( nrUnique, 
				format="d", big.mark=","), "\t", 
				as.percent( nrUnique, big.value=nReadsIn)))
		txt <- base::append( txt, base::paste( "\nN_Multi Reads:        \t", prettyNum( nrMulti, 
				format="d", big.mark=","), "\t", 
				as.percent( nrMulti, big.value=nReadsIn)))
		txt <- base::append( txt, "\n")
		return(txt)
	}

	if (verbose) {
		cat( verboseOutputDivider)
	}
	txtSummary <- summarize()
	cat( txtSummary)
	myTime <- elapsedProcTime( startT, proc.time(), N=nReadsIn)
	cat( "\n\nTiming Stats:\n")
	print( myTime)

	summary.path <- file.path( resultsPath, "summary")
	if ( ! file.exists( summary.path)) dir.create( summary.path, recursive=T)
	summaryFile <- file.path( summary.path, paste( sampleID,"splice.Summary.txt", sep="."))
	conTxt <- file( summaryFile, open="w")
	writeLines( txtSummary, con=conTxt, sep="")
	close( conTxt)

	gc()

	return( list( "RawReads"=nReadsIn, "UniqueReads"=nrUnique, "MultiReads"= nrMulti, 
			"NoHitReads"=nrNohit, "Time"=myTime))
}

