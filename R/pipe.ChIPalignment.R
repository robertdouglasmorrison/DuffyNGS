# pipe.ChIPalignment.R

# do the A to Z alignment pipeline on a sample

`pipe.ChIPalignment` <- function( sampleID=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt",
				verbose=TRUE) {

	if (verbose) {
		cat( verboseOutputDivider)
		curHost <- system( command="hostname", intern=TRUE)
		cat( "\nStarting  'ChIP Alignment Pipeline' on Sample:     ", sampleID,
			"\nSample Data Type:  \tChIP-seq",
			"\n\nHost Machine:      \t", curHost, "\n", R.version.string,
			"\nStart Date/Time:    \t", date(), "\n", sep="")
	}

	target <- getAndSetTarget( optionsFile, verbose=T)
	results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)

	# file(s) to process comes from annotation and options...
	rawFastq <- getRawFastqFileNames( sampleID, annotationFile, optionsFile, verbose=verbose)
	inputFastqFiles <- rawFastq$files
	asMatePairs <- rawFastq$asMatePairs

	gc()
	startT <- proc.time()

	nReadsIn <- NULL
	nNoHit <- 0

	# Part 1:  ribo clearing:  never do it
	doRiboClearing <- FALSE
	if ( ! doRiboClearing) {
		nRiboU <- nRiboM <- nRiboMoreK <- 0
		notRiboFiles <- inputFastqFiles
	}
	nRibo <- ( nRiboU + nRiboM)

	# Part 2:   genomic -- always do
	nGenomicU <- nGenomicM <- 0
	genomicInFiles <- notRiboFiles
	notGenomicFiles <- getNotGenomicFastqFileNames( sampleID, asMatePairs)
	status2 <- pipe.GenomicAlign( inputFastqFile=genomicInFiles, sampleID, optionsFile=optionsFile, 
				asMatePairs=asMatePairs, verbose=verbose, rawReadCount=nReadsIn)

	# now we really know how many reads we started with...
	nReadsInGenomic <- status2$RawReads
	if ( is.null( nReadsIn)) {
		nReadsIn <- nReadsInGenomic
	}
	nGenomicU <- status2$UniqueReads
	nGenomicM <-  status2$MultiReads
	nGenomic <- (nGenomicU + nGenomicM)
	nNotRibo <- nReadsIn - nRibo
	nNotGenomic <- nNotRibo - nGenomic

	# Part 3:  splices:  never do it
	nNotSplice <- nNotGenomic
	doSplices <- FALSE
	if ( ! doSplices) {
		nSpliceU <- nSpliceM <- 0
		notSpliceFiles <- notGenomicFiles
	}
	nSplice <- ( nSpliceU + nSpliceM)
	nNotSplice <- nNotGenomic - nSplice

	# final no-hits
	nNoHit <- nNotSplice
	noHitFile <- notSpliceFiles
	finalNoHitFile <- paste( sampleID,"noHits.fastq.gz", sep=".")
	finalNoHitFile <- file.path( results.path, "fastq", finalNoHitFile)

	if (asMatePairs) {
		# 2 .gz files in the current folder...
		# discordant reads went to both alignments and noHits...
		# so we need to strip out discordant and combine into one final no hits file
		noHitAns <- pipe.RemoveDiscordantAlignsFromNoHits( sampleID, file=noHitFile, mode="pair")
		nNoHit <- noHitAns$NoHitReads
		quickFileLineCountRecord( finalNoHitFile, sampleID, lineCount=nNoHit*4, readCount=nNoHit)
		file.delete( noHitFile)
	} else {
		file.rename( noHitFile, finalNoHitFile)
		quickFileLineCountRecord( finalNoHitFile, sampleID, lineCount=nNoHit*4, readCount=nNoHit)
	}

	# convert any BAM files that need it
	pipe.ConvertAllBAMs( sampleID, annotationFile=annotationFile, optionsFile=optionsFile, 
				rawReadCount=nReadsIn, dataType="ChIP-seq", verbose=verbose)

	# now do the cleanup...
	pipe.FileCleanup( sampleID, optionsFile=optionsFile, verbose=verbose)

	myTime <- elapsedProcTime( startT, proc.time(), N=nReadsIn)

	# local function to make the results summary
	`makeAlignSummary` <- function() {
	out <- vector()
	out <- base::append( out, base::paste( "\n\nInput File:           \t", inputFastqFiles))
	out <- base::append( out, base::paste( "\nN_Raw Reads:          \t", 
			prettyNum( nReadsIn, width=12, format="d", big.mark=","),"\n"))

	out <- base::append( out, base::paste( "\nNoHits File:          \t", finalNoHitFile))
	out <- base::append( out, base::paste( "\nN_NoHit Reads:        \t", 
			prettyNum( nNoHit, width=12, format="d", big.mark=","), "\t", 
			as.percent( nNoHit, big.value=nReadsIn),"\n"))

	out <- base::append( out, base::paste( "\nN_Unique Genomic:     \t", 
			prettyNum( nGenomicU, width=12, format="d", big.mark=","), "\t", 
			as.percent( nGenomicU, big.value=nReadsIn)))
	out <- base::append( out, base::paste( "\nN_Multi Genomic:      \t", 
			prettyNum( nGenomicM, width=12, format="d", big.mark=","), "\t", 
			as.percent( nGenomicM, big.value=nReadsIn)))
	out <- base::append( out, base::paste( "\nAll Genomic Reads:    \t", 
			prettyNum( nGenomic, width=12, format="d", big.mark=","), "\t", 
			as.percent( nGenomic, big.value=nReadsIn),"\n"))

	return(out)
	}  # end of local results summary


	cat( verboseOutputDivider)
	cat( "\n\nFinished 'ChIP Alignment Pipeline' on Sample:     ", sampleID, "\n")
	cat( "\nPipeline:             \t Results:")

	alignSummary <- makeAlignSummary()
	cat( alignSummary)
	summary.path <- file.path( results.path, "summary")
	if ( ! file.exists( summary.path)) dir.create( summary.path, recursive=T)
	summaryFile <- file.path( summary.path, paste( sampleID, "pipeline.Summary.txt", sep="."))
	writeLines( alignSummary, con=summaryFile, sep="")

	cat( "\n\nTiming Stats: \n")
	print( myTime)
	gc()

	out <- list( "nReadsIn"=nReadsIn, "nNoHit"=nNoHit, "nRibo"=nRibo, "nGenomic"=nGenomic, "nSplice"=nSplice)

	return( out)
}

