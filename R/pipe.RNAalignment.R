# pipe.RNAalignment.R

# do the A to Z alignment pipeline on a sample

`pipe.RNAalignment` <- function( sampleID=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt",
				verbose=TRUE) {

	if (verbose) {
		cat( verboseOutputDivider)
		curHost <- system( command="hostname", intern=TRUE)
		cat( "\nStarting  'RNA Alignment Pipeline' on Sample:     ", sampleID,
			"\nSample Data Type:  \tRNA-seq",
			"\n\nHost Machine:      \t", curHost, "\n", R.version.string,
			"\nStart Date/Time:    \t", date(), "\n", sep="")
	}

	target <- getAndSetTarget( optionsFile, sampleID=sampleID, annotationFile=annotationFile)
	results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)

	# file(s) to process comes from annotation and options...
	rawFastq <- getRawFastqFileNames( sampleID, annotationFile, optionsFile, verbose=verbose)
	inputFastqFiles <- rawFastq$files
	asMatePairs <- rawFastq$asMatePairs

	gc()
	startT <- proc.time()

	nReadsIn <- NULL
	nNoHit <- 0

	# Part 1:  ribo clearing
	# either bypasss or do it, and gather the counts and filenames as we go
	doRiboClearing <- getOptionTrue( optionsFile, "RiboIndex")

	# mated paired alignment is incompatible with ribo clearing.  BOWTIE requires both mates align to the same
	# contig for the 2 mates to be removed from the set of reads.  So often, a gene is "ribo cleared", but many 
	# reads still align there at genomic step...
	if (doRiboClearing && asMatePairs) {
		cat( "\n\n-------------------------")
		cat( "\nBOWTIE paired end alignment mode incompatible with Ribo Clearing..")
		forcePairs <- getOptionTrue( optionsFile, "forcePairedEnd", notfound=FALSE, verbose=T)
		if ( !forcePairs) {
			cat( "\nEither set option 'forcePairedEnd TRUE' to override.")
			cat( "\nOr set 'StrandSpecific' to FALSE in the annotation file.")
			cat( "\n\nPipeline halting now.  ")
			cat( "\n-------------------------", "\n")
			stop()
		} else {
			cat( "\nForce Paired End override:  keeping paired end mode anyway..")
			cat( "\n-------------------------", "\n")
		}
	}

	if ( ! doRiboClearing) {
		nRiboU <- nRiboM <- nRiboMoreK <- 0
		notRiboFiles <- inputFastqFiles
	} else {
		status1 <- pipe.RiboClear( inputFastqFiles, sampleID, optionsFile=optionsFile, 
					asMatePairs=asMatePairs, verbose=verbose, rawReadCount=nReadsIn)
		nReadsIn <- status1$RawReads
		nNotRibo <- status1$NoHitReads
		nRiboU <- status1$UniqueReads
		nRiboM <-  status1$MultiReads
		notRiboFiles <- getNotRiboFastqFileNames( sampleID, asMatePairs)
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
	nNotGenomic <- status2$NoHitReads

	if ( ! exists( "nNotRibo")) nNotRibo <- nReadsIn - nRibo

	# Part 3:  splices
	# either bypasss or do it, and gather the counts and filenames as we go
	# the splices are not done as pairs, no matter what...  but do them separately if needed
	nNotSplice <- nNotGenomic
	doSplices <- getOptionTrue( optionsFile, "SpliceIndex")
	if ( ! doSplices) {
		nSpliceU <- nSpliceM <- 0
		notSpliceFiles <- notGenomicFiles
	} else {
		spliceInFiles <- notGenomicFiles
		notSpliceFiles <- getNotSpliceFastqFileNames( sampleID, asMatePairs=asMatePairs)
		status3 <- pipe.SpliceAlign( inputFastqFile=spliceInFiles, sampleID, 
				annotationFile=annotationFile, optionsFile=optionsFile,
				asMatePairs=asMatePairs, verbose=verbose, rawReadCount=nReadsIn)
		nSpliceU <- status3$UniqueReads
		nSpliceM <-  status3$MultiReads
		nNotSplice <- status3$NoHitReads
	}
	nSplice <- ( nSpliceU + nSpliceM)

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
				rawReadCount=nReadsIn, dataType="RNA-seq", verbose=verbose)

	# now do the cleanup...
	pipe.FileCleanup( sampleID, optionsFile=optionsFile, verbose=verbose)

	myTime <- elapsedProcTime( startT, proc.time(), N=nReadsIn)

	# local function to make the results summary
	`makeAlignSummary` <- function() {
	out <- vector()
	out <- base::append( out, "\n")
	out <- base::append( out, base::paste( "\nInput File:           \t", inputFastqFiles))
	out <- base::append( out, base::paste( "\nN_Raw Reads:          \t", 
			prettyNum( nReadsIn, width=12, format="d", big.mark=","),"\n"))

	out <- base::append( out, base::paste( "\nNoHits File:          \t", finalNoHitFile))
	out <- base::append( out, base::paste( "\nN_NoHit Reads:        \t", 
			prettyNum( nNoHit, width=12, format="d", big.mark=","), "\t", 
			as.percent( nNoHit, big.value=nReadsIn),"\n"))

	out <- base::append( out, base::paste( "\nN_Unique Ribo:        \t", 
			prettyNum( nRiboU, width=12, format="d", big.mark=","), "\t", 
			as.percent( nRiboU, big.value=nReadsIn)))
	out <- base::append( out, base::paste( "\nN_Multi Ribo:         \t", 
			prettyNum( nRiboM, width=12, format="d", big.mark=","), "\t", 
			as.percent( nRiboM, big.value=nReadsIn)))
	out <- base::append( out, base::paste( "\nAll Ribo Reads:       \t", 
			prettyNum( nRibo, width=12, format="d", big.mark=","), "\t", 
			as.percent( nRibo, big.value=nReadsIn),"\n"))

	out <- base::append( out, base::paste( "\nN_Unique Genomic:     \t", 
			prettyNum( nGenomicU, width=12, format="d", big.mark=","), "\t", 
			as.percent( nGenomicU, big.value=nReadsIn)))
	out <- base::append( out, base::paste( "\nN_Multi Genomic:      \t", 
			prettyNum( nGenomicM, width=12, format="d", big.mark=","), "\t", 
			as.percent( nGenomicM, big.value=nReadsIn)))
	out <- base::append( out, base::paste( "\nAll Genomic Reads:    \t", 
			prettyNum( nGenomic, width=12, format="d", big.mark=","), "\t", 
			as.percent( nGenomic, big.value=nReadsIn),"\n"))

	out <- base::append( out, base::paste( "\nN_Unique Splice:      \t", 
			prettyNum( nSpliceU, width=12, format="d", big.mark=","), "\t", 
			as.percent( nSpliceU, big.value=nReadsIn)))
	out <- base::append( out, base::paste( "\nN_Multi Splice:       \t", 
			prettyNum( nSpliceM, width=12, format="d", big.mark=","), "\t", 
			as.percent( nSpliceM, big.value=nReadsIn)))
	out <- base::append( out, base::paste( "\nAll Splice Reads:     \t", 
			prettyNum( nSplice, width=12, format="d", big.mark=","), "\t", 
			as.percent( nSplice, big.value=nReadsIn),"\n"))
	return(out)
	}  # end of local results summary


	cat( verboseOutputDivider)
	cat( "\n\nFinished 'RNA Alignment Pipeline' on Sample:     ", sampleID, "\n")
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

