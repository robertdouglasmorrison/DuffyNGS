# pipe_RiboClear.R

# run a sample's  .fastq file thru the riboClearing pipeline

`pipe.RiboClear` <-
function( inputFastqFile, sampleID, annotationFile="Annotation.txt", optionsFile="Options.txt",
		asMatePairs=FALSE, verbose=TRUE, rawReadCount=NULL)
{

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nStarting pipe 'RiboClear' on Sample:     ", sampleID, 
			"\n\nInput Fastq file:  ", inputFastqFile,"\n")
	}
	
	startT <- proc.time()
	gc()

	# prep the filenames and remove any old copy...
	finalHit <- paste( sampleID, "ribo.bam", sep=".")
	finalHit2 <- paste( sampleID, "ribo.converted.bam", sep=".")
	finalNohit <- getNotRiboFastqFileNames( sampleID, asMatePairs)
	file.delete( c( finalHit, finalHit2, finalNohit))

	# make sure there are reads...
	nReadsIn <- getFileLineCount( inputFastqFile[1], sampleID, verbose=TRUE, mode="quick") / 4
	if ( nReadsIn < 1) stop( paste( "pipe.RiboCLear:  missing or empty input file:  ", inputFastqFile[1]))

	optT <- readOptionsTable( optionsFile)
	alignIndexPath <- getOptionValue( optT, "bowtie2Index.path")
	alignIndex <- getOptionValue( optT, "RiboIndex")
	alignPolicy <- getOptionValue( optT, "RiboAlignmentPolicy")
	riboMap <- getOptionValue( optT, "RiboMap")
	riboMap <- file.path( alignIndexPath, riboMap)
	maxMultiHits <- as.numeric( getOptionValue( optT, "maxMultiHits", notfound="10"))

	maxFastqReads <- as.numeric( getOptionValue( optT, "maxFastqReads", notfound=""))
	if ( is.na( maxFastqReads)) maxFastqReads <- NULL

	resultsPath <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	finalHit <- file.path( resultsPath, "riboClear", finalHit)
	finalHit2 <- file.path( resultsPath, "riboClear", finalHit2)
	# also delete any final results copy
	file.delete( c( finalHit, finalHit2))

	# if doing pairs, the no hit prefix looks like the unpaired format
	nohitPrefix <- finalNohit
	if ( asMatePairs) nohitPrefix <- getNotRiboFastqFileNames( sampleID, asMatePairs=FALSE)

	# we do it all in one pass now
	nrHit <- nrNohit <- nrForward <- 0
	ans <- fastqToBAM( inputFastqFile, finalHit, k=maxMultiHits, sampleID=sampleID,
			noHitsFile=nohitPrefix, optionsFile=optionsFile, 
			alignIndex=alignIndex, alignPolicy=alignPolicy, 
			asMatePairs=asMatePairs, maxReads=maxFastqReads, verbose=verbose)

	nReadsIn <- ans$RawReads
	nrUnique <- ans$UniqueReads
	nrMulti <- ans$MultiReads
	nrNohit <- ans$NoHitReads
	nrAligned <- nrUnique + nrMulti

	if ( nReadsIn > 0) quickFileLineCountRecord( inputFastqFile, sampleID, lineCount=nReadsIn*4, readCount=nReadsIn)
	if ( nrNohit >= 0) quickFileLineCountRecord( finalNohit, sampleID, lineCount=nrNohit*4, readCount=nrNohit)
	if ( nrAligned >= 0) quickFileLineCountRecord( finalHit, sampleID, lineCount=nrAligned)

	# at this point we now know how many reads we started with...
	if ( is.null( rawReadCount)) rawReadCount <- nReadsIn

	# in the case there were no ribo alignments, redelete that BAM file
	if ( nrAligned < 1) file.delete( finalHit)


	# local function to print stats
	makeRiboSummary <- function() {

		out <- vector( mode="character")
		out <- base::append( out, "\n")
		out <- base::append( out, base::paste( "\nInput File:        \t", inputFastqFile))
		out <- base::append( out, base::paste( "\nRaw Reads:         \t", prettyNum( nReadsIn, format="d", 
				big.mark=",", width=12)))
		out <- base::append( out, base::paste( "\n\nAlignment Policy:  \t", alignPolicy))
		out <- base::append( out, base::paste( "\n\nNot Ribo File:      \t", finalNohit))
		out <- base::append( out, base::paste( "\nN_NotRibo Reads:   \t", prettyNum( nrNohit, format="d", 
				big.mark=",", width=12), "\t", as.percent( nrNohit, big.value=nReadsIn)))
		out <- base::append( out, base::paste( "\n\nHit Ribo File:  \t", finalHit))
		out <- base::append( out, base::paste( "\nN_Unique Reads:    \t", prettyNum( nrUnique, format="d", 
				big.mark=",", width=12), "\t", as.percent( nrUnique, big.value=nReadsIn)))
		out <- base::append( out, base::paste( "\nN_Multi Reads:     \t", prettyNum( nrMulti, format="d", 
				big.mark=",", width=12), "\t", as.percent( nrMulti, big.value=nReadsIn)))

		out <- base::append( out, base::paste( "\n\nTotal Cleared Reads:\t", prettyNum( nrAligned, format="d", 
				big.mark=",", width=12), "\t", as.percent( nrAligned, big.value=nReadsIn)))
		return( out)
	}  # end local function...

	myTime <- elapsedProcTime( startT, proc.time(), N=nReadsIn)
	cat( "\nFinished pipe 'RiboClear':    \t Results: ")
	outText <- makeRiboSummary()
	cat( outText)

	# the more detailed summary must wait till after genomic alignment is done
	allSummary <- outText
	summaryFile <- file.path( resultsPath, "summary", paste( sampleID, "ribo.Summary.txt", sep="."))
	writeLines( allSummary, con=summaryFile, sep="")

	cat( "\n\nTiming Stats: \n")
	print( myTime)
	gc()

	return( list( "RawReads"=nReadsIn, "UniqueReads"=nrUnique, "MultiReads"= nrMulti, "NoHitReads"=nrNohit, 
			"Time"=myTime))
}

