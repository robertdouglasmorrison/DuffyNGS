# pipe.USR.R

# standalone version of building Unique Short Reads data structure, if need to
# re-assess without a full pipeline alignment...


`pipe.USR` <- function( sampleID, annotationFile="Annotation.txt",
		optionsFile="Options.txt", resultsPath=NULL, 
		mode=c("normal", "QuickQC"), 
		trim5=as.numeric( getOptionValue( optionsFile, "trim5", notfound=0)), 
		trim3=as.numeric( getOptionValue( optionsFile, "trim3", notfound=0)),
		nUSRkeep=1000000, verbose=TRUE) {

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\n\nBuilding 'USR' object for Sample:     ", sampleID, "\n")
	}

	mode <- match.arg( mode)

	if( is.null(resultsPath)) {
		resultsPath <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=verbose)
		if ( mode == "QuickQC") resultsPath <- "./QuickQC"
	}

	gc()
	
	sample <- sampleID

	# get the file of noHit reads
	infile <- paste( sample, "noHits", "fastq", sep=".") 
	if ( mode == "QuickQC") infile <- paste( sample, "QuickQC.noHits", "fastq", sep=".") 
	infile <- file.path( resultsPath, "fastq", infile)
	infile <- allowCompressedFileName( infile)
	if ( !file.exists( infile) && mode == "normal") {
		notgenomicfile <- paste( sample, "not.genomic", "fastq", sep=".") 
		notgenomicfile <- file.path( resultsPath, "fastq", notgenomicfile)
		cat( "\nFile not found:  ", infile, "\nTrying: ", notgenomicfile)
		infile <- allowCompressedFileName( notgenomicfile)
	}
	if ( ! file.exists( infile)) {
		stop( paste( "\nUnable to find 'NoHits' file:   ", infile))
	}
	nNoHits <- getFileLineCount( infile, sample) / 4


	# build that USR data structure
	ansUSR <- USR_setup( infile, sample, resultsPath, Nkeep=nUSRkeep, trim5=trim5, trim3=trim3)
	USRfile <- ansUSR$USR_File

	# now look for adapter tails
	ans <- noHit.adapterTails( sample, annotationFile, optionsFile, nRawReads=NULL)
	adapterCounts <- ans$adapterCounts

	# stash adapter tail info with the USR, upate and then we're done with USR data
	USR_AdapterHits <<- adapterCounts
	saveUSRcontext( USRfile)
	USR_cleanup()

	#}  # end of each pair of the sampleID

	gc()

	return( ansUSR)
}
