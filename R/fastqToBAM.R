# fastqToBAM.R

# turn a .fastq file into an aligned .bam file, by calling Bowtie2

`fastqToBAM` <- function( inputFastqFile, outputFile=sub( "(fastq|fq|gz)$", "bam", inputFastqFile[1]),  
		k=NULL, sampleID="", optionsFile="Options.txt",  annotationFile="Annotation.txt", noHitsFile=NULL, 
		alignIndex=getOptionValue( optionsFile, "GenomicIndex"), index.path=NULL,
		alignPolicy=getOptionValue( optionsFile, "GenomicAlignmentPolicy", verbose=F), 
		maxReads=NULL, skipReads=NULL, asMatePairs=FALSE, keepUnaligned=TRUE,
		wait=TRUE, verbose=FALSE, label=sampleID) {

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nCalling fastqToBAM:     ", label, "\n")
	}
	gc()


	# we can be given a vector of 2 fastq filenames...
	if ( outputFile == inputFastqFile[1]) stop( "FastqToBAM:  output file is same name as input file")

	# remove any old file copies...
	file.delete( outputFile)
	file.delete( noHitsFile)

	# make sure the input file is readable
	nReads <- getFileLineCount( inputFastqFile[1], sampleID=sampleID, mode="quick", verbose=verbose) / 4
	if ( nReads == 0) {
		warning( paste( "fastqToBAM:  input file contains no reads to align: ", inputFastqFile[1]))
		return( list( "RawReads"=0, "UniqueReads"=0, "MultiReads"=0,
			"NoHitReads"=0, "Time"=0))
	}

	if ( asMatePairs && !is.null(k)) {
		if ( k > 1) {
			cat( "\nMate Pair alignment not compatible with K > 1    Ignoring K\n")
			k <- NULL
		}
	}

	# turn all the arguments into one Bowtie command line
	metricsFile <- paste( sampleID, "Bowtie2.AlignMetrics.txt", sep=".")
	cmd <- buildBowtie2CommandLine( inputFastqFile=inputFastqFile, outputFile=outputFile, 
			optionsFile=optionsFile, annotationFile=annotationFile, metricsFile=metricsFile,
			k=k, noHitsFile=noHitsFile, 
			alignIndex=alignIndex, index.path=index.path, alignPolicy=alignPolicy, 
			maxReads=maxReads, skipReads=skipReads, asMatePairs=asMatePairs, 
			keepUnaligned=keepUnaligned, sampleID=sampleID, verbose=verbose)

	# Run the alignment program
	bowtieTiming <- callBowtie2( cmd, verbose=verbose)

	# retrieve the read counts
	ans <- getBowtie2AlignMetrics( metricsFile, pairedEnd=asMatePairs, bowtieTiming=bowtieTiming, verbose=verbose)
	nReads <- ans$RawReads
	nUnique <- ans$UniqueReads
	nMulti <- ans$MultiReads
	nNoHit <- ans$NoHitReads

	if (verbose) {
		cat( "\nfastqToBAM     \t Results: \n")
		cat( "\nRaw Reads:       \t", formatC( nReads, width=12, big.mark=",", format="d"), "\t", 
				basename(inputFastqFile))
		cat( "\nUniqueHit Reads: \t", formatC( nUnique, width=12, big.mark=",", format="d"), "\t", outputFile)
		cat( "\nMultiHit Reads:  \t", formatC( nMulti, width=12, big.mark=",", format="d"), "\t", outputFile)
		cat( "\nNoHit Reads:     \t", formatC( nNoHit, width=12, big.mark=",", format="d"), "\t", noHitsFile)
		cat( "\n", verboseOutputDivider)
	}
	gc()

	return( list( "RawReads"=nReads, "UniqueReads"=nUnique,
			"MultiReads"=nMulti, "NoHitReads"=nNoHit, "Time"=ans$Time))
}

