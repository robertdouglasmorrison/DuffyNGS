# buildBowtie2CommandLine.R


`buildBowtie2CommandLine` <- 
function( inputFastqFile, outputFile=sub( "fastq$", "bam", inputFastqFile[1]),  
		optionsFile="Options.txt", metricsFile="bowtie2Metrics.txt",
		k=NULL, asMatePairs=FALSE, noHitsFile=NULL, 
		alignIndex=getOptionValue( optionsFile, "GenomicIndex", verbose=F), index.path=NULL,
		alignPolicy=getOptionValue( optionsFile, "GenomicAlignmentPolicy", verbose=F),
		maxReads=NULL, skipReads=NULL, keepUnaligned=FALSE, quiet=FALSE, 
		verbose=TRUE, debug=FALSE) {

	optT <- readOptionsTable( optionsFile)

	mypaste <- base::paste

	# force the look at options file to initialize bowtie.
	bowtie2Par.defaults( optionsFile, verbose=verbose)

	## build the Bowtie alignment program's command line
	progName <- getOptionValue( optT, "bowtie2Program", verbose=verbose)
	if (debug) progName <- mypaste( progName, "-debug", sep="")
	out <- progName

	# next is "input options"
	out <- mypaste( out, getOptionValue( optT, "bowtie2InputOptions", notfound="", verbose=verbose))

	# next is any hard trimming
	trim5 <- as.integer( getOptionValue( optT, "trim5", notfound="0", verbose=verbose))
	trim3 <- as.integer( getOptionValue( optT, "trim3", notfound="0", verbose=verbose))
	if ( trim5 > 0) out <- mypaste( out, " --trim5", trim5)
	if ( trim3 > 0) out <- mypaste( out, " --trim3", trim3)

	# next is explicit counts of max reads or skipping reads
	if ( !is.null( skipReads)) out <- mypaste( out, " --skip", as.integer( skipReads))
	if ( !is.null( maxReads)) out <- mypaste( out, " --qupto", as.integer( maxReads))
	# next is explicit alignment policy
	out <- mypaste( out, alignPolicy)

	# next is "scoring options"
	out <- mypaste( out, getOptionValue( optT, "bowtie2ScoringOptions", notfound="", verbose=verbose))

	# explicit K value for 'maximum alignments' for each read
	if ( ! is.null(k) && k > 1) {
		out <- mypaste( out,  " -k", k)
	}

	if ( quiet) out <- mypaste( out, " --quiet")

	# next is "pair options"
	out <- mypaste( out, getOptionValue( optT, "bowtie2PairOptions", notfound="", verbose=verbose))

	# add in files to catch the non-aligned reads
	if ( ! is.null( noHitsFile)) {
		if ( regexpr( "\\.gz$", noHitsFile) < 1) noHitsFile <- paste( noHitsFile, "gz", sep=".")
		unFlag <- " --un-gz "
		if (asMatePairs) unFlag <- " --un-conc-gz "
		out <- mypaste( out, unFlag, noHitsFile)

		# if catching unaligned, do not send them to the SAM file too
		out <- mypaste( out, "  --no-unal ")
	} else {
		# when not catching the unaligned, allow explicit throw-away of unaligned
		if ( ! keepUnaligned) {
			out <- mypaste( out, "  --no-unal ")
		}
	}
	
	# next is threads
	nCores <- as.integer( getOptionValue( optT, "nCores", notfound="4", verbose=verbose))
	out <- mypaste( out, " --threads", nCores)

	# next is "other options"
	out <- mypaste( out, getOptionValue( optT, "bowtie2OtherOptions", notfound="", verbose=verbose))

	# the index
	if ( is.null( index.path)) index.path <- getOptionValue( optT, "bowtie2Index.path", notfound=".", verbose=verbose)
	myIndexFile <- file.path( index.path, alignIndex)

	# test to see that a file is really there...
	if ( ! file.exists(  mypaste( myIndexFile, ".1.bt2", sep=""))) {
	    # new long indexex are possible
	    if ( ! file.exists(  mypaste( myIndexFile, ".1.bt2l", sep=""))) {
		stop( paste( "buildBowtie2CommandLine:  Error:  bowtie index Path and/or File not found. \n",
			"Filename as given:  ", myIndexFile))
	    }
	}
	out <- mypaste( out, " -x", myIndexFile)
	
	# there MAY be pipes, etc. at the front...
	zipCmd <- headCmd <- ""

	# the input file
	if ( (length( inputFastqFile) == 2) && asMatePairs) {
		inputFileTerm <- paste( " -1", inputFastqFile[1], " -2", inputFastqFile[2])
	} else {
		inputFileTerm <- paste( " -U",  paste(inputFastqFile,collapse=","))
	}
	
	# OK, with all pipes built if we needed any, we now know what to use for the input 
	out <- mypaste( out, inputFileTerm)

	# redirect stderr to catch the reported read counts, etc.
	out <- mypaste( out, " 2> ", metricsFile, " ")

	# lastly, the output destination
	outputTerm <- paste ( " | samtools view -bS -o", outputFile, " - ")
	out <- mypaste( out, outputTerm, sep="  ")

	# pre-pend all the special bits...
	out <- mypaste( zipCmd, headCmd, out)

	return( out)
}


`buildBowtie2BuildCommandLine` <- 
function( inputFastaFile, outputIndexFile, optionsFile="Options.txt", verbose=TRUE, debug=FALSE) {

	optT <- readOptionsTable( optionsFile)

	# force the look at options file to initialize bowtie.
	bowtie2Par.defaults( optionsFile, verbose=verbose)

	## build the Bowtie alignment program's command line
	progName <- getOptionValue( optT, "bowtie2Program", verbose=verbose)
	progName <- base::paste( progName, "-build", sep="")
	if (debug) progName <- base::paste( progName, "-debug", sep="")
	out <- progName

	out <- base::paste( out, " -f ", paste(inputFastaFile,collapse=","), " ")

	out <- base::paste( out, outputIndexFile, sep=" ")
	
	return( out)
}
