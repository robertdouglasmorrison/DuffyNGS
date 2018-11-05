# buildBowtieCommandLine.R


`buildBowtie1CommandLine` <- 
function( inputFastqFile, outputAlignFile=sub( "fastq$", "align", inputFastqFile),  
		m=1, k=NULL, 
		noHitsFile=NULL, multiHitsFile=NULL, optionsFile="Options.txt", 
		alignIndex=getOptionValue( optionsFile, "GenomicIndex"),
		alignPolicy=getOptionValue( optionsFile, "GenomicAlignmentPolicy"),
		maxReads=NULL, verbose=TRUE, debug=FALSE) {

	optT <- readOptionsTable( optionsFile)

	# force the look at options file to initialize bowtie.
	bowtiePar.defaults( optionsFile, verbose=verbose)

	## build the Bowtie alignment program's command line
	progName <- getOptionValue( optT, "bowtieProgram", verbose=verbose)
	if (debug) progName <- base::paste( progName, "-debug", sep="")
	out <- progName

	# any hard and fast rules that our package assumes go first...
	out <- base::paste( out, " -B 1 ")

	# next is "input options"
	out <- base::paste( out, getOptionValue( optT, "bowtieInputOptions", notfound="", verbose=verbose))

	# next is "reporting options"
	out <- base::paste( out, getOptionValue( optT, "bowtieReportingOptions", notfound="", verbose=verbose))

	# next is "performance options"
	out <- base::paste( out, getOptionValue( optT, "bowtiePerformanceOptions", notfound="", verbose=verbose))

	# next is explicit alignment policy
	out <- base::paste( out, alignPolicy)

	# explicit M value for 'maximum alignments' for each read
	if ( is.null(k)) {
		kValue <- as.numeric( getOptionValue( optionsFile, "maxMultiReport", notfound=m, verbose=verbose))
	} else {
		kValue <- k
	}
	if ( ! is.null(m)) {
		if ( kValue > m) kValue <- m
		out <- base::paste( out,  " -m", m, " -k", kValue, " ", sep=" ")
	} else {
		out <- base::paste( out,  " -k", kValue, " ", sep=" ")
	}

	# add in files to catch the non-aligned reads
	if ( ! is.null( noHitsFile)) {
		out <- base::paste( out,  bowtiePar("NoHitFileOption"), noHitsFile, " ", sep=" ")
	}
	if ( ! is.null( multiHitsFile)) {
		out <- base::paste( out,  bowtiePar("MultiHitFileOption"), multiHitsFile, " ", sep=" ")
	}
	
	# the index
	myIndexFile <- file.path( getOptionValue( optT, "bowtieIndex.path", notfound=".", verbose=verbose), alignIndex)

	# test to see that a file is really there...
	if ( ! file.exists(  base::paste( myIndexFile, ".1.ebwt", sep=""))) {
		stop( paste( "buildBowtieCommandLine:  Error:  bowtie index Path and/or File not found. \n",
			"Filename as given:  ", myIndexFile))
	}
	out <- base::paste( out, myIndexFile, " ", sep="  ")
	
	# the input file
	inputFileTerm <- inputFastqFile
	if ( regexpr( "gz$", inputFastqFile) > 0 ) {
		GUNZIP <- TRUE
		cat( "\n\nUsing 'gunzip' to extract reads...\n")
		zipCmd <- base::paste( "gunzip -c ", inputFastqFile, "  | ")
		inputFileTerm <- " - "
	} else {
		GUNZIP <- FALSE
		zipCmd <- ""
	}

	if ( ! is.null( maxReads)) {
		cat( "\nLimiting Input to the first ", prettyNum( maxReads, big.mark=","), " reads.\n")
		# remember 4 lines per read
		maxReads <- as.integer( maxReads * 4)
		if (GUNZIP) {
			headCmd <- base::paste( "  head -n", maxReads, "| ", sep="  ")
		} else {
			headCmd <- base::paste( "head -n", maxReads, inputFastqFile, "| ", sep="  ")
		}
		inputFileTerm <- " - "
	} else {
		headCmd <- ""
	}
	# OK, with all pipes built if we needed any, we now know what to use for the input 
	out <- base::paste( out, inputFileTerm, " ", sep="  ")

	# lastly, the output file
	out <- base::paste( out, outputAlignFile, sep="  ")

	# pre-pend all the special bits...
	out <- base::paste( zipCmd, headCmd, out, sep="  ")

	return( out)
}


`buildBowtie1BuildCommandLine` <- 
function( inputFastaFile, outputIndexFile, optionsFile="Options.txt", verbose=TRUE, debug=FALSE) {

	optT <- readOptionsTable( optionsFile)

	# force the look at options file to initialize bowtie.
	bowtiePar.defaults( optionsFile, verbose=verbose)

	## build the Bowtie alignment program's command line
	progName <- getOptionValue( optT, "bowtieProgram", verbose=verbose)
	progName <- base::paste( progName, "-build", sep="")
	if (debug) progName <- base::paste( progName, "-debug", sep="")
	out <- progName

	out <- base::paste( out, " -f ")

	out <- base::paste( out, inputFastaFile, outputIndexFile, sep="  ")
	
	return( out)
}
