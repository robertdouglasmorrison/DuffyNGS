# callBowtie2.R

# pass a Unix command line(s) to the System to run the Bowtie2 alignment program
`callBowtie2` <- function( bowtie2CommandLine, wait=TRUE, verbose=FALSE)
{
	if (verbose) {
		cat( "\n\nCalling Bowtie2: \t", date());
		cat( "\nBowtie2 Version:   \t", getBowtie2Version());
		cat( "\nCommand Line:\n", bowtie2CommandLine, "\n\n")
	}
	
	## Run the alignment program
	timeIn <- proc.time()
	ans = catch.system( command=bowtie2CommandLine, wait=wait);
	timeOut <- proc.time()
	return( list( "timeIn"=timeIn, "timeOut"=timeOut))
}


`callBowtie2Build` <- function( bowtie2CommandLine, wait=TRUE, verbose=FALSE)
{
	if (verbose) {
		cat( "\n\nCalling Bowtie2-build: \t", date());
		cat( "\nCommand Line:\n", bowtie2CommandLine, "\n\n")
	}
	
	## Run the alignment program
	timeIn <- proc.time()
	ans = catch.system( command= paste( bowtie2CommandLine, "  2>&1 "), intern=TRUE, wait=wait);
	timeOut <- proc.time()

	if (verbose) {
		cat( ans, collapse="\n");
		cat( "\nBowtie2-build done:    \t", date(), "\n")
	}
}

