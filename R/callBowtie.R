# callBowtie.R

# pass a Unix command line(s) to the System to run the Bowtie alignment program
`callBowtie` <- function( bowtieCommandLine, wait=TRUE, verbose=FALSE)
{
	if (verbose) {
		cat( "\n\nCalling Bowtie: \t", date());
		cat( "\nBowtie Version:   \t", getBowtieVersion());
		cat( "\nCommand Line:\n", bowtieCommandLine, "\n\n")
	}
	
	## Run the alignment program
	timeIn <- proc.time()
	ans = catch.system( command=bowtieCommandLine, wait=wait);
	timeOut <- proc.time()
	return( list( "timeIn"=timeIn, "timeOut"=timeOut))
}


`callBowtieBuild` <- function( bowtieCommandLine, wait=TRUE, verbose=FALSE)
{
	if (verbose) {
		cat( "\n\nCalling Bowtie-build: \t", date());
		cat( "\nCommand Line:\n", bowtieCommandLine, "\n\n")
	}
	
	## Run the alignment program
	timeIn <- proc.time()
	ans = catch.system( command= paste( bowtieCommandLine, "  2>&1 "), intern=TRUE, wait=wait);
	timeOut <- proc.time()

	if (verbose) {
		cat( ans, collapse="\n");
		cat( "\nBowtie-build done:    \t", date(), "\n")
	}
}

