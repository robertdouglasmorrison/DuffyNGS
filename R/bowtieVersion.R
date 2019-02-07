# bowtieVersion.R


`checkBowtieVersion` <- function() {

	# if no one has initialized any Bowtie actions, do that first
	bp <- bowtiePar()
	if ( ! length(bp)) {
		if ( file.exists( "./Options.txt")) {
			bowtiePar.defaults( optionsFile="./Options.txt")
		} else {
			bowtiePar.defaults()
		}
	}

	# see if this program is readable and executable
	myProgram <- bowtiePar( "Program")
	if ( ! (( file.access( myProgram, mode=0) == 0) &&  
		( file.access( myProgram, mode=1) == 0))) {
			cat( "\nCan not execute Bowtie version test.  Executable may be missing or not permissioned \n")
			return()
	}

	versionCallCommand <- paste( myProgram, "  --version")
	versionCall <- system( command=versionCallCommand, intern=TRUE)
	terms <- strsplit( versionCall[1], split=" ", fixed=TRUE)[[1]]

	callSuccess <- (terms[2] == "version")
	if ( ! callSuccess) {
		stop( paste( "checkBowtieVersion:  Error:  bowtie Program executable version detection failed. \n",
			"Filename as given:  ", bowtiePar("Program"), "\n",
			"Result: ", versionCall))
	}

	# defaults are based on 'current verison'...  if older, override to whats expected...
	thisVersion <- terms[3]
	assign( "CurrentVersion", value=thisVersion, envir=BowtieEnv)
	return( thisVersion)
}


`getBowtieVersion` <- function() return( get( "CurrentVersion", envir=BowtieEnv))

