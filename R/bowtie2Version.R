# bowtie2Version.R


`checkBowtie2Version` <- function() {

	# see if this program is readable and executable
	if ( ! (( file.access( bowtie2Par( "Program"), mode=0) == 0) &&  
		( file.access( bowtie2Par( "Program"), mode=1) == 0))) {
			cat( "\nCan not execute Bowtie2 version test.  Executable may be missing or not permissioned \n")
			return()
	}

	versionCallCommand <- paste( bowtie2Par("Program"), "  --version")
	versionCall <- system( command=versionCallCommand, intern=TRUE)
	terms <- strsplit( versionCall[1], split=" ", fixed=TRUE)[[1]]

	callSuccess <- (terms[2] == "version")
	if ( ! callSuccess) {
		stop( paste( "checkBowtie2Version:  Error:  bowtie2 Program executable version detection failed. \n",
			"Filename as given:  ", bowtie2Par("Program"), "\n",
			"Result: ", versionCall))
	}

	# defaults are based on 'current verison'...  if older, override to whats expected...
	thisVersion <- terms[3]
	assign( "CurrentVersion", value=thisVersion, envir=Bowtie2Env)
	#assign( "NoHitFileOption", value=" --un ", envir=Bowtie2Env)
	#assign( "MultiHitFileOption", value=" --max ", envir=Bowtie2Env)
	return( thisVersion)
}


`getBowtie2Version` <- function() return( get( "CurrentVersion", envir=Bowtie2Env))
