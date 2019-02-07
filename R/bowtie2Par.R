# bowtie2Par.R

# a collection of routines to work with Bowtie2...


# get and set Bowtie2 default parameters
bowtie2Par <- function( ... ) {

	args <- list( ... )
	if ( length( args) < 1) {
		return( as.list( Bowtie2Env))
	}

	if ( is.null( names(args))) {
		out <- list()
		inEnv <- names( as.list( Bowtie2Env))
		for ( j in 1:length(args)) {
			argName <- args[[ j]]
			out <- base::append( out, Bowtie2Env[[ argName ]] )
			if( !( argName %in% inEnv)) warning( paste( "'", argName, "'",  "  is not a bowtie2Par()  parameter"))
		}
		names( out) <- base::unlist( args)
		return( if ( length(args) > 1) out else base::unlist(out))
	}

	for ( j in 1:length(args)) {
		nam <- names( args)[j]
		val <- args[j][[1]]
		names(val) <- names(nam) <- NULL
		assign( nam[[1]], val[[1]], envir=Bowtie2Env)
	}

	return()
}
		

bowtie2Par.defaults <- function( optionsFile="", verbose=TRUE) {

	if ( optionsFile == "") {
		if (verbose) cat( "\nSetting Bowtie2 Parameters to default values...")
		myProgram <- Sys.which( "bowtie2")
		programEnvVar <- Sys.getenv( "BOWTIE2_PROGRAM", unset=NA)
		if ( !is.na( programEnvVar)) myProgram <- programEnvVar

		myPath <- "~/NGS/Bowtie2Indexes/"
		pathEnvVar <- Sys.getenv( "BOWTIE2_INDEX_PATH", unset=NA)
		if ( !is.na( pathEnvVar)) myPath <- pathEnvVar
		
		myIndex <- "HsPf.genomic_idx"
		indexEnvVar <- Sys.getenv( "BOWTIE2_INDEX_FILE", unset=NA)
		if ( !is.na( indexEnvVar)) myIndex <- indexEnvVar
		
	} else {
		if ( verbose) cat( "\nSetting Bowtie2 Parameters from Options File:  \t", optionsFile)
		optT <- readOptionsTable( optionsFile)
		myProgram <- getOptionValue( optT, "bowtie2Program", verbose=verbose)
		myPath <- getOptionValue( optT, "bowtie2Index.path", verbose=verbose)
		myIndex <- getOptionValue( optT, "GenomicIndex", verbose=verbose)
	}

	if ( myProgram == "") {
		cat( "\nWarning:  failed to find/set Bowtie2 program.")
		cat( "\nTried:    Sys.which('bowtie2')  and .Sys.getenv('BOWTIE2_PROGRAM')")
	}

	assign( "Program", myProgram, envir=Bowtie2Env)
	assign( "IndexPath", myPath, envir=Bowtie2Env)
	assign( "GenomicIndex", myIndex, envir=Bowtie2Env)

	# based on the version, we may tweak some settings...
	if ( ! exists( "CurrentVersion", envir=Bowtie2Env)) {
		ver <- checkBowtie2Version()
		return( ver)
	} else {
		return( get( "CurrentVersion", envir=Bowtie2Env))
	}
}

