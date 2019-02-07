# bowtiePar.R

# a collection of routines to work with Bowtie2...


# get and set Bowtie2 default parameters
bowtiePar <- function( ... ) {

	args <- list( ... )
	if ( length( args) < 1) {
		return( as.list( BowtieEnv))
	}

	if ( is.null( names(args))) {
		out <- list()
		inEnv <- names( as.list( BowtieEnv))
		for ( j in 1:length(args)) {
			argName <- args[[ j]]
			out <- base::append( out, BowtieEnv[[ argName ]] )
			if( !( argName %in% inEnv)) warning( paste( "'", argName, "'",  "  is not a bowtiePar()  parameter"))
		}
		names( out) <- base::unlist( args)
		return( if ( length(args) > 1) out else base::unlist(out))
	}

	for ( j in 1:length(args)) {
		nam <- names( args)[j]
		val <- args[j][[1]]
		names(val) <- names(nam) <- NULL
		assign( nam[[1]], val[[1]], envir=BowtieEnv)
	}

	return()
}
		

bowtiePar.defaults <- function( optionsFile="", verbose=TRUE) {

	if ( optionsFile == "") {
		if (verbose) cat( "\nSetting Bowtie Parameters to default values...")
		myProgram <- Sys.which( "bowtie")
		programEnvVar <- Sys.getenv( "BOWTIE_PROGRAM", unset=NA)
		if ( !is.na( programEnvVar)) myProgram <- programEnvVar

		myPath <- "~/NGS/BowtieIndexes/"
		pathEnvVar <- Sys.getenv( "BOWTIE_INDEX_PATH", unset=NA)
		if ( !is.na( pathEnvVar)) myPath <- pathEnvVar
		
		myIndex <- "HsPf.genomic_idx"
		indexEnvVar <- Sys.getenv( "BOWTIE_INDEX_FILE", unset=NA)
		if ( !is.na( indexEnvVar)) myIndex <- indexEnvVar
		
	} else {
		if ( verbose) cat( "\nSetting Bowtie Parameters from Options File:  \t", optionsFile)
		optT <- readOptionsTable( optionsFile)
		myProgram <- getOptionValue( optT, "bowtieProgram", verbose=verbose)
		myPath <- getOptionValue( optT, "bowtieIndex.path", verbose=verbose)
		myIndex <- getOptionValue( optT, "GenomicIndex", verbose=verbose)
	}

	if ( myProgram == "") {
		cat( "\nWarning:  failed to find/set Bowtie program.")
		cat( "\nTried:    Sys.which('bowtie')  and .Sys.getenv('BOWTIE_PROGRAM')")
	}
	assign( "Program", myProgram, envir=BowtieEnv)
	assign( "IndexPath", myPath, envir=BowtieEnv)
	assign( "GenomicIndex", myIndex, envir=BowtieEnv)

	# based on the version, we may tweak some settings...
	if ( ! exists( "CurrentVersion", envir=BowtieEnv)) {
		ver <- checkBowtieVersion()
		return( ver)
	} else {
		return( get( "CurrentVersion", envir=BowtieEnv))
	}
}

