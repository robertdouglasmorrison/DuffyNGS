# expressionUnits.R -- decide what units for gene expression data to use  RPKM, TPM, etc


# get the units for expression from the Options file
`expressionUnits` <- function( optionsFile="Options.txt", verbose=TRUE) {

	knownUnits <- c( "RPKM", "READS", "TPM")
	defaultUnits <- "RPKM"

	# default has always be RPKM.  Allow reading from the Options file
	if ( ! file.exists( optionsFile)) {
		if ( verbose) {
			cat( "\nOptions file not found:  ", optionsFile)
			cat( "\nUsing default expression units:  ", defaultUnits)
		}
		return( defaultUnits)
	}

	units <- getOptionValue( optionsFile, arg="expression.units", notfound="RPKM", verbose=verbose)

	# force it to be something we know
	units <- toupper( units)
	if ( ! (units %in% knownUnits)) {
		cat( "\nInvalid expression units: ", units, "\tMust be one of: ", knownUnits)
		return( defaultUnits)
	}

	return( units)
}


# create the text string of the expression column name we will use for expression
`getExpressionUnitsColumn` <- function( optionsFile="Options.txt", useMultiHits=TRUE, verbose=TRUE) {

	units <- expressionUnits( optionsFile, verbose=verbose)
	suffix <- if (useMultiHits) "M" else "U"
	return( paste( units, suffix, sep="_"))
}

