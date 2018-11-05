# loadMultipleWIGs.R

# gather the compressed WIG datasets for multiple samples

`loadMultipleWIGs` <- function( sampleSet=NULL, speciesID=NULL, results.path=".", verbose=F) {

	wigPath <- file.path( results.path, "wig")

	## set up for this species
	if ( ! is.null( speciesID)) setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()

	# load the data...
	WIGlist <- vector( mode="list", length=length(sampleSet))
	if (verbose) cat( "\nLoading wiggle datasets...")
	
	for( i in 1:length(sampleSet)) {
		file <- paste( sampleSet[i], prefix, "WIG.rda", sep=".")
		file <- file.path( wigPath, file)
		if ( ! file.exists( file)) {
			cat( "\nWIG file not found: ", file)
			WIGlist[[ i]] <- NULL
		} else {
			if (verbose) cat( "  ", sampleSet[i])
			load( file)
			WIGlist[[i]] <- wiggles
		}
	}
	names( WIGlist) <- sampleSet
	if (verbose) cat( "  Done.\n")

	return( WIGlist)
}
