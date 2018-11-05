`rosettaSigma` <- function( nReads) {

	# sigma based on Rosetta error model of Poisson Distribution by John Castle
	# originally term 2 was (0.05 x nReads)^2, but that is underweighting sigma for large counts
	return( sqrt( nReads + (0.125*nReads)^2 + 25 ))
}

