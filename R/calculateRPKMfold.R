`calculateRPKMfold` <-
function( rpkmA, rpkmB, minRPKM=1, clipFold=10.0) {

	# try to turn the fold to something sensible when the math goes funny...
	if ( is.na( rpkmA) || is.na( rpkmB)) return( NA)

	# use some minimum detectable RPKM
	minSeen <- min( rpkmA, rpkmB, minRPKM)
	if ( minSeen < minRPKM) {
		offset <- minRPKM - minSeen
		rpkmA <- rpkmA + offset
		rpkmB <- rpkmB + offset
	}

	rpkmFold <- log2( rpkmA/rpkmB)
	if ( rpkmFold > clipFold) rpkmFold <- clipFold
	if ( rpkmFold < (-clipFold)) rpkmFold <- (-clipFold)

	return( rpkmFold)
}

