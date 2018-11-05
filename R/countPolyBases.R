# countPolyBases.R 


`countPolyBases` <- function() {

	# see how many USR sequences are poly nucleotides...

	# the colors are integers 1..5 for ACGTN

	npoly <- rep( 0, times=5)
	bases <- c("A","C","G","T","N")
	names( npoly) <- paste( "poly", bases, sep="")
	
	if ( ! exists( "USR_seq")) stop( "Run 'USR_setup()' first...")

	minNpct <- 0.25
	minPct <- 0.8

	cat( "\n")
	cat( "\nCounting bases...")
	USR_colors <<- lapply( USR_seq, USR_BaseFactorsFunc)

	# now there are '.' as N's too, colorNumber is 7
	N <- length(USR_colors)
	pcts <- matrix( 0, nrow=5, ncol=N)
	i <- 0
	sapply( USR_colors, function(x) {
			nN <- sum( x == 5 | x == 7)
			nA <- sum( x == 1)
			nC <- sum( x == 2)
			nG <- sum( x == 3)
			nT <- sum( x == 4)
			i <<- i + 1
			pcts[ ,i] <<- c( nA, nC, nG, nT, nN) / length(x)
		})

	npoly[5] <- sum( pcts[ 5, ] >= minNpct)
	# these get counted when the USR object is made.. should be zero here
	#cat( "\nPoly N: ", npoly[5])

	for ( j in 1:4) {
		npoly[j] <- sum( pcts[ j, ] >= minPct)
		cat( "   Poly ", bases[j],": ", npoly[j], sep="")
	}

	return( list("nAs"=npoly[1], "nCs"=npoly[2], "nGs"=npoly[3], "nTs"=npoly[4], "nNs"=npoly[5]))
}
