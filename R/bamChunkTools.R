# bamChunkTools.R -- helper functions that operate on a "chunk" from a BAM file


`is.bamChunk` <- function( chunk) {

	if ( class(chunk) != "bamRange") stop( "argument 1 must be a 'BAM chunk' object")

	# these functions assume the reads are in the raw order coming off the alignment tool

	# so try to catch any sorted bamRange object as invalid
	Ntest <- min( 1000, size(chunk))
	if (Ntest < 1) return( FALSE)

	# method 1:  the alignments should not be in order
	refIDs <- refID( chunk)[1:Ntest]
	pos <- position( chunk)[1:Ntest]
	ord <- order( refIDs, pos)
	if ( any( ord != 1:Ntest)) return( TRUE)
	return( FALSE)
}


`weight.align` <- function( chunk) {

	if ( ! is.bamChunk( chunk)) return(NULL)
	N <- size( chunk)
	if ( N == 1) return(1)

	# we can turn where the primary aligns are into weights, by using the distance between them
	primary <- which( ! secondaryAlign( chunk))
	pts <- if ( primary[1] == 1) c( primary, N+1) else c( 1, primary, N+1)
	difs <- diff( pts)
	cnts <- rep( difs, times=difs)
	return( rep.int( 1,N) / cnts)
}


`factor.align` <- function( chunk) {

	if ( ! is.bamChunk( chunk)) return(NULL)
	N <- size( chunk)

	# we can turn where the primary aligns are into factors, by using the distance between them
	primary <- which( ! secondaryAlign( chunk))
	pts <- if ( primary[1] == 1) c( primary, N+1) else c( 1, primary, N+1)
	difs <- diff( pts)
	levs <- 1:length(difs)
	fac <- rep( levs, times=difs)
	levels(fac) <- levs
	class(fac) <- "factor"
	return( fac)
}


`which.unique.align` <- function( chunk) {

	if ( ! is.bamChunk( chunk)) return(NULL)
	return( which( weight.align( chunk) == 1.0))
}


`which.multi.align` <- function( chunk) {

	if ( ! is.bamChunk( chunk)) return(NULL)
	return( which( weight.align( chunk) < 1.0))
}


`which.multi.read` <- function( chunk) {

	if ( ! is.bamChunk( chunk)) return(NULL)
	isMulti <- which( weight.align( chunk) < 1.0)
	isPrimary <- which( ! secondaryAlign( chunk))
	return( intersect( isMulti, isPrimary))
}
