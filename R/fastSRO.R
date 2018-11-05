# fastSRO.R

fastSRO <- function( s1, s2, mismatch=(-6.0), gap=(-10.0)) {

	ans <- .C( "fastShortReadLocal", 
			as.integer(s1), 
			as.integer(length(s1)), 
			as.integer(s2), 
			as.integer(length(s2)), 
			as.double( mismatch),
			as.double( gap),
			score=double(1), offset=integer(1), gaps=integer(4))
	return( list( "score"=ans$score[1], "offset"=ans$offset[1], "gaps"=ans$gaps[1:4]))
}

