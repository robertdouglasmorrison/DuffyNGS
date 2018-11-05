# adapterTailSearch.R 

#  use the reads stored in the current USR (unique short read) objects
#  and test them for hits to adapter tails


`countAdapterHits` <- function(  adapters, minOverlap=20, logFreq= 50000, Nquit=1000000) {

	# the typical score is 2 per base, so turn length & overlap fraction into a score
	#len <- base::nchar( USR_seq[1])
	myMinScore <- 2 * minOverlap

	# count them
	ans <- adapterTailSearch( adapters, goodScore=myMinScore, 
			logFreq=logFreq, Nquit=Nquit)

	return( ans)
}



# wrapper to fast "Short Read Overlap" C code...
`fastSRO` <- function( s1, s2, mismatch=(-6.0), gap=(-10.0)) {

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


`fastSROscore` <- function( s1, s2, mismatch=(-6.0), gap=(-10.0)) {

	ans <- .C( "fastShortReadLocal", 
			as.integer(s1), 
			as.integer(length(s1)), 
			as.integer(s2), 
			as.integer(length(s2)), 
			as.double( mismatch),
			as.double( gap),
			score=double(1), offset=integer(1), gaps=integer(4))
	return( ans$score[1])
}


`adapterTailSearch` <- function( adapters=list( 
			"rc(FwdAdapt)"=myReverseComplement("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"),
			"rc(RevAdapt)"=myReverseComplement("CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT")),
			goodScore=24, seq=NULL, logFreq=50000, Nquit=1000000) {


	tails <- adapters
	nTails <- length( tails)

	# allow a one time pass if a sequence was given
	if ( ! is.null(seq)) {
		alignScore <- vector()
		for( j in 1:length( tails)) {
			alignScore[j] <- pairwiseAlignment( seq, tails[[j]], type="overlap", 
					gapOpen=(-15), scoreOnly=TRUE)
			cat( "\n", names(tails)[j], "\t", alignScore[j])
		}
		cat( "\n")
		names(alignScore) <- names(tails)
		return( alignScore)
	}

	cat( "\n\nSearching reads for Adapter Tail Hits\n")
	for( j in 1:nTails) cat( "\n\t", names(tails)[j], "\t", as.character( tails[[j]]))
	cat( "\n")

	cSum <- cumsum( USR_counts)
	names( cSum) <- NULL
	bigTotal <- sum( USR_counts)

	GAP_OPEN_COST <- (-30)
	MISMATCH_COST <- (-8)
	tailCounts <- rep( 0, time=length(tails))
	
	# turn the targets to color numbers
	tailColors <- lapply( tails, function(x) USR_BaseFactorsFunc( as.character(x), needSplit=TRUE))

	nToDo <- min( USR_nTotal, Nquit)

	# visit them all, in chunks
	for ( i in seq( 1, nToDo, by=logFreq)) {

		thisFrom <- i
		thisTo <- min( (i+logFreq-1), nToDo)
		thisWho <- thisFrom:thisTo
		thisChunk <- lapply( USR_seq[ thisWho], FUN=USR_BaseFactorsFunc) 

		for ( j in 1:nTails) {
			theseScores <- sapply( thisChunk, FUN=fastSROscore, tailColors[[j]], 
					mismatch=MISMATCH_COST, gap=GAP_OPEN_COST)
			yesHits <- which( theseScores >= goodScore)
			thisTotal <- sum( USR_counts[ thisWho[ yesHits]])
			tailCounts[j] <- tailCounts[j] + thisTotal
		}

		# watch...
		thisBigCount <- cSum[ thisTo]
		pct <- tailCounts * 100 / thisBigCount
		cat( "\n(nTested, %All,", 
			paste( rep( c("N","%"),each=nTails), rep( names(tails), times=2), sep="", collapse=", "),
			") ",
			prettyNum( thisBigCount, big.mark=","), as.percent( thisBigCount, big.value=bigTotal),
			tailCounts, formatC( pct, digits=2, format="f"))
	}

	names(tailCounts) <- names(tails)
	names(pct) <- names(tails)
	cat( "\n\nFinal Adapter Tail Hit Count:  \n")
	print( tailCounts)
	cat( "\nFinal Adapter Tail Hit Percentage:  \n")
	print( pct)

	return( list( "Adapters"=tails, "nTails"=tailCounts, "nTested"=thisBigCount, 
			"PctTails"=pct))
}

