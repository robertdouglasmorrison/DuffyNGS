# CR.R - Consensus Reads



# constants for the DNA bases
N_BASES <- 6
BASES <- c( "A", "C", "G", "T", "N", "-")
BASE_T <- 4
BASE_N <- 5
BASE_GAP <- 6
BASE_COLORS <- c( "red", "blue", "orange", "green", "brown", "black")
GAP_OPEN_COST <- (-30)
MISMATCH_COST <- (-12)


autoRunCR <- function( nBest=5, maxCR=4000, nIterations=1000, maxTime=1000, maxCycles=10, 
			ratePerCycle=NULL, pause=0, makePlots=TRUE,
			contextFile="CRTcontext.rda", pngPath=".", label="") {

	if ( ! exists( "USR_nTotal")) stop( "No USR data is loaded...")

	myReadLen <- base::nchar( USR_seq[1])
	greatScore <- curScore <- round( myReadLen * 1.6)
	if ( is.null(ratePerCycle)) ratePerCycle <- max( 1, round(greatScore*0.025))
	okScore <- round( myReadLen * 1.2)
	deltaTime <- 0
	nCycles <- 0
	whoToPlot <- vector()

	while ( deltaTime < maxTime && curScore >= okScore && nCycles < maxCycles) {
	
		firstPass <- (deltaTime == 0)
		ans <- runCR( curScore, maxCR=maxCR, restart=firstPass, iterations=nIterations, 
				contextFile=contextFile, max.time=maxTime)
		nCycles <- nCycles + 1

		cat( "\nCycle: ", nCycles, "\tN_Combines: ", ans$Combos, "\tElased time: ", 
				round( ans$Time), "seconds.")
		if ( ans$Combos < nIterations) cat("\nEarly exit of 'runCR'...")
		deltaTime <- ans$Time

		cat( "\nTop CR's by length: ")
		whoLen <- CRplotter( seconds=NULL, plotOrder="length", N=nBest, doPlot=FALSE)
		cat( "   ", whoLen)
		cat( "\nTop CR's by read count: ")
		whoCnt <- CRplotter( seconds=NULL, plotOrder="count", N=nBest, doPlot=FALSE)
		cat( "   ", whoCnt)
		whoToPlot <- base::sort( base::union( whoLen, whoCnt))

		# save this set of the best consensus reads every time thru
		CRT_best <<- whoToPlot
		saveCRTcontext( contextFile, verbose=FALSE)

		CRplotter( seconds=pause, plotOrder="", doPlot=makePlots, whoToPlot=whoToPlot, 
				asPNG=T, pngPath=pngPath, label=label)

		# make it a bit easier, and repeat
		curScore <- ans$CurrentScore - ratePerCycle
	}

	# save this set of the best consensus reads
	CRT_best <<- whoToPlot

	return( whoToPlot)
}
	

# 
# run the Consensus Reads finder
#
runCR <- function( goodScore=40, iterations=1000, restart=TRUE, maxCR=4000, contextFile="CRTcontext.rda",
			max.time=1000) {

	if ( restart) {
		newCRT( maxCR)
		cat( "\n")
	}
	nHit <- nNoLuck <- 0
	newHit <- FALSE
	grow <- TRUE
	startTime <- proc.time()
	
	# score for combining are based on what we call 'good enough'
	CRT_combineCutoff <<- goodScore * 1.00
	CRT_reAlignCutoff <<- goodScore / 2
	CRT_pairAlignCutoff <<- goodScore / 2
	stallLimit <- iterations / 2

	if ( ! exists( "USR_nTotal")) stop( "No USR data is loaded...")

	cat( "\n\nSearching for Consensus Reads...  Good Score=", goodScore, "\n")

	grow <- ( USR_nUsed < USR_nTotal)
	while( TRUE) {

		# grow the set of Consensus Reads by one
		# if room, add one; otherwise loosen match criteria
		

		if (grow) {
		    if ( CRT_nCR < CRT_Size && USR_nUsed < USR_nTotal) {
			addOneUSRtoCRT()
		    } else {
			CRT_combineCutoff <<- CRT_combineCutoff - 1
			if ( CRT_combineCutoff < goodScore) break
		    }
		}

		if ( newHit && nHit %% 100 == 0) {
			cat("\n(nUSR, nCR, nPairs, Score, nCombo, %All)", USR_nUsed, 
					CRT_nCR, CRT_nPairs, CRT_combineCutoff, nHit, 
					as.percent( sum(USR_counts[1:USR_nUsed]), big.value=USR_GrandTotal)) 
			# allow a loosening on the fly
			if ( USR_nUsed > 1000 &&  USR_nUsed/CRT_nCR < 1.5) {
				CRT_combineCutoff <<- CRT_combineCutoff - 1
			}
		}

		if ( nHit >= iterations) break
		if ( nNoLuck >= stallLimit) break
		
		# find the best now
		bestHit <- getBestCRTalignment()

		if ( bestHit$score >= CRT_combineCutoff) {
			grow <- FALSE
			who <- combineAndReduceCRT( bestHit)
			if ( who  <= 0) next
			newHit <- TRUE
			nHit <- nHit + 1
			nNoLuck <- round( nNoLuck / 2)
			dtime <- elapsedProcTime( startTime, proc.time(), N=nHit, what="")
			if ( dtime$Raw.Seconds >= max.time) break
		} else {
			newHit <- FALSE
			grow <- TRUE
			nNoLuck <- nNoLuck + 1
			if ( nNoLuck %% 100 == 0) cat( "  Stall", nNoLuck)
		}

	}
	saveCRTcontext( file=contextFile)

	stopTime <- proc.time()
	dtime <- elapsedProcTime( startTime, stopTime, N=nHit, what="combinations")
	out <- list( "Combos"=nHit, "Time"=dtime$Raw.Seconds, "CurrentScore"=CRT_combineCutoff)
	return( out)
}



#
# high level self test
#
testAll <- function() {

	# not much in the USR needs testing...
	ok <- TRUE
	if ( USR_nUsed < 0 || USR_nUsed > USR_nTotal) ok <- FALSE
	if ( ! ok) cat( "\nUSR test failed")

	# most everything is in the CRT
	if ( testCRT() == FALSE) ok <- FALSE
	if ( ! ok) cat( "\nCRT test failed")
	return( ok)
}


#
#	CR - Consensus Read
#
#	Growable Meta-read that encompasses all Unique Short Reads with 'good' 
#	between-read alignment
#
newCR <- function( str="", cnt=1, USRid=1, CRTid=1) {
	
	# make a new CR object from a character string of DNA
	emptyCR <- list( "len"=0)
	str <- toupper( str[1] )
	len <- base::nchar( str)
	if ( len == 0) return( emptyCR)
	bases <- strsplit( str, split="")[[1]]
	badBases <- bases[ which(  ! ( bases %in% BASES))]
	if ( length( badBases) > 0) {
		cat( "\nnewCR:  non-DNA bases: ", badBases)
		return( emptyCR)
	}
	dna <- DNAString( str)
	depth <- rep( 1, times=len)
	scores <- rep( 100, times=len)
	bcm <- newBCM( bases, cnt)
	thisID <- USRid
	thisPosition <- 1
	thisGaps <- thisColors <- vector( mode="list", length=1)
	thisGaps[[1]] <- c(0,0)
	if ( USRid > 0) {
		thisColors[[1]] <- getUSRentry( USRid)$colors
	} else {
		thisColors[[1]] <- USR_BaseFactorsFunc( bases, needSplit=FALSE)
	}
	names( thisGaps)[1] <- names( thisColors)[1] <- USRid

	return( list( "len"=len, "bases"=bases, "scores"=scores, "dna"=dna, "BCM"=bcm, "depth"=depth, 
			"baseIDs"=USR_BaseFactorsFunc( bases, needSplit=FALSE), 
			"USRid"=thisID, "USRposition"=thisPosition, "USRgaps"=thisGaps, 
			"USRcolors"=thisColors, "CRTid"=CRTid, "totReads"=cnt ))
}


testCR <- function( cr) {
	ok <- TRUE
	if ( ! ("len" %in% names(cr))) return( FALSE)
	if ( is.na( cr$len)) ok <- FALSE
	if ( cr$len < 0) ok <- FALSE
	if ( ! ok) cat( "\nCR length test failed")
	if ( cr$len == 0) {
		cat("\ntesting an unused CR...")
		return(FALSE)
	}
	if ( ! all( names( cr) == c( "len", "bases", "scores", "dna", "BCM", "depth", "baseIDs",
			"USRid", "USRposition", "USRgaps", "USRcolors", "CRTid", "totReads" ))) ok <- FALSE
	if ( ! ok) cat( "\nCR list names test failed")
	if ( length( cr$bases) != cr$len) ok <- FALSE
	if ( base::nchar( cr$dna) != cr$len) ok <- FALSE
	if ( ncol( cr$BCM) != cr$len) ok <- FALSE
	if ( ! ok) cat( "\nCR seq length test failed")
	if ( max( cr$depth) > length( cr$USRid)) ok <- FALSE
	if ( any( cr$USRposition < 1 | cr$USRposition > cr$len)) ok <- FALSE
	if ( ! ok) cat( "\nCR USRpositions test failed")
	for ( i in 1:length(cr$USRid)) {
		thisPos <- cr$USRposition[i]
		thiscolorset <- cr$USRcolors[[i]]
		if ( (thisPos + length(thiscolorset) - 1) > cr$len) ok <- FALSE
	}
	if ( ! ok) cat( "\nCR USRcolor lengths test failed")
	if ( ! testBCM( cr$BCM)) ok <- FALSE
	if ( ! ok) cat( "\nCR BCM test failed")
	ans <- baseConsensusFromBCM( cr$BCM)
	if ( any( ans$bases != cr$bases)) ok <- FALSE
	if ( any( ans$scores != cr$scores)) ok <- FALSE
	if ( ! ok) cat( "\nCR consensus seq test failed")
	if ( ! ok) cat( " \t CRid= ", cr$CRTid)
	# looks good!!
	return( ok)
}


revCompCR <- function( cr) {

	# turn this into its reverse complement
	newcr <- cr
	newcr$bases <- strsplit( myReverseComplement( paste( cr$bases, collapse="")), split="")[[1]]
	newcr$scores <- rev( cr$scores)
	newcr$dna <- reverseComplement( cr$dna)
	newcr$BCM <- revCompBCM( cr$BCM)
	newcr$depth <- rev( cr$depth)
	newcr$baseIDs <- USR_BaseFactorsFunc( newcr$bases, needSplit=FALSE) 
	for ( k in 1:length( newcr$USRcolors)) {
		newcr$USRcolors[[k]] <- revCompCRcolors( newcr$USRcolors[[k]])
	}

	return( newcr)
}


revCompCRcolors <- function(x) {

	# do reverse complement on the the color points 1:4 (A:T)
	newx <- revx <- rev(x)
	newx[ revx %in% 1:6] <- c( 4,3,2,1,5,6)[ match( revx, 1:6)]
	return( newx)
}


gapShiftCR <- function( cr, gapDesc) {

	# apply a gap shift to a CR because an alignment says it needs one.
	if ( gapDesc[1] <= 1) return( cr)
	if ( gapDesc[2] == 0) return( cr)

	gapStart <- gapDesc[1]
	gapLength <- gapDesc[2]

	oldlen <- cr$len
	oldbcm <- cr$BCM

	# change the vector of bases
	newbases <- c( cr$bases[1:(gapStart-1)], rep( "-", times=gapLength), cr$bases[ gapStart:oldlen])
	newlen <- length( newbases)
	if ( newlen != (oldlen+gapLength)) {
		cat( "\n\ngapShift Dump:\n(gapS, gapL, oldL,newL,CRid)", gapStart, gapLength,
			oldlen,newlen, cr$CRTid,"\n")
		stop("bad math...")
	}

	# change the BCM matrix
	newbcm <- newBCM( newbases, cr$totReads)
	newbcm[ , 1:(gapStart-1)] <- oldbcm[ , 1:(gapStart-1)]
	newbcm[ , (gapStart+gapLength):ncol(newbcm)] <- oldbcm[ , gapStart:oldlen]

	# rebuild the scores and DNA string
	ans <- baseConsensusFromBCM( newbcm)
	newdna <- DNAString( base::paste( newbases, collapse=""))
	newscores <- ans$scores

	# adjust the depth vector
	newdepth <- c( cr$depth[1:(gapStart-1)], rep( max(cr$depth[gapStart:(gapStart+gapLength-1)]), 
			times=gapLength), cr$depth[ gapStart:oldlen])

	# adjust the gaps list and base colors and positions
	newposition <- cr$USRposition
	newgaplist <- cr$USRgaps
	newcolorslist <- cr$USRcolors
	for ( i in 1:length( cr$USRid)) {
		# the start position of each USR may change...
		if ( newposition[i] >= gapStart) newposition[i] <- newposition[i] + gapLength

		# if no gaps yet, this is first, else extend the set...
		tmpgap <- newgaplist[[i]]
		if ( tmpgap[1] == 0) {
			tmpgap <- c( gapStart, gapLength)
		} else {
			tmpgap <- base::append( tmpgap, c( gapStart, gapLength))
		}
		newgaplist[[i]] <- tmpgap

		# each color vector may be in a different relative location for this gap
		thisNewColors <- oldColors <- newcolorslist[[i]]
		lenOC <- length( oldColors)
		myGapStart <- gapStart - newposition[i] + 1
		if ( (myGapStart > 1) && (myGapStart <= lenOC)) {
			thisNewColors <- c( oldColors[1:(myGapStart-1)], rep( BASE_GAP, times=gapLength), 
						oldColors[ myGapStart:lenOC])
		}
		newcolorslist[[i]] <- thisNewColors
	}

	return( list( "len"=newlen, "bases"=newbases, "scores"=newscores, "dna"=newdna, 
			"BCM"=newbcm, "depth"=newdepth, 
			"baseIDs"=USR_BaseFactorsFunc( newbases, needSplit=FALSE), 
			"USRid"=cr$USRid, "USRposition"=newposition, "USRgaps"=newgaplist, 
			"USRcolors"=newcolorslist, "CRTid"=cr$CRTid, "totReads"=cr$totReads ))
}


combine2CR <- function( cr1, cr2, offset=0, gaps=NULL) {
	
	# join 2 CRs that have been pre-selected for overlap into 1:
			
	# watch for gapping found by alignment...
	if ( ! is.null( gaps)) {
		if ( gaps[1,1] > 1) cr1 <- gapShiftCR( cr1, gaps[ 1, ])
		if ( gaps[2,1] > 1) cr2 <- gapShiftCR( cr2, gaps[ 2, ])
	}
	
	# offset is how many characters from head of CR1 to head of CR2,
	# in the 'after combine' state.   positive offset means CR1 starts first...
	# turn that into combined size and where the two go in the result
	
	len1 <- cr1$len;  len2 <- cr2$len
	beg1 <- beg2 <- 1
	end1 <- len1;  end2 <- len2;
	if ( offset > 0) {
		beg2 <- offset + 1
		end2 <- offset + len2
	}
	if ( offset < 0) {
		beg1 <- 1 - offset
		end1 <- len1 - offset
	}
	# now we know the new size, and all 'end points in new bigger read' are known
	newlen <- max( c(end1, end2))
	
	# make the Base Counts Matrix
	newbcm <- blankBCM( newlen)
	newbcm[ , beg1:end1] <- cr1$BCM
	newbcm[ , beg2:end2] <- newbcm[ , beg2:end2] + cr2$BCM
	
	# re-calculate the new consensus Bases
	ans <- baseConsensusFromBCM( newbcm)
	bases <- ans$bases
	scores <- ans$scores
	dna <- DNAString( base::paste( bases, collapse=""))
	len <- length( bases)
	
	# update the positions of where each USR is now in the new CR
	if ( any( cr1$USRid %in% cr2$USRid)) {
		cat( "\ncombine2CR: some USR IDs seen twice:\n")
		cat( "cr1: ", base::sort(cr1$USRid), " cr2: ", base::sort(cr1$USRid))
	}
	thisID <- c( cr1$USRid, cr2$USRid)
	thisPosition <- c( ( cr1$USRposition + beg1 - 1), ( cr2$USRposition + beg2 - 1))
	thisCRTid <- min( c( cr1$CRTid, cr2$CRTid))
	thisTotReads <- cr1$totReads + cr2$totReads
	
	# update the depth
	newdepth <- tmp1 <- tmp2 <- rep( 0, times=newlen)
	newdepth[ beg1:end1] <- tmp1[ beg1:end1] <- cr1$depth
	newdepth[ beg2:end2] <- tmp2[ beg2:end2] <- cr2$depth
	overlap <- max( c( beg1, beg2)) : min(c( end1, end2))
	newdepth[ overlap] <- tmp1[ overlap] + tmp2[ overlap]

	newgapslist <- base::append( cr1$USRgaps, cr2$USRgaps)
	newcolorslist <- base::append( cr1$USRcolors, cr2$USRcolors)
	
	newcr <- list( "len"=newlen, "bases"=bases, "scores"=scores, "dna"=dna, "BCM"=newbcm, 
			"depth"=newdepth, "baseIDs"=USR_BaseFactorsFunc( bases, needSplit=FALSE), 
				"USRid"=thisID, "USRposition"=thisPosition, 
				"USRgaps"=newgapslist, "USRcolors"=newcolorslist, "CRTid"=thisCRTid, 
				"totReads"=thisTotReads)

	return( newcr)
}


phredScoreCR <- function( cr) {

	# turn the consensus scores into Phred style scores
	phredScores <- 40.0 - (log( 10001 - (cr$scores * 100), 10) * 10)
	return( as.integer( round( phredScores)))
}


calc2CRalignment <- function( cr1, cr2) {

	ansNew <-  fastSRO( cr1$baseIDs, cr2$baseIDs, mismatch=MISMATCH_COST, gap=GAP_OPEN_COST)
	return( list( "ansNew"=ansNew, "ansOld"=NULL, "score"=ansNew$score))
}


getAlignmentOffsetAndGap <- function( align) {

	# turn the fast linear gaps into a matrix
	gapsM <- t(matrix( align$rawalignment$ansNew$gaps[1:4], 2 ,2))
	ansNew <- list( "offset"=align$rawalignment$ansNew$offset, "gaps"=gapsM)
	return( list( "ansNew"=ansNew, "ansOld"=NULL))
}


plotOneCRT <- function( i, label="") {

	# plot one CR in the CRT, with some labels...
	par( mfcol=c(1,1))
	readplotCR( CRT_List[[ i]], label=label)
	return()
}


plotOneCRTpair <- function( i, label="") {

	# plot one CR in the CRT, with some labels...
	par( mfcol=c(2,1))
	id1 <- CRT_pairTable_cr1[i]
	id2 <- CRT_pairTable_cr2[i]
	label <- paste( "alignScore:",round( CRT_pairTable_score[i]))
	readplotCR( CRT_List[[ id1]], label=label)
	readplotCR( CRT_List[[ id2]], label=label)
	par( mfcol=c(1,1))
	return()
}

  
plotTwoCR <- function( cr1, cr2, label="") {

	# plot 2 CR, with some labels...
	par( mfcol=c(2,1))
	readplotCR( cr1, label=label)
	readplotCR( cr2, label=label)
	par( mfcol=c(1,1))
	return()
}



readplotCR <- function( cr, label="") {
	
	# draw the reads in a 'pretty' way
	starts <- cr$USRposition
	ids <- cr$USRid
	visitOrd <- base::order( starts)
	savpar <- par("mai") 
	par("mai"=c( 0.25, savpar[2:4]))

	# we need a map of where on the plotted reads already cover
	freeSpots <- matrix( TRUE, nrow=round( max(cr$depth)*1.6), ncol=(cr$len+1))

	maintxt <- paste( "CR", cr$CRTid, "--   ", label,"\n N_USRs: ", length(cr$USRid), 
			"     N_Reads: ", formatC( as.integer(cr$totReads), format="d", big.mark=","),
			"     Percent of 'No Hits' file: ", 
			as.percent( cr$totReads, big.value=USR_GrandTotal, digits=3))
	plot( 1,1, type="n", xlim=c( 1, cr$len + 1), ylim=c(0.5, max(cr$depth)*1.1+1), xaxt="n", 
		main=maintxt, xlab="Consensus Sequence", ylab="Distinct USR Depth")

	bigY <- 1
	addBB <- (length( ids) < 20)
	for( i in 1:length( ids)) {
		thisID <- ids[ visitOrd[i]]
		xbeg <- starts[ visitOrd[i]]
		thisColors <- cr$USRcolors[[visitOrd[i]]]
		thisLen <- length( thisColors)
		# get where theres room to draw it, then mark that space as used
		ybeg <- findFreeSpot( freeSpots, xbeg, thisLen)
		freeSpots[ ybeg, xbeg:(xbeg+thisLen)] <- FALSE
		if ( ybeg > bigY) bigY <- ybeg
		
		drawOneRead( xbeg, ybeg, thisLen, thisColors, bbox=addBB)
	}
	mycex <- 1.0 * sqrt(40/cr$len)
	text( x=1:cr$len, y=1.0, label=cr$bases, cex=mycex, pos=1)
	text( x=1:cr$len, y=bigY+1.0, label=cr$bases, cex=mycex, pos=3)
	dev.flush()

	par("mai"=savpar)
	return()
}



#
#	BCM - Base Counts Matrix
#
#
# 	counts every occurrance of every ACGTN in this growing CR.  Used to
#	calculate the consensus sequence and the score for each base
#
newBCM <- function( bases, cnt) {
	
	# make a new Base Count Matrix
	bcm <- blankBCM( length( bases))
	if ( cnt < 1) {
		warning( paste( "newBCM: bad Reads Count: ", cnt))
		cnt <- 1
	}
	for ( i in 1:N_BASES) bcm[ i, (bases == BASES[i])] <- cnt
	return( bcm)
}


blankBCM <- function( len) {
	bcm <- matrix( 0, nrow=N_BASES, ncol=len)
	rownames(bcm) <- BASES
	return( bcm)
}


testBCM <- function( bcm) {
	if( ! all( rownames( bcm) == BASES)) return( FALSE)
	sums <- apply( bcm, MARGIN=2, sum)
	maxs <- apply( bcm, MARGIN=2, max)
	best <- apply( bcm, MARGIN=2, which.max)
	if ( any( sums == 0)) return( FALSE)
	if ( any( is.na( best))) return( FALSE)
	if ( any( best < 1)) return( FALSE)
	if ( any( best > N_BASES)) return( FALSE)
	scores <- 100 * maxs / sums
	if ( any( scores > 100)) return( FALSE)
	if ( any( scores < (100/N_BASES))) return( FALSE)
	# looks good!!
	return( TRUE)
}


revCompBCM <- function( BCM) {

	# do the 'reverse complement' of these base counts
	n <- ncol(BCM)
	# reverse the columns
	newbcm <- BCM[ , n:1]
	# flip the AT and CG rows
	newbcm <- newbcm[ c(4,3,2,1,5,6), ]
	rownames(newbcm) <- BASES
	return( newbcm)
}


baseConsensusFromBCM <- function( bcm) {
	
	sums <- apply( bcm, MARGIN=2, sum)
	maxs <- apply( bcm, MARGIN=2, max)
	best <- apply( bcm, MARGIN=2, which.max)
	scores <- 100 * maxs / sums
	bases <- BASES[ best]
	return( list( "bases"=bases, "scores"=scores))
}


basePercentsFromBCM <- function( bcm, scale=1) {
	
	sums <- apply( bcm, MARGIN=2, sum)
	pcts <- t(  t(bcm) / (sums/scale))
	return( pcts)
}



#
#	CRT -- Consensus Reads Table
#
#	the main data structure for maintaining the set of all CRs and how
#	well they score against each other
#
#
newCRT <- function( n) {
	
	# allocate space for N CRs
	crt <- vector( mode="list", length=n)

	# and a table of pairs of CRs, with their current scores
	maxPairs <- n * 100
	zeros <- rep( 0, times=maxPairs)
	aligns <- vector( mode="list", length=maxPairs)

	CRT_nCR <<- 0
	CRT_Size <<- n
	CRT_freeCR <<- 1:n
	CRT_List <<- crt

	CRT_nPairs <<- 0
	CRT_Size_Pairs <<- maxPairs
	CRT_freePairs <<- 1:maxPairs
	CRT_pairTable_cr1 <<- zeros
	CRT_pairTable_cr2 <<- zeros
	CRT_pairTable_score <<- zeros
	CRT_pairTable_align <<- aligns

	CRT_best <<- vector()

	CRT_reAlignCutoff <<- CRT_pairAlignCutoff <<- 0

	# when this gets reset, we must reset a few bits of the USR too
	USR_nUsed <<- 0

	return()
}


saveCRTcontext <- function( file="CRTcontext.rda", selfTest=FALSE, verbose=TRUE) {
	if (verbose) cat( "\nSaving CRT context...\n")
	save( CRT_Size, CRT_Size_Pairs, CRT_nCR, CRT_nPairs, 
		USR_GrandTotal, USR_nTotal, USR_nUsed, 
		CRT_List, CRT_freeCR, CRT_freePairs,
		CRT_pairTable_cr1, CRT_pairTable_cr2, 
		CRT_pairTable_score, CRT_pairTable_align, 
		CRT_best, file=file)

	if (selfTest) {
		if ( testAll() == FALSE) cat( "\nSelf-test failed...\n")
	}
	return()
}


CR_cleanup <- function() {

	# remove objects from global workspace that we used...
	rm( list=ls( pattern="^CRT_", envir=.GlobalEnv), envir=.GlobalEnv)
	if( exists( "CRplotted")) rm( CRplotted, envir=.GlobalEnv)
	if( exists( "USR_nUsed")) rm( USR_nUsed, envir=.GlobalEnv)
}


addOneUSRtoCRT <- function() {
	
	# get the next read not yet in the pool, allowing for invalid CRs to be skipped
	# and add it to the growing group

	newSlot <- CRT_freeCR[ (CRT_nCR + 1)]

	repeat {
		id <- USR_nUsed + 1
		thisUSR <- getUSRentry( id)
		USR_nUsed <<- id

		# skip over N's because they don't assemble
		if ( regexpr( "N", thisUSR$seq, fixed=T) > 0) next

		# make it its own 'small' consensus read
		cr <- newCR( thisUSR$seq, thisUSR$count, USRid=id, CRTid=newSlot)

		# make sure its valid
		if ( cr$len > 0) break
	}

	# add it to the CRT
	# (no need to update the free list, because it's location didn't change...)
	CRT_List[[ newSlot]] <<- cr
	CRT_nCR <<- CRT_nCR + 1

	# this new CR now gets a pair entry to everyone already in the CRT
	addPairSetsToCRT( newSlot)

	# and we're done
	return()
}


freeOneCRT <- function( id) {

	if ( id < 1 || id > CRT_Size) stop( paste( "freeOneCRT: invalid 'id': ", id))
	if ( CRT_nCR < 1) return()

	# push this location back onto the free stack
	lastInUse <- CRT_freeCR[ CRT_nCR]
	freeSpot <- CRT_nCR
	CRT_freeCR[ freeSpot] <<- id
	CRT_nCR <<- CRT_nCR - 1

	if ( lastInUse != id) {
		whereWasID <- base::match( id, CRT_freeCR[1:CRT_nCR])
		CRT_freeCR[whereWasID] <<- lastInUse
	}

	# erase what was there
	CRT_List[[ id]] <<- list( "len"=0)
	return()
}


inUseCRTs <- function() {

	# the first N of the free list is who is in use...
	if ( CRT_nCR < 1) {
		return( vector())
	}
	return( base::sort( CRT_freeCR[ 1:CRT_nCR]))
}


addPairSetsToCRT <- function( id) {

	# one new CR has been added to the CRT at position 'id'
	if ( CRT_nCR < 2) return()

	# it needs new pair entries with all existing CRs in the table.
	oldIDs <- setdiff( inUseCRTs(), id)
	newID <- id
	nNewPairs <- length( oldIDs)
	nOldPairs <- CRT_nPairs
	if ( (nNewPairs + nOldPairs) > (CRT_Size_Pairs)) stop( "CRT size exceeded...")

	# rather than blindly add them all, get the alignment first, and only keep pairings
	# that may eventually get joined...
	nnew <- 0
	greatScore <- CRT_combineCutoff * 1.1

	# let's look at both the new CR and its reverse complement, and see which is 'better'
	score1 <- score2 <- vector()
	ntry <- 0
	thisCR <- CRT_List[[newID]]
	rcCR <- revCompCR( thisCR)
	for ( i in oldIDs) {
		thisAlign1 <- calc2CRalignment( CRT_List[[i]], thisCR)
		thisAlign2 <- calc2CRalignment( CRT_List[[i]], rcCR)
		ntry <- ntry + 1
		score1[ntry] <- thisAlign1$score
		score2[ntry] <- thisAlign2$score
		if ( max( thisAlign1$score, thisAlign2$score) >= greatScore) break
	}
	# if the revComp gave a better match, use it instead
	if ( ntry > 0 && (max( score2) > max( score1))) {
		CRT_List[[ newID]] <<- thisCR <- rcCR
	}

	for ( i in oldIDs) {
		thisAlign <- calc2CRalignment( CRT_List[[i]], thisCR)
		if ( thisAlign$score < CRT_pairAlignCutoff) next
		nnew <- nnew + 1
		newSpot <- CRT_freePairs[ (CRT_nPairs + nnew)]
		CRT_pairTable_cr1[newSpot] <<- i
		CRT_pairTable_cr2[newSpot] <<- newID
		CRT_pairTable_score[newSpot] <<- thisAlign$score
		CRT_pairTable_align[[newSpot]] <<- thisAlign

		# if we get a rally high score, quit early, because these pairs will
		# be stripped out on the next pass
		if ( thisAlign$score >= greatScore) break
	}
	CRT_nPairs <<- nOldPairs + nnew
	return()
}


freeCRTpairs <- function( idSet) {

	# free up each pair in reverse order
	for( id in rev( idSet)) {
		if ( id < 1 || id > CRT_Size_Pairs) stop( "freeCRTpair: invalid 'id'")
		if ( CRT_nPairs < 1) return()

		lastInUse <- CRT_freePairs[ CRT_nPairs]
		freeSpot <- CRT_nPairs
		CRT_freePairs[ freeSpot] <<- id
		CRT_nPairs <<- CRT_nPairs - 1

		if ( lastInUse != id) {
			whereWasID <- base::match( id, CRT_freePairs[1:CRT_nPairs])
			CRT_freePairs[whereWasID] <<- lastInUse
		}

		# erase what was there
		CRT_pairTable_cr1[id] <<- 0
		CRT_pairTable_cr2[id] <<- 0
		CRT_pairTable_score[id] <<- 0
		CRT_pairTable_align[id] <<- list( "score"=0)
	}
	return()
}


inUseCRTpairs <- function() {

	# the first N of the free list is who is in use...
	if ( CRT_nPairs < 1) {
		return( vector())
	}
	return( base::sort( CRT_freePairs[ 1:CRT_nPairs]))
}


getBestCRTalignment <- function() {

	if ( CRT_nPairs < 1) return( list( "score"=0))

	myPairs <- inUseCRTpairs()
	best <- which.max( CRT_pairTable_score[ myPairs])
	best <- myPairs[ best]
	return( list( "pairID"=best, "score"=CRT_pairTable_score[best], 
		"rawalignment"=CRT_pairTable_align[[best]]))
}


combineAndReduceCRT <- function( align) {

	# this CRT pair has the best alignment score, and is to be combined into 1 new CR.
	thisID <- align$pairID
	if ( ! (thisID %in% inUseCRTpairs())) {
		stop( paste( "combineAndReduceCRT:  pairTable ID is not in use: ", thisID))
	}

	# the offset for combining is taken from the aligment details
	ans <- getAlignmentOffsetAndGap( align)
	offset <- ans$ansNew$offset
	gaps <- ans$ansNew$gaps

	# get the two old entries in the CRlist
	id1 <- CRT_pairTable_cr1[ thisID]
	id2 <- CRT_pairTable_cr2[ thisID]
	keeper <- min( c( id1,id2))
	dropper <- max( c( id1,id2))
	cr1 <- CRT_List[[ id1]]
	cr2 <- CRT_List[[ id2]]

	if ( ! is.null( gaps)) {
	if ((gaps[1,1] >= cr1$len) || (gaps[2,1] >= cr2$len)) {
		cat( "\n\nBad Gap: (cr1,cr2) ", id1, id2,"\n")
		print( gaps)
		print( align)
		saveCRTcontext()
		stop()
	}}

	# make the new bigger CR and stuff it back into the 'higher' spot in CRlist
	newcr <- combine2CR( cr1, cr2, offset=offset, gaps=gaps)
	if ( newcr$CRTid != keeper) {
		cat( "\n\nBad CRT pointer returned from 'combine2CR'", id1, id2, keeper, 
			dropper, newcr$CRTid,"\n")
		stop()
	}
	CRT_List[[ keeper]] <<- newcr
	freeOneCRT( dropper)


	# eliminate all pairs that use the CR being dropped...
	who <- inUseCRTpairs()
	pairsToDrop <- which( CRT_pairTable_cr1[who] == dropper | CRT_pairTable_cr2[who] == dropper )
	pairsToDrop <- who[ pairsToDrop]
	freeCRTpairs( pairsToDrop)

	# almost done...  both the CRlist and pairTable are now up to date for the smaller number of entries.
	# now, rescore the pairTable entries that the new bigger CR is party to...
	who <- inUseCRTpairs()
	reVisits <- which( (CRT_pairTable_cr1[who] == keeper | CRT_pairTable_cr2[who] == keeper))
	reVisits <- who[ reVisits]

	# make a growing list of who else to drop and a shrinking list of who else to add...
	newDrops <- vector()
	allOtherCRids <- setdiff( inUseCRTs(), keeper)

	for ( j in reVisits) {
		id1 <- CRT_pairTable_cr1[j]
		id2 <- CRT_pairTable_cr2[j]
		allOtherCRids <- setdiff( allOtherCRids, c(id1,id2))
		thisAlign <- calc2CRalignment( CRT_List[[ id1]], CRT_List[[ id2]])
		if ( thisAlign$score < CRT_pairAlignCutoff) {
			newDrops <- base::append( newDrops, j)
			next
		}
		CRT_pairTable_score[j] <<- thisAlign$score
		CRT_pairTable_align[[j]] <<- thisAlign
	}
	# now re-drop these newly 'no good' pairs
	if ( length( newDrops) > 0) {
		freeCRTpairs( newDrops)
	} # end 'newDrops'

	# if there are any CR left in the list of 'others', they need new pairs from our keeper...
	if ( length( allOtherCRids) > 0) {
		nnew <- 0
		for ( j in 1:length( allOtherCRids)) {
			id1 <- min( c( keeper, allOtherCRids[j]))
			id2 <- max( c( keeper, allOtherCRids[j]))
			thisAlign <- calc2CRalignment( CRT_List[[id1]], CRT_List[[id2]])
			if ( thisAlign$score < CRT_pairAlignCutoff) next
			nnew <- nnew + 1
			newSpot <- CRT_freePairs[ (CRT_nPairs + nnew)]
			CRT_pairTable_cr1[newSpot] <<- id1
			CRT_pairTable_cr2[newSpot] <<- id2
			CRT_pairTable_score[newSpot] <<- thisAlign$score
			CRT_pairTable_align[[newSpot]] <<- thisAlign
		}
		CRT_nPairs <<- CRT_nPairs + nnew
	}

	# really done...
	return( keeper)
}



testCRT <- function() {

	ok <- TRUE
	# basic sizes
	if ( CRT_nCR < 0 || CRT_nCR > CRT_Size) ok <- FALSE
	if ( CRT_nPairs < 0 || CRT_nPairs > (CRT_Size * CRT_Size)) ok <- FALSE
	if ( ! ok) cat( "\nCRT size test failed")
	if ( CRT_nCR == 0) return( ok)

	# all the CRs are ok?
	myIDs <- inUseCRTs()
	tests <- sapply( CRT_List[ myIDs], FUN=testCR)
	if ( any( tests == FALSE)) {
		ok <- FALSE
		cat( "\nbad tests: ", myIDs[ which( tests == FALSE)])
	}
	if ( ! ok) cat( "\nCRT test of all CRs failed")
	
	idsUsed <- CRT_freeCR[ 1:CRT_nCR]
	idsFree <- CRT_freeCR[ (CRT_nCR+1):CRT_Size]
	if ( length( overlap <- intersect( idsUsed, idsFree)) > 0) {
		ok <- FALSE
		cat( "\nSome CRT ids in both 'used' and 'free' lists")
		cat( "\nN=", length(overlap), "who=", overlap)
	}
	if ( length( overlap <- base::union( idsUsed, idsFree)) != CRT_Size ) {
		ok <- FALSE
		cat( "\nSome CRT ids not in either 'used' or 'free' lists")
		missing <- setdiff( 1:CRT_Size, overlap)
		cat( "\nN=", length(missing), "who=", missing)
	}

	# are all the CRTids inside the CRs ok?
	ids <- sapply( CRT_List[ myIDs], FUN=function(x) return( x$CRTid))
	if ( length( who <- which( ids != myIDs)) > 0) {
		ok <- FALSE
		cat( "\nBad CRTids in CRs: ", who)
	}
	if ( ! ok) cat( "\nCRT test of all CRs embedded CRTids failed")

	# all the pair ptrs to CRs ok?
	if ( CRT_nPairs > 0) {
	myPairs <- inUseCRTpairs()
	allptrs <- base::union( unique.default( CRT_pairTable_cr1[ myPairs]), unique.default( CRT_pairTable_cr2[ myPairs])) 
	if ( ! all( allptrs %in% myIDs)) ok <- FALSE
	if ( ! ok) cat( "\nCRT pair pointers failed")

	notPaired <- setdiff( myIDs, allptrs)
	if ( length( notPaired) > 0) {
		cat( "\nInfo:  Some CRs are in no pairs: ", notPaired)
	}}

	# try to prove that no USR is a member of more that one CR
	# and that every USR up to now is in some one CR... 'No longer True' since polyNs get skipped
	allUSRids <- base::unlist( lapply( CRT_List[myIDs], FUN=function(x) {return( x$USRid)}))
	allUSRids <- base::sort( allUSRids)
	expect <- 1:USR_nUsed
	if ( length(dups <- which( duplicated( allUSRids))) > 0) {
		ok <- FALSE
		cat( "\nDuplicate USR ids: ", allUSRids[dups],"\n")
	}
	if ( ! ok) cat( "\nCRT CR_to_USR pointers failed")
	return( ok)
}



findFreeSpot <- function( isFree, x, len) {

	# find the lowest Y that has room at [X:X=len]

	xlast <- x + len
	for (iy in 1:nrow(isFree)) {
		if ( all( isFree[ iy, x:xlast])) return( iy)
	}
	return( NA)
}


drawOneRead <- function( x, y, len, colorVec, bbox=TRUE) {

	# draw a colored bar at the given spot...
	xl <- seq( x-0.5, (x+len-1.5), by=1.0)
	xr <- xl + 1
	yu <- y + 0.85
	# we may need to leave out the black edges...
	if (bbox) {
		rect( xl, y, xr, yu, col=BASE_COLORS[colorVec])
		rect( xl[1], y, xr[len], yu, col=NULL, lwd=2)
	} else {
		rect( xl, y, xr, yu, col=BASE_COLORS[colorVec], border=NA)
	}
	return()
}



CRplotter <- function( seconds=5, plotOrder="index", N=NULL, label="", asPNG=FALSE, 
			doPlot=TRUE, whoToPlot=NULL, pngPath=".") {
	
	# turning this into being able to do it in 2 parts, find the 'whos' only,
	# and then draw them on a later pass
	if ( is.null(whoToPlot)) {

	who <- 1:CRT_nCR
	len <- sapply( CRT_List[ who], function(x) return( x$len))
	who <- who[ len > 0]
	len <- len[ len > 0]
	nreads <<- sapply( CRT_List[ who], function(x) return( x$totReads))

	if ( plotOrder == "length") {
		ord <- base::order( len, decreasing=TRUE)
		who <- who[ ord]
		len <- len[ ord]
		nreads <<- nreads[ ord]
	}
	if ( plotOrder == "count") {
		ord <- base::order( nreads, decreasing=TRUE)
		who <- who[ ord]
		len <- len[ ord]
		nreads <<- nreads[ ord]
	}

	CRplotted <<- vector()

	if ( ! is.null(N)) {
		who <- who[ 1:N]
		len <- len[ 1:N]
		nreads <- nreads[ 1:N]
	}

	# OK, now we know who...
	whoToPlot <- who
	}

	if ( ! doPlot) return( whoToPlot) 

	# now do those plots...
	for( i in whoToPlot) {
		plotOneCRT( i, label=label)

		if ( asPNG) {
			f <- file.path( pngPath, paste( "cr_",i,".png", sep=""))
			png( f, width=1000, height=700, bg="white")
			plotOneCRT( i, label=label)
			dev.off()
		}

		CRplotted <<- base::append( CRplotted, i)

		# perhaps pause
		if ( ! is.null(seconds)) {
			if ( seconds > 0) {
				Sys.sleep( seconds)
			} else if (seconds < 0) {
				locator(1)
			}
		}
	}
	return( CRplotted)
}

