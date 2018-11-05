# CP.R - Consensus Peptides

require( Biostrings)

# constants for the Amino Acid residues
N_RESIDUES <- 26
RESIDUES <- c( "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K",
		"M", "F", "P", "S", "T", "W", "Y", "V", "B", "J", "Z", "X", "?", "-")
RESIDUE_GAP <- 26
RESIDUE_COLORS <- c( rainbow( 20, end=0.80), rep.int('gray60',3), rep.int('brown',2), 'black')
GAP_OPEN_COST <- (-12)
MISMATCH_COST <- (-6)



# 
# run the Consensus Peptides finder
#
runCP <- function( goodScore=40, iterations=1000, restart=TRUE, maxCP=4000, contextFile="CPTcontext.rda",
			max.time=1000) {

	if ( restart) {
		newCPT( maxCP)
		cat( "\n")
	}
	nHit <- nNoLuck <- 0
	newHit <- FALSE
	grow <- TRUE
	startTime <- proc.time()
	
	# score for combining are residued on what we call 'good enough'
	CPT_combineCutoff <<- goodScore * 1.00
	CPT_reAlignCutoff <<- goodScore / 2
	CPT_pairAlignCutoff <<- goodScore / 2
	stallLimit <- iterations / 2

	if ( ! exists( "USP_nTotal")) stop( "No USP data is loaded...")

	cat( "\n\nSearching for Consensus Peptides...  Good Score=", goodScore, "\n")

	grow <- ( USP_nUsed < USP_nTotal)
	while( TRUE) {

		# grow the set of Consensus Peptides by one
		# if room, add one; otherwise loosen match criteria
		

		if (grow) {
		    if ( CPT_nCP < CPT_Size && USP_nUsed < USP_nTotal) {
			addOneUSPtoCPT()
		    } else {
			CPT_combineCutoff <<- CPT_combineCutoff - 1
			if ( CPT_combineCutoff < goodScore) break
		    }
		}

		if ( newHit && nHit %% 100 == 0) {
			cat("\n(nUSP, nCP, nPairs, Score, nCombo, %All)", USP_nUsed, 
					CPT_nCP, CPT_nPairs, CPT_combineCutoff, nHit, 
					as.percent( sum(USP_counts[1:USP_nUsed]), big.value=USP_GrandTotal)) 
			# allow a loosening on the fly
			if ( USP_nUsed > 1000 &&  USP_nUsed/CPT_nCP < 1.5) {
				CPT_combineCutoff <<- CPT_combineCutoff - 1
			}
		}

		if ( nHit >= iterations) break
		if ( nNoLuck >= stallLimit) break
		
		# find the best now
		bestHit <- getBestCPTalignment()

		if ( bestHit$score >= CPT_combineCutoff) {
			grow <- FALSE
			who <- combineAndReduceCPT( bestHit)
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
	saveCPTcontext( file=contextFile)

	stopTime <- proc.time()
	dtime <- elapsedProcTime( startTime, stopTime, N=nHit, what="combinations")
	out <- list( "Combos"=nHit, "Time"=dtime$Raw.Seconds, "CurrentScore"=CPT_combineCutoff)
	return( out)
}



#
#	CP - Consensus Peptide
#
#	Growable Meta-Peptide that encompasses all Unique Short Peptides with 'good' 
#	between-peptide alignment
#
newCP <- function( str="", cnt=1, USPid=1, CPTid=1) {
	
	# make a new CP object from a character string of amino acids
	emptyCP <- list( "len"=0)
	str <- toupper( str[1] )
	len <- nchar( str)
	if ( len == 0) return( emptyCP)
	residues <- strsplit( str, split="")[[1]]
	badResidues <- residues[ which(  ! ( residues %in% RESIDUES))]
	if ( length( badResidues) > 0) {
		cat( "\nnewCP:  non-AA residues: ", badResidues)
		return( emptyCP)
	}
	aa <- AAString( str)
	depth <- rep( 1, times=len)
	scores <- rep( 100, times=len)
	bcm <- newRCM( residues, cnt)
	thisID <- USPid
	thisPosition <- 1
	thisGaps <- thisColors <- vector( mode="list", length=1)
	thisGaps[[1]] <- c(0,0)
	if ( USPid > 0) {
		thisColors[[1]] <- getUSPentry( USPid)$colors
	} else {
		thisColors[[1]] <- USP_ResidueFactorsFunc( residues, needSplit=FALSE)
	}
	names( thisGaps)[1] <- names( thisColors)[1] <- USPid

	return( list( "len"=len, "residues"=residues, "scores"=scores, "aa"=aa, "RCM"=bcm, "depth"=depth, 
			"residueIDs"=USP_ResidueFactorsFunc( residues, needSplit=FALSE), 
			"USPid"=thisID, "USPposition"=thisPosition, "USPgaps"=thisGaps, 
			"USPcolors"=thisColors, "CPTid"=CPTid, "totPeptides"=cnt ))
}


testCP <- function( cp) {
	ok <- TRUE
	if ( ! ("len" %in% names(cp))) return( FALSE)
	if ( is.na( cp$len)) ok <- FALSE
	if ( cp$len < 0) ok <- FALSE
	if ( ! ok) cat( "\nCP length test failed")
	if ( cp$len == 0) {
		cat("\ntesting an unused CP...")
		return(FALSE)
	}
	if ( ! all( names( cp) == c( "len", "residues", "scores", "aa", "RCM", "depth", "residueIDs",
			"USPid", "USPposition", "USPgaps", "USPcolors", "CPTid", "totPeptides" ))) ok <- FALSE
	if ( ! ok) cat( "\nCP list names test failed")
	if ( length( cp$residues) != cp$len) ok <- FALSE
	if ( nchar( cp$aa) != cp$len) ok <- FALSE
	if ( ncol( cp$RCM) != cp$len) ok <- FALSE
	if ( ! ok) cat( "\nCP seq length test failed")
	if ( max( cp$depth) > length( cp$USPid)) ok <- FALSE
	if ( any( cp$USPposition < 1 | cp$USPposition > cp$len)) ok <- FALSE
	if ( ! ok) cat( "\nCP USPpositions test failed")
	for ( i in 1:length(cp$USPid)) {
		thisPos <- cp$USPposition[i]
		thiscolorset <- cp$USPcolors[[i]]
		if ( (thisPos + length(thiscolorset) - 1) > cp$len) ok <- FALSE
	}
	if ( ! ok) cat( "\nCP USPcolor lengths test failed")
	if ( ! testRCM( cp$RCM)) ok <- FALSE
	if ( ! ok) cat( "\nCP RCM test failed")
	ans <- residueConsensusFromRCM( cp$RCM)
	if ( any( ans$residues != cp$residues)) ok <- FALSE
	if ( any( ans$scores != cp$scores)) ok <- FALSE
	if ( ! ok) cat( "\nCP consensus seq test failed")
	if ( ! ok) cat( " \t CPid= ", cp$CPTid)
	# looks good!!
	return( ok)
}


gapShiftCP <- function( cp, gapDesc) {

	# apply a gap shift to a CP because an alignment says it needs one.
	if ( gapDesc[1] <= 1) return( cp)
	if ( gapDesc[2] == 0) return( cp)

	gapStart <- gapDesc[1]
	gapLength <- gapDesc[2]

	oldlen <- cp$len
	oldbcm <- cp$RCM

	# change the vector of residues
	newresidues <- c( cp$residues[1:(gapStart-1)], rep( "-", times=gapLength), cp$residues[ gapStart:oldlen])
	newlen <- length( newresidues)
	if ( newlen != (oldlen+gapLength)) {
		cat( "\n\ngapShift Dump:\n(gapS, gapL, oldL,newL,CPid)", gapStart, gapLength,
			oldlen,newlen, cp$CPTid,"\n")
		stop("bad math...")
	}

	# change the RCM matrix
	newbcm <- newRCM( newresidues, cp$totPeptides)
	newbcm[ , 1:(gapStart-1)] <- oldbcm[ , 1:(gapStart-1)]
	newbcm[ , (gapStart+gapLength):ncol(newbcm)] <- oldbcm[ , gapStart:oldlen]

	# rebuild the scores and AA string
	ans <- residueConsensusFromRCM( newbcm)
	newaa <- AAString( paste( newresidues, collapse=""))
	newscores <- ans$scores

	# adjust the depth vector
	newdepth <- c( cp$depth[1:(gapStart-1)], rep( max(cp$depth[gapStart:(gapStart+gapLength-1)]), 
			times=gapLength), cp$depth[ gapStart:oldlen])

	# adjust the gaps list and residue colors and positions
	newposition <- cp$USPposition
	newgaplist <- cp$USPgaps
	newcolorslist <- cp$USPcolors
	for ( i in 1:length( cp$USPid)) {
		# the start position of each USP may change...
		if ( newposition[i] >= gapStart) newposition[i] <- newposition[i] + gapLength

		# if no gaps yet, this is first, else extend the set...
		tmpgap <- newgaplist[[i]]
		if ( tmpgap[1] == 0) {
			tmpgap <- c( gapStart, gapLength)
		} else {
			tmpgap <- c( tmpgap, c( gapStart, gapLength))
		}
		newgaplist[[i]] <- tmpgap

		# each color vector may be in a different relative location for this gap
		thisNewColors <- oldColors <- newcolorslist[[i]]
		lenOC <- length( oldColors)
		myGapStart <- gapStart - newposition[i] + 1
		if ( (myGapStart > 1) && (myGapStart <= lenOC)) {
			thisNewColors <- c( oldColors[1:(myGapStart-1)], rep( RESIDUE_GAP, times=gapLength), 
						oldColors[ myGapStart:lenOC])
		}
		newcolorslist[[i]] <- thisNewColors
	}

	return( list( "len"=newlen, "residues"=newresidues, "scores"=newscores, "aa"=newaa, 
			"RCM"=newbcm, "depth"=newdepth, 
			"residueIDs"=USP_ResidueFactorsFunc( newresidues, needSplit=FALSE), 
			"USPid"=cp$USPid, "USPposition"=newposition, "USPgaps"=newgaplist, 
			"USPcolors"=newcolorslist, "CPTid"=cp$CPTid, "totPeptides"=cp$totPeptides ))
}


combine2CP <- function( cp1, cp2, offset=0, gaps=NULL) {
	
	# join 2 CPs that have been pre-selected for overlap into 1:
			
	# watch for gapping found by alignment...
	if ( ! is.null( gaps)) {
		if ( gaps[1,1] > 1) cp1 <- gapShiftCP( cp1, gaps[ 1, ])
		if ( gaps[2,1] > 1) cp2 <- gapShiftCP( cp2, gaps[ 2, ])
	}
	
	# offset is how many characters from head of CP1 to head of CP2,
	# in the 'after combine' state.   positive offset means CP1 starts first...
	# turn that into combined size and where the two go in the result
	
	len1 <- cp1$len;  len2 <- cp2$len
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
	
	# make the Residue Counts Matrix
	newbcm <- blankRCM( newlen)
	newbcm[ , beg1:end1] <- cp1$RCM
	newbcm[ , beg2:end2] <- newbcm[ , beg2:end2] + cp2$RCM
	
	# re-calculate the new consensus Residues
	ans <- residueConsensusFromRCM( newbcm)
	residues <- ans$residues
	scores <- ans$scores
	aa <- AAString( paste( residues, collapse=""))
	len <- length( residues)
	
	# update the positions of where each USP is now in the new CP
	if ( any( cp1$USPid %in% cp2$USPid)) {
		cat( "\ncombine2CP: some USP IDs seen twice:\n")
		cat( "cp1: ", sort(cp1$USPid), " cp2: ", sort(cp1$USPid))
	}
	thisID <- c( cp1$USPid, cp2$USPid)
	thisPosition <- c( ( cp1$USPposition + beg1 - 1), ( cp2$USPposition + beg2 - 1))
	thisCPTid <- min( c( cp1$CPTid, cp2$CPTid))
	thisTotPeptides <- cp1$totPeptides + cp2$totPeptides
	
	# update the depth
	newdepth <- tmp1 <- tmp2 <- rep( 0, times=newlen)
	newdepth[ beg1:end1] <- tmp1[ beg1:end1] <- cp1$depth
	newdepth[ beg2:end2] <- tmp2[ beg2:end2] <- cp2$depth
	overlap <- max( c( beg1, beg2)) : min(c( end1, end2))
	newdepth[ overlap] <- tmp1[ overlap] + tmp2[ overlap]

	newgapslist <- c( cp1$USPgaps, cp2$USPgaps)
	newcolorslist <- c( cp1$USPcolors, cp2$USPcolors)
	
	newcp <- list( "len"=newlen, "residues"=residues, "scores"=scores, "aa"=aa, "RCM"=newbcm, 
			"depth"=newdepth, "residueIDs"=USP_ResidueFactorsFunc( residues, needSplit=FALSE), 
				"USPid"=thisID, "USPposition"=thisPosition, 
				"USPgaps"=newgapslist, "USPcolors"=newcolorslist, "CPTid"=thisCPTid, 
				"totPeptides"=thisTotPeptides)

	return( newcp)
}


phredScoreCP <- function( cp) {

	# turn the consensus scores into Phred style scores
	phredScores <- 40.0 - (log( 10001 - (cp$scores * 100), 10) * 10)
	return( as.integer( round( phredScores)))
}


calc2CPalignment <- function( cp1, cp2) {

	ansNew <-  fastSRO( cp1$residueIDs, cp2$residueIDs, mismatch=MISMATCH_COST, gap=GAP_OPEN_COST)
	return( list( "ansNew"=ansNew, "ansOld"=NULL, "score"=ansNew$score))
}


getAlignmentOffsetAndGap <- function( align) {

	# turn the fast linear gaps into a matrix
	gapsM <- t(matrix( align$rawalignment$ansNew$gaps[1:4], 2 ,2))
	ansNew <- list( "offset"=align$rawalignment$ansNew$offset, "gaps"=gapsM)
	return( list( "ansNew"=ansNew, "ansOld"=NULL))
}


plotOneCPT <- function( i, label="") {

	# plot one CP in the CPT, with some labels...
	par( mfcol=c(1,1))
	readplotCP( CPT_List[[ i]], label=label)
	return()
}


plotOneCPTpair <- function( i, label="") {

	# plot one CP in the CPT, with some labels...
	par( mfcol=c(2,1))
	id1 <- CPT_pairTable_cp1[i]
	id2 <- CPT_pairTable_cp2[i]
	label <- paste( "alignScore:",round( CPT_pairTable_score[i]))
	readplotCP( CPT_List[[ id1]], label=label)
	readplotCP( CPT_List[[ id2]], label=label)
	par( mfcol=c(1,1))
	return()
}

  
plotTwoCP <- function( cp1, cp2, label="") {

	# plot 2 CP, with some labels...
	par( mfcol=c(2,1))
	readplotCP( cp1, label=label)
	readplotCP( cp2, label=label)
	par( mfcol=c(1,1))
	return()
}



readplotCP <- function( cp, label="") {
	
	# draw the reads in a 'pretty' way
	starts <- cp$USPposition
	ids <- cp$USPid
	visitOrd <- order( starts)
	savpar <- par("mai") 
	par("mai"=c( 0.25, savpar[2:4]))

	# we need a map of where on the plotted reads already cover
	freeSpots <- matrix( TRUE, nrow=round( max(cp$depth)*1.6), ncol=(cp$len+1))

	maintxt <- paste( "CP", cp$CPTid, "--   ", label,"\n N_USPs: ", length(cp$USPid), 
			"     N_Peptides: ", formatC( as.integer(cp$totPeptides), format="d", big.mark=","),
			"     Percent of 'No Hits' file: ", 
			as.percent( cp$totPeptides, big.value=USP_GrandTotal, digits=3))
	plot( 1,1, type="n", xlim=c( 1, cp$len + 1), ylim=c(0.5, max(cp$depth)*1.1+1), xaxt="n", 
		main=maintxt, xlab="Consensus Sequence", ylab="Distinct USP Depth")

	bigY <- 1
	addBB <- (length( ids) < 20)
	for( i in 1:length( ids)) {
		thisID <- ids[ visitOrd[i]]
		xbeg <- starts[ visitOrd[i]]
		thisColors <- cp$USPcolors[[visitOrd[i]]]
		thisLen <- length( thisColors)
		# get where theres room to draw it, then mark that space as used
		ybeg <- findFreeSpot( freeSpots, xbeg, thisLen)
		freeSpots[ ybeg, xbeg:(xbeg+thisLen)] <- FALSE
		if ( ybeg > bigY) bigY <- ybeg
		
		drawOnePeptide( xbeg, ybeg, thisLen, thisColors, bbox=addBB)
	}
	mycex <- 1.0 * sqrt(40/cp$len)
	text( x=1:cp$len, y=1.0, label=cp$residues, cex=mycex, pos=1)
	text( x=1:cp$len, y=bigY+1.0, label=cp$residues, cex=mycex, pos=3)

	par("mai"=savpar)
	return()
}



#
#	RCM - Residue Counts Matrix
#
#
# 	counts every occurrance of every ACGTN in this growing CP.  Used to
#	calculate the consensus sequence and the score for each residue
#
newRCM <- function( residues, cnt) {
	
	# make a new Residue Count Matrix
	bcm <- blankRCM( length( residues))
	if ( cnt < 1) {
		warning( paste( "newRCM: bad Peptides Count: ", cnt))
		cnt <- 1
	}
	for ( i in 1:N_RESIDUES) bcm[ i, (residues == RESIDUES[i])] <- cnt
	return( bcm)
}


blankRCM <- function( len) {
	bcm <- matrix( 0, nrow=N_RESIDUES, ncol=len)
	rownames(bcm) <- RESIDUES
	return( bcm)
}


testRCM <- function( bcm) {
	if( ! all( rownames( bcm) == RESIDUES)) return( FALSE)
	sums <- apply( bcm, MARGIN=2, sum)
	maxs <- apply( bcm, MARGIN=2, max)
	best <- apply( bcm, MARGIN=2, which.max)
	if ( any( sums == 0)) return( FALSE)
	if ( any( is.na( best))) return( FALSE)
	if ( any( best < 1)) return( FALSE)
	if ( any( best > N_RESIDUES)) return( FALSE)
	scores <- 100 * maxs / sums
	if ( any( scores > 100)) return( FALSE)
	if ( any( scores < (100/N_RESIDUES))) return( FALSE)
	# looks good!!
	return( TRUE)
}


residueConsensusFromRCM <- function( bcm) {
	
	sums <- apply( bcm, MARGIN=2, sum)
	maxs <- apply( bcm, MARGIN=2, max)
	best <- apply( bcm, MARGIN=2, which.max)
	scores <- 100 * maxs / sums
	residues <- RESIDUES[ best]
	return( list( "residues"=residues, "scores"=scores))
}


residuePercentsFromRCM <- function( bcm, scale=1) {
	
	sums <- apply( bcm, MARGIN=2, sum)
	pcts <- t(  t(bcm) / (sums/scale))
	return( pcts)
}



#
#	CPT -- Consensus Peptides Table
#
#	the main data structure for maintaining the set of all CPs and how
#	well they score against each other
#
#
newCPT <- function( n) {
	
	# allocate space for N CPs
	cpt <- vector( mode="list", length=n)

	# and a table of pairs of CPs, with their current scores
	maxPairs <- n * 100
	zeros <- rep( 0, times=maxPairs)
	aligns <- vector( mode="list", length=maxPairs)

	CPT_nCP <<- 0
	CPT_Size <<- n
	CPT_freeCP <<- 1:n
	CPT_List <<- cpt

	CPT_nPairs <<- 0
	CPT_Size_Pairs <<- maxPairs
	CPT_freePairs <<- 1:maxPairs
	CPT_pairTable_cp1 <<- zeros
	CPT_pairTable_cp2 <<- zeros
	CPT_pairTable_score <<- zeros
	CPT_pairTable_align <<- aligns

	CPT_best <<- vector()

	CPT_reAlignCutoff <<- CPT_pairAlignCutoff <<- 0

	# when this gets reset, we must reset a few bits of the USP too
	USP_nUsed <<- 0

	return()
}


saveCPTcontext <- function( file="CPTcontext.rda", selfTest=FALSE) {
	cat( "\nSaving CPT context...\n")
	save( CPT_Size, CPT_Size_Pairs, CPT_nCP, CPT_nPairs, 
		USP_GrandTotal, USP_nTotal, USP_nUsed, 
		CPT_List, CPT_freeCP, CPT_freePairs,
		CPT_pairTable_cp1, CPT_pairTable_cp2, 
		CPT_pairTable_score, CPT_pairTable_align, 
		CPT_best, file=file)

	if (selfTest) {
		if ( testAll() == FALSE) cat( "\nSelf-test failed...\n")
	}
	return()
}


CP_cleanup <- function() {

	# remove objects from global workspace that we used...
	rm( list=ls( pattern="^CPT_", envir=.GlobalEnv), envir=.GlobalEnv)
	if( exists( "CPplotted")) rm( CPplotted, envir=.GlobalEnv)
	if( exists( "USP_nUsed")) rm( USP_nUsed, envir=.GlobalEnv)
}


addOneUSPtoCPT <- function() {
	
	# get the next read not yet in the pool, allowing for invalid CPs to be skipped
	# and add it to the growing group

	newSlot <- CPT_freeCP[ (CPT_nCP + 1)]

	repeat {
		id <- USP_nUsed + 1
		thisUSP <- getUSPentry( id)
		USP_nUsed <<- id

		# skip over N's because they don't assemble
		if ( regexpr( "N", thisUSP$seq, fixed=T) > 0) next

		# make it its own 'small' consensus read
		cp <- newCP( thisUSP$seq, thisUSP$count, USPid=id, CPTid=newSlot)

		# make sure its valid
		if ( cp$len > 0) break
	}

	# add it to the CPT
	# (no need to update the free list, because it's location didn't change...)
	CPT_List[[ newSlot]] <<- cp
	CPT_nCP <<- CPT_nCP + 1

	# this new CP now gets a pair entry to everyone already in the CPT
	addPairSetsToCPT( newSlot)

	# and we're done
	return()
}


freeOneCPT <- function( id) {

	if ( id < 1 || id > CPT_Size) stop( paste( "freeOneCPT: invalid 'id': ", id))
	if ( CPT_nCP < 1) return()

	# push this location back onto the free stack
	lastInUse <- CPT_freeCP[ CPT_nCP]
	freeSpot <- CPT_nCP
	CPT_freeCP[ freeSpot] <<- id
	CPT_nCP <<- CPT_nCP - 1

	if ( lastInUse != id) {
		whereWasID <- match( id, CPT_freeCP[1:CPT_nCP])
		CPT_freeCP[whereWasID] <<- lastInUse
	}

	# erase what was there
	CPT_List[[ id]] <<- list( "len"=0)
	return()
}


inUseCPTs <- function() {

	# the first N of the free list is who is in use...
	if ( CPT_nCP < 1) {
		return( vector())
	}
	return( sort( CPT_freeCP[ 1:CPT_nCP]))
}


addPairSetsToCPT <- function( id) {

	# one new CP has been added to the CPT at position 'id'
	if ( CPT_nCP < 2) return()

	# it needs new pair entries with all existing CPs in the table.
	oldIDs <- setdiff( inUseCPTs(), id)
	newID <- id
	nNewPairs <- length( oldIDs)
	nOldPairs <- CPT_nPairs
	if ( (nNewPairs + nOldPairs) > (CPT_Size_Pairs)) stop( "CPT size exceeded...")

	# rather than blindly add them all, get the alignment first, and only keep pairings
	# that may eventually get joined...
	nnew <- 0
	greatScore <- CPT_combineCutoff * 1.1

	thisCP <- CPT_List[[newID]]
	for ( i in oldIDs) {
		thisAlign <- calc2CPalignment( CPT_List[[i]], thisCP)
		if ( thisAlign$score < CPT_pairAlignCutoff) next
		nnew <- nnew + 1
		newSpot <- CPT_freePairs[ (CPT_nPairs + nnew)]
		CPT_pairTable_cp1[newSpot] <<- i
		CPT_pairTable_cp2[newSpot] <<- newID
		CPT_pairTable_score[newSpot] <<- thisAlign$score
		CPT_pairTable_align[[newSpot]] <<- thisAlign

		# if we get a rally high score, quit early, because these pairs will
		# be stripped out on the next pass
		if ( thisAlign$score >= greatScore) break
	}
	CPT_nPairs <<- nOldPairs + nnew
	return()
}


freeCPTpairs <- function( idSet) {

	# free up each pair in reverse order
	for( id in rev( idSet)) {
		if ( id < 1 || id > CPT_Size_Pairs) stop( "freeCPTpair: invalid 'id'")
		if ( CPT_nPairs < 1) return()

		lastInUse <- CPT_freePairs[ CPT_nPairs]
		freeSpot <- CPT_nPairs
		CPT_freePairs[ freeSpot] <<- id
		CPT_nPairs <<- CPT_nPairs - 1

		if ( lastInUse != id) {
			whereWasID <- match( id, CPT_freePairs[1:CPT_nPairs])
			CPT_freePairs[whereWasID] <<- lastInUse
		}

		# erase what was there
		CPT_pairTable_cp1[id] <<- 0
		CPT_pairTable_cp2[id] <<- 0
		CPT_pairTable_score[id] <<- 0
		CPT_pairTable_align[id] <<- list( "score"=0)
	}
	return()
}


inUseCPTpairs <- function() {

	# the first N of the free list is who is in use...
	if ( CPT_nPairs < 1) {
		return( vector())
	}
	return( sort( CPT_freePairs[ 1:CPT_nPairs]))
}


getBestCPTalignment <- function() {

	if ( CPT_nPairs < 1) return( list( "score"=0))

	myPairs <- inUseCPTpairs()
	best <- which.max( CPT_pairTable_score[ myPairs])
	best <- myPairs[ best]
	return( list( "pairID"=best, "score"=CPT_pairTable_score[best], 
		"rawalignment"=CPT_pairTable_align[[best]]))
}


combineAndReduceCPT <- function( align) {

	# this CPT pair has the best alignment score, and is to be combined into 1 new CP.
	thisID <- align$pairID
	if ( ! (thisID %in% inUseCPTpairs())) {
		stop( paste( "combineAndReduceCPT:  pairTable ID is not in use: ", thisID))
	}

	# the offset for combining is taken from the aligment details
	ans <- getAlignmentOffsetAndGap( align)
	offset <- ans$ansNew$offset
	gaps <- ans$ansNew$gaps

	# get the two old entries in the CPlist
	id1 <- CPT_pairTable_cp1[ thisID]
	id2 <- CPT_pairTable_cp2[ thisID]
	keeper <- min( c( id1,id2))
	dropper <- max( c( id1,id2))
	cp1 <- CPT_List[[ id1]]
	cp2 <- CPT_List[[ id2]]

	if ( ! is.null( gaps)) {
	if ((gaps[1,1] >= cp1$len) || (gaps[2,1] >= cp2$len)) {
		cat( "\n\nBad Gap: (cp1,cp2) ", id1, id2,"\n")
		print( gaps)
		print( align)
		saveCPTcontext()
		stop()
	}}

	# make the new bigger CP and stuff it back into the 'higher' spot in CPlist
	newcp <- combine2CP( cp1, cp2, offset=offset, gaps=gaps)
	if ( newcp$CPTid != keeper) {
		cat( "\n\nBad CPT pointer returned from 'combine2CP'", id1, id2, keeper, 
			dropper, newcp$CPTid,"\n")
		stop()
	}
	CPT_List[[ keeper]] <<- newcp
	freeOneCPT( dropper)


	# eliminate all pairs that use the CP being dropped...
	who <- inUseCPTpairs()
	pairsToDrop <- which( CPT_pairTable_cp1[who] == dropper | CPT_pairTable_cp2[who] == dropper )
	pairsToDrop <- who[ pairsToDrop]
	freeCPTpairs( pairsToDrop)

	# almost done...  both the CPlist and pairTable are now up to date for the smaller number of entries.
	# now, rescore the pairTable entries that the new bigger CP is party to...
	who <- inUseCPTpairs()
	reVisits <- which( (CPT_pairTable_cp1[who] == keeper | CPT_pairTable_cp2[who] == keeper))
	reVisits <- who[ reVisits]

	# make a growing list of who else to drop and a shrinking list of who else to add...
	newDrops <- vector()
	allOtherCPids <- setdiff( inUseCPTs(), keeper)

	for ( j in reVisits) {
		id1 <- CPT_pairTable_cp1[j]
		id2 <- CPT_pairTable_cp2[j]
		allOtherCPids <- setdiff( allOtherCPids, c(id1,id2))
		thisAlign <- calc2CPalignment( CPT_List[[ id1]], CPT_List[[ id2]])
		if ( thisAlign$score < CPT_pairAlignCutoff) {
			newDrops <- c( newDrops, j)
			next
		}
		CPT_pairTable_score[j] <<- thisAlign$score
		CPT_pairTable_align[[j]] <<- thisAlign
	}
	# now re-drop these newly 'no good' pairs
	if ( length( newDrops) > 0) {
		freeCPTpairs( newDrops)
	} # end 'newDrops'

	# if there are any CP left in the list of 'others', they need new pairs from our keeper...
	if ( length( allOtherCPids) > 0) {
		nnew <- 0
		for ( j in 1:length( allOtherCPids)) {
			id1 <- min( c( keeper, allOtherCPids[j]))
			id2 <- max( c( keeper, allOtherCPids[j]))
			thisAlign <- calc2CPalignment( CPT_List[[id1]], CPT_List[[id2]])
			if ( thisAlign$score < CPT_pairAlignCutoff) next
			nnew <- nnew + 1
			newSpot <- CPT_freePairs[ (CPT_nPairs + nnew)]
			CPT_pairTable_cp1[newSpot] <<- id1
			CPT_pairTable_cp2[newSpot] <<- id2
			CPT_pairTable_score[newSpot] <<- thisAlign$score
			CPT_pairTable_align[[newSpot]] <<- thisAlign
		}
		CPT_nPairs <<- CPT_nPairs + nnew
	}

	# really done...
	return( keeper)
}



testCPT <- function() {

	ok <- TRUE
	# basic sizes
	if ( CPT_nCP < 0 || CPT_nCP > CPT_Size) ok <- FALSE
	if ( CPT_nPairs < 0 || CPT_nPairs > (CPT_Size * CPT_Size)) ok <- FALSE
	if ( ! ok) cat( "\nCPT size test failed")
	if ( CPT_nCP == 0) return( ok)

	# all the CPs are ok?
	myIDs <- inUseCPTs()
	tests <- sapply( CPT_List[ myIDs], FUN=testCP)
	if ( any( tests == FALSE)) {
		ok <- FALSE
		cat( "\nbad tests: ", myIDs[ which( tests == FALSE)])
	}
	if ( ! ok) cat( "\nCPT test of all CPs failed")
	
	idsUsed <- CPT_freeCP[ 1:CPT_nCP]
	idsFree <- CPT_freeCP[ (CPT_nCP+1):CPT_Size]
	if ( length( overlap <- intersect( idsUsed, idsFree)) > 0) {
		ok <- FALSE
		cat( "\nSome CPT ids in both 'used' and 'free' lists")
		cat( "\nN=", length(overlap), "who=", overlap)
	}
	if ( length( overlap <- union( idsUsed, idsFree)) != CPT_Size ) {
		ok <- FALSE
		cat( "\nSome CPT ids not in either 'used' or 'free' lists")
		missing <- setdiff( 1:CPT_Size, overlap)
		cat( "\nN=", length(missing), "who=", missing)
	}

	# are all the CPTids inside the CPs ok?
	ids <- sapply( CPT_List[ myIDs], FUN=function(x) return( x$CPTid))
	if ( length( who <- which( ids != myIDs)) > 0) {
		ok <- FALSE
		cat( "\nBad CPTids in CPs: ", who)
	}
	if ( ! ok) cat( "\nCPT test of all CPs embedded CPTids failed")

	# all the pair ptrs to CPs ok?
	if ( CPT_nPairs > 0) {
	myPairs <- inUseCPTpairs()
	allptrs <- union( unique.default( CPT_pairTable_cp1[ myPairs]), unique.default( CPT_pairTable_cp2[ myPairs])) 
	if ( ! all( allptrs %in% myIDs)) ok <- FALSE
	if ( ! ok) cat( "\nCPT pair pointers failed")

	notPaired <- setdiff( myIDs, allptrs)
	if ( length( notPaired) > 0) {
		cat( "\nInfo:  Some CPs are in no pairs: ", notPaired)
	}}

	# try to prove that no USP is a member of more that one CP
	# and that every USP up to now is in some one CP... 'No longer True' since polyNs get skipped
	allUSPids <- unlist( lapply( CPT_List[myIDs], FUN=function(x) {return( x$USPid)}))
	allUSPids <- sort( allUSPids)
	expect <- 1:USP_nUsed
	if ( length(dups <- which( duplicated( allUSPids))) > 0) {
		ok <- FALSE
		cat( "\nDuplicate USP ids: ", allUSPids[dups],"\n")
	}
	if ( ! ok) cat( "\nCPT CP_to_USP pointers failed")
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


drawOnePeptide <- function( x, y, len, colorVec, bbox=TRUE) {

	# draw a colored bar at the given spot...
	xl <- seq( x-0.5, (x+len-1.5), by=1.0)
	xr <- xl + 1
	yu <- y + 0.85
	# we may need to leave out the black edges...
	if (bbox) {
		rect( xl, y, xr, yu, col=RESIDUE_COLORS[colorVec])
		rect( xl[1], y, xr[len], yu, col=NULL, lwd=2)
	} else {
		rect( xl, y, xr, yu, col=RESIDUE_COLORS[colorVec], border=NA)
	}
	return()
}



CPplotter <- function( seconds=5, plotOrder="index", N=NULL, label="", asPNG=FALSE, 
			doPlot=TRUE, whoToPlot=NULL, pngPath=".") {
	
	# turning this into being able to do it in 2 parts, find the 'whos' only,
	# and then draw them on a later pass
	if ( is.null(whoToPlot)) {

	who <- 1:CPT_nCP
	len <- sapply( CPT_List[ who], function(x) return( x$len))
	who <- who[ len > 0]
	len <- len[ len > 0]
	nreads <<- sapply( CPT_List[ who], function(x) return( x$totPeptides))

	if ( plotOrder == "length") {
		ord <- order( len, decreasing=TRUE)
		who <- who[ ord]
		len <- len[ ord]
		nreads <<- nreads[ ord]
	}
	if ( plotOrder == "count") {
		ord <- order( nreads, decreasing=TRUE)
		who <- who[ ord]
		len <- len[ ord]
		nreads <<- nreads[ ord]
	}

	CPplotted <<- vector()

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
		plotOneCPT( i, label=label)

		if ( asPNG) {
			f <- file.path( pngPath, paste( "cp_",i,".png", sep=""))
			png( f, width=1000, height=700, bg="white")
			plotOneCPT( i, label=label)
			dev.off()
		}

		CPplotted <<- c( CPplotted, i)

		# perhaps pause
		if ( ! is.null(seconds)) {
			if ( seconds > 0) {
				Sys.sleep( seconds)
			} else if (seconds < 0) {
				locator(1)
			}
		}
	}
	return( CPplotted)
}



#	USP	- Unique Short Peptides
#	The universe of all unique short peptides, to be searched for Consensus Peptides.


`USP_setup` <- function( file, sampleID, resultsPath=".", reload=FALSE, dropPolyX=TRUE) {

	filePath <- file.path( resultsPath, "USP")
	if ( ! file.exists( filePath)) dir.create( filePath, showWarnings=FALSE)
	USPfile <- file.path( filePath, paste( sampleID, "USP.rda", sep="."))

	if ( (!file.exists( USPfile)) || reload) {
		# create the universe of Unique Short Peptides
		cat( "\n\nCreating USP dataset for:  ", sampleID)
		newUSP( file, sep=",", dropPolyX=dropPolyX)
		saveUSPcontext( USPfile)
	} else {
		cat( "\n\nLoading existing USP dataset for:  ", sampleID)
		load( USPfile, envir=.GlobalEnv)
		cat( "   N_USP: ", USP_nTotal)
	}

	return( list( "nUSP"=USP_nTotal, "nPeptides"=USP_GrandTotal, "nPolyX"=USP_nDropPolyX, "USP_File"=USPfile))
}


`USP_cleanup` <- function() {

	# clean up
	if ( ! exists( "USP_seq")) return()
	rm( USP_seq, USP_nTotal, USP_counts, USP_GrandTotal, envir=.GlobalEnv)
	rm( list=c("USP_nDropPolyX", "USP_pctCoverage", "USP_AdapterHits"),  envir=.GlobalEnv)
	if ( exists("USP_colors")) rm( USP_colors, envir=.GlobalEnv)
}


newUSP <- function( file, peptideColumn="Peptide", sep="\t", dropPolyX=TRUE) {

	ans <- buildUSPfromFile( file, peptideColumn=peptideColumn, sep=sep, dropPolyX=dropPolyX)
	if ( length(ans$seq) < 2) stop( "newUSP:  not enough short peptides to proceed.")
	
	# for every unique peptide, get its repeat count and factoring of its residues
	seq <- toupper( ans$seq)
	cnts <- ans$counts

	# order from 'most frequent down'
	ord <- base::order( cnts, decreasing=TRUE)

	# save to global storage for speed
	USP_nTotal <<- ans$N
	USP_seq <<- seq[ord]
	USP_counts <<- cnts[ord]
	USP_GrandTotal <<- sum( USP_counts)
	USP_nDropPolyX <<- ans$nPolyX
	USP_pctCoverage <<- 1.0
	USP_AdapterHits <<- NULL
	return()
}


getUSPentry <- function( id) {
	if ( id < 1 | id > USP_nTotal) stop( paste( "getUSPentry:  invalid USP id numer: ", id))
	return( list( "seq"=USP_seq[id], "count"=USP_counts[id], "colors"=USP_ResidueFactorsFunc(USP_seq[id])))
}


writeUSPasFasta <- function( fileout, N=NULL) {

	if ( is.null(N)) N <- length( USP_counts)
	cat( "\nWriting ", N, "Unique Short Peptides to file: ", fileout)
	nams <- base::paste( "USP_", 1:N, "::", USP_counts[1:N], sep="")
	txt <- base::paste( ">", nams, "\n", USP_seq[1:N], sep="")
	writeLines( txt, con=fileout)
	cat( "  Done.\n")
	return()
}


buildUSPfromFile <- function( file, peptideColumn="Peptide", sep="\t", dropPolyX=TRUE) {
	
	raw <- read.delim( file, sep=sep, as.is=T)
	NR <- nrow(raw)
	nPolyX <- 0

	thisSet <- raw[[ peptideColumn]]
	allSeqs <- base::table( thisSet)

	if ( dropPolyX) {
		nXs <- sapply( names(allSeqs), FUN=function(x) {
			 	chr <- strsplit( x, split="")[[1]]
				return( ( sum( chr == "X") + sum( chr == "?")))
			})
		chrLen <- base::nchar( names(allSeqs)[1])
		drops <- which( nXs > (chrLen * 0.2))
		if ( length(drops) > 0) {
			allSeqs <- allSeqs[ -drops]
			class( allSeqs) <- "table"
			cat( "  PolyX: ", length(drops), "  ")
		}
		nPolyX <- length(drops)
	}

	return( list( "N"=length(allSeqs), "seq"=names( allSeqs), "counts"=allSeqs, "nPolyX"=nPolyX))
}


saveUSPcontext <- function( fileout) {
	# cat( "\nSaving USP context...\n")
	save( USP_nTotal, USP_GrandTotal, USP_seq, USP_counts, USP_nDropPolyX, USP_pctCoverage, 
			USP_AdapterHits, file=fileout)
	cat( "\nWrote USP file: ", fileout, "\n")
	return()
}


calc2USPalignment <- function( id1, id2) {

	GAP_OPEN_COST <- (-30)
	MISMATCH_COST <- (-12)
	return( pairwiseAlignment( USP_seq[id1], USP_seq[id2], type="overlap", gapOpening=GAP_OPEN_COST))
}


`USP_ResidueFactorsFunc` <- function(residues, needSplit=TRUE) {

	# turn a character string of DNA into a vector of base IDs
	if ( needSplit) residues <- strsplit(residues, split="")[[1]]

	# this may be faster...
	residues[ residues == "."] <- "N"
	return( base::match( residues, RESIDUES, nomatch=0))
}

