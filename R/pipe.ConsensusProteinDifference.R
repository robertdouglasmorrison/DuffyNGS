# pipe.ConsensusProteinDifference -- contrast proteins from 2 groups to find what's different

`pipe.ConsensusProteinDifference` <- function( sampleIDset, geneID, geneName=geneID, 
						annotationFile="Annotation.txt", optionsFile="Options.txt",
						groupColumn="Group", results.path=NULL, substitutionMatrix=NULL,
						doMSA=TRUE) {

	# set up the distance metric for assessing difference between sequences
	GAP <<- "-"
	if ( is.null(substitutionMatrix)) {
		require(Biostrings)
		if (version$major == "4" && as.numeric( version$minor) >= 4) require( pwalign)
		data( PAM70)
		substitutionMatrix <- PAM70
	}
	distMatrix <- asDistanceMatrix( substitutionMatrix)

	# the universe of possible AA calls is now fixed
	AA_LEVELS <<- sort( unique( c( GAP, colnames(distMatrix))))

	# get the samples and their group membership
	annT <- readAnnotationTable( annotationFile)
	annT <- subset( annT, SampleID %in% sampleIDset)
	sampleSet <- annT$SampleID
	nSamples <- length( sampleSet)
	if (is.null( groupColumn)) {
		groupSet <- rep.int( "Self", nSamples)
		myGroups <- "Self"
		nGroups <- 1
	} else {
		groupSet <- checkGroupNames( annT[[ groupColumn]])
		myGroups <- sort( unique( groupSet))
		nGroups <- length(myGroups)
	}
	if ( nGroups > 2) {
		cat("\n'ConsensusProteinDifference()' expects at most 2 groups")
		cat( "\nGiven:  ", myGroups)
		stop()
	}
	groupCompareName <- paste( myGroups, collapse=".v.")

	# confirm the needed files are found for all samples
	if ( is.null(results.path)) results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	consensus.path <- file.path( results.path, "ConsensusProteins")
	consensusDetailsFiles <- paste( sampleSet, geneName, "ConsensusAnswer.rda", sep=".")
	consensusDetailsFiles <- file.path( consensus.path, sampleSet, consensusDetailsFiles)
	found <- file.exists( consensusDetailsFiles)
	if ( any( ! found)) {
		cat( "\nSome sample consensus details files not found: \n", consensusDetailsFiles[!found])
		stop()
	}
	fastaFiles <- sub( "Answer.rda", "AA.fasta", consensusDetailsFiles)
	found <- file.exists( fastaFiles)
	if ( any( ! found)) {
		cat( "\nSome sample consensus AA files not found: \n", fastaFiles[!found])
		stop()
	}

	# step 1: make a FASTA of just these AA files
	cat( "\nBuilding FASTA from ", nSamples, "samples")
	seq <- vector()
	for ( i in 1:nSamples) {
		fa <- loadFasta( fastaFiles[i], verbose=F)
		seq[i] <- fa$seq[1]
	}
	myFA <- as.Fasta( desc=sampleSet, seq=seq)
	fastaInFile <- paste( "MAFFT", geneName, groupCompareName, "fasta", sep=".")
	writeFasta( myFA, fastaInFile, line.width=100)

	# step 2:  run the MAFFT tool to get a global MSA consensus
	alnOutFile <- paste( "MAFFT", geneName, "aln", sep=".")
	if (doMSA || ! file.exists(alnOutFile)) {
		cat( "\nCalling MAFFT to do MSA alignment..")
		mafftAns <- mafft( fastaInFile, alnOutFile, mode="global", outfmt="clustal", verbose=FALSE)
	} else {
		cat( "\nReading Previous MAFFT result..")
		mafftAns <- readALN( alnOutFile, verbose=F)
	}
	bigAAmatrix <- mafftAns$alignment
	rownames(bigAAmatrix) <- sampleSet
	colnames(bigAAmatrix) <- apply( bigAAmatrix, 2, function(x) {
				chars <- sort( table( x), decreasing=T)
				return( names(chars)[1])
		})
	N_AA <- ncol( bigAAmatrix)
	cat( "\nMSA results: ", N_AA, "|", dim( bigAAmatrix))

	# step 3:  get each sample's details, and perhaps expand to conform to the MSA consensus
	cat( "\nLoading consensus from each isolate..")
	consensusCalls <- consensusWts <- vector( length=nSamples, mode="list")
	for ( i in 1:nSamples) {
		
		# get the 3 items stored for this sample:  "AA_Calls, AA_Weights, Construct"
		cat( "  ", sampleSet[i])
		consensusAns <- NULL
		load( consensusDetailsFiles[i])
		AA_Calls <- consensusAns$AA_Calls
		AA_Weights <- consensusAns$AA_Weights
		Construct <- consensusAns$Construct
		faSeq <- myFA$seq[i]

		# expand them to fit the MSA
		ans <- expandToFitMSA( AA_Calls, AA_Weights, Construct, faSeq, bigAAmatrix[i, ],
					column.names=colnames(bigAAmatrix))
		consensusCalls[[i]] <- ans$AA_Calls
		consensusWts[[i]] <- ans$AA_Weights
		myDimC <- dim( ans$AA_Calls)
		myDimW <- dim( ans$AA_Weights)
		if ( i == 1) {
			myMsize <- myDimC
			if ( ! all( myMsize == myDimW)) stop( "Weight Matrix size error")
		} else {
			if ( ! all( myMsize == myDimC)) stop( "Call Matrix size error")
			if ( ! all( myMsize == myDimW)) stop( "Weight Matrix size error")
		}
	}
	names(consensusCalls) <- names(consensusWts) <- sampleSet
	
	# step 4.1:  convert AA counts to AA percentages
	cat( "\nConvert AA counts to percentages..")
	consensusPcts <- AAcountsToPercents( consensusCalls, consensusWts)
	names( consensusPcts) <- sampleSet

	# step 4.2:  average the percent AA distributions for each group
	finalPcts <- vector( length=nGroups, mode="list")
	names(finalPcts) <- myGroups
	nSamplesPerGroup <- vector()
	for ( i in 1:nGroups) {
		inGrp <- which( groupSet == myGroups[i])
		nSamplesPerGroup[i] <- length(inGrp)
		if ( length(inGrp) == 1) {
			finalPcts[[i]] <- consensusPcts[[ inGrp[1] ]]
		} else {
			cat( "\nAveraging group percentages: ", myGroups[i], "  ")
			finalPcts[[i]] <- averageIsolatePercents( consensusPcts[ inGrp])
		}
	}

	# step 5.1:  find the difference between groups at every loci
	if ( nGroups > 1) {
		cat( "\nMeasure group differences..")
		finalAns <- measureSequenceDifference( finalPcts, distMatrix, nSamplesPerGroup)

		# step 5.2:   make a distance matrix over all
		cat( "\nMeasure all pairwise distances..")
		distAns <- measureSequenceDistances( consensusPcts, distMatrix, groupNames=groupSet)
	
		# step 6:  summarize and report
		out <- list( "GroupDifference"=finalAns, "Distance"=distAns, "GeneID"=geneID, "GeneName"=geneName)

	} else {
		cat( "\nMeasure group variation..")
		finalAns <- measureSequenceVariation( finalPcts, distMatrix, nSamplesPerGroup)
		# step 6:  summarize and report
		out <- list( "GroupVariation"=finalAns, "GeneID"=geneID, "GeneName"=geneName)
	}

	# lastly, clean up
	out
}


expandToFitMSA <- function( AA_Calls, AA_Weights, Construct, seq, bigAAvec, column.names) {

	# given the consensus for one isolate, and the MSA gapped consensus for this same isolate
	N_AA_final <- length( bigAAvec)

	# this isolate is in matrix and string form, force them to be consistent
	N_AA_in_str <- length( Construct)
	seqLen <- nchar( seq)
	if ( N_AA_in_str < ncol(AA_Calls)) {
		AA_Calls <- AA_Calls[ ,1:N_AA_in_str]
		AA_Weights <- AA_Weights[ ,1:N_AA_in_str]
	}
	if ( seqLen < ncol(AA_Calls)) {
		AA_Calls <- AA_Calls[ ,1:seqLen]
		AA_Weights <- AA_Weights[ ,1:seqLen]
	}

	# OK to expand this isolate
	N_AA_in <- ncol( AA_Calls)
	NR <- nrow(AA_Calls)

	# the length of the MSA dominates -- we will either pad or clip ours as needed
	callsOut <- matrix( "", nrow=NR, ncol=N_AA_final)
	wtsOut <- matrix( 0, nrow=NR, ncol=N_AA_final)

	# step along, either moving our data or inserting gaps that have the weight of whereever we are in the sequence
	curX <- 1
	for ( i in 1:N_AA_final) {
	
		# what AA did the consensus MSA have for our isolate?
		thisAA <- bigAAvec[i]
		# match or mismatch:  move my data
		if (thisAA != GAP) {
			callsOut[ ,i] <- AA_Calls[ ,curX]
			wtsOut[ ,i] <- AA_Weights[ ,curX]
			curX <- curX + 1
			if ( curX > N_AA_in) break
			next
		}
		# the MSA says insert a gap, give it the same weight as we have for whereever we are now in our isolate
		hasCalls <- which( AA_Weights[ ,curX] > 0)
		callsOut[ hasCalls, i] <- GAP
		wtsOut[ hasCalls, i] <- AA_Weights[ hasCalls, curX]
	}
	colnames(callsOut) <- colnames(wtsOut) <- column.names

	# if we broke out before the end of the MSA, then this isolate has no data to add
	# thus nothing else needs to happen...

	out <- list( "AA_Calls"=callsOut, "AA_Weights"=wtsOut)
	return( out)
}


AAcountsToPercents <- function( consensusCalls, consensusWts) {

	# given AA calls and the counts for each, convert them to a percentage distribution of 
	# all possible AAs in standard fixed order
	N <- length(consensusCalls)
	out <- vector( length=N, mode="list")
	for ( i in 1:N) {
		cat( "  ", names(consensusCalls)[i])
		myCalls <- consensusCalls[[i]]
		myWts <- consensusWts[[i]]
		nr <- nrow(myCalls)
		nc <- ncol(myCalls)
		sml <- matrix( 0, nrow=length(AA_LEVELS), ncol=nc)
		rownames(sml) <- AA_LEVELS
		colnames(sml) <- colnames(myCalls)
		for ( j in 1:nc) {
			aa <- myCalls[ ,j]
			wt <- myWts[ ,j]
			allAA <- rep( aa, times=wt)
			if ( length(allAA) < 1) allAA <- GAP
			aaTbl <- table( factor( allAA, levels=AA_LEVELS))
			if (any(is.na(aaTbl))) {
				cat( "\nNA in AA table!:\n")
				print( aaTbl, exclude=NULL)
			}
			# convert counts to percents, round to integers on 0..100
			sml[ ,j] <- as.0.100.percent( aaTbl)
		}
		out[[i]] <- sml
	}
	out
}


averageIsolatePercents <- function( consensusPcts) {

	# given more than one isolate, all as percent AA distribution, merge/average them all
	N <- length(consensusPcts)

	# the result is exact same shape as one of the inputs
	out <- consensusPcts[[1]]
	NR <- nrow(out)
	NC <- ncol(out)
	for ( i in 1:NC) out[ ,i] <- 0
	
	# add in all the isolates
	for ( i in 1:N) {
		cat( "  ", names(consensusPcts)[i])
		myPcts <- consensusPcts[[i]]
		out <- out + myPcts
	}

	# rescale to 0..100 percents
	for ( i in 1:NC) {
		v <- out[ ,i]
		out[ ,i] <- as.0.100.percent(v)
	}
	out
}


measureSequenceDifference <- function( finalPcts, distMatrix, nSamplesPerGroup) {

	# given the 2 final AA percentage profiles, in identically sized tables along
	# the MSA sequence, find the 'difference distance' at every AA
	pctTbl1 <- finalPcts[[1]]
	pctTbl2 <- finalPcts[[2]]
	NC <- ncol( pctTbl1)
	NR <- nrow( pctTbl1)

	# storage for distance and p-value calls
	outD <- outP <- vector( length=NC)
	names(outD) <- names(outP) <- colnames(pctTbl1)
	outStr1 <- outStr2 <- rep.int("", NC)

	# we may want to scale to adjust for the number of samples
	nSamples <- sum( nSamplesPerGroup)

	# visit every column in both sets
	for ( i in 1:NC) {
		p1 <- pctTbl1[ ,i]
		p2 <- pctTbl2[ ,i]
		same <- (p1 == p2)
		nonZero1 <- (p1 > 0)
		nonZero2 <- (p2 > 0)
		nonZero <- (nonZero1 | nonZero2)

		# when identical percentages, nothing to evaluate
		if ( any( is.na(same))) {
			cat( "\nDebug:  NA in percent compares!:  AA=",i,"\n")
			print( p1)
			print( p2)
		}
		if ( all(same, na.rm=T)) {
			outD[i] <- 0
			outP[i] <- 1
		} else {
			# distance only needs the AA where they differ
			outD[i] <- calcOneSeqPctDistance( p1, p2, same, distMatrix)
			# p-value uses all AA with non-zero entries
			outP[i] <- calcOneSeqPctPvalue( p1, p2, nonZero, N=nSamples)
		}

		# make a character string that describes whats there
		ord <- order( p1[nonZero1], decreasing=T)
		who <- which( nonZero1)[ord] 
		outStr1[i] <- paste( AA_LEVELS[who], p1[who], sep=":", collapse=";")
		ord <- order( p2[nonZero2], decreasing=T)
		who <- which( nonZero2)[ord] 
		outStr2[i] <- paste( AA_LEVELS[who], p2[who], sep=":", collapse=";")
	}
	
	out <- data.frame( "Position"=1:NC, "ConsensusAA"=colnames(pctTbl1), "Distance"=outD, "Pvalue"=outP, 
			"String1"=outStr1, "String2"=outStr2, stringsAsFactors=F)
	GRP1_COL <<- 5
	GRP2_COL <<- 6
	colnames(out)[GRP1_COL:GRP2_COL] <- names(finalPcts)
	rownames(out) <- 1:nrow(out)
	out
}


measureSequenceVariation <- function( finalPcts, distMatrix, nSamplesPerGroup) {

	# given the 1 final AA percentage profile (not 2 groups given), 
	# the MSA sequence, find the 'variation' at every AA
	pctTbl1 <- finalPcts[[1]]
	NC <- ncol( pctTbl1)
	NR <- nrow( pctTbl1)

	# storage for distance and p-value calls
	outD <- outP <- vector( length=NC)
	names(outD) <- names(outP) <- colnames(pctTbl1)
	outStr1 <- rep.int("", NC)

	# we may want to scale to adjust for the number of samples
	nSamples <- sum( nSamplesPerGroup)

	# visit every column in both sets
	for ( i in 1:NC) {
		p1 <- pctTbl1[ ,i]
		nonZero1 <- (p1 > 0)
		useP1 <- p1[ nonZero1]
		if ( length(useP1) == 1) {
			outD[i] <- 0
			outP[i] <- 1
		} else {
			outD[i] <- 100 - mean( rep( useP1, times=useP1))
			outP[i] <- t.test( x=c(100, rep( useP1, times=useP1)), mu=max(useP1))$p.value
		}
		# make a character string that describes whats there
		ord <- order( useP1, decreasing=T)
		who <- which( nonZero1)[ord] 
		outStr1[i] <- paste( AA_LEVELS[who], p1[who], sep=":", collapse=";")
	}
	
	out <- data.frame( "Position"=1:NC, "ConsensusAA"=colnames(pctTbl1), "Distance"=outD, "Pvalue"=outP, 
			"String1"=outStr1, stringsAsFactors=F)
	GRP1_COL <<- 5
	colnames(out)[GRP1_COL] <- names(finalPcts)[1]
	rownames(out) <- 1:nrow(out)
	out
}


measureSequenceDistances <- function( allPcts, distMatrix, groupNames=NULL) {

	# given the full set of AA percentage profiles, in identically sized tables along
	# the MSA sequence, find the 'difference distance' at every AA
	N <- length( allPcts)
	labels <- names(allPcts)
	if ( ! is.null(groupNames)) labels <- paste( labels, groupNames, sep="_")
	NC <- ncol( allPcts[[1]])

	# storage for distance matrix to send back
	NOUT <- N * (N-1) / 2
	outDM <- rep.int( 0, NOUT)
	n <- 0

	# visit all pairs of isolates
	for ( i1 in 1:(N-1)) {
		pctTbl1 <- allPcts[[i1]]
		cat( "  ", labels[i1])
		for ( i2 in (i1+1):N) {
			pctTbl2 <- allPcts[[i2]]
			# visit every column in both sets
			outD <- rep.int( 0, NC)
			for ( i in 1:NC) {
				p1 <- pctTbl1[ ,i]
				p2 <- pctTbl2[ ,i]
				same <- (p1 == p2)
				# when identical percentages, nothing to evaluate
				if ( all(same)) next
				# distance only needs the AA where they differ
				outD[i] <- calcOneSeqPctDistance( p1, p2, same, distMatrix)
			}
			n <- n + 1
			outDM[n] <- sum( outD)
		}
	}

	# add the attributes to make this a distance object
	attributes(outDM) <- list( "Size"=N, "Labels"=labels, "Diag"=FALSE, "Upper"=FALSE)
	outDM
}


as.0.100.percent <- function( aaTbl) {
	
	# turn to intger percents, and gaurantee they sum to exactly 100
	aaPct <- round( aaTbl * 100 / sum(aaTbl))

	# slight chance that rounding causes non-100 total
	roundError <- 100 - sum(aaPct, na.rm=TRUE)
	if ( roundError != 0) {
		whoModify <- which.max( aaPct)
		aaPct[whoModify] <- aaPct[whoModify] + roundError
	}
	aaPct
}


calcOneSeqPctDistance <- function( p1, p2, same, distMatrix) {

	# only the AA that are not the same percentage go into the distance measure
	# we want the counts of how each AA differs
	who <- which( ! same)
	dif <- p1[who] - p2[who]

	# turn this into a vector of the AA that are greater in each isotate
	up1 <- which( dif > 0)
	up2 <- which( dif < 0)
	aaSet1 <- rep( AA_LEVELS[up1], times=dif[up1])
	aaSet2 <- rep( AA_LEVELS[up2], times=abs(dif[up2]))

	# use a string form of their combinations to know the counts of each pair
	aaKey <- paste( aaSet1, aaSet2, sep=":")
	aaKeyPatt <- unique( aaKey)
	totalDist <- 0
	for (k in aaKeyPatt) {
		who <- which( aaKey == k)
		N <- length(who)
		aa1 <- aaSet1[ who[1]]
		aa2 <- aaSet2[ who[1]]
		thisDist <- distMatrix[ aa1, aa2]
		totalDist <- totalDist + N*thisDist
	}
	# there are always 100 AA, so divide by that to get final distance
	return( totalDist / 100)
}


calcOneSeqPctPvalue <- function( p1, p2, nonZero, N=100) {

	# all AA that are above zero percent go into the proportion test
	m <- matrix( c( p1[nonZero], p2[nonZero]), nrow=sum(nonZero), ncol=2)

	# we may choose to scale the probabilities to account for how many samples...  not implemented yet...
	# 'N' is total number of samples

	ans <- prop.test( m)
	return( ans$p.value)
}


asDistanceMatrix <- function( substitutionMatrix) {

	N <- ncol(substitutionMatrix) + 1
	mOut <- mIn <- matrix( -5, nrow=N, ncol=N) # give the GAP a weight of -5
	mIn[ 1:(N-1), 1:(N-1)] <- substitutionMatrix
	mIn[N,N] <- 5
	colnames(mOut) <- rownames(mOut) <- c( colnames(substitutionMatrix), GAP)

	for ( i in 1:N) {
		v <- mIn[,i]
		vMax <- max( v)
		vNew <- -(v - vMax)
		mOut[,i] <- vNew
		mOut[i,] <- vNew
	}
	mOut
}


visualizeGroupDifference <- function( tbl, label="") {

	checkX11( "bg"='white', width=12, height=7)
	devSize <- dev.size( units="in")
	lettersWidth <- strwidth( paste(LETTERS,collapse=""), units="in")
	lettersHeight <- strheight( paste(LETTERS,collapse=""), units="in")

	# given a 'small segment' of the group difference result, try to plot it...
	N <- nrow(tbl)
	x <- as.numeric( tbl$Position)
	aa <- tbl$ConsensusAA
	grpNames <- colnames(tbl)[GRP1_COL:GRP2_COL]
	str1 <- tbl[,GRP1_COL]
	str2 <- tbl[,GRP2_COL]
	pvals <- tbl$Pvalue

	# distance of zero means both agree
	agrees <- which( tbl$Distance == 0)
	whoDif <- which( tbl$Distance > 0)

	# parse the AA profiles bck to AA letters and percent sizes
	aaCntList <- lapply( 1:N, function(i) {
				if ( i %in% agrees) return( list( "AA1"=aa[i], "PCT1"=100, "AA2"=aa[i], "PCT2"=100))
				terms1 <- strsplit( str1[i], split=";")[[1]]
				aaTerms1 <- strsplit( terms1, split=":")
				aa1 <- sapply( aaTerms1, `[`, 1)
				pct1 <- as.numeric( sapply( aaTerms1, `[`, 2))
				terms2 <- strsplit( str2[i], split=";")[[1]]
				aaTerms2 <- strsplit( terms2, split=":")
				aa2 <- sapply( aaTerms2, `[`, 1)
				pct2 <- as.numeric( sapply( aaTerms2, `[`, 2))
				return( list( "AA1"=aa1, "PCT1"=pct1, "AA2"=aa2, "PCT2"=pct2))
		})
	
	# letter coloring...
	AA_Color <- rep.int( 'black', length(AA_LEVELS))
	AA_Color[ AA_LEVELS %in% c("A","C","F","I","L","V","W","M")] <- 'blue'
	AA_Color[ AA_LEVELS %in% c("N","Q","S","T")] <- 'green'
	AA_Color[ AA_LEVELS %in% c("D","E")] <- 'magenta'
	AA_Color[ AA_LEVELS %in% c("K","R")] <- 'red'
	AA_Color[ AA_LEVELS %in% c("H")] <- 'deeppink'
	AA_Color[ AA_LEVELS %in% c("G")] <- 'orange'
	AA_Color[ AA_LEVELS %in% c("P")] <- 'yellow'
	AA_Color[ AA_LEVELS %in% c("Y")] <- 'turquoise'

	# scaling to fit the window:  total inches of window divided by how many inches of characters
	myCEXx<- (devSize[1] * 0.85) / ( N * lettersWidth / 26) 
	myCEXy<- (devSize[2] * 0.8) / lettersHeight
	myCEX <- min( myCEXx, myCEXy)
	yMag <- ((devSize[2] * 0.8) / (lettersHeight * myCEX * 2)) * 0.95

	plot( 1,1, type="n", xlim=range(x), ylim=c(-yMag,yMag), yaxt="n", ylab=NA, xlab="Amino Acid Position",
			main=paste( "Consensus Protein Difference:   ", label))
	lines( range(x)+c(-2,2), c(0,0), col='gray50', lwd=1, lty=1)

	for ( i in 1:N) {
		if ( i %in% agrees) {
			text( x[i], 0, aa[i], col=AA_Color[match(aa[i],AA_LEVELS)], cex=myCEX, font=2)
			next
		}
		# grab the AA and PCT details at this location
		tmp <- aaCntList[[i]]
		n1 <- length( tmp$AA1)
		n2 <- length( tmp$AA2)
		# special case if 'almost the same'...
		if ( min( tmp$PCT1[1], tmp$PCT2[1]) >= 95 && tmp$AA1[1] == tmp$AA2[1]) {
			aaNow <- tmp$AA1[1]
			text( x[i], 0, aaNow, col=AA_Color[match(aaNow,AA_LEVELS)], cex=myCEX, font=2)
			next
		}
		# group 1 goes up, group 2 goes down
		yNow <- 0.8
		for ( j in 1:n1) {
			aaNow <- tmp$AA1[j]
			pctNow <- sqrt( tmp$PCT1[j] / 100)
			colNow <- AA_Color[match(aaNow,AA_LEVELS)] 
			text( x[i], yNow, aaNow, col=colNow, cex=myCEX*pctNow, font=2)
			yNow <- yNow + pctNow + 0.1
		}
		yNow <- -0.8
		for ( j in 1:n2) {
			aaNow <- tmp$AA2[j]
			pctNow <- sqrt( tmp$PCT2[j] / 100)
			colNow <- AA_Color[match(aaNow,AA_LEVELS)] 
			text( x[i], yNow, aaNow, col=colNow, cex=myCEX*pctNow, font=2)
			yNow <- yNow - pctNow - 0.1
		}
		if ( pvals[i] <= 0.05) {
			nstar <- 1
			if ( pvals[i] <= 0.005) nstar <- 2
			if ( pvals[i] <= 0.0005) nstar <- 3
			if ( pvals[i] <= 0.00005) nstar <- 4
			dy <- yMag * 0.08
			for ( j in 1:nstar) { text( c(x[i],x[i]), c(yMag-dy,-(yMag-dy)), "*", col=1, cex=2, font=2); dy <- dy + (yMag*0.06) }
		}
	}

	text( rep.int(mean(x),2), c(yMag,-yMag), grpNames, col=1, cex=2, font=2)

	legend( "topleft", paste( "P <", c(0.05,0.005,0.0005,0.00005), " = ", c( "   *","  **", " ***", "****")), cex=1.2, bg='white')
	return(NULL)
}
