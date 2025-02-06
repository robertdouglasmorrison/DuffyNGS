# protienConstructPeptidePileups.R -- show evidence that raw read peptides land on a protein construct
#									and use it to extract the consensus best protein


`proteinConstructPeptidePileups` <- function( sampleID, geneName, constructFile, peptide.path=".",
					constructName=paste(sampleID,geneName,sep="_"), 
					txt.cex=0.25, maxNoHits=1000000, max.depth=60, pct.aligned.depth=0.75, 
					max.drawnPerSite=3, mode=c("normal", "realigned"), draw.box=FALSE, chunkSize=20000, 
					showFrameShiftPeptides=TRUE, extraGeneNames=NULL, intronMaskFasta=NULL, ...) {

	SAPPLY <- base::sapply
	WHICH <- base::which
	UNLIST <- base::unlist
	NCHAR <- base::nchar
	SUBSTR <- base::substr
	ORDER <- base::order
	TAPPLY <- base::tapply
	LEVELS <- base::levels

	# get the protein construct from a FASTA file
	if ( ! file.exists( constructFile)) {
		cat( "\nProtein construct FASTA file not found.")
		cat( "\nTried: ", constructFile)
		return(NULL)
	}
	fa <- loadFasta( constructFile, verbose=T)
	who <- intersect( grep ( sampleID, fa$desc, fixed=T), grep( geneName, fa$desc, fixed=T))
	if ( length(who) != 1) {
		who <- grep ( constructName, fa$desc, fixed=T)
		if ( length(who) != 1) {
			cat( "\nFailed to find a single protein construct in FASTA file")
			cat( "\nLooked for: ", constructName, " or ", sampleID, "\nConstruct names found:\n")
			print( fa$desc)
			return(NULL)
		}
	}
	construct <- fa$seq[ who[1]]
	nAA <- NCHAR( construct)
	seqAA <- strsplit( construct, split="")[[1]]
	cat( "\nConstruct: ", constructName, "\tN_AA: ", nAA)
	constructAAstr <- AAString( construct)

	mode <- match.arg( mode)

	# make sure we have the raw reads as peptides ready.  Small chance that gene name could have invalid filename characters
	geneFileName <- file.cleanSpecialCharactersFromFileName( geneName)
	constructFileName <- file.cleanSpecialCharactersFromFileName( constructName)
	
	alignedReadsFile <- file.path( peptide.path, paste( sampleID, geneFileName, "RawReadPeptides.txt", sep="."))
	nohitReadsFile <- file.path( peptide.path, paste( sampleID, "NoHits", "RawReadPeptides.txt", sep="."))
	if ( ! file.exists( alignedReadsFile)) {
		cat( "\nMissing Raw Reads Peptides files for", sampleID, geneName)
		return(NULL)
	}

	# if we were given any intron masking info, set up that data
	if ( ! is.null(intronMaskFasta)) {
		intronMaskInfo <- IntronMaskSetup( intronMaskFasta, refAA=construct, geneName=geneName, path=peptide.path)
	}
	
	# set up storage to determine the final consensus call for every AA, and the current depth on the pileup plot
	maxDepth <- round(max.depth)
	if ( pct.aligned.depth > 1 || pct.aligned.depth < 0.5) stop( "function argument 'pct.aligned.depth' must be in range 0.5 to 1.0")
	alignedReadsDepth <- round( maxDepth * pct.aligned.depth)
	pepFree <- matrix( TRUE, nrow=maxDepth, ncol=nAA+30)
	aaCalls <- matrix( "", nrow=maxDepth, ncol=nAA+30)
	aaWeights <- matrix( 0, nrow=maxDepth, ncol=nAA+30)

	# also store up the good scoring peptides
	outPep <- outCnt <- outStart <- outScore <- outSource <- vector()
	nPepOut <- 0
	pepOutfile <- file.path( peptide.path, paste( constructFileName, "GoodScorePeptides.txt", sep="."))

	# lastly, lets try to track alternate reading frames too, if possible
	constructRF2 <- constructRF3 <- NULL
	baseCallsFile <- file.path( peptide.path, paste( sampleID, geneFileName, "ConsensusBaseCalls.txt", sep="."))
	if ( mode == "realigned") {
		baseCallsFile <- file.path( peptide.path, paste( sampleID, geneFileName, "RealignedConsensusBaseCalls.txt", sep="."))
	}
	if ( file.exists( baseCallsFile)) {
		tmp <- read.delim( baseCallsFile, as.is=T)
		constructRF1 <- paste( tmp$Frame1, collapse="")
		constructRF2 <- paste( tmp$Frame2, collapse="")
		constructRF3 <- paste( tmp$Frame3, collapse="")
		# RF1 'should' match our construct
		nch <- min( nAA, NCHAR( constructRF1))
		if ( SUBSTR(constructRF1,1,nch) != SUBSTR(construct,1,nch)) {
			cat( "\nWarning:  constuct does not match ConsensusBaseCalls AA seq!!\n")
			tmpSeq <- tmp$Frame1
			tmpSeq <- tmpSeq[ tmpSeq != ""]
			bad <- WHICH( seqAA[1:nch] != tmpSeq[1:nch])
			print( table( seqAA[1:nch] == tmpSeq[1:nch]))
			print( head( data.frame( "Loc"=bad, "Fasta"=seqAA[bad], "BaseCalls"=tmpSeq[bad]), 22))
		}
	}

	if ( ! showFrameShiftPeptides) constructRF2 <- constructRF3 <- NULL

	# tiny chance of invalid DNA/AA calls existing in the consensus construct.  Catch to prevent breaking Biostrings compares
	constructAAstr <- AAString( gsub( "?", "A", construct, fixed=T))
	if ( ! is.null( constructRF2)) constructRF2 <- gsub( "?", "A", constructRF2, fixed=T)
	if ( ! is.null( constructRF3)) constructRF3 <- gsub( "?", "A", constructRF3, fixed=T)

	# local function to find first free spot on the plot to lay down a peptide
	findFreeDepth <- function( xfrom, xto, maxUseDepth=maxDepth) {

		xfrom <- max( xfrom, 1)
		xto <- min( xto, ncol(pepFree))
		for ( i in 1:maxDepth) {
			if ( all( pepFree[ i, xfrom:xto])) {
				# save some space for the No Hits reads too
				if (i > maxUseDepth) return(NA)
				pepFree[ i, xfrom:xto] <<- FALSE
				return( i)
			}
		}
		return(NA)
	}


	# set up to plot
	xlim <- c( -10, nAA+25)
	ylim <- c( -1, maxDepth+2)
	par( mai=c(0.8, 0.9, 0.8, 0.3))
	plot( 1,1, type='n', main=paste( "Protein Construct:   ", constructName, "    Gene: ", geneName, 
			"\nPileup of Illumina short reads translated to peptides"),
		xlab="Protein Construct Location  (amino acids)", ylab="Unique Peptides Depth",
		xlim=xlim, ylim=ylim, xaxs="i", xaxt="n")

	# put ticks and counts along top and bottom
	tickAts <- numberAts <- 1
	if (nAA > 10) tickAts <- c( 1, seq( 10, nAA, by=10))
	if (nAA > 50) numberAts <- c( 1, seq( 50, nAA, by=50))
	axis( side=1, at=tickAts, label=FALSE, tcl=-0.2)
	axis( side=1, at=numberAts, label=TRUE, tcl=-0.2, padj=-2.0, cex.axis=0.5)
	axis( side=3, at=tickAts, label=FALSE, tcl=-0.2)
	axis( side=3, at=numberAts, label=TRUE, tcl=-0.2, padj=3.3, cex.axis=txt.cex*2)

	# put the reference below and above
	lines( c(-1, nAA+2), c( 0,0), col=1, lwd=1)
	text( 1:nAA, rep.int(-1,nAA), seqAA, cex=txt.cex)
	lines( c(-1, nAA+2), rep.int(maxDepth+1,2), col=1, lwd=1)
	text( 1:nAA, rep.int(maxDepth+2,nAA), seqAA, cex=txt.cex)

	# if we were given any intron masking info, show it
	if ( ! is.null(intronMaskFasta)) IntronMaskDisplay( intronMaskInfo, cex=txt.cex)

	
	# explain the coloring
	mtext( "BLACK = 'Good' read aligned to Gene", side=3, line=0.5, adj=0, col=1, cex=txt.cex*1.5)
	mtext( " BLUE = 'NoHit' read failed genomic alignment", side=3, line=1, adj=0, col='dodgerblue', cex=txt.cex*1.5)
	mtext( " RED  =  Amino Acid differs from Construct", side=3, line=1.5, adj=0, col=2, cex=txt.cex*1.5)
	mtext( "PURPLE = Frame Shift +1 required to fit", side=3, line=1, adj=1, col='purple', cex=txt.cex*1.5)
	mtext( " GREEN = Frame Shift +2 required to fit", side=3, line=1.5, adj=1, col='green3', cex=txt.cex*1.5)

	# put domain boxes on?
	if ( geneName == "Var2csa") {
		yLocs <- c( -2.4, maxDepth + 3.4)
		rectHalfHt <- 0.6
		domAns <- findVar2csaDomains( aaSeq=construct, minScorePerAA=1.5, one.time=F, plot.Y=yLocs, 
				rect.halfHeight=rectHalfHt, rect.border='white',
				rect.lwd=0.5, rect.fill='gold', text.cex=txt.cex, text.X.repeat=35)
	}
	

	# local function to try to preserve gaps '-' in the peptides
	preserveLocalInsertions <- function( peps, paAns) {

		# given a set of full length original peptides, and the pairwiseAlignment object,
		# try to reinsert local gaps
		outPeps <- peps

		localPatts <- as.character( alignedPattern( paAns))
		hasGap <- grep( "-", localPatts, fixed=TRUE)
		for (i in hasGap) {
			withGaps <- localPatts[i]
			noGaps <- gsub( "-", "", withGaps)
			outPeps[i] <- sub( noGaps, withGaps, peps[i])
		}
		outPeps
	}


	# local function that finds and places peptides
	addPeptidePileups <- function( peps, wt=rep.int(1,length(peps)), chunk.size=10000, txt.col=1, 
					maxUseDepth=maxDepth, max.drawnPerSite=3, 
					source="aligned.read", maxReads=NULL) {
	
		# make sure we combine duplicates
		pepFac <- factor( peps)
		pepCnts <- TAPPLY( wt, pepFac, sum)
		peps <- LEVELS(pepFac)
		wt <- pepCnts
		rm( pepFac, pepCnts)

		# find high scoring peptides, a chunk at a time
		goodScore <- mean( NCHAR( peps)) * 2
		done <- 0
		repeat {
			from <- done + 1
			to <- min( from+chunk.size-1, length( peps))
			if ( to <= done) break
			if ( !is.null(maxReads) && done >= maxReads) break
			chunk <- peps[ from:to]
			chunkWts <- wt[ from:to]
			done <- to
			cat( "\rN_Pep: ", formatC(done, format="d",big.mark=","), "  Aligning..")

			# lets use 2 ways to find alignments, for both speed and to find multiple locations for a read
			# method 1:  perfect matches that see multiple sites
			chunkAAstr <- AAStringSet( chunk)
			mpAns <- matchPDict( chunkAAstr, constructAAstr)
			mpStarts <- startIndex( mpAns)
			mpLen <- SAPPLY( mpStarts, FUN=length)
			hit1 <- WHICH( mpLen > 0)
			if ( length(hit1)) {
				starts1 <- UNLIST( mpStarts[hit1])
				peptides1 <- rep( chunk[hit1], times=mpLen[hit1])
				wts1 <- rep( chunkWts[hit1], times=mpLen[hit1])
				scores1 <- NCHAR( peptides1) * 4
				chunk <- chunk[ -hit1]
				chunkWts <- chunkWts[ -hit1]
				chunkAAstr <- AAStringSet( chunk)
			} else {
				starts1 <- vector()
				peptides1 <- vector()
				wts1 <- vector()
				scores1 <- vector()
			}
			rm( mpAns)

			# method 2:  pairwise with gaps and mismatches, on whats left...
			scores <- pairwiseAlignment( chunkAAstr, constructAAstr, type='local', 
							scoreOnly=T, gapOpen=-6, gapExt=-3)
			hit2 <- WHICH( scores >= goodScore)
			cat( "  Hits (match, mismatch): ", length(hit1), length(hit2))
			if ( length(hit2)) {
				# for each high scoring hit, see where to place it
				paAns <- pairwiseAlignment( chunkAAstr[hit2], constructAAstr, type='local', 
							gapOpen=-6, gapExt=-3)
				starts2 <- start( subject( paAns))
				scores2 <- round( score( paAns))
				peptides2 <- preserveLocalInsertions( chunk[hit2], paAns)
				patternStart <- start( pattern(paAns))
				wts2 <- chunkWts[hit2]
				for ( i in WHICH( patternStart > 1)) {
					correctedStart <- starts2[i] - patternStart[i] + 1
					# trap falling off left side before zero
					if ( correctedStart < 1) {
						peptides2[i] <- SUBSTR( peptides2[i], patternStart[i], NCHAR(peptides2[i]))
					} else {
						starts2[i] <- correctedStart
					}
				}
				rm( paAns)
			} else {
				starts2 <- vector()
				peptides2 <- vector()
				wts2 <- vector()
				scores2 <- vector()
			}
			hits <- c( hit1, hit2)
			if ( length(hits) < 1) next
			myStart <- c( starts1, starts2)
			mySeq <- c( peptides1, peptides2)
			myWts <- c( wts1, wts2)
			myScore <- c( scores1, scores2)

			# let the high weight peptides be drawn first
			cat( "  plotting..")
			ord <- ORDER( myStart, -myWts)

			# let's not always draw every peptide.  If too many start at the same location, it's a sign
			# of low quality read data, etc.  Only keep the first (highest weight) K at each location.
			# But let location 1 have more.
			allStartLocs <- sort( unique( myStart))
			orderedStarts <- myStart[ ord]
			finalOrd <- vector()
			for ( k in allStartLocs) {
				thisNmax <- if (k == 1) max.drawnPerSite*2 else max.drawnPerSite
				thisOrdSet <- which( orderedStarts == k)
				if ( length(thisOrdSet) > thisNmax) length(thisOrdSet) <- thisNmax
				finalOrd <- c( finalOrd, ord[thisOrdSet])
			}
			ord <- finalOrd

			# now draw and remember
			for ( j in ord) {
				# save it for future assembly use...
				nPepOut <<- nPepOut + 1
				outPep[nPepOut] <<-  mySeq[j]
				outCnt[nPepOut] <<-  myWts[j]
				outStart[nPepOut] <<-  myStart[j]
				outScore[nPepOut] <<-  myScore[j] 
				outSource[nPepOut] <<-  source
				# see where to draw it, but skip after we fill the Y axis
				xLoc <- myStart[j]
				xEnd <- xLoc + NCHAR( mySeq[j]) - 1
				if ( xEnd - xLoc < 5) next
				useY <- findFreeDepth( xLoc, xEnd+1, maxUseDepth=maxUseDepth)
				if( is.na( useY)) next
				naa <- xEnd - xLoc + 1
				sameAA <- rep.int( FALSE, naa)
				myAAs <- strsplit( mySeq[j], split="")[[1]] 
				if ( length( myAAs) < 5) next
				xxEnd <- min( length(seqAA), xEnd)
				nOverlap <- xxEnd - xLoc + 1
				sameAA[1:nOverlap] <- (myAAs[1:nOverlap] == seqAA[ xLoc:xxEnd])
				text( xLoc:xEnd, rep.int( useY, (xEnd-xLoc+1)), myAAs, cex=txt.cex, 
						col=ifelse( sameAA, txt.col, 2))
				if (draw.box) rect( xLoc-0.6, useY-0.4, xEnd+0.6, useY+0.4, border=txt.col,
						lwd=0.1)
				# stuff this into our consensus calls matrix too
				myEnd <- min( xEnd, ncol(aaCalls))
				aaCalls[ useY, xLoc:myEnd ] <<- myAAs[ 1:(myEnd-xLoc+1)]
				aaWeights[ useY, xLoc:myEnd ] <<- myWts[j]
			}
			dev.flush()

			# if we have alternate reading frame sequences, 
			# then try to place the peptides not yet placed...
			if ( is.null( constructRF2)) next
			if ( regexpr( "align", source) < 1) next

			# OK, try the other 2 reading frames too
			# with new first pass of 'matchPDict()', the peptides in 'chunk' got reduced...
			# just the 'hit2' hits should be removed here...
			if ( ! length(chunk)) next
			notHit <- setdiff( 1:length(chunk), hit2)
			if ( ! length(notHit)) next
			scoresRF2 <- pairwiseAlignment( chunk[notHit], constructRF2, type='local', scoreOnly=T,
						gapOpen=-6, gapExt=-3)
			scoresRF3 <- pairwiseAlignment( chunk[notHit], constructRF3, type='local', scoreOnly=T,
						gapOpen=-6, gapExt=-3)
			hitsRF2 <- notHit[ WHICH( scoresRF2 >= goodScore & scoresRF2 > scoresRF3)]
			hitsRF3 <- notHit[ WHICH( scoresRF3 >= goodScore & scoresRF3 > scoresRF2)]

			cat( "\nFrame Shifted Petide Hits: \tRF2: ", length( hitsRF2), "\tRF3: ", length(hitsRF3))
			# for each high scoring hit, see where to place it
			hits <- myStart <- patternStart <- myScore <- mySeq <- myWts <- vector()
			pepColor <- pepShift <- vector()
			if ( length(hitsRF2)) {
				ans <- pairwiseAlignment( chunk[hitsRF2], constructRF2, type='local', gapOpen=-6, gapExt=-3)
				tmpStart <- start( subject( ans))
				tmpPatternStart <- start( pattern(ans))
				tmpScore <- round( score( ans))
				tmpSeq <- preserveLocalInsertions( chunk[hitsRF2], ans)
				tmpWts <- chunkWts[hitsRF2]
				for ( i in WHICH( tmpPatternStart > 1)) {
					correctedStart <- tmpStart[i] - tmpPatternStart[i] + 1
					# trap falling off left side before zero
					if ( correctedStart < 1) {
						tmpSeq[i] <- SUBSTR( tmpSeq[i], tmpPatternStart[i], NCHAR(tmpSeq[i]))
					} else {
						tmpStart[i] <- correctedStart
					}
				}
				hits <- hitsRF2
				myStart <- tmpStart
				patternStart <- tmpPatternStart
				myScore <- tmpScore
				mySeq <- tmpSeq
				myWts <- tmpWts
				pepColor <- rep.int( 'purple', length(hitsRF2))
				pepShift <- rep.int( 0.33, length(hitsRF2))
			}
			if ( length(hitsRF3)) {
				ans <- pairwiseAlignment( chunk[hitsRF3], constructRF3, type='local', gapOpen=-6, gapExt=-3)
				tmpStart <- start( subject( ans))
				tmpPatternStart <- start( pattern(ans))
				tmpScore <- round( score( ans))
				tmpSeq <- preserveLocalInsertions( chunk[hitsRF3], ans)
				tmpWts <- chunkWts[hitsRF3]
				for ( i in WHICH( tmpPatternStart > 1)) {
					correctedStart <- tmpStart[i] - tmpPatternStart[i] + 1
					# trap falling off left side before zero
					if ( correctedStart < 1) {
						tmpSeq[i] <- SUBSTR( tmpSeq[i], tmpPatternStart[i], NCHAR(tmpSeq[i]))
					} else {
						tmpStart[i] <- correctedStart
					}
				}
				hits <- c( hits, hitsRF3)
				myStart <- c( myStart, tmpStart)
				patternStart <- c( patternStart, tmpPatternStart)
				myScore <- c( myScore, tmpScore)
				mySeq <- c( mySeq, tmpSeq)
				myWts <- c( myWts, tmpWts)
				pepColor <- c( pepColor, rep.int( 'green3', length(hitsRF3)))
				pepShift <- c( pepShift, rep.int( 0.67, length(hitsRF3)))
			}
			
			# let the high weight peptides be drawn first
			ord <- ORDER( myStart, -myWts)
			for ( j in ord) {
				# see where to draw it
				xLoc <- myStart[j]
				xEnd <- xLoc + NCHAR( mySeq[j]) - 1
				if ( xEnd - xLoc < 5) next
				useY <- findFreeDepth( xLoc, xEnd+1, maxUseDepth=maxUseDepth)
				if( is.na( useY)) next
				naa <- xEnd - xLoc + 1
				if ( naa < 5) next
				myAAs <- strsplit( mySeq[j], split="")[[1]] 
				xxEnd <- min( length(seqAA), xEnd)
				text( (xLoc:xEnd)+pepShift[j], rep.int( useY, (xEnd-xLoc+1)), myAAs, cex=txt.cex, 
						col=pepColor[j])
			}
		}
		# finished adding this type of peptides, show them
		dev.flush()
	}


	# ready to go:
	# get the 'gene' peptides
	if ( file.exists( alignedReadsFile)) {
		cat( "\nReading Gene peptides file: ", alignedReadsFile)
		vgTbl <- read.delim( alignedReadsFile, as.is=T)
		cat( "\nDone.  \tN_Peptides:  ", nrow(vgTbl), "\n")
		# see if we were given any extra gene names to also try to draw, like from highly similar genes
		if ( ! is.null( extraGeneNames)) {
			extraGeneNames <- setdiff( extraGeneNames, geneName)
			for ( extraName in extraGeneNames) {
				extraReadsFile <- file.path( peptide.path, paste( sampleID, extraName, "RawReadPeptides.txt", sep="."))
				if ( ! file.exists( extraReadsFile)) next
				cat( "\nReading Extra Gene peptides file: ", extraReadsFile)
				vgTbl2 <- read.delim( extraReadsFile, as.is=T)
				cat( "\nDone.  \tN_Extra_Peptides:  ", nrow(vgTbl2), "\n")
				vgTbl <- rbind( vgTbl, vgTbl2)
			}
		}
		# we can cap how many we keep, based on protein length, since we can only pile up so deep...
		avgPepLen <- round( mean( nchar( vgTbl$Peptide)))
		maxPeptidesNeeded <- nAA * max.depth * 1.5
		minPeptideCount <- 0
		while ( nrow(vgTbl) > maxPeptidesNeeded) {
			minPeptideCount <- minPeptideCount + 1
			drops <- which( vgTbl$Count <= minPeptideCount)
			if ( length(drops)) {
				vgTbl <- vgTbl[ -drops, ]
				cat( "  Dropping Peptides with Count at/below", minPeptideCount, " N_Peptides now: ", nrow(vgTbl), "\n")
			}
		}
		addPeptidePileups( vgTbl$Peptide, wt=vgTbl$Count, chunk.size=chunkSize, 
				txt.col=1, maxUseDepth=alignedReadsDepth, 
				max.drawnPerSite=max.drawnPerSite, source="aligned.read")
	}

	# also try the 'NoHits'
	if ( file.exists( nohitReadsFile) && maxNoHits > 0) {
		cat( "\nReading NoHits peptides file: ", nohitReadsFile)
		nhTbl <- read.delim( nohitReadsFile, as.is=T, nrows=maxNoHits)
		cat( "\nDone.  \tN_Peptides:  ", nrow(nhTbl), "\n")
		addPeptidePileups( nhTbl$Peptide, wt=nhTbl$Count, chunk.size=100000, 
				txt.col='dodgerblue', maxUseDepth=maxDepth,
				max.drawnPerSite=max.drawnPerSite, source="nohit.read",
				maxReads=maxNoHits)
	}

	# all done
	pepTbl <- data.frame( "Peptide"=outPep, "Count"=outCnt, "AA_Start"=outStart, "Score"=outScore,
				"Source"=outSource, stringsAsFactors=FALSE)
	ord <- ORDER( pepTbl$AA_Start, -pepTbl$Count)
	pepTbl <- pepTbl[ ord, ]
	if (nrow(pepTbl)) rownames(pepTbl) <- 1:nrow(pepTbl)
	write.table( pepTbl, pepOutfile, sep="\t", quote=F, row.names=F)
	totalPeptides <- sum( pepTbl$Count, na.rm=T)

	# create final data object, and summarize
	out <- list( "AA_Calls"=aaCalls, "AA_Weights"=aaWeights, "Construct"=seqAA, "N_Peptides"=totalPeptides)
	cat( "\nConstruct: ", constructName, "\tN_GoodScoring_Peptides: ", 
		formatC( totalPeptides, format="d", big.mark=","), "\n")
	
	# if given intron masking, apply those corrections now to the matrix & weights
	if ( ! is.null(intronMaskFasta)) {
		out <- IntronMaskAdjustPileups( out, maskInfo=intronMaskInfo)
	}
	
	return( out)
}


