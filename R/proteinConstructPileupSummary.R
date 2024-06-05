# proteinConstructPileupSummary.R -- summarize the pileup results to give the best one protein 
#								and it's supporting details


proteinConstructPileupSummary <- function( constructSaveFile, sampleID, geneName="Var2csa", 
					constructName=paste(sampleID,geneName,sep="."), txt.cex=0.25,
					doPlot=TRUE, verbose=doPlot, summaryRange=NULL, intronMaskFasta=NULL) {

	NCHAR <- base::nchar

	# grab the Consensus Answer data object 
	consensusAns <- NULL
	load( constructSaveFile)
	if (is.null(consensusAns)) {
		cat( "\nError: unable to load consensus answer for summarizing.. ", sampleID, geneName)
		return(NULL)
	}
	
	# if we are given intron masking, do that now, silently
	if ( ! is.null( intronMaskFasta)) {
		peptide.path <- dirname(constructSaveFile)
		refAA <- paste( consensusAns$Construct, collapse="")
		intronMaskInfo <- IntronMaskSetup( intronMaskFasta, refAA=refAA, 
						geneName=geneName, path=peptide.path, verbose=TRUE)
		if ( ! is.null( intronMaskInfo)) {
			consensusAns <- IntronMaskAdjustPileups( consensusAns, maskInfo=intronMaskInfo)
		}
	}
	
	# now grab the components and evaluate
	aaCalls <- consensusAns$AA_Calls
	aaWeights <- consensusAns$AA_Weights
	seqAA <- consensusAns$Construct
	construct <- paste( seqAA, collapse="")

	nx <- ncol( aaCalls)
	ny <- nrow( aaCalls)

	pctAA <- bestAA <- strAApcts <- strAAcnts <- depth <- nAAcalls <- vector( length=nx)
	for ( i in 1:nx) {
		myAAs <- aaCalls[ , i]
		myWts <- aaWeights[ , i]
		drops <- which( myAAs == "")
		if ( length( drops)) {
			myAAs <- myAAs[ -drops]
			myWts <- myWts[ -drops]
		}
		if ( length(myAAs)) {
			bigSetAA <- rep( myAAs, times=myWts)
			myFreqs <- sort( table( bigSetAA), decreasing=T)
			bestAA[i] <- names( myFreqs)[1]
			allCnts <- as.numeric( myFreqs)
			allPcts <- allCnts * 100 / sum(allCnts)
			pctAA[i] <- allPcts[1]
			allPcts <- round( allPcts)
			names(allPcts) <- names(allCnts) <- names(myFreqs)
			keep <- which( allPcts > 0)
			allPcts <- allPcts[ keep]
			allCnts <- allCnts[ keep]
			nAAcalls[i] <- length(keep)
			strAApcts[i] <- paste( names(allPcts), as.vector(allPcts), sep=":", collapse="; ")
			strAAcnts[i] <- paste( names(allCnts), as.vector(allCnts), sep=":", collapse="; ")
			depth[i] <- sum(allCnts)
		} else {
			bestAA[i] <- "X"
			pctAA[i] <- 0
			strAApcts[i] <- strAAcnts[i] <- "X:0"
			depth[i] <- 0
			nAAcalls[i] <- 0
		}
	}
	nAA <- length( seqAA)
	len <- min( nAA, nx)
	same <- rep.int( FALSE, nx)
	same[ 1:len] <- (seqAA[1:len] == bestAA[1:len])
	domainTxt <- rep.int( "", length(strAApcts))

	if (doPlot) {
		plot( 1:nx, pctAA, main=paste( "Final Consensus:    ", constructName, "   Gene: ", geneName), 
				type="l", col=2, lwd=1, xlim=c( -10, nx+10), ylim=c( 0, 105), xlab="Amino Acid Location", 
				xaxs="i", ylab="Percent Consensus")

		# put ticks and counts along top 
		tickAts <- numberAts <- 1
		if (nAA > 10) tickAts <- c( 1, seq( 10, nAA, by=10))
		if (nAA > 50) numberAts <- c( 1, seq( 50, nAA, by=50))
		axis( side=3, at=tickAts, label=FALSE, tcl=-0.2)
		axis( side=3, at=numberAts, label=TRUE, tcl=-0.2, padj=3.3, cex.axis=txt.cex*2)
	
		# put domain boxes on?
		if ( geneName == "Var2csa") {
			yLocs <- c( 107)
			rectHalfHt <- 1.0
			domAns <- findVar2csaDomains( aaSeq=construct, minScorePerAA=2, one.time=F, plot.Y=yLocs, 
					rect.halfHeight=rectHalfHt, rect.border='white',
					rect.lwd=0.5, rect.fill='gold', text.cex=txt.cex, text.X.repeat=35)
			if ( nrow(domAns)) {
				for ( k in 1:nrow(domAns)) {
					from <- domAns$QUERY_START[k]
					to <- domAns$QUERY_STOP[k]
					domainTxt[ from:to] <- domAns$QUERY_DOMAIN_ID[k]
				}
			}
		}

		text( 1:nAA, rep.int(104,nAA), seqAA, cex=txt.cex, col=1)
		text( 1:nx, rep.int(102,nx), bestAA, cex=txt.cex, col=ifelse( same, 1,2))
	} else {
		# do the domain call no matter what
		if ( geneName == "Var2csa") {
			domAns <- findVar2csaDomains( aaSeq=construct, minScorePerAA=2, one.time=F, plot.Y=NULL)
			if ( nrow(domAns)) {
				for ( k in 1:nrow(domAns)) {
					from <- domAns$QUERY_START[k]
					to <- domAns$QUERY_STOP[k]
					domainTxt[ from:to] <- domAns$QUERY_DOMAIN_ID[k]
				}
			}
		}
	}
	

	# show the disagreements
	notSame <- which( ! same)
	notSame <- notSame[ notSame %in% 1:nAA]
	if ( length(notSame)) {
		badAAlocs <- data.frame( "Location"=notSame, "CurrentAA"=seqAA[notSame], "ConsensusAA"=bestAA[notSame],
					"PctConsensus"=as.percent(pctAA[notSame],big=100,digit=0), stringsAsFactors=F)
		if (verbose) {
			cat( "\nDiscrepency AAs: ", length(notSame), "\n")
			print( head( badAAlocs, 22))
		}
	} else {
		if (verbose) cat( "\nNo Discrepencies !!\n")
		badAAlocs <- data.frame()
	}

	# write the summary table
	nAAout <- min( nAA, nx)
	sumTbl <- data.frame( "Position"=1:nAAout, "ConsensusAA"=bestAA[1:nAAout], "ReadDepth"=depth[1:nAAout],
				"Percentages"=strAApcts[1:nAAout], "Counts"=strAAcnts[1:nAAout], 
				"Domain"=domainTxt[1:nAAout], stringsAsFactors=FALSE)
	sumFile <- file.path( dirname( constructSaveFile), paste( constructName, "ConsensusProteinSummary.txt", sep="."))
	write.table( sumTbl, sumFile, sep="\t", quote=F, row.names=F)

	toEval <- 1:nAAout
	if ( ! is.null( summaryRange)) toEval <- intersect( toEval, summaryRange[1]:summaryRange[2])
	pctAA <- pctAA[ toEval]
	nAAcalls <- nAAcalls[ toEval]

	finalPct <- round( mean( pctAA, na.rm=T), digits=2)
	if (verbose) cat( "\nConstruct: ", constructName, "\tFinal Percent Consensus: ", finalPct, "\n")

	pctPerfect <- round( sum( round(pctAA) == 100) * 100 / length(toEval), digits=2)
	avgNonperfectPercent <- round( mean( pctAA[ pctAA < 99.9], na.rm=T), digits=2)
	avgNvariants <- round( mean( nAAcalls[ pctAA < 99.9]), digits=3)

	return( invisible( list("Discrepencies"=badAAlocs, "FinalConsensusPercentage"=finalPct, "PctPerfect"=pctPerfect, 
			"AvgNonperfectPercent"=avgNonperfectPercent, "AvgVariants"=avgNvariants, "SummaryTable"=sumTbl)))
}
