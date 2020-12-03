# pipe.BarcodeSNPs.R

`pipe.BarcodeMotif` <- function( sampleID, annotationFile="Annotation.txt", optionsFile="Options.txt", 
				speciesID=getCurrentSpecies(), results.path=NULL, reload=FALSE, 
				nFDR=2000, makePlots=TRUE, verbose=TRUE) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	target <- getAndSetTarget( optionsFile, sampleID=sampleID, annotationFile=annotationFile)

	# set up for this species
	needSetup <- FALSE
	if ( speciesID != getCurrentSpecies()) {
		setCurrentSpecies(speciesID)
		needSetup <- TRUE
	}
	prefix <- getCurrentSpeciesFilePrefix()
	geneMap <- getCurrentGeneMap()

	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	} else {
		results.path <- results.path
	}
	genomeFastaFile <- getOptionValue( optT, "genomicFastaFile", verbose=FALSE)

	barcode.path <- file.path( results.path, "BarcodeMotifs", sampleID)
	if ( ! file.exists( barcode.path)) dir.create( barcode.path, recursive=TRUE)

	# set up if needed
	if ( ! exists( "BarcodeSNPs") || needSetup) setupBarcodeSNPs()

	cat( "\n\nProcessing Sample:  ", sampleID, "\n")

	barcodeDetailsFile <- file.path( barcode.path, paste( sampleID, prefix, "AllBarcodeSites.csv", sep="."))
	# make sure existing data matches the SNP set
	wrongSize <- TRUE
	if ( file.exists(barcodeDetailsFile)) {
		bamAns <- read.csv( barcodeDetailsFile, as.is=T)
		wrongSize <- ( nrow(bamAns) != nrow(BarcodeSNPs))
	}
	readBAM <- (reload || ! file.exists(barcodeDetailsFile) || wrongSize)

	if (readBAM) {
		# interogate this samples BAM file at all barcode marker sites
		bamfile <- file.path( results.path, "align", paste( sampleID, "genomic.sorted.bam", sep="."))
		bamAns <- getBarcodeFromBAMfile( bamfile, genomeFastaFile=genomeFastaFile)
		if ( is.null( bamAns)) return(NULL)
		write.table( bamAns, barcodeDetailsFile, sep=",", quote=T, row.names=F)
	} else {
		bamAns <- read.csv( barcodeDetailsFile, as.is=T)
		# it is sorted by Position, but force it to be sure
		ord <- order( bamAns$SEQ_ID, bamAns$POSITION) 
		bamAns <- bamAns[ ord, ]
	}

	barcodeMutationsFile <- file.path( barcode.path, paste( sampleID, prefix, "MutatedBarcodeSites.csv", sep="."))
	mutations <- subset( bamAns, IsMutant | (!IsClonal))
	write.table( mutations, barcodeMutationsFile, sep=",", quote=T, row.names=F)

	# see which 'lab line' barcode we are most like
	motif <- extractMotifFromBAMans( bamAns) 
	barcodeAns <- bestBarcodeFromMotif( motif)
	if ( is.null( barcodeAns)) {
		# exit early with what we have...
		return( list( "Motif"=motif, "Clonal_YesNo"=table( bamAns$IsClonal), 
				"Mutant_YesNo"=table( bamAns$IsMutant), "Avg_Pct_Clonal"=mean( bamAns$CalledPct,na.rm=T)))
	}

	myBarcode <- barcodeAns$details$BarcodeID[1]
	myDistance <- barcodeAns$details$EditDistance[1]

	# estimate a False Discovery Rate, and append it to the final results
	myFDR <- barcodeFDR( motif, barcodeAns, nSimulations=nFDR)
	fdrOut <- rep.int( NA, nrow( barcodeAns$details))
	fdrOut[1] <- myFDR
	barcodeAns$details$FDR <- fdrOut

	resultFile <- file.path( barcode.path, paste( sampleID, prefix, "FinalBarcodeCalls.csv", sep="."))
	write.table( barcodeAns$details, resultFile, sep=",", quote=T, row.names=F)

	# measure for any mixed infection
	if (FALSE) {
		minorBarcodes <- measureInfectionMixture( sampleID, barcode=barcodeAns, results.path=results.path, makePlots=makePlots)
		if (makePlots) {
			pngFile <- file.path( barcode.path, paste( sampleID, "MinorBarcode.Frequency.Plot.png", sep="."))
			dev.print( png, pngFile, width=1000, height=700)
		}
	}

	return( list( "Motif"=motif, "BarcodeID"=myBarcode, "EditDistance"=myDistance, "FDR"=myFDR))
}


`setupBarcodeSNPs` <- function() {

	BarcodeSNPs <<- loadBarcodeSNPs()
	# BarcodeMotifs <<- makeAllBarcodeMotifs()
}

	
# load the reference set of all SNPs that are barcode markers
`loadBarcodeSNPs` <- function() {

	prefix <- getCurrentSpeciesFilePrefix()
	f <- paste( prefix, "Barcode.SNPs", sep=".")
	data( list=list(f), package="DuffyTools", envir=environment())

	out <- snpTbl
	out
}


extractMotifFromBAMans <- function( bamAns) {

	# make sure we respect any 'Excluded' barcode markers
	keep <- 1:nrow(bamAns)
	if ( "Exclude" %in% colnames(BarcodeSNPs)) {
		keep <- which( BarcodeSNPs$Exclude == FALSE)
	}
	motif <- paste( bamAns$CalledBase[keep], collapse="")
	return( motif)
}


bestBarcodeFromMotif <- function( motif, verbose=T) {

	# we have a fixed length string, to compare to all known motifs
	# regardless of 'N' from missing data, a simple edit distance is fine

	# for now the 'known' motifs are extra columns in the SNP matrix
	motifColumns <- grep( "^Motif", colnames(BarcodeSNPs))
	if ( ! length( motifColumns)) return(NULL)
	N <- length( motifColumns)
	editDist <- rep.int( NA, N)

	myMotifV <- strsplit( motif, split="")[[1]]

	# disount 'N's from no coverage, so they don't raise the distance
	# there 'could' be N's in the reference motifs too...
	isN <- which( myMotifV == "N")

	for ( i in 1:N) {
		k <- motifColumns[i]
		thatMotifV <- BarcodeSNPs[[k]]
		thatN <- which( thatMotifV == "N")
		editDist[i] <- sum( myMotifV != thatMotifV, na.rm=T)
		nN <- length( setdiff( isN, thatN))
		if ( nN > 0) editDist[i] <- editDist[i] - nN
	}
	names( editDist) <- sub( "^Motif_", "", colnames(BarcodeSNPs)[motifColumns])
	editDist <- sort( editDist, decreasing=F)

	bestDist <- editDist[1]
	whoBest <- which( editDist == bestDist)
	if ( length( whoBest) > 1 && verbose) {
		cat( "\nNon-Unique Barcode:  Can't distinquish: \n")
		print( editDist[whoBest])
	}
	bestBarcode <- names(editDist)[whoBest]

	# show it?
	details <- data.frame( "BarcodeID"=names(editDist), "EditDistance"=as.numeric(editDist),
				stringsAsFactors=F)
	rownames(details) <- 1:nrow(details)
	if (verbose) {
		cat( "\nTop Barcode Matches:\n")
		print( head( details, 4))
	}

	out <- list( "barcode"=bestBarcode, "distance"=bestDist, "details"=details)
	return(out)
}


barcodeFDR <- function( motif, barcodeAns=NULL, nSimulations=2000) {

	if ( nSimulations < 1) return(NA)

	# given a sample motif, accumulate random mutations to see
	# how often other barcodes might get called by chance
	if ( is.null( barcodeAns)) barcodeAns <- bestBarcodeFromMotif( motif, verbose=F)
	bestBarcode <- barcodeAns$barcode[1]

	# set up to do the random permutations
	motifLen <- nchar(motif)
	randomWhere <- round( runif( nSimulations, 0.51, motifLen+0.49))
	randomBase <- c("A","C","G","T")[ round( runif( nSimulations, 0.51, 4.49))]
	motifNow <- motif
	nBetterHits <- 0
	cat( "\nCalculating Barcode FDR..\n")
	for ( i in 1:nSimulations) {
		j <- randomWhere[i]
		ch <- randomBase[i]
		substr(motifNow,j,j) <- ch
		# see who is best now
		ansNow <- bestBarcodeFromMotif( motifNow, verbose=F)
		ansDF <- ansNow$details
		where <- which( ansDF$BarcodeID == bestBarcode)[1]
		myDist <- ansDF$EditDistance[where]
		if ( i %% 20 == 0 || where > 1) cat( "\r", i, ansDF$BarcodeID[1], ansDF$EditDistance[1])
		if ( where > 1) {
			nBetterHits <- sum( ansDF$EditDistance[1:where] < myDist)
			if (nBetterHits) break
		}
	}
	# we either broke out early due to better, or never found better
	nTrial <- i
	# count up how often did someone have a better scoring barcode
	fdr <- round( nBetterHits / i, digits=4)
	return( fdr)
}


# extract from a BAM file all the base calls at the Barcode sites
getBarcodeFromBAMfile <- function( bamfile, genomeFastaFile, verbose=T) {

	if ( ! file.exists(bamfile)) {
		cat( "BAM file not found: ", bamfile)
		return( NULL)
	}
	# we will build a N by 4 matrix of depth counts
	N <- nrow( BarcodeSNPs)
	depthM <- matrix( 0, nrow=N, ncol=4)
	colnames(depthM) <- c( "A","C","G","T")
	rownames(depthM) <- paste( BarcodeSNPs$SEQ_ID, BarcodeSNPs$POSITION, sep="::")
	motifV <- rep.int( "N", N)
	if (verbose) cat( "\nExtracting Base Calls from", N, "Barcode marker sites..\n")
	for ( i in 1:N) {
		thisSeqID <- BarcodeSNPs$SEQ_ID[i]
		thisPos <- BarcodeSNPs$POSITION[i]
		thisRef <- BarcodeSNPs$REF_ALLELE[i]
		showBAMmessage <- (i == 1 && verbose)
		bamAns <- BAM.mpileup( bamfile, seqID=thisSeqID, fastaFile=genomeFastaFile, start=thisPos, stop=thisPos, 
					min.depth=1, summarize=TRUE, verbose=F)
		# a total fail is different than an empty result
		if ( is.null( bamAns)) return(NULL)
		# check for no depth coverage
		if ( ! nrow(bamAns)) next
		motifV[i] <- bamAns$CALL_BASE[1]

		# a deleted base comes back as an empty character call.  
		# Set it to be a '-' gap character to keep the size of the final notif correct
		if ( motifV[i] %in% c( "", " ", "*")) motifV[i] <- "-"

		# also a slight chance of having an insertion called.  Just keep the first base
		if ( nchar( motifV[i]) > 1) motifV[i] <- substr( motifV[i], 1,1)

		# also disect the full set of base calls
		myBaseTbl <- MPU.callStringToTable( bamAns$BASE_TABLE[1])[[1]]
		whoRef <- which( names(myBaseTbl) == ",")
		if ( length(whoRef)) names(myBaseTbl)[whoRef[1]] <- thisRef
		whereM <- match( colnames(depthM), names(myBaseTbl), nomatch=0)
		depthM[ i, whereM > 0] <- myBaseTbl[ whereM]
		if (verbose && i %% 50 == 0) cat( "\r", i, thisSeqID, thisPos, motifV[i], bamAns$BASE_TABLE[1])
	}
	if (verbose) cat( "\nDone.")

	# interogate the depth matrix to get all the details we need
	colnames(depthM) <- paste( colnames(depthM), "Depth", sep="_")
	ttlDepth <- apply( depthM, MARGIN=1, sum)
	pctM <- depthM
	colnames(pctM) <- sub( "Depth", "Pct", colnames(pctM))
	for ( i in 1:nrow(pctM)) pctM[ i, ] <- round( depthM[ i, ] * 100 / ttlDepth[i], digits=1)
	maxPct <- apply( pctM, MARGIN=1, max)

	# seems to be a small chance that the called base is empty
	isNoCall <- which( motifV == "")
	if ( length(isNoCall)) {
		whomax <- apply( pctM[ isNoCall, , drop=F], 1, which.max)
		motifV[ isNoCall] <- colnames(pctM)[ whomax]
	}
	isClonal <- ( maxPct > 90.0)
	isMutate <- ( BarcodeSNPs$REF_ALLELE != motifV)
	# the SNP table has Major/Minor calls, and a reference.  Turn this into what is the 'expected' variant
	barcodeBase <- vector( length=N)
	for ( i in 1:N) {
		barcodeBase[i] <- setdiff( c( BarcodeSNPs$MAJOR_ALLELE[i], BarcodeSNPs$MINOR_ALLELE[i]), BarcodeSNPs$REF_ALLELE[i])[1]
	}
	barcodeBase[ is.na(barcodeBase)] <- ""
	expectedMutate <- (barcodeBase == motifV)

	isMutate[ motifV == "N"] <- NA
	expectedMutate[ motifV == "N"] <- NA

	out <- cbind( BarcodeSNPs[,1:6], "CalledBase"=motifV, "CalledPct"=as.numeric(maxPct), "IsMutant"=isMutate, 
			"IsBarcode"=expectedMutate, "IsClonal"=isClonal, depthM, pctM, 
			stringsAsFactors=F)
	# it is sorted by Position, but force it to be sure
	ord <- order( out$SEQ_ID, out$POSITION)
	out <- out[ ord, ]
	rownames(out) <- 1:N
	return( out)
}


calcInfectionMixure <- function( sampleID, results.path="./results") {

	# set up if needed
	if ( ! exists( "BarcodeSNPs")) setupBarcodeSNPs()

	# read in the details of base counts
	barcode.path <- file.path( results.path, "BarcodeSNPs", sampleID)
	barcodeDetailsFile <- file.path( barcode.path, paste( sampleID, "AllBarcodeSites.csv", sep="."))
	bamAns <- read.csv( barcodeDetailsFile, as.is=T)

	# we will inspect the base depth counts
	cntsM <- as.matrix( bamAns[ , grep("_Depth",colnames(bamAns))])
	bigCnt <- apply( cntsM, 1, max, na.rm=T)
	totalCnt <- apply( cntsM, 1, sum, na.rm=T)

	barcodeFac <- factor( bamAns$BarcodeID)
	nSites <- tapply( 1:nrow(bamAns), barcodeFac, length)

	# now measure the minor allele frequencies 'p' as M/D (minor count / total count)
	pCnts <- tapply( 1:nrow(bamAns), barcodeFac, function(x) {
			myD <- sum( bigCnt[x], na.rm=T)
			myT <- sum( totalCnt[x], na.rm=T)
			myM <- myT - myD
			return( myM / myT)
		})

	# also measure the minor allele frequencies 'p' as the average of the percentages
	pPcts <- tapply( 1:nrow(bamAns), barcodeFac, function(x) {
			myD <- bigCnt[x]
			myT <- totalCnt[x]
			myM <- myT - myD
			myPcts <- myM / myT
			# ignore any divide by zeros
			myPcts <- myPcts[ myT > 0]
			return( mean( myPcts, na.rm=T))
		})

	# sort from most mixed to least
	pCnts <- round( pCnts, digits=4)
	pPcts <- round( pPcts, digits=4)
	pFinal <- round( (pCnts + pPcts) / 2, digits=4)

	out <- data.frame( "BarcodeID"=levels(barcodeFac), "N_Sites"=nSites, "MinorFreq"=pFinal, 
			"MinorFreq_ByCounts"=pCnts, "MinorFreq_ByPct"=pPcts,
			stringsAsFactors=F)
	ord <- order( pFinal, decreasing=T)
	out <- out[ ord, ]
	rownames(out) <- 1:nrow(out)
	return( out)
}


measureInfectionMixture <- function( sampleID, barcode=NULL, results.path="./results", F2.cutoff=0.05, makePlots=TRUE) {

	# first step is to calculate all the minor allele frequencies
	ans1 <- calcInfectionMixure( sampleID,  results.path=results.path)
	outfile <- file.path( results.path, BarcodeResultsFolder, paste( sampleID, "MinorBarcode.Frequency.csv", sep="."))
	write.table( ans1, outfile, sep=",", quote=T, row.names=F)

	N <- nrow(ans1)
	x <- 1:N
	y <- ans1$MinorFreq

	# show it visually
	if (makePlots) {
		nams <- ans1$BarcodeID
		nams <- sub( "Barcode", "L", nams)
		yMax <- max( y*1.1, 0.20, na.rm=T)
		yMin <- (yMax/10) * -1
		yTextShift <- diff( range( c( yMin, yMax))) * 0.06
		mainText <- paste( "Minor Allele Frequencies for all Barcode Sets:  ", sampleID)

		par( mai=c(1,1,0.8,0.4))
		plot( x, y, type='p', pch=21, bg=1, cex=0.9, ylim=c(yMin,yMax), 
				xlab="Barcode Sets", ylab="Minor Allele Frequency", 
				main=mainText, yaxt='n')
		axis( side=2, at=seq( 0, yMax, by=0.05))
	
		# add some lines that represent 'random minor noise'
		lines( c(0,N+1), c( 0.01, 0.01), lty=2, lwd=1, col='brown')
		text( N-6, 0.01, "Random Sequencing Error Rate (1%)", pos=3, cex=0.65, col='brown')
		lines( c(0,N+1), c( F2.cutoff, F2.cutoff), lty=2, lwd=1, col='blue')
		text( N-6, F2.cutoff, "Second Minor Barcode Threshold (5%)", pos=3, cex=0.7, col='blue')

		# add labels
		drawVert <- which( y <= 0.02)
		if( length( drawVert)) text( x[drawVert]-0.25, y[drawVert]-yTextShift, padStringWidth(nams[drawVert],pad='left'), 
						srt=90, pos=1, cex=0.7)
		drawHorz <- which( y > 0.02)
		if( length( drawHorz)) text( x[drawHorz], y[drawHorz], nams[drawHorz], pos=4, cex=0.9)

		# show which are the called barcode?
		if ( ! is.null( barcode)) {
			bestCalls <- barcode$barcode
			who <- which( ans1$BarcodeID %in% bestCalls)
			if ( length(who)) points( x[who], y[who], pch=21, cex=1.3, bg='red')
			# also flag other top contenders
			detailDF <- barcode$details
			#whereRef <- match( "H37Rv", detailDF$BarcodeID, nomatch=0)
			whereRef <- match( "Barcode 4.9", detailDF$BarcodeID, nomatch=0)
			if ( whereRef > 1) {
				others <- setdiff( detailDF$Barcode[1:(whereRef-1)], bestCalls)
				if ( length( others)) {
					who <- which( ans1$BarcodeID %in% others)
					if ( length(who)) points( x[who], y[who], pch=21, cex=1.1, bg='orange')
					legend( 'topright', c("Called Barcode", "Other Contenders"), pch=21, pt.bg=c('red','orange'), bg='white')
				} else {
					legend( 'topright', c("Called Barcode"), pch=21, pt.bg=c('red'), bg='white')
				}
			}
			shortCalls <- sub( "Barcode", "L", bestCalls)
			legend( 'topleft', paste( "Barcode: ", paste(shortCalls,collapse=", ")), bg='white')
		}
	}

	# calc the F2 and F47 metrics
	F2 <- round( mean( y[1:2], na.rm=T), digits=3)
	F47 <- round( mean( y[(N-46):N], na.rm=T), digits=3)
	infCall <- "'Unmixed'"
	if ( F47 > 0.02) infCall <- "'Contamination'"
	if ( F2 > F2.cutoff && F47 <= 0.02) infCall <- "'Mixed Barcodes'"
	callString <- c( paste( "'F2' metric =    ", F2), paste("'F47' metric =   ", F47), paste("Infection call =", infCall))

	if (makePlots) {
		legend( 'top', callString, bg='white')
		dev.flush()
	}

	# say it to standard out too...
	cat( "\nFinal Call:    ", sampleID, "\n")
	cat( callString, sep="\n")
	cat( "\n")

	# if it is mixed, send back the barcodes that were too high
	otherBarcodes <- ans1$BarcodeID[ which( y > F2.cutoff)]
	if ( length(otherBarcodes)) {
		return( otherBarcodes)
	} else {
		return( NULL)
	}
}


# grab the barcode calls previously done
extractBarcodeCalls <- function( sampleIDset, results.path="./results") {

	barcode.path <- file.path( results.path, BarcodeResultsFolder)

	# read each answer and grab the first entry
	outSID <- outLin <- outDist <- outFDR <- outMixed <- outMinorLin <- outMinorPct <- vector()
	for( i in 1:length(sampleIDset)) {
		 s <- sampleIDset[i]
		f <- file.path( barcode.path, paste( s, "FinalBarcodeCalls.csv", sep="."))
		if ( ! file.exists(f)) {
			cat( "\nNo Barcode final call file found.  Tried: ", f)
			next
		}
		tbl <- read.csv( f, as.is=T)
		outSID[i] <- s
		outLin[i] <- bestBarcode <- tbl$BarcodeID[1]
		outDist[i] <- tbl$EditDistance[1]
		if ( "FDR" %in% colnames(tbl)) outFDR[i] <- tbl$FDR[1]

		# also grab the minor barcode info
		f2 <- file.path( barcode.path, paste( s, "MinorBarcode.Frequency.csv", sep="."))
		if ( ! file.exists(f2)) next
		tbl2 <- read.csv( f2, as.is=T)
		# don't let the major barcode call come up here
		majorRow <- which( tbl2$BarcodeID == bestBarcode)
		if ( length(majorRow)) tbl2 <- tbl2[ -majorRow, ]
		minorPct <- tbl2$MinorFreq[1]
		if ( minorPct >= 0.05) {
			outMixed[i] <- "Mixed"
			outMinorPct[i] <- minorPct
			outMinorLin[i] <- tbl2$BarcodeID[1]
		} else {
			outMixed[i] <- "Clonal"
			outMinorPct[i] <- 0
			outMinorLin[i] <- ""
		}
	}
	out <- data.frame( "SampleID"=outSID, "BarcodeID"=outLin, "EditDistance"=outDist, 
				"FDR"=outFDR, "Population"=outMixed, "MinorID"=outMinorLin, 
					"MinorPct"=round( outMinorPct * 100, digits=1),
					stringsAsFactors=F)
	out
}

