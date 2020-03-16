# plotSNP.R -- plot the base call pileups at a genomic location


plotSNP <- function( position, seqID, sampleID, bamfile, vcffile, fastaFile, geneID=NULL, gmap=NULL,
			tailWidth=50, cex.text=1.0, cex.legend=1.0, mode=c("single","multi"), 
			label=sampleID, forceYmax=NULL, linewidth.factor=1, SNPwidth.factor=1.0,
			max.depth=10000, show.legends=c( "all", "gene", "none"), verbose=TRUE, ...) {

	mode <- match.arg( mode)
	show.legends <- match.arg( show.legends)

	# various scaling based on range and counts... and the size of the current window
	figInch <- par("fin")
	xFig <- figInch[1] - 1
	xLo <- max( 1, round( position-tailWidth))
	xHi <- min( subset( getCurrentSeqMap(), SEQ_ID == seqID)$LENGTH, round( position+tailWidth))
	xLim <- c( xLo, xHi)
	xWide <- base::diff( xLim)
	lineWidth <- (600/xWide) * (xFig/10) * linewidth.factor
	baseCex <- (60/xWide) * (xFig/10)
	if ( baseCex > 1.5) baseCex <- 1.5
	mainCex <- 1

	# control the plot boundary explicitly
	savMAI <- par( "mai")
	par( "mai"=c( 0.45,0.96,0.8,0.3))
	if ( mode == "multi") {
		par( "mai"=c( 0.42,0.66,0.44,0.15))
	}

	# get the gene annotation for this region
	geneLocQualifier <- ""
	if ( is.null(gmap)) {
		#gmap <- subset( getCurrentGeneMap(), SEQ_ID == seqID & position >= POSITION & position <= END)
		gmap <- subset( getCurrentGeneMap(), SEQ_ID == seqID)
		if ( ! nrow(gmap)) {
			cat( "\nError:  No SEQ_ID in current species has the name: ", seqID)
			return()
		}
	}
	if ( ! is.null( geneID)) {
		if ( geneID %in% gmap$GENE_ID) {
			gmap <- subset( gmap, GENE_ID == geneID)
		}
		if ( nrow(gmap) > 1) {
			gmap <- subset( gmap, shortGeneName(gmap$GENE_ID,keep=1) == shortGeneName(geneID,keep=1))
		}
	}
	if ( nrow( gmap) < 1) {
		#cat( "\nNo gene found at:  ", seqID, position, geneID, "  Skipping plot...")
		#return()
		# find whatever row is closest
		fullGeneMap <- subset( getCurrentGeneMap(), SEQ_ID == seqID)
		dStarts <- abs( position - fullGeneMap$POSITION)
		dStops <- abs( position - fullGeneMap$END)
		bestStart <- which.min( dStarts)
		bestStop <- which.min( dStops)
		gmap <- fullGeneMap[ bestStart, ]
		if ( dStops[bestStop] < dStarts[bestStart]) {
			gmap <- fullGeneMap[ bestStop, ]
			geneLocQualifier <- "Intergenic after "
		} else {
			gmap <- fullGeneMap[ bestStart, ]
			geneLocQualifier <- "Intergenic before "
		}
	}
	if ( nrow( gmap) > 1) {
		bestOne <- which.min( gmap$N_EXON_BASES)
		gmap <- gmap[ bestOne, ]
	}
	geneID <- gmap$GENE_ID[1]
	geneName <- shortGeneName( geneID, keep=1)
	geneStrand <- gmap$STRAND[1]
	gProd <- gmap$PRODUCT[1]
	cmap <- subset( getCurrentCdsMap(), GENE_ID %in% gmap$GENE_ID)

	# get the bases to plot, and load that portion of the PILEUPS
	loadKnownSNPtable(seqID, verbose=verbose)
	SNP_curMPU <<- BAM.mpileup( bamfile, seqID, fastaFile, start=xLo, stop=xHi, summarize.calls=TRUE,
					max.depth=max.depth, verbose=FALSE)
	if ( is.null( SNP_curMPU)) return(NULL)
	SNP_curVCF <<- if ( file.exists(vcffile)) read.delim( vcffile, as.is=T) else data.frame()
	genomicStr <- getFastaSeqFromFilePath( fastaFile, seqID)
	curGenomeDNA <<- strsplit( as.character(genomicStr), split="")[[1]]
	SNP_curGeneMapSubset <- subset.data.frame( getCurrentGeneMap(), SEQ_ID == seqID)
	SNP_curCdsMapSubset <- subset.data.frame( getCurrentCdsMap(), SEQ_ID == seqID)

	hasStuffToPlot <- ( nrow(SNP_curMPU) > 0)

	# we may find details for a SNP at the center position, but pre-make in case we don't
	snpCenter <- data.frame()

	if ( hasStuffToPlot) {

	# get the flip counts and total bases for each nucleotide with any reads
	snpCenter <- subset.data.frame( SNP_curMPU, POSITION == position)
	xLocs <- SNP_curMPU$POSITION

	# turn the base calls into an explicit matrix in order: genomic,A,C,G,T,N, indel
	flips <- MPU.callStringsToMatrix( SNP_curMPU$BASE_TABLE)

	# the Indels are tougher, we need to get the actual bases for both the reference and the indels and figure if
	# its an insertion or deletion
	indelDetails <- MPU.getIndelDetails( flips, SNP_curMPU)
	isIndel <- isDeletions <- vector()
	if ( indelDetails$nIndels > 0) {
		isIndel <- indelDetails$who
		#snpTopBase[isIndel] <- indelDetails$bases
		# indels that are the start of deletions are actually the genomic base, followed by 1+ '*' deletions
		# we want to visualize these as genomic at the first location, not as an indel, so fudge those counts
		isDeletion <- which( substr( indelDetails$indelText, 1, 2) == ",-")
		# those deletions we need to fix are the intersect of Indel and Deletion
		isDeletion <- intersect( isIndel, isDeletion)
		if ( length( isDeletion)) {
			delRows <- isDeletion
			# extract the first numeric
			delStrings <- indelDetails$indelText[ isDeletion]
			delCountString <- sub( "^,\\-[ACGT]+=", "", delStrings)
			delCountString <- sub( ";.+","",delCountString)
			delCount <- as.integer( delCountString)
			delCount[is.na(delCount)] <- 0
			# now tweak those counts
			gCntsNow <- flips[ delRows, 1] + delCount
			indelCntsNow <- flips[ delRows, 7] - delCount
			indelCntsNow[ indelCntsNow < 0] <- 0
			flips[ delRows, 1] <- gCntsNow
			flips[ delRows, 7] <- indelCntsNow
			# lastly, flag these deletion starts as looking like the genome
			indelDetails$bases[ which( isIndel %in% isDeletion)] <- ","
		}
	}

	# now tabulate the cummulative sums to know what to draw
	csums <- apply( flips, MARGIN=1, cumsum)
	csums <- t( csums)
	rownames( flips) <- rownames( csums) <- xLocs

	# get the majority base at each SNP spot, if not a SNP it will be ',' (matches the reference
	snpTopBase <- apply( flips, MARGIN=1, function(x) colnames(flips)[ which.max(x)])
	names( snpTopBase) <- xLocs
	
	# now stuff in the indel calls
	if ( length(isIndel)) {
		snpTopBase[isIndel] <- indelDetails$bases
	}

	# build the set of colored bars to be drawn
	allX <- allY <- allColor <- vector()
	isSNPbars <- vector()
	for ( i in 1:nrow(csums)) {
		# make a mask to keep the right colors
		mask <- rep( TRUE, times=7)
		# first part is to only draw when the height is higher than previous color
		mask[ 2:7] <- ( csums[ i, 2:7] > csums[ i, 1:6])
		# second part is to not draw the black if 'too small'
		mask[1] <- ( csums[ i, 1] > ( 0.02 * csums[ i, 7]))
		allX <- base::append( allX, rep( xLocs[i], times=7)[mask])
		allY <- base::append( allY, (csums[ i, 1:7])[mask])
		allColor <- base::append( allColor, c("gray50","red","dodgerblue","orange","green","brown","lightpink")[mask])
		# try to keep track of which are SNPs that deserve a re-draw (black less than 80%)
		smallSNPbars <- rep( FALSE, times=7)
		if ( csums[ i, 1] < (csums[ i, 7] * 0.90)) smallSNPbars <- mask
		isSNPbars <- append( isSNPbars, smallSNPbars[mask])
	}
	# sort them into tallest to shortest at every x
	ord <- base::order( allX, -allY)
	allX <- allX[ ord]
	allY <- allY[ ord]
	allColor <- allColor[ ord]
	isSNPbars <- isSNPbars[ ord]

	# now we can set the Y scaling.  Then stretch the lower bound and set where the annotations go
	yLim <- range( allY)
	if ( yLim[2] < 10) yLim[2] <- 10
	yFig <- figInch[2]
	yAnnoScale <- ((xFig+2.5)/xWide) * (8/yFig)
	if ( yAnnoScale < 0.10) yAnnoScale <- 0.10
	if ( yAnnoScale > 0.60) yAnnoScale <- 0.60
	yLim[2] <- yLim[2] * 1.15
	if ( !is.null(forceYmax)) yLim[2] <- forceYmax
	yLim[1] <-  (-(yLim[2] * yAnnoScale))
	extraYscale <- sqrt( 0.25 / yAnnoScale)
	if ( extraYscale > 0.85) extraYscale <- 0.85
	genTextOffset <- yLim[1] * 0.2 * extraYscale
	snpTextOffset <- yLim[1] * 0.45 * extraYscale
	genAAoffset <- yLim[1] * 0.95 * extraYscale
	snpAAoffset <- yLim[1] * 1.15 * extraYscale
	aaBoxOffset <- abs( snpAAoffset - genAAoffset) * 0.6
	if ( !is.null(forceYmax) && mode == "multi") yLim[1] <- 0

	mainText <- paste( "SNP Evidence:          Sample = ", sampleID, 
			paste( "\nChrom=", seqID, "       Loc=", position, "        Strand= '", 
			geneStrand, "'" ), "\nGeneID=", paste(geneLocQualifier,geneName), "       ", gProd) 
	yLab <- "Raw Reads per Base"
	if (mode == "multi" || xWide > 300) {
		mainText <- paste( "Sample:  ", sampleID, "      GeneID:  ", paste(geneLocQualifier,geneName), "\n", gProd)
	}

	if ( xFig < 7) {
		mainCex <- sqrt( 7/xFig)
		mainText <- paste( "Sample = ", sampleID, "      GeneID:  ", paste(geneLocQualifier,geneName), "\n", gProd) 
	}
	if ( xFig < 5.5) {
		mainCex <- 1.5
		mainText <- paste( label, "     ", paste(geneLocQualifier,geneName))
		par( "mai"=c( 0.40,0.36,0.44,0.10))
		yLab <- NA

	}
	if ( xFig < 3) {
		mainCex <- 1.5
		mainText <- label
	}

	# if user expanded the text size for AA and Bases, expand the bar width too
	if ( cex.text > 1) lineWidth <- lineWidth * (cex.text * 0.9)

	plot( allX, allY, xlim=xLim, ylim=yLim, main=mainText, cex.main=mainCex, xlab=NA,
		ylab=yLab, type="h", lwd=lineWidth, lend=2, col=allColor, 
		xaxt="n", yaxt="n", font.lab=2, ...)
	if ( xFig > 5) {
		axis( side=1, font=2)
	} else {
		ats <- axTicks( side=1)
		ats <- ats[ seq( 2, length(ats), by=2)]
		axis( side=1, at=ats, font=2)
	}
	ats <- axTicks( side=2)
	ats <- ats[ ats >= 0]
	if ( yFig < 5) {
		if ( length(ats) > 3) ats <- ats[ seq( 1, length(ats), by=2)]
	}
	axis( side=2, at=ats, font=2)

	# draw the SNP parts again a bit wider, when the plot is big
	if ( xWide > 275 && SNPwidth.factor > 1.0) {
		thisLWD <- lineWidth * SNPwidth.factor
		lines( allX[isSNPbars], allY[isSNPbars], type="h", lwd=thisLWD,
			col=allColor[isSNPbars])
	}

	# draw the arrow
	whoCtr <- base::which( allX == position)
	if ( show.legends == "all" && length( whoCtr) > 0) {
		myYs <- allY[whoCtr]
		best <- which.max( myYs)
		gap <- (yLim[2]-yLim[1]) * 0.05
		# get the color for the biggest component
		snpCenter <- flips[ rownames(flips) == position, , drop=FALSE]
		cnts <- as.vector( snpCenter[ 1, ])
		myColor <- c( "gray50", "red", "dodgerblue", "orange", "green", "brown", "lightpink")[ which.max( cnts)]
		arrows( position, (yLim[2]*1.1), position, (myYs[best]+gap), col=myColor, lwd=2, lty=1, 
			length=0.12) 
	}

	} else {  # no bases to plot

		snpTopBase <- ","
		yLim <- c( 0, 10)
		yLim[2] <- yLim[2] * 1.15
		if ( !is.null(forceYmax)) yLim[2] <- forceYmax
		yLim[1] <-  (-(yLim[2]/8))
		genTextOffset <- yLim[1] * 0.2
		snpTextOffset <- yLim[1] * 0.45
		genAAoffset <- yLim[1] * 0.95
		snpAAoffset <- yLim[1] * 1.15
		aaBoxOffset <- abs( snpAAoffset - genAAoffset) * 0.6
		if ( !is.null(forceYmax) && mode == "multi") yLim[1] <- 0

		mainText <- paste( "SNP Evidence:        Sample = ", sampleID, 
				"\nGeneID=", geneName, "       ", gProd)
		yLab <- "Raw Reads per Base"
		if (mode == "multi" || xWide > 300) {
			mainText <- paste( "Sample:  ", sampleID, "      GeneID:  ", geneName, "\n", gProd)
		}
		if ( xFig < 7) {
			mainCex <- sqrt( 7/xFig)
			mainText <- paste( "Sample = ", sampleID, "      GeneID:  ", geneName, "\n", gProd) 
		}
		if ( xFig < 5.5) {
			mainCex <- 1.5
			mainText <- paste( label, "     ", geneName)
			par( "mai"=c( 0.40,0.36,0.44,0.10))
			yLab <- NA
		}
		if ( xFig < 3) {
			mainCex <- 1.5
			mainText <- label
		}

		plot( 1, 1, xlim=xLim, ylim=yLim, main=mainText, cex.main=mainCex, xlab=NA, 
			ylab=yLab, type="n", font.axis=2, font.lab=2, ...) 
	}

	# draw the genome Base letters
	allBases <- xLim[1] : xLim[2]
	genomeSNPtext <- genomeBaseText <- curGenomeDNA[ allBases]
	text( allBases, rep(genTextOffset,times=length(allBases)), genomeBaseText, cex=baseCex*cex.text)

	# see what bases are different from genomic, and draw the different AA
	whoSNP <- base::which( snpTopBase != ",")
	if ( length( whoSNP) > 0) {
		myLocs <- as.integer( names(snpTopBase)[whoSNP])
		where <- base::match( myLocs, allBases)
		genomeSNPtext[ where] <- snpTopBase[ whoSNP]
	}

	# draw the protein amino acid letters. with and without SNP effect
	ans <- convertGenomicBasesToCodingAA( seqID, position=xLim[1], end=xLim[2], 
			strand=geneStrand, dnaQuery=genomeSNPtext, genomeDNA=curGenomeDNA,
			geneMap=SNP_curGeneMapSubset, cdsMap=SNP_curCdsMapSubset, geneID=geneID)
	

	genomeAminoText <- ans$genomic
	text( allBases, rep(genAAoffset,times=length(allBases)), genomeAminoText, cex=baseCex*cex.text)
	if ( nrow(cmap) > 0) {
		for ( j in 1:nrow(cmap)) rect( cmap$POSITION[j] - 1, genAAoffset-aaBoxOffset, cmap$END[j] + 1, 
						genAAoffset+aaBoxOffset, border=1, lwd=1)
	}
	if ( length( whoSNP) > 0) {
		aminoSNPtext <- ans$query
		whoDiff <- base::which( genomeAminoText != aminoSNPtext)
		if ( length( whoDiff) > 0) {
			# try to color the SNP'ed AA
			snpColors <- c("red","dodgerblue","orange","green","brown")[ base::match( genomeSNPtext, 
					c("A","C","G","T","N"))] 
			points( allBases[whoDiff], rep(snpAAoffset,times=length(whoDiff)), pch=22, 
					bg="white",  col=snpColors[whoDiff], cex=baseCex*3.4*cex.text)
			text( allBases[whoDiff], rep(snpAAoffset,times=length(whoDiff)), 
					aminoSNPtext[whoDiff], font=2, cex=baseCex*cex.text)
		}
	}

	# top legends...
	legendCex <- cex.legend * 0.9
	if ( xFig < 9) legendCex <- legendCex * sqrt( xFig/9)

	if( show.legends == "all" && nrow(snpCenter) > 0 && xFig > 4) {
		legendText <- paste( c("Genome", "Flip_A  ", "Flip_C  ", "Flip_G  ", "Flip_T   ", "Flip_N  ", "Indels   "), 
					" =", snpCenter[ 1, ])
		legend( "topleft", legendText, col=c("gray50","red","dodgerblue","orange","green", "brown", "lightpink"),
			lwd=6, cex=legendCex, bg="white")

		# grab the Phred scaled likelihood from the VCF
		if ( nrow( SNP_curVCF) > 0) {
			vcfCenter <- subset( SNP_curVCF, POSITION == position)
		} else {
			vcfCenter <- data.frame()
		}
		if ( nrow( vcfCenter) > 0) {
			plTerm <- vcfCenter[ 1, 10]
			score <- VCF.Score( plTerm)
		} else {
			score <- 1
		}
		legendText <- paste( "SNP Score:  ", score)
		legend( "topright", legendText, cex=legendCex, bg="white")
	}
	if( show.legends == "gene" && xFig > 4) {
		legendText <- c("Genome", "Flip_A  ", "Flip_C  ", "Flip_G  ", "Flip_T   ", "Flip_N  ", "Indels   ")
		legend( "topleft", legendText, col=c("gray50","red","dodgerblue","orange","green", "brown", "lightpink"),
			lwd=6, cex=legendCex, bg="white")
	}

	# draw in known SNPs if we have any
	didLowerLegend <- FALSE
	if ( ! is.null( SNP_curSNPtable) && (nrow(SNP_curSNPtable) > 0)) {
	    #snps <- subset.data.frame( SNP_curSNPtable, POSITION >= xLim[1] & POSITION <= xLim[2])
	    snps <- which( SNP_curSNPtable$POSITION >= xLim[1] & SNP_curSNPtable$POSITION <= xLim[2])
	    nSNP <- length( snps) 
	    #if ( nrow( snps) > 0) {
	    if ( nSNP > 0) {
	        snpPOSITION <- SNP_curSNPtable$POSITION[snps]
	        snpMAJOR_ALLELE <- SNP_curSNPtable$MAJOR_ALLELE[snps]
	        snpMINOR_ALLELE <- SNP_curSNPtable$MINOR_ALLELE[snps]
		lines( xLim, rep( genTextOffset*1.4, 2))
		who <- base::match( snpMAJOR_ALLELE, c("A","C","G","T","N"))
		snpColors <- c("red","dodgerblue","orange","green","brown")[who]
		points( snpPOSITION, rep( snpTextOffset, times=nSNP), pch=22, bg=snpColors, 
				col=snpColors, cex=baseCex*3.1*cex.text)
		text( snpPOSITION, rep( snpTextOffset, times=nSNP), labels=snpMAJOR_ALLELE, 
				col=1, font=2, cex=baseCex*cex.text)

		# minor alleles
		minorTerms <- strsplit( snpMINOR_ALLELE, split=", *")
		minors <- base::sapply( minorTerms, function(x) return( names( 
				base::sort( base::table(x), decreasing=TRUE))[1]))
		who <- base::match( minors, c("A","C","G","T","N"))
		snpColors <- c("red","dodgerblue","orange","green","brown")[who]
		points( snpPOSITION, rep( snpTextOffset*1.5, times=nSNP), pch=22, bg=snpColors, 
				col=snpColors, cex=baseCex*3.1*cex.text)
		text( snpPOSITION, rep( snpTextOffset*1.5, times=nSNP), labels=minors, col=1, 
				font=2, cex=baseCex*cex.text)
		if ( show.legends %in% c("all","gene") && xFig > 4) {
			legend( "bottomleft", c( "Major Allele", "Minor Allele", "Amino Acid", "SNP AA"), 
				cex=legendCex*1.0, bg="white")
			didLowerLegend <- TRUE
		}
	    }
	}
	if ( ! didLowerLegend && xFig > 4 && show.legends %in% c("all","gene") ) {
		legend( "bottomleft", c( "Amino Acid", "SNP AA"), cex=legendCex*1.0, bg="white")
	}
	dev.flush()

	par( "mai"=savMAI)

	return()
}



multiSample.plotSNP <- function( position, seqID, sampleSet, bamfileSet, vcffileSet, groupSet, geneID=NULL, gmap=NULL,
			tailWidth=30, cex.text=1.0, cex.legend=1.0, forceYmax=NULL, labelSet="",
			show.legends=c("all", "gene", "none"), mf=NULL, verbose=TRUE, ...) {

	# given a set of samples, draw that one location in all samples
	nSamples <- length( sampleSet)
	if ( nSamples < 2) stop( "Need at least 2 smaples for 'multiSample.plotSNP()'")
	seqIDset <- rep( seqID, length.out=length(sampleSet))
	labelSet <- rep( labelSet, length.out=length(sampleSet))

	if ( is.null(mf)) {
		mf <- c(1,1)
		if ( nSamples > 1) mf <- c(1,2)
		if ( nSamples > 2) mf <- c(2,2)
		if ( nSamples > 4) mf <- c(2,3)
		if ( nSamples > 6) mf <- c(3,3)
		if ( nSamples > 9) mf <- c(3,4)
		if ( nSamples > 12) mf <- c(4,4)
		if ( nSamples > 16) mf <- c(4,5)
		if ( nSamples > 20) mf <- c(5,5)
		if ( nSamples > 25) mf <- c(5,6)
		if ( nSamples > 30) mf <- c(6,6)
		if ( nSamples > 36) mf <- c(6,7)
		if ( nSamples > 42) mf <- c(6,8)
		if ( nSamples > 48) mf <- c(6,9)
		if ( nSamples > 54) mf <- c(7,9)
		if ( nSamples > 63) mf <- c(7,10)
		if ( nSamples > 70) mf <- c(7,11)
		if ( nSamples > 77) mf <- c(8,11)
		if ( nSamples > 88) mf <- c(8,12)
	}
	par( mfrow=mf)

	show.legends <- match.arg( show.legends)

	for ( k in 1:nSamples) {
		sample <- sampleSet[k]
		bamfile <- bamfileSet[k]
		vcffile <- vcffileSet[k]
		# when few enough plots, the label will be made from sample/gene/etc.
		# when too many, the label may get over simplified...
		label <- groupSet[k]
		if (nSamples > 2 && all(labelSet == "")) {
			# combine the sampleID and the group into the label if its not already there
			if ( regexpr( groupSet[k], sample) < 1) {
				label <- paste( sample, groupSet[k], sep=":  ")
			} else {
				label <- sample
			}
		} else {
			label <- labelSet[k]
		}

		plotSNP( position, seqID, sample, bamfile, vcffile, tailWidth=tailWidth, geneID=geneID, 
				gmap=gmap, cex.text=cex.text, cex.legend=cex.legend, mode="multi", label=label, 
				forceYmax=forceYmax, show.legends=show.legends, verbose=verbose, ...)
	}
	par( mfrow=c(1,1))
}


multiPosition.plotSNP <- function( positionSet, seqIDset, sample, bamfile, vcffile, labelSet, geneID=NULL, gmap=NULL,
			tailWidth=30, cex.text=1.0, cex.legend=1.0, forceYmax=NULL, 
			show.legends=c("all", "gene", "none"), mf=NULL, verbose=TRUE, ...) {

	# given one sample, draw multiple locations
	nPlots <- length( positionSet)
	if ( nPlots < 2) stop( "Need at least 2 positions for 'multiPosition.plotSNP()'")
	seqIDset <- rep( seqIDset, length.out=length(positionSet))
	labelSet <- rep( labelSet, length.out=length(positionSet))

	if ( is.null(mf)) {
		mf <- c(1,1)
		if ( nPlots > 1) mf <- c(1,2)
		if ( nPlots > 2) mf <- c(2,2)
		if ( nPlots > 4) mf <- c(2,3)
		if ( nPlots > 6) mf <- c(3,3)
		if ( nPlots > 9) mf <- c(3,4)
		if ( nPlots > 12) mf <- c(4,4)
		if ( nPlots > 16) mf <- c(4,5)
		if ( nPlots > 20) mf <- c(5,5)
		if ( nPlots > 25) mf <- c(5,6)
		if ( nPlots > 30) mf <- c(6,6)
		if ( nPlots > 36) mf <- c(6,7)
		if ( nPlots > 42) mf <- c(6,8)
		if ( nPlots > 48) mf <- c(6,9)
		if ( nPlots > 54) mf <- c(7,9)
		if ( nPlots > 63) mf <- c(7,10)
		if ( nPlots > 70) mf <- c(7,11)
		if ( nPlots > 77) mf <- c(8,11)
		if ( nPlots > 88) mf <- c(8,12)
	}
	par( mfrow=mf)

	show.legends <- match.arg( show.legends)

	for ( k in 1:nPlots) {
		position <- positionSet[k]
		seqID <- seqIDset[k]
		label <- labelSet[k]

		plotSNP( position, seqID, sample, bamfile, vcffile, tailWidth=tailWidth, geneID=geneID, 
				gmap=gmap, cex.text=cex.text, cex.legend=cex.legend, mode="multi", label=label, 
				forceYmax=forceYmax, show.legends=show.legends, verbose=verbose, ...)
	}
	par( mfrow=c(1,1))
}

