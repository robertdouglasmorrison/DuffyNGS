# plotWIG.R

# routines to plot a gene or region of one or two WIG compressed objects

`plotWIGgene` <- function( WIG, gene, tailWidth=1000, plotType=c("boxes", "lines", "segments", "addedLines", "none"), 
			altGeneMap=NULL, label="", useLog=FALSE, col=c(4,2,1), lwd=2, addStrands=FALSE, 
			forceYmax=NULL, minYmax=10, showDetectability=TRUE, legend.cex=1, triggerWarning=NULL) {

	plotType <- match.arg( plotType)

	# check species
	if ( is.null( altGeneMap)) {
		if ( ! verifyWIGspecies( WIG, seqID=NULL, geneID=gene)) return()
	}

	curSpecies <- getCurrentSpecies()
	seqMap <- getCurrentSeqMap()
	geneMap <- getCurrentGeneMap()
	exonMap <- getCurrentExonMap()
	cdsMap <- getCurrentCdsMap()

	# this may be used for detectability....
	WB_setBinSizeBySpecies( curSpecies)

	# allow an alternate gene map
	if ( ! is.null( altGeneMap)) {
		geneMap <- altGeneMap
		exonMap <- altGeneMap
		cdsMap <- altGeneMap
	}

	# use the gene if given...allow partial matches
	gmap <- subset.data.frame( geneMap, GENE_ID==gene)
	if ( nrow(gmap) < 1) {
		gmap <- subset.data.frame( geneMap, GENE_ID %in% grep( paste( "^",gene,sep=""), 
				geneMap$GENE_ID, value=TRUE))
		if ( nrow( gmap) < 1) {
			cat( "\nNo gene of that name: ", gene)
			return()
		}
		gene <- gmap$GENE_ID[1]
	}
	# force exactly one row...
	gmap <- as.data.frame( gmap[ 1:1, ])
	seqid <- gmap$SEQ_ID
	if ( !(seqid %in% seqMap$SEQ_ID)) {
		cat( "\nGeneID and/or SeqID not in current Species.\n")
		return()
	}

	# turn the gene location and tails into the base limits
	minBase <- 1
	maxBase <- subset.data.frame( seqMap, SEQ_ID==seqid)$LENGTH
	leftBase <- gmap$POSITION - tailWidth
	rightBase <- gmap$END + tailWidth
	if ( leftBase < minBase) leftBase <- minBase
	if ( rightBase > maxBase) rightBase <- maxBase

	# scale based on the one gene of interest
	leftBaseScaling <- gmap$POSITION
	rightBaseScaling <- gmap$END

	# get the actual data
	wigdata <- WIG_getWigglesOneSeq( WIG, seqid)
	hasData <- ( nrow(wigdata$Plus) > 0) || ( nrow(wigdata$Minus) > 0)

	if ( hasData) {
		yData <- list( "Plus"=baseDepthTableSubset( wigdata$Plus, leftBase, rightBase),
			"Minus"=baseDepthTableSubset( wigdata$Minus, leftBase, rightBase),
			"PlusUnique"=baseDepthTableSubset( wigdata$PlusUnique, leftBase, rightBase),
			"MinusUnique"=baseDepthTableSubset( wigdata$MinusUnique, leftBase, rightBase))

		if ( addStrands) {
			yData$Plus <- merge.baseDepthTables( yData$Plus, yData$Minus)
			yData$Minus <- emptyBaseDepthTable()
			yData$PlusUnique <- merge.baseDepthTables( yData$PlusUnique, yData$MinusUnique)
			yData$MinusUnique <- emptyBaseDepthTable()
		}

		# pick sensible plot limits, given the data and user given limits
		bigY <- max( c( minYmax, 
				baseDepthTableSubset( yData$Plus, leftBaseScaling, rightBaseScaling)$DEPTH, 
				baseDepthTableSubset( yData$Minus, leftBaseScaling, rightBaseScaling)$DEPTH)
			)

		# certain plots may need more room -- a summed plot may be up to 2x tall
		if ( regexpr( "added", plotType) > 0) bigY <- bigY * 2

	} else {  # no data to draw...
		bigY <- minYmax
	}
	if ( !is.null( forceYmax)) bigY <- forceYmax
	xlim <- c(leftBase, rightBase)

	if ( useLog) {
		useBigY <- log2( bigY)
		smallestY <- 0.25
	} else {
		useBigY <- bigY
		smallestY <- useBigY / 800
	}

	ylim <- c( 0, useBigY*1.1)

	# now stretch the lower end of the plot to fit the annotations
	ylim[1] <- ylim[1] - (( ylim[2]-ylim[1]) * 0.20)
	annoSize <- (ylim[2]-ylim[1]) * 0.04

	dataType <- WIG$Info$DataType
	if ( is.null( dataType)) dataType <- "RNA-seq"
	ylabel <- "Raw Read Depth"
	pltTxt <- paste( dataType, label, sep=":  ") 
	if ( ! is.null( WIG$Info$DataType) && WIG$Info$DataType == "Peptides") {
		pltTxt <- paste( "Proteomics:  ",label) 
		ylabel <- "Raw Peptide Count"
	}
	mainText <- paste( pltTxt, "      ", shortGeneName(gene), "\n", gmap$PRODUCT)
	if ( gene2OrigID(gene) != gene) {
		mainText <- paste( pltTxt, "      ", shortGeneName(gene), "  (", gene2OrigID(gene), ")\n", 
				gmap$PRODUCT, sep="")
	}
	if ( useLog) ylabel <- paste( ylabel, "   (log scale)")

	plot( x=xlim, y=c(0,0), main=mainText, type="l", xlim=xlim, ylim=ylim, xaxs="r",
		xlab=paste( "Chromosome: ", gmap$SEQ_ID[1], "     Nucleotide Location"), 
		ylab=ylabel, yaxt="n", font.axis=2, font.lab=2)


	if ( useLog) {
		ticks <- rep( c( 1,2,5), times=7) * rep( c( 1,10,100,1000,10000,100000,1000000), each=3)
		ticks <- c( 3,10,30,100,300,1000,3000,10000,30000,100000,3000000,1000000)
		useTicks <- ticks[ ticks < bigY]
		if ( length( useTicks) < 4) {
			useTicks <- ticks[1:4]
		} else {
			useTicks <- rev( useTicks)[ seq( 1, length(useTicks), by=2)]
		}
		at <- log2( useTicks)
		if ( hasData) {
			yData <- log2Scale( yData)
			WIG_Draw( yData, plotType=plotType, col=col, lwd=lwd, smallestDepth=smallestY,
					show=if (dataType == "ChIP-seq") "multiOnly" else "both")
		}
		WIG_Ticks( xlim, at, labels=as.integer(useTicks))

	} else {
		ticks <- seq( from=0, to=1000000, length.out=5)[2:5]
		yBreaks <- c(  800000, 600000, 400000, 200000, 160000, 120000, 80000, 60000, 40000, 20000, 
				16000, 12000, 8000, 6000, 4000, 2000, 1600, 1200, 800, 600, 400, 200, 
				160, 120, 80, 60, 40, 20, 16, 12, 8, 4)
		for ( ycut in yBreaks) {
			if ( bigY < ycut) ticks <- seq( from=0, to=ycut, length.out=5)[2:5]
		}
		# draw the wanted style...
		if ( hasData) {
			WIG_Draw( yData, plotType=plotType, col=col, lwd=lwd, smallestDepth=smallestY,
					show=if (dataType == "ChIP-seq") "multiOnly" else "both")
		}
		WIG_Ticks( xlim, ticks, labels=ticks)
	}

	annotateWIGplot( seqID=seqid, gene=gene, xlim=xlim, annoSize, speciesID=curSpecies, 
			altGeneMap=altGeneMap, showDetectability=showDetectability)

	if ( exists( "infoContent")) plotInfoContent( seqid, lwd=2, legendLocation="topright")

	if ( legend.cex > 0) {
		if ( showDetectability) {
			legend( x="topleft", legend=c("(+) Strand Reads", " (-) Strand Reads", "DNA Uniqueness"), 
				col=c(col[1],col[2],3), text.col=c(col[1],col[2],3), lwd=c(5,5,1), 
				bg="white", cex=legend.cex)
		} else {
			legend( x="topleft", legend=c("(+) Strand Reads", " (-) Strand Reads"), 
				col=c(col[1],col[2]), text.col=c(col[1],col[2]), lwd=c(5,5), 
				bg="white", cex=legend.cex)
		}
	}

	if ( ! is.null( triggerWarning)) {
		xmid <- mean( xlim)
		ymid <- mean( ylim)
		text( xmid, ymid, triggerWarning, font=2, cex=1.2)
	}

	return()
}


`plotWIGregion` <- function( WIG, seqid, position=1, end=NULL, 
			plotType=c("lines", "boxes", "segments", "addedLines", "none"), label="", 
			minYmax=40, forceYmax=NULL, showDetectability=TRUE, useLog=FALSE,
			legend.cex=1, lwd=2, col=c(4,2,1)) {

	plotType <- match.arg( plotType)

	# check species
	if ( ! verifyWIGspecies( WIG, seqID=seqid, geneID=NULL)) return()

	curSpecies <- getCurrentSpecies()
	seqMap <- getCurrentSeqMap()
	geneMap <- getCurrentGeneMap()
	exonMap <- getCurrentExonMap()
	cdsMap <- getCurrentCdsMap()
	WB_setBinSizeBySpecies( curSpecies)

	smap <- subset.data.frame( seqMap, SEQ_ID==seqid)
	if ( nrow(smap) < 1) {
		cat( "\nNo SeqID with name: ", seqid, " in the current Species.\n")
		return()
	}

	# force exactly one row...
	smap <- as.data.frame( smap[ 1:1, ])

	# turn the gene location and tails into the base and bin limits
	minBase <- 1
	maxBase <- smap$LENGTH
	leftBase <- position
	rightBase <- end
	if ( is.null( leftBase)) leftBase <- minBase
	if ( leftBase < minBase) leftBase <- minBase
	if ( is.null( rightBase)) rightBase <- maxBase
	if ( rightBase > maxBase) rightBase <- maxBase

	# get the actual data
	wigdata <- WIG_getWigglesOneSeq( WIG, seqid)
	hasData <- ( nrow(wigdata$Plus) > 0) || ( nrow(wigdata$Minus) > 0)

	if ( hasData) {
		yData <- list( "Plus"=baseDepthTableSubset( wigdata$Plus, leftBase, rightBase),
			"Minus"=baseDepthTableSubset( wigdata$Minus, leftBase, rightBase),
			"PlusUnique"=baseDepthTableSubset( wigdata$PlusUnique, leftBase, rightBase),
			"MinusUnique"=baseDepthTableSubset( wigdata$MinusUnique, leftBase, rightBase))

		# pick sensible plot limits, given the data and user given limits
		bigY <- max( c( minYmax, yData$Plus$DEPTH, yData$Minus$DEPTH))
		# certain plots may need more room -- a summed plot may be up to 2x tall
		if ( regexpr( "added", plotType) > 0) bigY <- bigY * 2

	} else {  # no data to draw...
		bigY <- minYmax
	}
	if ( !is.null( forceYmax)) bigY <- forceYmax
	xlim <- c(leftBase, rightBase)

	if ( useLog) {
		useBigY <- log2( bigY)
		smallestY <- 0.25
	} else {
		useBigY <- bigY
		smallestY <- useBigY / 800
	}

	ylim <- c( 0, useBigY*1.1)

	# now stretch the lower end of the plot to fit the annotations
	ylim[1] <- ylim[1] - (( ylim[2]-ylim[1]) * 0.20)
	annoSize <- (ylim[2]-ylim[1]) * 0.04

	dataType <- WIG$Info$DataType
	if ( is.null( dataType)) dataType <- "RNA-seq"
	ylabel <- "Raw Read Depth"
	if ( useLog) ylabel <- paste( ylabel, "   (log scale)")

	pltTxt <- paste( "RNA_Seq:  ",label) 
	if ( ! is.null( dataType)) pltTxt <- paste( WIG$Info$DataType, ":   ",label, sep="") 
	if ( ! is.null( dataType) && WIG$Info$DataType == "Peptides") {
		pltTxt <- paste( "Proteomics:  ",label) 
		ylabel <- "Raw Peptide Count"
	}
	pltTxt <- paste( pltTxt, "\n", seqid, ":     ", formatC( leftBase, format="d", big.mark=","),
			"  to  ", formatC( rightBase, format="d", big.mark=","))

	plot( x=xlim, y=c(0,0), main=pltTxt, type="l", xlim=xlim, ylim=ylim, xaxs="r", xaxt="n",
		xlab=paste( "Chromosome: ", seqid, "     Nucleotide Location"), ylab=ylabel, yaxt="n",
		font.axis=2, font.lab=2)
	
	xAt <- pretty( xlim)
	axis( side=1, at=xAt, labels=formatC( as.integer( xAt), format="d", big.mark=","), font=2)

	if ( useLog) {
		ticks <- rep( c( 1,2,5), times=7) * rep( c( 1,10,100,1000,10000,100000,1000000), each=3)
		ticks <- c( 3,10,30,100,300,1000,3000,10000,30000,100000,3000000,1000000)
		useTicks <- ticks[ ticks < bigY]
		if ( length( useTicks) < 4) {
			useTicks <- ticks[1:4]
		} else {
			useTicks <- rev( useTicks)[ seq( 1, length(useTicks), by=2)]
		}
		at <- log2( useTicks)
		if ( hasData) {
			yData <- log2Scale( yData)
			WIG_Draw( yData, plotType=plotType, col=col, lwd=lwd, smallestDepth=smallestY,
					show=if (dataType == "ChIP-seq") "multiOnly" else "both")
		}
		WIG_Ticks( xlim, at, labels=as.integer(useTicks))

	} else {
		ticks <- seq( from=0, to=1000000, length.out=5)[2:5]
		yBreaks <- c(  800000, 600000, 400000, 200000, 160000, 120000, 80000, 60000, 40000, 20000, 
				16000, 12000, 8000, 6000, 4000, 2000, 1600, 1200, 800, 600, 400, 200, 
				160, 120, 80, 60, 40, 20, 16, 12, 8, 4)
		for ( ycut in yBreaks) {
			if ( bigY < ycut) ticks <- seq( from=0, to=ycut, length.out=5)[2:5]
		}
		# draw the wanted style...
		if ( hasData) {
			WIG_Draw( yData, plotType=plotType, col=col, lwd=lwd, smallestDepth=smallestY,
					show=if (dataType == "ChIP-seq") "multiOnly" else "both")
		}
		WIG_Ticks( xlim, ticks, labels=ticks)
	}

	annotateWIGplot( seqID=seqid, gene=NULL, xlim=xlim, annoSize, speciesID=curSpecies, 
			showDetectability=showDetectability)

	if ( exists( "infoContent")) plotInfoContent( seqid, lwd=2, legendLocation="topright")

	if ( legend.cex > 0) {
	if ( showDetectability) {
		legend( x="topleft", legend=c("(+) Strand Reads", " (-) Strand Reads", "DNA Uniqueness"), 
			col=c(col[1],col[2],3), text.col=c(col[1],col[2],3), lwd=c(5,5,1), 
			bg="white", cex=legend.cex)
	} else {
		legend( x="topleft", legend=c("(+) Strand Reads", " (-) Strand Reads"), 
			col=c(col[1],col[2]), text.col=c(col[1],col[2]), lwd=c(5,5), 
			bg="white", cex=legend.cex)
	}
	}

	return()
}



`plotMultiWIGgene` <- function( WIGlist, colors=1:length(WIGlist),  gene, addStrands=TRUE, tailWidth=1000, lwd=2,
			plotType=c( "boxes", "lines", "segments", "none"), label="", forceYmax=NULL, minYmax=10, 
			altGeneMap=NULL, showDetectability=TRUE, useLog=FALSE, legend.cex=1.0,
			PLOT.FUN=NULL) {

	if ( ! is.null(PLOT.FUN) && ! is.function(PLOT.FUN)) return()

	plotType <- match.arg( plotType)

	# check species and extract total reads for every dataset
	nWIG <- length( WIGlist)
	if ( nWIG < 2) {
		cat( "Need at least 2 WIG datasets for 'plotMultiWIGgene()'")
		return()
	}
	mySpecies <- WIGlist[[1]]$Info$Species
	dataType <- WIGlist[[1]]$Info$DataType
	if ( is.null( dataType)) dataType <- "RNA-seq"

	totalReadsPerWIG <- vector( length=nWIG)
	medianDepthPerWIG <- vector( length=nWIG)
	MIN_TOTAL_READS_SCALING <- 100000
	MIN_MEDIAN_READS_SCALING <- 10
	if ( dataType == "DNA-seq") {
	    for( i in 1:nWIG) {
		thisInfo <- WIGlist[[i]]$Info
		if ( mySpecies != thisInfo$Species) {
			cat( "\nGiven species: ", mySpecies)
			cat( "\nFound species: ", thisInfo$Species)
			stop( "All WIG datasets are not the same speicesID")
		}
		thisDepthM <- thisInfo$MedianDepth
		# don't let very low counts turn into very high scaling...
		medianDepthPerWIG[i] <- median( pmax( apply( thisDepthM, MARGIN=1, sum, na.rm=T), MIN_MEDIAN_READS_SCALING))
	    }
	} else {
	    for( i in 1:nWIG) {
		thisInfo <- WIGlist[[i]]$Info
		if ( mySpecies != thisInfo$Species) {
			cat( "\nGiven species: ", mySpecies)
			cat( "\nFound species: ", thisInfo$Species)
			stop( "All WIG datasets are not the same speicesID")
		}
		thisCountsM <- thisInfo$RawReads
		# don't let very low counts turn into very high scaling...
		totalReadsPerWIG[i] <- max( sum( thisCountsM), MIN_TOTAL_READS_SCALING)
	    }
	}

	# set up for this speciesID
	curSpecies <- getCurrentSpecies()
	if ( mySpecies != curSpecies) {
		setCurrentSpecies( mySpecies)
		curSpecies <- getCurrentSpecies()
	}
	seqMap <- getCurrentSeqMap()
	geneMap <- getCurrentGeneMap()
	exonMap <- getCurrentExonMap()
	cdsMap <- getCurrentCdsMap()
	WB_setBinSizeBySpecies( curSpecies)

	# allow alternate gene map
	if ( ! is.null( altGeneMap)) {
		geneMap <- altGeneMap
	}

	# use the gene if given...allow partial matches
	gmap <- subset.data.frame( geneMap, GENE_ID==gene)
	if ( nrow(gmap) < 1) {
		gmap <- subset.data.frame( geneMap, GENE_ID %in% grep( paste( "^",gene,sep=""), 
					geneMap$GENE_ID, value=TRUE))
		if ( nrow( gmap) < 1) {
			cat( "\nNo gene of that name: ", gene)
			return()
		}
		gene <- gmap$GENE_ID[1]
	}
	# force exactly one row...
	gmap <- as.data.frame( gmap[ 1:1, ])
	seqid <- gmap$SEQ_ID
	if ( !(seqid %in% seqMap$SEQ_ID)) {
		cat( "\nGeneID and/or SeqID not in current Species.\n")
		return()
	}

	# turn the gene location and tails into the base and bin limits
	minBase <- 1
	maxBase <- subset.data.frame( seqMap, SEQ_ID==seqid)$LENGTH
	leftBase <- gmap$POSITION - tailWidth
	rightBase <- gmap$END + tailWidth
	if ( leftBase < minBase) leftBase <- minBase
	if ( rightBase > maxBase) rightBase <- maxBase

	# scale based on the one gene of interest
	if ( dataType == "RNA-seq") {
		leftBaseScaling <- gmap$POSITION
		rightBaseScaling <- gmap$END
	} else {
		leftBaseScaling <- leftBase
		rightBaseScaling <- rightBase
	}
	
	# get all the actual data we need, and normalize for reads per file...
	if ( dataType == "DNA-seq") {
		avgDepthPerWIG <- mean( medianDepthPerWIG)
	} else {
		avgReadsPerWIG <- mean( totalReadsPerWIG)
	}
	yDataList <- vector( mode="list", length=nWIG)
	yMax <- 0
	netReadsPerSample <- vector( length=nWIG)
	for ( i in 1:nWIG) {
		oneWIG <- WIGlist[[i]]	
		wigdata <- WIG_getWigglesOneSeq( oneWIG, seqid)
		if ( dataType == "DNA-seq") {
			myScale <- avgDepthPerWIG / medianDepthPerWIG[i]
		} else {
			myScale <- avgReadsPerWIG / totalReadsPerWIG[i]
		}

		# get the subset in view
		tmpPos <- baseDepthTableSubset( wigdata$Plus, leftBase, rightBase)
		tmpNeg <- baseDepthTableSubset( wigdata$Minus, leftBase, rightBase)
		tmpPos$DEPTH <- tmpPos$DEPTH * myScale
		tmpNeg$DEPTH <- tmpNeg$DEPTH * myScale

		if ( addStrands) {
			tmpPos <- merge.baseDepthTables( tmpPos, tmpNeg)
			tmpNeg <- emptyBaseDepthTable()
		}
		# some decisions are based on just the gene of interest without the flanking regions
		genePos <- baseDepthTableSubset( tmpPos, leftBaseScaling, rightBaseScaling)
		geneNeg <- baseDepthTableSubset( tmpNeg, leftBaseScaling, rightBaseScaling)
		yMax <- max( yMax, genePos$DEPTH, geneNeg$DEPTH)
		netReadsPerSample[i] <- sum( genePos$DEPTH) + sum( geneNeg$DEPTH)
		#clip so they stay in view...
		if ( !is.null( forceYmax)) {
			tmpPos$DEPTH[ tmpPos$DEPTH > forceYmax] <- forceYmax
			tmpNeg$DEPTH[ tmpNeg$DEPTH > forceYmax] <- forceYmax
		}
		# OK, its ready for later drawing...
		yDataList[[i]] <- list( "Plus"=tmpPos, "Minus"=tmpNeg)
	}

	# pick sensible plot limits, given the data and user given limits
	bigY <- max( c( minYmax, yMax))
	if ( !is.null( forceYmax)) bigY <- forceYmax
	xlim <- c(leftBase, rightBase)
	if ( useLog) {
		useBigY <- log2( bigY)
		smallestY <- 0.25
	} else {
		useBigY <- bigY
		smallestY <- useBigY / 800
	}
	ylim <- c( 0, useBigY*1.1)

	# now stretch the lower end of the plot to fit the annotations
	annoBase <- ylim[1] <- ylim[1] - (( ylim[2]-ylim[1]) * 0.20)
	annoSize <- (ylim[2]-ylim[1]) * 0.04

	ylabel <- "Normalized Read Depth"
	pltTxt <- paste( "RNA_Seq:  ",label) 
	if ( ! is.null( WIGlist[[1]]$Info$DataType) && WIGlist[[1]]$Info$DataType == "DNA-seq") {
		pltTxt <- paste( "DNA_seq:  ",label) 
		ylabel <- "Median Normalized Read Count"
	}
	if ( ! is.null( WIGlist[[1]]$Info$DataType) && WIGlist[[1]]$Info$DataType == "Peptides") {
		pltTxt <- paste( "Proteomics:  ",label) 
		ylabel <- "Normalized Peptide Count"
	}
	if ( dataType == "ChIP-seq") {
		pltTxt <- paste( "ChIP-seq:  ",label) 
	}
	if ( dataType == "RIP-seq") {
		pltTxt <- paste( "RIP-seq:  ",label) 
	}
	mainText <- paste( pltTxt, "\n", shortGeneName(gene), " -  ", gmap$PRODUCT)
	if ( gene2OrigID(gene) != gene) {
		mainText <- paste( pltTxt, "\n", shortGeneName(gene), "  (", gene2OrigID(gene), ")  - ", 
				gmap$PRODUCT, sep="")
	}
	if ( useLog) ylabel <- paste( ylabel, "   (log scale)")

	plot( x=xlim, y=c(0,0), main=mainText, type="l", xlim=xlim, ylim=ylim, xaxs="r",
		xlab=paste( "Chromosome: ", gmap$SEQ_ID[1], "     Nucleotide Location"), ylab=ylabel, yaxt="n",
		font.axis=2, font.lab=2)

	if ( useLog) {
		ticks <- rep( c( 1,2,5), times=7) * rep( c( 1,10,100,1000,10000,100000,1000000), each=3)
		ticks <- c( 3,10,30,100,300,1000,3000,10000,30000,100000,3000000,1000000)
		useTicks <- ticks[ ticks < bigY]
		if ( length( useTicks) < 4) {
			useTicks <- ticks[1:4]
		} else {
			useTicks <- rev( useTicks)[ seq( 1, length(useTicks), by=2)]
		}
		at <- log2( useTicks)
	} else {
		ticks <- seq( from=0, to=1000000, length.out=5)[2:5]
		yBreaks <- c(  800000, 600000, 400000, 200000, 160000, 120000, 80000, 60000, 40000, 20000, 
				16000, 12000, 8000, 6000, 4000, 2000, 1600, 1200, 800, 600, 400, 200, 
				160, 120, 80, 60, 40, 20, 16, 12, 8, 4)
		for ( ycut in yBreaks) {
			if ( bigY < ycut) ticks <- seq( from=0, to=ycut, length.out=5)[2:5]
		}
		at <- useTicks <- ticks
	}


	# draw the wanted style...in the order from tallest to shortest...
	#  can now draw them all at once, but gather in order anyway...
	drawOrder <- base::order( netReadsPerSample, decreasing=TRUE)

	# adjust any colors that are exact matches to differentiate better
	colors <- adjustColorSet( colors)

	for ( i in 1:nWIG) {
		who <- drawOrder[i]
		yData <- yDataList[[who]]
		if ( useLog) {
			yData <- log2Scale( yData, uniqueToo=FALSE)
		}
		myColors <- rep( colors[who], times=3)
		WIG_Draw( yData, col=myColors, plotType=plotType, lwd=lwd, smallestDepth=smallestY,
				show="multiOnly")
	}

	WIG_Ticks( xlim, at, labels=as.integer(useTicks))

	annotateWIGplot( seqID=seqid, gene=gene, xlim=xlim, annoSize, speciesID=curSpecies, 
			altGeneMap=altGeneMap, showDetectability=showDetectability)

	if ( exists( "infoContent")) plotInfoContent( seqid, lwd=2, legendLocation="topright")

	if ( legend.cex > 0) {
	if ( plotType != "boxes") {
		legend( x="topleft", legend=names( WIGlist), col=colors, lwd=4, bg="white", cex=legend.cex)
	} else {
		legend( x="topleft", legend=names( WIGlist), fill=colors, bg="white", cex=legend.cex)
	}
	}

	return()
}


`plotMultiWIGregion` <- function( WIGlist, colors=1:length(WIGlist),  seqid, position=1, end=NULL, 
			addStrands=TRUE, tailWidth=1000, lwd=2,
			plotType=c( "boxes", "lines", "segments", "none"), label="", forceYmax=NULL, minYmax=10, 
			altGeneMap=NULL, showDetectability=TRUE, useLog=FALSE, legend.cex=1.0,
			PLOT.FUN=NULL) {

	if ( ! is.null(PLOT.FUN) && ! is.function(PLOT.FUN)) return()

	plotType <- match.arg( plotType)

	# check species and extract total reads for every dataset
	nWIG <- length( WIGlist)
	if ( nWIG < 2) {
		cat( "Need at least 2 WIG datasets for 'plotMultiWIGgene()'")
		return()
	}
	mySpecies <- WIGlist[[1]]$Info$Species
	dataType <- WIGlist[[1]]$Info$DataType
	if ( is.null( dataType)) dataType <- "RNA-seq"

	totalReadsPerWIG <- vector( length=nWIG)
	for( i in 1:nWIG) {
		thisInfo <- WIGlist[[i]]$Info
		if ( mySpecies != thisInfo$Species) {
			cat( "\nGiven species: ", mySpecies)
			cat( "\nFound species: ", thisInfo$Species)
			stop( "All WIG datasets are not the same speicesID")
		}
		thisCountsM <- thisInfo$RawReads
		# don't let very low counts turn into very high scaling...
		totalReadsPerWIG[i] <- max( sum( thisCountsM), 1000000)
	}

	# set up for this speciesID
	curSpecies <- getCurrentSpecies()
	if ( mySpecies != curSpecies) {
		setCurrentSpecies( mySpecies)
		curSpecies <- getCurrentSpecies()
	}
	seqMap <- getCurrentSeqMap()
	geneMap <- getCurrentGeneMap()
	exonMap <- getCurrentExonMap()
	cdsMap <- getCurrentCdsMap()
	WB_setBinSizeBySpecies( curSpecies)

	# allow alternate gene map
	if ( ! is.null( altGeneMap)) {
		geneMap <- altGeneMap
	}

	smap <- subset.data.frame( seqMap, SEQ_ID==seqid)
	if ( nrow(smap) < 1) {
		cat( "\nNo SeqID with name: ", seqid, " in the current Species.\n")
		return()
	}

	# force exactly one row...
	smap <- as.data.frame( smap[ 1:1, ])

	# turn the gene location and tails into the base and bin limits
	minBase <- 1
	maxBase <- smap$LENGTH
	leftBase <- position
	rightBase <- end
	if ( is.null( leftBase)) leftBase <- minBase
	if ( leftBase < minBase) leftBase <- minBase
	if ( is.null( rightBase)) rightBase <- maxBase
	if ( rightBase > maxBase) rightBase <- maxBase

	# use the gene if given...allow partial matches
	gmap <- subset.data.frame( geneMap, SEQ_ID == seqid & POSITION < rightBase & END < leftBase)

	# scale based on the region given
	leftBaseScaling <- leftBase
	rightBaseScaling <- rightBase
	
	# get all the actual data we need, and normalize for reads per file...
	avgReadsPerWIG <- mean( totalReadsPerWIG)
	yDataList <- vector( mode="list", length=nWIG)
	yMax <- 0
	netReadsPerSample <- vector( length=nWIG)
	for ( i in 1:nWIG) {
		oneWIG <- WIGlist[[i]]	
		wigdata <- WIG_getWigglesOneSeq( oneWIG, seqid)
		myScale <- avgReadsPerWIG / totalReadsPerWIG[i]

		# get the subset in view
		tmpPos <- baseDepthTableSubset( wigdata$Plus, leftBase, rightBase)
		tmpNeg <- baseDepthTableSubset( wigdata$Minus, leftBase, rightBase)
		tmpPos$DEPTH <- tmpPos$DEPTH * myScale
		tmpNeg$DEPTH <- tmpNeg$DEPTH * myScale

		if ( addStrands) {
			tmpPos <- merge.baseDepthTables( tmpPos, tmpNeg)
			tmpNeg <- emptyBaseDepthTable()
		}
		# some decisions are based on just the gene of interest without the flanking regions
		genePos <- baseDepthTableSubset( tmpPos, leftBaseScaling, rightBaseScaling)
		geneNeg <- baseDepthTableSubset( tmpNeg, leftBaseScaling, rightBaseScaling)
		yMax <- max( yMax, genePos$DEPTH, geneNeg$DEPTH)
		netReadsPerSample[i] <- sum( genePos$DEPTH) + sum( geneNeg$DEPTH)
		#clip so they stay in view...
		if ( !is.null( forceYmax)) {
			tmpPos$DEPTH[ tmpPos$DEPTH > forceYmax] <- forceYmax
			tmpNeg$DEPTH[ tmpNeg$DEPTH > forceYmax] <- forceYmax
		}
		# OK, its ready for later drawing...
		yDataList[[i]] <- list( "Plus"=tmpPos, "Minus"=tmpNeg)
	}

	# pick sensible plot limits, given the data and user given limits
	bigY <- max( c( minYmax, yMax))
	if ( !is.null( forceYmax)) bigY <- forceYmax
	xlim <- c(leftBase, rightBase)
	if ( useLog) {
		useBigY <- log2( bigY)
		smallestY <- 0.25
	} else {
		useBigY <- bigY
		smallestY <- useBigY / 800
	}
	ylim <- c( 0, useBigY*1.1)

	# now stretch the lower end of the plot to fit the annotations
	annoBase <- ylim[1] <- ylim[1] - (( ylim[2]-ylim[1]) * 0.20)
	annoSize <- (ylim[2]-ylim[1]) * 0.04

	ylabel <- "Normalized Read Depth"
	if ( useLog) ylabel <- paste( ylabel, "   (log scale)")

	pltTxt <- paste( "RNA_Seq:  ",label) 
	if ( dataType == "ChIP-seq") {
		pltTxt <- paste( "ChIP-seq:  ",label) 
	}
	if ( dataType == "RIP-seq") {
		pltTxt <- paste( "RIP-seq:  ",label) 
	}
	if ( ! is.null( dataType)) pltTxt <- paste( dataType, ":   ",label, sep="") 
	if ( ! is.null( dataType) && dataType == "Peptides") {
		pltTxt <- paste( "Proteomics:  ",label) 
		ylabel <- "Raw Peptide Count"
	}
	pltTxt <- paste( pltTxt, "\n", seqid, ":     ", formatC( leftBase, format="d", big.mark=","),
			"  to  ", formatC( rightBase, format="d", big.mark=","))
	mainText <- pltTxt

	plot( x=xlim, y=c(0,0), main=mainText, type="l", xlim=xlim, ylim=ylim, xaxs="r",
		xlab=paste( "Chromosome: ", seqid, "     Nucleotide Location"), ylab=ylabel, yaxt="n",
		font.axis=2, font.lab=2)

	if ( useLog) {
		ticks <- rep( c( 1,2,5), times=7) * rep( c( 1,10,100,1000,10000,100000,1000000), each=3)
		ticks <- c( 3,10,30,100,300,1000,3000,10000,30000,100000,3000000,1000000)
		useTicks <- ticks[ ticks < bigY]
		if ( length( useTicks) < 4) {
			useTicks <- ticks[1:4]
		} else {
			useTicks <- rev( useTicks)[ seq( 1, length(useTicks), by=2)]
		}
		at <- log2( useTicks)
	} else {
		ticks <- seq( from=0, to=1000000, length.out=5)[2:5]
		yBreaks <- c(  800000, 600000, 400000, 200000, 160000, 120000, 80000, 60000, 40000, 20000, 
				16000, 12000, 8000, 6000, 4000, 2000, 1600, 1200, 800, 600, 400, 200, 
				160, 120, 80, 60, 40, 20, 16, 12, 8, 4)
		for ( ycut in yBreaks) {
			if ( bigY < ycut) ticks <- seq( from=0, to=ycut, length.out=5)[2:5]
		}
		at <- useTicks <- ticks
	}


	# draw the wanted style...in the order from tallest to shortest...
	#  can now draw them all at once, but gather in order anyway...
	drawOrder <- base::order( netReadsPerSample, decreasing=TRUE)

	# adjust any colors that are exact matches to differentiate better
	colors <- adjustColorSet( colors)

	for ( i in 1:nWIG) {
		who <- drawOrder[i]
		yData <- yDataList[[who]]
		if ( useLog) {
			yData <- log2Scale( yData, uniqueToo=FALSE)
		}
		myColors <- rep( colors[who], times=3)
		WIG_Draw( yData, col=myColors, plotType=plotType, lwd=lwd, smallestDepth=smallestY,
				show="multiOnly")
	}

	WIG_Ticks( xlim, at, labels=as.integer(useTicks))

	annotateWIGplot( seqID=seqid, gene=NULL, xlim=xlim, annoSize, speciesID=curSpecies, 
			altGeneMap=altGeneMap, showDetectability=showDetectability)

	if ( exists( "infoContent")) plotInfoContent( seqid, lwd=2, legendLocation="topright")

	if ( legend.cex > 0) {
	if ( plotType != "boxes") {
		legend( x="topleft", legend=names( WIGlist), col=colors, lwd=4, bg="white", cex=legend.cex)
	} else {
		legend( x="topleft", legend=names( WIGlist), fill=colors, bg="white", cex=legend.cex)
	}
	}

	return()
}


`log2Scale` <- function( yData, uniqueToo=TRUE) {

	# fix the very low intensity data, because every '1 or smaller' would be zero or negative
	tmp <- yData$Plus$DEPTH
	yData$Plus$DEPTH <- ifelse( tmp < 2, tmp/2, log2(tmp))
	tmp <- yData$Minus$DEPTH
	yData$Minus$DEPTH <- ifelse( tmp < 2, tmp/2, log2(tmp))
	if ( "Combo" %in% names(yData)) {
		tmp <- yData$Combo$DEPTH
		yData$Combo$DEPTH <- ifelse( tmp < 2, tmp/2, log2(tmp))
	}
	if ( ! uniqueToo) return( yData)

	tmp <- yData$PlusUnique$DEPTH
	yData$PlusUnique$DEPTH <- ifelse( tmp < 2, tmp/2, log2(tmp))
	tmp <- yData$MinusUnique$DEPTH
	yData$MinusUnique$DEPTH <- ifelse( tmp < 2, tmp/2, log2(tmp))

	return( yData)
}


`yDataScale` <- function( yData, scaleFactor=1, uniqueToo=TRUE) {

	tmp <- yData$Plus$DEPTH
	yData$Plus$DEPTH <- tmp * scaleFactor
	tmp <- yData$Minus$DEPTH
	yData$Minus$DEPTH <- tmp * scaleFactor
	if ( "Combo" %in% names(yData)) {
		tmp <- yData$Combo$DEPTH
		yData$Combo$DEPTH <- tmp * scaleFactor
	}
	if ( ! uniqueToo) return( yData)

	tmp <- yData$PlusUnique$DEPTH
	yData$PlusUnique$DEPTH <- tmp * scaleFactor
	tmp <- yData$MinusUnique$DEPTH
	yData$MinusUnique$DEPTH <- tmp * scaleFactor

	return( yData)
}


`WIG_Ticks` <- function( xlim, ticks, labels) {

	xextra <- diff( xlim) * 0.5
	xlim2 <- xlim + c( -xextra, xextra)
	for( i in 1:length( ticks)) {
		y <- rep( ticks[i], times=2)
		lines( x=xlim2, y=y, lty=3)
		text( x=xlim, y=y, label=rep( labels[i], times=2), cex=0.95, pos=3, font=2)
	}
}


`WIG_Draw` <- function( yData, col=c(4,2,1), plotType="boxes", lwd=2, smallestDepth=0,
			show=c("both", "multiOnly")) {

	if ( plotType == "none") return()
	userUnits <- par("usr")
	deltaX <- userUnits[2] - userUnits[1]
	DXperPixel <- deltaX / 1000
	show <- match.arg( show)

	# turn the location and counts of wiggle bars into line segments on the plot
	if ( plotType == "lines") {
		WIG_Lines( yData, col, lwd=lwd)
		return()
	}
	if ( plotType == "segments") {
		WIG_Segments( yData, col, lwd=lwd)
		return()
	}
	if ( plotType == "addedLines") {
		WIG_AddedLines( yData, c( col, "magenta"), lwd=lwd)
		return()
	}

	# default is filled boxes
	if ( DXperPixel > 100) {
		WIG_Strokes( yData, col, smallestDepth=smallestDepth, lwd=lwd)
	} else {
		WIG_Bars( yData, col, DXperPixel, smallestDepth=smallestDepth, show=show)
	}
}


`WIG_Strokes` <- function( yData, col=c(4,2,1), smallestDepth=0, lwd=2) {

	# we know there's lots of bases per pixel, so try to reduce these to less rows
	ypos <- yData$Plus
	if ( nrow(ypos) > 0) {
		xpix <- as.integer( round( ypos$START / 100) * 100)
		xpixFac <- factor(xpix)
		xpixDepth <- tapply( ypos$DEPTH, xpixFac, max)
		xpix <- as.integer( levels(xpixFac))
		ypos <- data.frame( "START"=xpix, "DEPTH"=xpixDepth)
	}

	yneg <- yData$Minus
	if ( nrow(yneg) > 0) {
		xpix <- as.integer( round( yneg$START / 100) * 100)
		xpixFac <- factor(xpix)
		xpixDepth <- tapply( yneg$DEPTH, xpixFac, max)
		xpix <- as.integer( levels(xpixFac))
		yneg <- data.frame( "START"=xpix, "DEPTH"=xpixDepth)
	}

	ycombo <- NULL
	if ( "Combo" %in% names(yData)) {
	  ycombo <- yData$Combo
	  if ( nrow(ycombo) > 0) {
		xpix <- as.integer( round( ycombo$START / 100) * 100)
		xpixFac <- factor(xpix)
		xpixDepth <- tapply( ycombo$DEPTH, xpixFac, max)
		xpix <- as.integer( levels(xpixFac))
		ycombo <- data.frame( "START"=xpix, "DEPTH"=xpixDepth)
	  }
	}
	# at each line, decide who to draw first...
	ypos$COLOR <- rep( col[1], times=nrow( ypos))
	yneg$COLOR <- rep( col[2], times=nrow( yneg))
	yDF <- rbind( ypos, yneg)
	if ( ! is.null( ycombo)) {
		ycombo$COLOR <- rep( col[3], times=nrow( ycombo))
		yDF <- rbind( yDF, ycombo)
	}
	if ( nrow(yDF) < 1) return()

	if ( smallestDepth > 0) {
		drops <- which( yDF$DEPTH < smallestDepth)
		if ( length(drops) > 0) yDF <- yDF[ -drops, ]
		if ( nrow(yDF) < 1) return()
	}

	ord <- base::order( yDF$DEPTH, decreasing=TRUE)
	yDF <- yDF[ ord, ]
	tooSmall <- max( yDF$DEPTH) * 0.001

	# try to do without a For loop
	todo <- which( yDF$DEPTH > tooSmall)
	if (length(todo) < 1) return()
	xx <- yDF$START[ todo]
	yy <- yDF$DEPTH[ todo]
	cols <- yDF$COLOR[ todo]
	lines( xx, yy, type="h", lty=1, col=cols, lwd=lwd)
	return()
}


`WIG_Bars` <- function( yData, col=c(4,2,1), DXperPixel=1, mode=c("solid","shaded"), 
			smallestDepth=0, show=c("both", "multiOnly")) {

	mode <- match.arg( mode)
	show <- match.arg( show)

	# if we are doing a single sample, then do the multi's first and shade them differently
	shadeMultis <- FALSE
	if ( show == "both") shadeMultis <- TRUE
	if ( shadeMultis && mode == "solid") {

		# let's not draw any shading that will get completely covered by the unique wiggles
		visibleYdata <- visibleMultiWIGs( yData)
		WIG_Bars( visibleYdata, col=col, DXperPixel=DXperPixel, mode="shaded")

		# now relabel the uniques to act like normal
		yData <- list( "Plus"=yData$PlusUnique, "Minus"=yData$MinusUnique)
		mode <- "solid"
	}

	# extract what we need to know where to draw
	yDF <- data.frame()
	ypos <- yData$Plus
	if ( ! is.null(ypos)) {
		ypos$COLOR <- rep( col[1], times=nrow( ypos))
		yDF <- ypos
	}
	yneg <- yData$Minus
	if ( ! is.null(yneg)) {
		yneg$COLOR <- rep( col[2], times=nrow( yneg))
		yDF <- rbind( yDF, yneg)
	}
	ycombo <- yData$Combo
	if ( ! is.null( ycombo)) {
		ycombo$COLOR <- rep( col[3], times=nrow( ycombo))
		yDF <- rbind( yDF, ycombo)
	}

	# special mode for multi-WIG plots that allow lots of samples together...
	nDF <- length(yData)
	if ( is.null( ycombo) && nDF > 2) {
		# there are more sets of Plus/Minus here, add them in too...
		for (j in 3:nDF) {
			ymore <- yData[[j]]
			ymore$COLOR <- rep( col[j], times=nrow( ymore))
			yDF <- rbind( yDF, ymore)
		}
	}

	if ( nrow(yDF) < 1) return()

	if ( smallestDepth > 0) {
		drops <- which( yDF$DEPTH < smallestDepth)
		if ( length(drops) > 0) yDF <- yDF[ -drops, ]
		if ( nrow(yDF) < 1) return()
	}

	# we will draw tallest to shortest
	ord <- base::order( yDF$DEPTH, decreasing=TRUE)
	yDF <- yDF[ ord, ]

	# turn the X's and Y's into box corners
	x1 <- yDF$START - 0.5
	x2 <- yDF$STOP + 0.5
	tooThin <- which((x2-x1) < DXperPixel)
	if ( length(tooThin) > 0) {
		DXhalf <- DXperPixel * 0.5
		x1[tooThin] <- x1[tooThin] - DXhalf
		x2[tooThin] <- x2[tooThin] + DXhalf
	}
	y <- yDF$DEPTH
	zero <- rep( 0, times=nrow(yDF))

	allCol <- yDF$COLOR
	if ( mode == "shaded") {
		#rect( x1, zero, x2, y, density=20, col=allCol, border=NA)
		rect( x1, zero, x2, y, col=adjustColor(allCol,0.65), border=NA)
	} else { 
		rect( x1, zero, x2, y, col=allCol, border=allCol)
	}
	return()
}


`WIG_Lines` <- function( yData, col=c( 4,2,1), lwd=2) {

	# just draw one line along X for each strand
	if ( "Combo" %in% names(yData)) {
		drawOneWIGline( yData$Combo, col[3], lwd=lwd)
	}
	drawOneWIGline( yData$Plus, col[1], lwd=lwd)
	drawOneWIGline( yData$Minus, col[2], lwd=lwd)
}


`WIG_Segments` <- function( yData, col=c( 4,2,1), lwd=3) {

	# just draw one line along X for each strand
	if ( "Combo" %in% names(yData)) {
		drawOneWIGsegment( yData$Combo, col[3], lwd=lwd)
	}
	drawOneWIGsegment( yData$Plus, col[1], lwd=lwd)
	drawOneWIGsegment( yData$Minus, col[2], lwd=lwd)
}


`WIG_AddedLines` <- function( yData, col=c( 4,2,"magenta"), lwd=2) {

	# just draw one line along X for each strand
	drawOneWIGline( yData$Plus, col[1], lwd=lwd)
	drawOneWIGline( yData$Minus, col[2], lwd=lwd)
	summedData <- baseDepthVectorToTable( mergeIntegerTables( baseDepthTableToVector( yData$Plus),
					baseDepthTableToVector( yData$Minus)))
	drawOneWIGline( summedData, col[3], lwd=lwd)
}


`drawOneWIGline` <- function( mydf, col=1, lwd=2) {

	if ( nrow(mydf) < 1) return()
	yval <- mydf$DEPTH
	xstart <- as.integer(mydf$START)
	xstop <- as.integer(mydf$STOP)

	# turn these into a vector with NA where we do not draw..
	# extend both edges by 1 so can catch the first UP, last DOWN events
	minX <- as.integer( round( min( xstart)))
	maxX <- as.integer( round( max( xstop)))
	allX <- minX : maxX
	N <- length( allX)
	allY <- rep( 0, times=N)
	allY[ xstart - minX + 1] <- yval

	needFill <- which( xstart != xstop)
	for( j in needFill) {
		beg <- xstart[j] - minX + 1
		end <- xstop[j] - minX + 1
		allY[ beg : end] <- yval[j]
	}
	
	# no line where its all zero
	if ( N > 2) {
		flatZero <- ( (allY[1:(N-2)] == 0) & (allY[2:(N-1)] == 0) & (allY[3:N] == 0))
		allY[ (which(flatZero)+1)] <- NA
	}

	# put explicit zeros at the ends to bring lines back to base?
	if (FORCE_ZERO_EDGES <- FALSE) {
		xx <- c( minX, allX, maxX)
		yy <- c( 0, allY, 0)
		lines( x=xx, y=yy, type="l", lwd=lwd, col=col)
	} else {
		lines( x=allX, y=allY, type="l", lwd=lwd, col=col)
	}
}


`drawOneWIGsegment` <- function( mydf, col=1, lwd=3) {

	if ( nrow(mydf) < 1) return()
	yval <- mydf$DEPTH
	xstart <- mydf$START
	xstop <- mydf$STOP

	# draw the little level segments
	dashedLine( xb=xstart, xe=xstop, y=yval, lwd=lwd, col=col)
}


dashedLine <- function( xb, xe, y, ...) {

	# turn the begin/end X coordinates into a dashed line...  "lines" lets the NA show
	# where gaps are...

	# turn the y's into y1,y1,NA,y2,y2,NA, ...
	ynew <- rep( y, each=3)
	npts <- length( ynew)
	ynew[ seq( 3, npts, by=3)] <- NA
	# the X's get xb1,xe1,NA,xb2,xe2,NA, etc..
	# the X's get xb1,xe1,NA,xb2,xe2,NA, etc..
	xnew <- ynew
	xnew[ seq( 1, (npts-2), by=3)] <- xb
	xnew[ seq( 2, (npts-1), by=3)] <- xe
	
	# draw that line
	lines( xnew, ynew, type="l", ...)
}


`annotateWIGplot` <- function( seqID,  gene=NULL, xlim, annoSize, speciesID, offset=0, 
		altGeneMap=NULL, showDetectability=TRUE) {

	# given the plot region in bases, ( not bins)...
	leftBase <- xlim[1]
	rightBase <- xlim[2]
	seqid <- seqID
	geneMap <- getCurrentGeneMap()
	exonMap <- getCurrentExonMap()
	cdsMap <- getCurrentCdsMap()

	# get all the genes we cover
	gmap <- subset.data.frame( geneMap, ((SEQ_ID == seqid) & (POSITION <= rightBase) & 
				(END >= leftBase) & REAL_G == TRUE))
	if ( nrow(gmap) == 0 && !showDetectability) return()
	emap <- subset.data.frame( exonMap, ((SEQ_ID == seqid) & (POSITION <= rightBase) & (END >= leftBase)))
	cmap <- subset.data.frame( cdsMap, ((SEQ_ID == seqid) & (POSITION <= rightBase) & (END >= leftBase)))

	# there are some things that don't apply when the # of bases gets to high
	dx <- rightBase - leftBase
	fastMode <- ( dx > 1000000 && nrow(gmap) > 100)

	# try to trap the rare case of more than one gene in the window, but one of them overlaps all others
	if ( !fastMode && sum( gmap$REAL_G) > 1) {
		allsizes <- base::pmin( rightBase,gmap$END) - base::pmax( leftBase,gmap$POSITION)
		strnd <- ifelse( is.na(gmap$STRAND), " ", gmap$STRAND)
		real <- gmap$REAL_G
		isOversize <- which( (allsizes >= (rightBase-leftBase-10)) & real)
		if ( length(isOversize) > 0) {
			drops <- vector()
			for (k in isOversize) {
				nSameStrand <- sum( strnd == strnd[k])
				if ( nSameStrand > 1) drops <- c( drops, k)
			}
			if ( length(drops) > 0) {
				gmap <- gmap[ -drops, ]
				emap <- subset.data.frame( emap, GENE_ID %in% gmap$GENE_ID)
				cmap <- subset.data.frame( cmap, GENE_ID %in% gmap$GENE_ID)
			}
		}
	}

	# set up where things get drawn
	gap <- annoSize * 0.15
	exonYloPlus <- offset - annoSize*2 + gap
	exonYhiPlus <- offset - annoSize - gap
	exonYloMinus <- offset - annoSize*4 + gap
	exonYhiMinus <- offset - annoSize*3 - gap
	varYtoggleSize <- gap * 2

	# label and outline box for each 'real' gene
	# don't include 'locus' genes
	isLOCUS <- grep( '@', gmap$GENE_ID, fixed=T)
	if ( length(isLOCUS)) gmap <- gmap[ -isLOCUS, ]

	# some features based on how many in the window
	nRealGenes <- sum( gmap$REAL_G)
	if ( !fastMode && nrow( gmap) > 0) {
	   	ggmap <- subset.data.frame( gmap, REAL_G)
		xleft <- ggmap$POSITION
		xright <- ggmap$END
		strand <- ifelse( is.na(ggmap$STRAND), " ", ggmap$STRAND)
		strandPN <- (strand %in% c( "+", " "))
		exonYlo <- ifelse( strandPN, exonYloPlus, exonYloMinus)
		exonYhi <- ifelse( strandPN, exonYhiPlus, exonYhiMinus)
		exonYtext <- ifelse( strandPN, exonYhiPlus-gap, exonYloMinus+gap)
		exonTextPos <- ifelse( strandPN, 3, 1)
		rect( xleft, exonYlo, xright, exonYhi, border=1, lwd=1)

		textX <- (xleft+xright)/2
		textX <- pmax( textX, leftBase)
		textX <- pmin( textX, rightBase)
		mycex <- 1.1 - (nRealGenes*0.03)
		if ( mycex < 0.2) mycex <- 0.2
		if ( nRealGenes < 30 || !is.null(gene)) {
			geneTxt <- ggmap$NAME
			text( x=textX, y=exonYtext, labels=paste( strand, geneTxt, 
					strand, sep="  "), pos=exonTextPos, cex=mycex, font=2)
		}
	}

	# exons get filled box
	# as of summer 2025, most genome annotations are good enough, to support drawing CDS regions instead of exons
	# if there are no coding CDS entries, show the exons instead
	if ( ! nrow(cmap)) cmap <- emap
	if ( nrow( cmap) > 0) {
	   	if ( nrow( cmap) > 1) {
			allsizes <- base::pmin( rightBase,cmap$END) - base::pmax( leftBase,cmap$POSITION)
			isOversize <- which( allsizes >= (rightBase-leftBase-10))
			if ( (length(isOversize) > 0) && !is.null(gene) && !(gene %in% cmap$GENE_ID)) {
				if ( length( isOversize) < nrow( cmap)) {
					cmap <- cmap[ -isOversize, ]
				} else {
					cmap <- cmap[ -isOVersize[ which.max(allsizes)], ]
				}
			}
	   	}
		xleft <- cmap$POSITION
		xright <- cmap$END
		strand <- ifelse( is.na(cmap$STRAND), " ", cmap$STRAND)
		exonYlo <- ifelse( strand %in% c( "+", " "), exonYloPlus, exonYloMinus)
		exonYhi <- ifelse( strand %in% c( "+", " "), exonYhiPlus, exonYhiMinus)
		rect( xleft, exonYlo, xright, exonYhi, col=1)
	}

	# go back and white out a tiny spot on the gene edges to distinguish consecutive genes
	npixels <- (rightBase - leftBase) / 2000
	npx2 <- npixels * 2
	npb2 <- npixels / 2
	gmap <- subset.data.frame( gmap, REAL_G == TRUE)
	if ( !fastMode && nrow(gmap) > 0 && nrow(gmap) < 30) {
	   for( i in 1:nrow( gmap)) {
		xleft <- gmap$POSITION[i]
		xright <- gmap$END[i]
		strand <- if ( is.na( gmap$STRAND[i])) " " else gmap$STRAND[i]
		if ( strand %in% c( "+", " ")) {
			exonYlo <- exonYloPlus;  exonYhi <- exonYhiPlus;
		} else {
			exonYlo <- exonYloMinus;  exonYhi <- exonYhiMinus;
		}
		# do we need to white the left edge?
		if ( i > 1) {
			if ( (gmap$STRAND[i] == gmap$STRAND[i-1]) && (gmap$POSITION[i] - gmap$END[i-1] < npx2)) {
				rect( xleft-npb2, exonYlo, xleft+npb2, exonYhi, border="white", col="white")
			}
		}
		# do we need to white the right edge?
		if ( i < nrow(gmap)) {
			if ( (gmap$STRAND[i] == gmap$STRAND[i+1]) && (gmap$POSITION[i+1] - gmap$END[i] < npx2)) {
				rect( xright-npb2, exonYlo, xright+npb2, exonYhi, border="white", col="white")
			}
		}
	   }
	}

	# the alt gene map...
	if ( ! is.null( altGeneMap)) {
		# get all the alt genes we cover
		gmap <- subset.data.frame( altGeneMap, ((SEQ_ID == seqid) & (POSITION <= rightBase) & 
				(END >= leftBase)))
		if ( nrow( gmap) > 0) {
			for( i in 1:nrow( gmap)) {
				xleft <- gmap$POSITION[i]
				xright <- gmap$END[i]
				strand <- if ( is.na( gmap$STRAND[i])) " " else gmap$STRAND[i]
				if ( strand %in% c( "+", " ")) {
					exonYlo <- exonYloPlus;  exonYhi <- exonYhiPlus;
					exonYtext <- exonYhiPlus - gap;   exonTextPos <- 3;
				} else {
					exonYlo <- exonYloMinus;  exonYhi <- exonYhiMinus;
					exonYtext <- exonYloMinus + gap;   exonTextPos <- 1;
				}
				rect( xleft, exonYlo, xright, exonYhi, col="magenta", border=1, lwd=1)
				textX <- (xleft+xright)/2
				if (textX < leftBase) textX <- leftBase
				if (textX > rightBase) textX <- rightBase
				mycex <- 0.9 - (nrow(gmap)*0.03)
				if ( mycex < 0.2) mycex <- 0.2
				if ( nRealGenes < 12 || !is.null(gene)) {
					geneTxt <- gmap$NAME[i]
					text( x=textX, y=exonYtext, labels=geneTxt, pos=exonTextPos, 
						cex=mycex, col="magenta", font=2)
				}
			}
		}
	}

	# var gene domains
	dmap <- getVargeneDomains( gmap$GENE_ID)
	if ( nrow( dmap)) {
		strands <- ifelse( is.na( dmap$STRAND), " ", dmap$STRAND)
		vgname <- ""
		for ( i in 1:nrow( dmap)) {
			xleft <- dmap$DNA_START[i]
			xright <- dmap$DNA_STOP[i]
			strand <- strands[i]
			if ( dmap$GENE_ID[i] != vgname) {
				if ( strand %in% c( "+", " ")) {
					varYlo <- exonYloMinus - gap;  varYhi <- exonYhiMinus - gap;
				} else {
					varYlo <- exonYloPlus + gap;  varYhi <- exonYhiPlus + gap;
				}
				varYcenter <- mean( c( varYlo, varYhi))
				toggle <- 1
				vgname <- dmap$GENE_ID[i]
			}
			ytoggle <- varYtoggleSize * toggle
			lines( x=c(xleft,xright), y=rep( varYcenter + ytoggle, 2), lty=1, lwd=8, 
					col="goldenrod", lend=1)
			textX <- (xleft+xright)/2
			domText <- dmap$DOMAIN_ID[i]
			yOffText <- varYcenter - ytoggle
			yOffDC <- varYcenter + ytoggle*3.5
			if ( regexpr( "_", domText) > 0) {
				domText <- sub( "_", "\n", domText)
				yOffText <- varYcenter - ytoggle*1.8
			}
			text( x=textX, y=yOffText, labels=domText, pos=NULL, cex=mycex*0.88, font=2)
			myDC <- dmap$CASSETTE[i]
			if (myDC != "") {
				lines( x=c(xleft,xright), y=rep( varYcenter + ytoggle*2.1, 2), lty=1, lwd=6, 
						col="magenta", lend=1)
			text( x=textX, y=yOffDC, labels=myDC, pos=NULL, cex=mycex*0.88, font=2)
			}
			toggle <- ( - toggle)
		}
	} 

	# label the bins with how many hits they 'could' get
	mySpecies <- getCurrentSpecies()

	# there are a few species that will NEVER have detectability data...
	if ( mySpecies %in% c( "VarGenes", "VSA", "JOS")) return()

	WB_setSpeciesDetectability( )
	if ( showDetectability) {
	  if ( (!exists( "curDetectSelf")) || (curDetectSelf$Species != mySpecies) ) {
		curDetectSelf <<- WB_getCurrentSelfDetectability( )
		curDetectOther <<- WB_getCurrentOtherDetectability( )
	  }
	} 

	if ( showDetectability && !is.null( curDetectSelf) && !is.null( curDetectSelf$Unique)) {
	    binsetptr <- WB_getBinSetPtrFromSeqID( curDetectSelf$Unique, seqid)
	    # if we fail to find the detectability, flag it as not being done back up in the caller routine and quit
	    if ( binsetptr < 1) {
	    	showDetectability <<- FALSE
		return()
	    }
	    thisBinset <- curDetectSelf$Unique$BinData[[ binsetptr]]

	    # setup for one 'other species too'
	    hasOther <- FALSE
	    if ( !is.null( curDetectOther) && !is.null( curDetectOther$Detectable)) {
	    	otherbinsetptr <- WB_getBinSetPtrFromSeqID( curDetectOther$Unique, seqid)
	    	otherBinset <- curDetectOther$Detectable$BinData[[ binsetptr]]
	    	if ( otherbinsetptr >= 1) {
		    hasOther <- TRUE
	        }
	    }

	    firstBin <- WB_base2bin( leftBase)
	    lastBin <- WB_base2bin( rightBase)
	    ylow <- exonYhiMinus + gap
	    ymax <- exonYloPlus - gap
	    dy <- ymax - ylow
	    asbox <- ((lastBin - firstBin) < 500)

	    # try to draw them without a for loop
	    mybins <-  firstBin:lastBin
	    xl <- WB_bin2base( mybins, "left")
	    xr <- WB_bin2base( mybins, "right")
	    yl <- rep( (offset + ylow), times=length(xl))
	    thisPcts <- thisBinset[ mybins, WB_UNIQUE_PCT]
	    yh <- yl + (dy * thisPcts)
	    toDo <- which( thisPcts > 0)
	    if ( length(toDo) > 0) {
	    if ( asbox) {
		rect( xl[toDo], yl[toDo], xr[toDo], yh[toDo], border=3, density=NULL, col=NA)
	    } else {
	    	# make segments with NA's
		Ntodo <- length(toDo)  *3
		xx <- rep( NA, times=Ntodo)
		yy <- rep( NA, times=Ntodo)
		ats <- seq( 1, Ntodo, by=3)
		xx[ ats] <- xl[toDo]
		xx[ ats+1] <- xl[toDo]
		yy[ ats] <- yl[toDo]
		yy[ ats+1] <- yh[toDo]
		lines( xx, yy, col=3)
	    }}
	    if ( hasOther) {
	    	otherPcts <- otherBinset[ mybins, WB_UNIQUE_PCT]
	    	yh <- yl + (dy * otherPcts)
	    	toDo <- which( otherPcts > 0)
	    	if ( length(toDo) > 0) {
	    	if ( asbox) {
			rect( xl[toDo], yl[toDo], xr[toDo], yh[toDo], border=2, density=NULL, col=NA)
	    	} else {
	    		# make segments with NA's
			Ntodo <- length(toDo)  *3
			xx <- rep( NA, times=Ntodo)
			yy <- rep( NA, times=Ntodo)
			ats <- seq( 1, Ntodo, by=3)
			xx[ ats] <- xl[toDo]
			xx[ ats+1] <- xl[toDo]
			yy[ ats] <- yl[toDo]
			yy[ ats+1] <- yh[toDo]
			lines( xx, yy, col=2)
	    	}}
	    }

	}
}



`verifyWIGspecies` <- function( WIG, seqID=NULL, geneID=NULL) {

	allSeqs <- rownames( WIG$Info$RawReads)
	if ( ! is.null( seqID)) {
		if ( ! (seqID %in% allSeqs)) {
			cat( "\nWIG data does not contain seqID: ", seqID)
			return(FALSE)
		}
		if ( ! (seqID %in% getCurrentSeqMap()$SEQ_ID)) {
			cat( "\nCurrent Species does not contain seqID: ", seqID, getCurrentSpecies(), "\n")
			return(FALSE)
		}
	}
	if ( ! is.null( geneID)) {
		gmap <- getCurrentGeneMap()
		gHits <- grep( geneID, gmap$GENE_ID, fixed=TRUE)
		if ( ! ( length( gHits) > 0)) {
			cat( "\nCurrent Species does not contain geneID: ", geneID, getCurrentSpecies(), "\n")
			return(FALSE)
		}
		seqID <- gmap$SEQ_ID[ gHits[1]]
		if ( ! (seqID %in% allSeqs)) {
			cat( "\nWIG data does not contain seqID: ", seqID, "\n")
			return(FALSE)
		}
	}
	if ( ! ( getSpeciesFromSeqID( allSeqs[1]) %in% getCurrentSpecies())) {
		cat( "\nCurrent WIG dataset does not match current species\n")
		return(FALSE)
	}
	if ( ! ( getSpeciesFromSeqID( allSeqs[1]) %in% getCurrentTargetSpecies())) {
		cat( "Warning:  current WIG species not in current target species set.", 
			"\nBin Detectability effected.  use 'setCurrentTarget()' to fix...\n")
	}

	return( TRUE)
}


`makeAllWIGgenePlots` <- function( WIG, geneSet=NULL, path="png", fileExtension="png", tailWidth=2000, 
			plotType=c( "boxes", "lines", "segments"), label="", minYmax=5, useLog=FALSE, 
			pause=NULL, geneNameColumn="GENE_ID", altGeneMap=NULL, PLOT.FUN=NULL, ...) {

	if ( ! is.null(PLOT.FUN) && ! is.function(PLOT.FUN)) return()

	cat( "\n")
	geneMap <- getCurrentGeneMap()
	if ( ! is.null( altGeneMap)) geneMap <- altGeneMap

	if ( is.null( geneSet)) {
		geneSet <- geneMap$GENE_ID
	}

	if ( ! file.exists( path)) dir.create( path, recursive=TRUE, showWarnings=FALSE)

	# order by seqID for speedier WIG access
	if ( ! is.null( WIG)) {
		ptrs <- match( geneSet, geneMap$GENE_ID)
		ord <- order( ptrs)
		geneSet <- geneSet[ ord]
	}

	plotType <- match.arg( plotType)

	for ( i in 1:length( geneSet)) {
		g <- geneSet[i]

		gplotname <- g
		if ( geneNameColumn != "GENE_ID") {
			where <- match( g, geneMap$GENE_ID, nomatch=0)
			if ( where > 0) gplotname <- geneMap[[geneNameColumn]][ where]
		}

		# allow draw to screen too
		if ( ! is.null(pause)) {
			if ( is.null( PLOT.FUN)) {
				plotWIGgene( WIG, gene=g, tailWidth=tailWidth, plotType=plotType, label=label, 
					minYmax=minYmax, useLog=useLog, altGeneMap=altGeneMap, ...)
			} else {
				PLOT.FUN( g)
			}
			if ( pause > 0) Sys.sleep( pause)
			dev.flush()
		}

		file <- paste( gplotname, fileExtension, sep=".")
		# make sure its a safe name
		file <- file.cleanSpecialCharactersFromFileName( file)
		f <- file.path( path, file)
		# plot type decided by new options.txt based setttings, here is now generic
		openPlot( f, bg="white")
		if ( is.null( PLOT.FUN)) {
			plotWIGgene( WIG, gene=g, tailWidth=tailWidth, plotType=plotType, label=label, minYmax=minYmax,
				useLog=useLog, altGeneMap=altGeneMap, ...)
		} else {
			PLOT.FUN( g)
		}
		dev.flush()
		dev.off()
		cat( "\r", i, "\t", g, "    ")
	}
	if ( exists("curDetectSelf", envir=.GlobalEnv)) rm( curDetectSelf, envir=.GlobalEnv)
	if ( exists("curDetectOther", envir=.GlobalEnv)) rm( curDetectOther, envir=.GlobalEnv)
}


`makeAllMultiWIGgenePlots` <- function( WIGlist, colors, geneSet=NULL, path="png", fileExtension="png", 
					tailWidth=2000, plotType=c( "boxes", "lines", "segments", "addedLines"), 
					label="", forceYmax=NULL, minYmax=5, useLog=FALSE, pause=NULL,
					geneNameColumn="GENE_ID", altGeneMap=NULL, PLOT.FUN=NULL, ...) {

	if ( ! is.null(PLOT.FUN) && ! is.function(PLOT.FUN)) return()

	cat( "\n")
	geneMap <- getCurrentGeneMap()
	if ( ! is.null( altGeneMap)) geneMap <- altGeneMap

	if ( is.null( geneSet)) {
		geneSet <- geneMap$GENE_ID
	}

	plotType <- match.arg( plotType)

	if ( ! file.exists( path)) dir.create( path, recursive=TRUE, showWarnings=FALSE)

	# order by seqID for speedier WIG access
	if ( ! is.null( WIGlist)) {
		ptrs <- match( geneSet, geneMap$GENE_ID)
		ord <- order( ptrs)
		geneSet <- geneSet[ ord]
	}

	for ( i in 1:length(geneSet)) {
		g <- geneSet[i]

		gplotname <- g
		if ( geneNameColumn != "GENE_ID") {
			where <- match( g, geneMap$GENE_ID, nomatch=0)
			if ( where > 0) gplotname <- geneMap[[geneNameColumn]][ where]
		}

		# allow draw to screen too
		if ( ! is.null(pause)) {
			if ( ! is.null( PLOT.FUN)) {
				PLOT.FUN(g)
			} else {
				plotMultiWIGgene( WIGlist, colors, g, plotType=plotType, tailWidth=tailWidth, 
					forceYmax=forceYmax, minYmax=minYmax, label=label, useLog=useLog,
					altGeneMap=altGeneMap, ...)
			}
			if ( pause > 0) Sys.sleep( pause)
			dev.flush()
		}

		file <- paste( gplotname, fileExtension, sep=".")
		# make sure its a safe name
		file <- file.cleanSpecialCharactersFromFileName( file)
		f <- file.path( path, file)
		# plot type decided by new options.txt based setttings, here is now generic
		openPlot( f, bg="white")
		if ( ! is.null( PLOT.FUN)) {
			PLOT.FUN(g)
		} else {
			plotMultiWIGgene( WIGlist, colors, g, plotType=plotType, tailWidth=tailWidth, forceYmax=forceYmax, 
					minYmax=minYmax, label=label, useLog=useLog,
					altGeneMap=altGeneMap, ...)
		}
		dev.flush()
		dev.off()
		cat( "\r", i, "\t", g, "    ")
	}
	cat( "\nDone.\n")
	if ( exists("curDetectSelf", envir=.GlobalEnv)) rm( curDetectSelf, envir=.GlobalEnv)
	if ( exists("curDetectOther", envir=.GlobalEnv)) rm( curDetectOther, envir=.GlobalEnv)
}


`visibleMultiWIGs` <- function( yData) {

	# flatten out the base depth tables, and then find the parts that are higher in the multis
	if ( is.null( yData$PlusUnique) || nrow( yData$PlusUnique) < 1) {
		outMultiPlus <- yData$Plus
	} else {
		uniq <- baseDepthTableToVector( yData$PlusUnique)
		multi <- baseDepthTableToVector( yData$Plus)
		outMultiPlus <- yData$Plus
		# find the locations they have in common
		locs <- intersect( names(uniq), names(multi))
		if ( length(locs) > 0) {
			whU <- match( locs, names(uniq))
			whM <- match( locs, names(multi))
			# places of equal height will not be seen
			drops <- which( uniq[ whU] >= multi[whM])
			if ( length( drops) > 0) {
				newmulti <- multi[ -whM[drops]]
				outMultiPlus <- baseDepthVectorToTable( newmulti)
			}
		}
	}

	# do the same for minus strand
	if ( is.null( yData$MinusUnique) || nrow( yData$MinusUnique) < 1) {
		outMultiMinus <- yData$Minus
	} else {
		uniq <- baseDepthTableToVector( yData$MinusUnique)
		multi <- baseDepthTableToVector( yData$Minus)
		outMultiMinus <- yData$Minus
		locs <- intersect( names(uniq), names(multi))
		if ( length(locs) > 0) {
			whU <- match( locs, names(uniq))
			whM <- match( locs, names(multi))
			drops <- which( uniq[ whU] >= multi[whM])
			if ( length( drops) > 0) {
				newmulti <- multi[ -whM[drops]]
				outMultiMinus <- baseDepthVectorToTable( newmulti)
			}
		}
	}

	out <- list( "Plus"=outMultiPlus, "Minus"=outMultiMinus)
	return(out)
}


`plotChIPpeak` <- function( WIG, chipTbl, chipRow=1, tailWidth=200, plotType=c("lines", "boxes"), 
			label="", useLog=FALSE, forceYmax=NULL, minYmax=10, legend.cex=1.15,
			scaledReadCount=NULL) {

	plotType <- match.arg( plotType)

	presentationMode <- TRUE
	cex.lab <- 1.0
	if (presentationMode) cex.lab=1.5


	# check species
	#if ( ! verifyWIGspecies( WIG, seqID=NULL, geneID=gene)) return()

	curSpecies <- getCurrentSpecies()
	seqMap <- getCurrentSeqMap()
	geneMap <- getCurrentGeneMap()
	exonMap <- getCurrentExonMap()

	# estimate a gene if we need
	if ( ! "GENE_ID" %in% colnames(chipTbl)) {
		gmap <- subset( geneMap, REAL_G == TRUE & POSITION < (chipTbl$Cstop[chipRow] + tailWidth) &
				END > (chipTbl$Cstart[chipRow] - tailWidth))
		# force exactly one row...
		best <- which.min( abs( chipTbl$Ccenter[chipRow] - (gmap$POSITION + gmap$END)/2))
		gmap <- as.data.frame( gmap[ best, ])
		gene <- gmap$GENE_ID[1]
	} else {
		gene <- chipTbl$GENE_ID[chipRow]
		gmap <- subset( geneMap, GENE_ID == gene)
	}
	seqid <- gmap$SEQ_ID
	if ( !(seqid %in% seqMap$SEQ_ID)) {
		cat( "\nGeneID and/or SeqID not in current Species.\n")
		return()
	}

	# turn the peak location and tails into the base limits
	minBase <- 1
	maxBase <- subset.data.frame( seqMap, SEQ_ID==seqid)$LENGTH
	leftBase <- chipTbl$Fstart[chipRow] - tailWidth
	if ( is.na( leftBase)) leftBase <- chipTbl$Cstart[chipRow] - tailWidth
	if ( is.na( leftBase)) leftBase <- chipTbl$Rstart[chipRow] - tailWidth*2
	rightBase <- chipTbl$Rstop[chipRow] + tailWidth
	if ( is.na( rightBase)) rightBase <- chipTbl$Cstop[chipRow] + tailWidth
	if ( is.na( rightBase)) rightBase <- chipTbl$Fstop[chipRow] + tailWidth*2
	if ( any( is.na( c( leftBase, rightBase)))) {
		cat( "\nNo non-NA peak edges.. Skipping")
		return()
	}
	if ( any( is.infinite( c( leftBase, rightBase)))) {
		cat( "\nNo non-NA peak edges.. Skipping")
		return()
	}
	if ( is.na( rightBase)) rightBase <- chipTbl$Fstop[chipRow] + tailWidth*2
	if ( leftBase < minBase) leftBase <- minBase
	if ( rightBase > maxBase) rightBase <- maxBase

	# scale based on the one peak of interest
	leftBaseScaling <- chipTbl$Cstart[chipRow]
	if ( is.na( leftBaseScaling)) leftBaseScaling <- leftBase
	rightBaseScaling <- chipTbl$Cstop[chipRow]
	if ( is.na( rightBaseScaling)) rightBaseScaling <- rightBase

	# get the actual data
	wigdata <- WIG_getWigglesOneSeq( WIG, seqid)
	hasData <- ( nrow(wigdata$Plus) > 0) || ( nrow(wigdata$Minus) > 0)
	totalReads <- WIG$Info$TotalReads
	scaleFactor <- 1
	if ( ! is.null( scaledReadCount)) scaleFactor <- scaledReadCount / totalReads

	if ( hasData) {
		wigP <- baseDepthTableSubset( wigdata$Plus, leftBase, rightBase)
		wigM <- baseDepthTableSubset( wigdata$Minus, leftBase, rightBase)
		wigC <- merge.baseDepthTables( wigP, wigM)
		yData <- list( "Plus"=wigP, "Minus"=wigM, "Combo"=wigC)
		if ( ! is.null( scaledReadCount)) yData <- yDataScale( yData, scaleFactor=scaleFactor)

		# pick sensible plot limits, given the data and user given limits
		bigY <- max( c( minYmax, baseDepthTableSubset( wigC, leftBaseScaling, 
				rightBaseScaling)$DEPTH * scaleFactor)) 

		# certain plots may need more room -- a summed plot may be up to 2x tall
		if ( regexpr( "added", plotType) > 0) bigY <- bigY * 2

	} else {  # no data to draw...
		bigY <- minYmax
	}
	if ( !is.null( forceYmax)) bigY <- forceYmax
	dx5 <- (rightBase - leftBase) * 0.05
	xlim <- c(leftBase+dx5, rightBase-dx5)

	if ( useLog) {
		useBigY <- log2( bigY)
		smallestY <- 0.25
	} else {
		useBigY <- bigY
		smallestY <- useBigY / 800
	}

	ylim <- c( 0, useBigY*1.1)

	# now stretch the lower end of the plot to fit the annotations
	ylim[1] <- ylim[1] - (( ylim[2]-ylim[1]) * 0.20)
	annoSize <- (ylim[2]-ylim[1]) * 0.04

	xlabel <- paste( "Chromosome: ", gmap$SEQ_ID[1], "     Nucleotide Location") 
	ylabel <- "ChIP Read Depth"
	if ( ! is.null(scaledReadCount)) ylabel <- "ChIP Read Depth after 'Read Count Scaling'"
	pltTxt <- paste( "ChIP-Seq:  ",label) 
	pltTxt <- paste( pltTxt, "\n", shortGeneName(gene), " -  ", gmap$PRODUCT)
	if ( useLog) ylabel <- paste( ylabel, "   (log scale)")

	if (presentationMode) {
		xlabel <- NA
		ylabel <- NA
		pltTxt <- shortGeneName( gene)
	}

	plot( x=xlim, y=c(0,0), main=pltTxt, type="l", xlim=xlim, ylim=ylim, xaxs="r",
		xlab=xlabel, ylab=ylabel, yaxt="n", font.axis=2, font.lab=2, cex.lab=cex.lab, cex.axis=cex.lab)


	if ( useLog) {
		ticks <- rep( c( 1,2,5), times=7) * rep( c( 1,10,100,1000,10000,100000,1000000), each=3)
		ticks <- c( 3,10,30,100,300,1000,3000,10000,30000,100000,3000000,1000000)
		useTicks <- ticks[ ticks < bigY]
		if ( length( useTicks) < 4) {
			useTicks <- ticks[1:4]
		} else {
			useTicks <- rev( useTicks)[ seq( 1, length(useTicks), by=2)]
		}
		at <- log2( useTicks)
		if ( hasData) {
			yData <- log2Scale( yData, uniqueToo=FALSE)

			# kyle says don't draw the combo peak
			if ( presentationMode) yData <- list( "Plus"=yData$Plus, "Minus"=yData$Minus)
			WIG_Draw( yData, plotType=plotType, lwd=3, smallestDepth=smallestY)
		}
		if ( ! presentationMode) {
			WIG_Ticks( xlim, at, labels=as.integer(useTicks))
		} else {
			axis( side=2, at=at, labels=as.integer(useTicks))
		}

	} else {
		ticks <- seq( from=0, to=1000000, length.out=5)[2:5]
		yBreaks <- c(  800000, 600000, 400000, 200000, 160000, 120000, 80000, 60000, 40000, 20000, 
				16000, 12000, 8000, 6000, 4000, 2000, 1600, 1200, 800, 600, 400, 200, 
				160, 120, 80, 60, 40, 20, 16, 12, 8, 4)
		for ( ycut in yBreaks) {
			if ( bigY < ycut) ticks <- seq( from=0, to=ycut, length.out=5)[2:5]
		}
		at <- useTicks <- ticks
		# draw the wanted style...
		if ( hasData) {
			if ( presentationMode) yData <- list( "Plus"=yData$Plus, "Minus"=yData$Minus)
			WIG_Draw( yData, plotType=plotType, lwd=3, smallestDepth=smallestY)
		}
		if( !presentationMode) {
			WIG_Ticks( xlim, ticks, labels=ticks)
		} else {
			axis( side=2, at=ticks, labels=ticks)
		}
	}

	annotateWIGplot( seqID=seqid, gene=gene, xlim=xlim, annoSize, speciesID=curSpecies, 
			altGeneMap=NULL, showDetectability=FALSE)

	if (! presentationMode && legend.cex > 0) {
		legend( x="topleft", legend=c("(+) Fwd Strand", " (-) Rev Strand", "Combined"), 
			col=c(4,2,1), text.col=c(4,2,1), lwd=c(3,3,3), bg="white", cex=legend.cex)
		legend( x="topright", legend=paste( c( "Peak Rank:        ", "Volume (VPM): ", "Score:     ", "P Value:  "), 
			c( chipRow, as.integer( chipTbl$VPM[chipRow]), 
			formatC(chipTbl$Score[chipRow], format="f", digits=3),
			formatC(chipTbl$P_Value[chipRow], format="f", digits=4))),
			bg="white", cex=legend.cex)
	}

	# draw all the models in the window
	toShow <- which( chipTbl$Cstart < rightBase & chipTbl$Cstop > leftBase)
	ltBlk <- adjustColor( 1, 0.55)
	ltBlue <- adjustColor( 4, 0.55)
	ltRed <- adjustColor( 2, 0.55)

	functSet <- list( gaussian, gumbel, lorentzian)
	functNames <- c( "gaussian", "gumbel", "lorentz")

	for ( ipk in toShow) {
		if ( ipk == chipRow) {
			LWD <- 3
			LTY <- 1
		} else {
			LWD <- 2
			LTY <- 3
		}
		if( ! any( is.na( chipTbl[ ipk, c( "Cstart", "Cstop")]))) {
			FUN <- functSet[[ match( chipTbl$Ctype[ipk], functNames) ]]
			comboX <- chipTbl$Cstart[ipk] : chipTbl$Cstop[ipk]
			comboY <- FUN( comboX, chipTbl$Ccenter[ipk], chipTbl$Cwidth[ipk], 
					chipTbl$Cheight[ipk], chipTbl$Cfloor[ipk])
			lines( comboX, comboY, col=ltBlk, lty=LTY, lwd=LWD)
		}

		if( ! any( is.na( chipTbl[ ipk, c( "Fstart", "Fstop")])) &&
			( all( is.finite( as.numeric(chipTbl[ ipk, c( "Fstart", "Fstop")]))))) {
			FUN <- functSet[[ match( chipTbl$Ftype[ipk], functNames) ]]
			plusX <- chipTbl$Fstart[ipk] : chipTbl$Fstop[ipk]
			plusY <- FUN( plusX, chipTbl$Fcenter[ipk], chipTbl$Fwidth[ipk], chipTbl$Fheight[ipk], 
					chipTbl$Ffloor[ipk])
			lines( plusX, plusY, col=ltBlue, lty=LTY, lwd=LWD)
			lines( c( rep( chipTbl$Fstart[ipk], 2), rep( chipTbl$Fstop[ipk], 2)), 
					c( plusY[1], rep( chipTbl$Ffloor[ipk],2), plusY[length(plusY)]), 
					col=ltBlue, lty=LTY, lwd=LWD)
		}

		if( ! any( is.na( chipTbl[ ipk, c( "Rstart", "Rstop")])) &&
			( all( is.finite( as.numeric(chipTbl[ ipk, c( "Rstart", "Rstop")]))))) {
			FUN <- functSet[[ match( chipTbl$Rtype[ipk], functNames) ]]
			minusX <- chipTbl$Rstart[ipk] : chipTbl$Rstop[ipk]
			minusY <- FUN( minusX, chipTbl$Rcenter[ipk], chipTbl$Rwidth[ipk], chipTbl$Rheight[ipk], 
					chipTbl$Rfloor[ipk])
			lines( minusX, minusY, col=ltRed, lty=LTY, lwd=LWD)
			lines( c( rep( chipTbl$Rstart[ipk], 2), rep( chipTbl$Rstop[ipk], 2)), 
					c( minusY[1], rep( chipTbl$Rfloor[ipk],2), minusY[length(minusY)]),
					col=ltRed, lty=LTY, lwd=LWD)
		}
	}

	return()
}


`plotRIPpeak` <- function( WIG, ripTbl, ripRow=1, tailWidth=1200, plotType=c("lines", "boxes"), 
			label="", useLog=FALSE, forceYmax=NULL, minYmax=10, legend.cex=1.15,
			scaledReadCount=NULL) {

	plotType <- match.arg( plotType)

	presentationMode <- FALSE
	cex.lab <- 1.0
	if (presentationMode) cex.lab=1.5


	# check species
	#if ( ! verifyWIGspecies( WIG, seqID=NULL, geneID=gene)) return()

	curSpecies <- getCurrentSpecies()
	seqMap <- getCurrentSeqMap()
	geneMap <- getCurrentGeneMap()
	exonMap <- getCurrentExonMap()

	# estimate a gene if we need
	if ( ! "GENE_ID" %in% colnames(ripTbl)) {
		gmap <- subset( geneMap, REAL_G == TRUE & POSITION < (ripTbl$Stop[ripRow] + tailWidth) &
				END > (ripTbl$Start[ripRow] - tailWidth))
		# force exactly one row...
		best <- which.min( abs( ripTbl$Center[ripRow] - (gmap$POSITION + gmap$END)/2))
		gmap <- as.data.frame( gmap[ best, ])
		gene <- gmap$GENE_ID[1]
	} else {
		gene <- ripTbl$GENE_ID[ripRow]
		gmap <- subset( geneMap, GENE_ID == gene)
	}
	seqid <- gmap$SEQ_ID
	if ( !(seqid %in% seqMap$SEQ_ID)) {
		cat( "\nGeneID and/or SeqID not in current Species.\n")
		return()
	}

	# turn the peak location and tails into the base limits
	minBase <- 1
	maxBase <- subset.data.frame( seqMap, SEQ_ID==seqid)$LENGTH
	leftBase <- ripTbl$Start[ripRow] - tailWidth
	if ( is.na( leftBase)) leftBase <- ripTbl$Center[ripRow] - tailWidth
	rightBase <- ripTbl$Stop[ripRow] + tailWidth
	if ( is.na( rightBase)) rightBase <- ripTbl$Center[ripRow] + tailWidth
	if ( any( is.na( c( leftBase, rightBase)))) {
		cat( "\nNo non-NA peak edges.. Skipping")
		return()
	}
	if ( any( is.infinite( c( leftBase, rightBase)))) {
		cat( "\nNo non-NA peak edges.. Skipping")
		return()
	}
	if ( leftBase < minBase) leftBase <- minBase
	if ( rightBase > maxBase) rightBase <- maxBase

	# scale based on the one peak of interest
	leftBaseScaling <- ripTbl$Start[ripRow]
	if ( is.na( leftBaseScaling)) leftBaseScaling <- leftBase
	rightBaseScaling <- ripTbl$Stop[ripRow]
	if ( is.na( rightBaseScaling)) rightBaseScaling <- rightBase

	# get the actual data
	wigdata <- WIG_getWigglesOneSeq( WIG, seqid)
	hasData <- ( nrow(wigdata$Plus) > 0) || ( nrow(wigdata$Minus) > 0)
	totalReads <- WIG$Info$TotalReads
	scaleFactor <- 1
	if ( ! is.null( scaledReadCount)) scaleFactor <- scaledReadCount / totalReads

	if ( hasData) {
		wigP <- baseDepthTableSubset( wigdata$Plus, leftBase, rightBase)
		wigM <- baseDepthTableSubset( wigdata$Minus, leftBase, rightBase)
		yData <- list( "Plus"=wigP, "Minus"=wigM)
		if ( ! is.null( scaledReadCount)) yData <- yDataScale( yData, scaleFactor=scaleFactor)

		# pick sensible plot limits, given the data and user given limits
		bigY <- max( c( minYmax, baseDepthTableSubset(wigP,leftBaseScaling,rightBaseScaling)$DEPTH * scaleFactor, 
				baseDepthTableSubset(wigM,leftBaseScaling,rightBaseScaling)$DEPTH * scaleFactor)) 
	} else {  # no data to draw...
		bigY <- minYmax
	}
	if ( !is.null( forceYmax)) bigY <- forceYmax
	dx5 <- (rightBase - leftBase) * 0.05
	xlim <- c(leftBase+dx5, rightBase-dx5)

	if ( useLog) {
		useBigY <- log2( bigY)
		smallestY <- 0.25
	} else {
		useBigY <- bigY
		smallestY <- useBigY / 800
	}

	ylim <- c( 0, useBigY*1.1)

	# now stretch the lower end of the plot to fit the annotations
	ylim[1] <- ylim[1] - (( ylim[2]-ylim[1]) * 0.20)
	annoSize <- (ylim[2]-ylim[1]) * 0.04

	xlabel <- paste( "Chromosome: ", gmap$SEQ_ID[1], "     Nucleotide Location") 
	ylabel <- "RIP Read Depth"
	if ( ! is.null(scaledReadCount)) ylabel <- "RIP Read Depth after 'Read Count Scaling'"
	pltTxt <- paste( "RIP-Seq:  ",label) 
	pltTxt <- paste( pltTxt, "\n", shortGeneName(gene), " -  ", gmap$PRODUCT)
	if ( useLog) ylabel <- paste( ylabel, "   (log scale)")

	if (presentationMode) {
		xlabel <- NA
		ylabel <- NA
		pltTxt <- shortGeneName( gene)
	}

	plot( x=xlim, y=c(0,0), main=pltTxt, type="l", xlim=xlim, ylim=ylim, xaxs="r",
		xlab=xlabel, ylab=ylabel, yaxt="n", font.axis=2, font.lab=2, cex.lab=cex.lab, cex.axis=cex.lab)


	if ( useLog) {
		ticks <- rep( c( 1,2,5), times=7) * rep( c( 1,10,100,1000,10000,100000,1000000), each=3)
		ticks <- c( 3,10,30,100,300,1000,3000,10000,30000,100000,3000000,1000000)
		useTicks <- ticks[ ticks < bigY]
		if ( length( useTicks) < 4) {
			useTicks <- ticks[1:4]
		} else {
			useTicks <- rev( useTicks)[ seq( 1, length(useTicks), by=2)]
		}
		at <- log2( useTicks)
		if ( hasData) {
			yData <- log2Scale( yData, uniqueToo=FALSE)

			# kyle says don't draw the combo peak
			if ( presentationMode) yData <- list( "Plus"=yData$Plus, "Minus"=yData$Minus)
			WIG_Draw( yData, plotType=plotType, lwd=3, smallestDepth=smallestY)
		}
		if ( ! presentationMode) {
			WIG_Ticks( xlim, at, labels=as.integer(useTicks))
		} else {
			axis( side=2, at=at, labels=as.integer(useTicks))
		}

	} else {
		ticks <- seq( from=0, to=1000000, length.out=5)[2:5]
		yBreaks <- c(  800000, 600000, 400000, 200000, 160000, 120000, 80000, 60000, 40000, 20000, 
				16000, 12000, 8000, 6000, 4000, 2000, 1600, 1200, 800, 600, 400, 200, 
				160, 120, 80, 60, 40, 20, 16, 12, 8, 4)
		for ( ycut in yBreaks) {
			if ( bigY < ycut) ticks <- seq( from=0, to=ycut, length.out=5)[2:5]
		}
		at <- useTicks <- ticks
		# draw the wanted style...
		if ( hasData) {
			if ( presentationMode) yData <- list( "Plus"=yData$Plus, "Minus"=yData$Minus)
			WIG_Draw( yData, plotType=plotType, lwd=3, smallestDepth=smallestY)
		}
		if( !presentationMode) {
			WIG_Ticks( xlim, ticks, labels=ticks)
		} else {
			axis( side=2, at=ticks, labels=ticks)
		}
	}

	annotateWIGplot( seqID=seqid, gene=gene, xlim=xlim, annoSize, speciesID=curSpecies, 
			altGeneMap=NULL, showDetectability=FALSE)

	if (! presentationMode && legend.cex > 0) {
		legend( x="topleft", legend=c("(+) Fwd Strand", " (-) Rev Strand"), 
			col=c(4,2), text.col=c(4,2), lwd=c(3,3), bg="white", cex=legend.cex)
		legend( x="topright", legend=paste( c( "Peak Rank:        ", "Volume (VPM): ", "Score:     ", "P Value:  "), 
			c( ripRow, as.integer( ripTbl$VPM[ripRow]), 
			formatC(ripTbl$Score[ripRow], format="f", digits=3),
			formatC(ripTbl$P_Value[ripRow], format="f", digits=4))),
			bg="white", cex=legend.cex)
	}

	# draw all the models in the window
	toShow <- which( ripTbl$Start < rightBase & ripTbl$Stop > leftBase)
	ltBlue <- adjustColor( 4, 0.55)
	ltRed <- adjustColor( 2, 0.55)

	functSet <- list( gaussian, gumbel, lorentzian, pulse)
	functNames <- c( "gaussian", "gumbel", "lorentz", "pulse")

	for ( ipk in toShow) {
		if ( ipk == ripRow) {
			LWD <- 3
			LTY <- 1
		} else {
			LWD <- 2
			LTY <- 3
		}
		myCol <- if ( ripTbl$Strand[ipk] == "Plus") ltBlue else ltRed

		if( ! any( is.na( ripTbl[ ipk, c( "Start", "Stop")]))) {
			FUN <- functSet[[ match( ripTbl$Type[ipk], functNames) ]]
			# extend X a bit to get PULSE legs to show
			comboX <- (ripTbl$Start[ipk] - 25) : (ripTbl$Stop[ipk] + 25)
			comboY <- FUN( comboX, ripTbl$Center[ipk], ripTbl$Width[ipk], 
					ripTbl$Height[ipk], ripTbl$Floor[ipk])
			lines( comboX, comboY, col=myCol, lty=LTY, lwd=LWD)
		}
	}

	return()
}

