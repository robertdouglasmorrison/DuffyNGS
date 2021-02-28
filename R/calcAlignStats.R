# calcAlignStats.R

`calcAlignStats` <- function( filein, sampleID, alignPhase=c("Genomic","RiboClear","Splice"), 
			what="SGBIDMAN", statsPath="AlignStats", reload=TRUE, totalReads=NULL, 
			chunkSize=200000, maxReads=NULL, plot=TRUE, pause=0, banner="",
			verbose=TRUE) {

	alignPhase <- match.arg( alignPhase)
	if ( is.null( what)) what <- "SGBIDMAN"

	# file and path for stats results and plots
	statsFile <- paste( sampleID, alignPhase, "alignStatsData.rda", sep=".")
	if ( ! file.exists( statsPath)) dir.create( statsPath, recursive=TRUE, showWarnings=TRUE)
	statsFile <- file.path( statsPath, statsFile)

	# when doing a subset of the data, adjust what we report
	if ( ! is.null( maxReads)) totalReads <- maxReads

	# set up all results for by-buffer progress
	checkX11( width=9, height=7, xpos=20, ypos=20, bg='white');

	# allowing this file to append to existing or start new results
	if ( reload) {
		ansS <- ansG <- data.frame()
		ansB <- ansI <- ansM <- ansA <- ansN <- NULL
		nReads <- nAligns <- 0
		ReadLength <- NULL
		fileset <- vector()
		fileset[1] <- filein
	} else {
		load( statsFile)
		ansS <- alignStatsData$seqStatistics
		ansG <- alignStatsData$geneStatistics
		ansM <- alignStatsData$marStatistics
		ansA <- alignStatsData$alignStatistics
		ansB <- alignStatsData$baseStatistics
		ansI <- alignStatsData$indelStatistics
		ansN <- alignStatsData$insertStatistics
		nReads <- alignStatsData$nReads
		nAligns <- alignStatsData$nAligns
		ReadLength <- alignStatsData$ReadLength
		fileset <- alignStatsData$filename
		fileset <- c( fileset, filein)
	}

	# open that BAM file
	reader <- bamReader( filein)
	refData <- getRefData(reader)

	hasMore <- TRUE
	while ( hasMore) {
		if ( ! is.null( maxReads) && nReads >= maxReads) break
		gc()
		cat( "\nReadBAM..")
		chunk <- getNextChunk( reader, n=chunkSize, alignedOnly=T, primaryOnly=F)
		nNow <- size( chunk)
		if ( nNow < 1) break
		if ( nNow < chunkSize) hasMore <- FALSE

		isPrimary <- which( !secondaryAlign(chunk))
		nNowReads <- length( isPrimary)
		nAligns <- nAligns + nNow;
		nReads <- nReads + nNowReads;
		cat( "  N_Aligns:", prettyNum( as.integer( nAligns), big.mark=","))
		cat( "  N_Reads:", prettyNum( as.integer( nReads), big.mark=","))

		# get all the fields we need/want from the BAM data and then release that memory
		refIDs <- refID( chunk)
		if ( alignPhase != "Splice") {
			geneIDs <- getTag( chunk, BAMTAG_GENEID)
		}
		# see how long the reads are if we don't know yet
		if ( is.null( ReadLength)) {
			seqs <- alignSeq(chunk)
			ReadLength <- round( mean( nchar( seqs)))
			rm( seqs)
		}
		readFac <- factor.align( chunk)
		scores <- as.integer( getTag( chunk, tag=BAMTAG_ALIGNSCORE))
		mismatchs <- NULL
		mismatchs <- mismatchData( chunk, readUnits=TRUE)
		inserts <- insertSize(chunk)

		# at this point, the 'chunk' is done...
		rm( chunk)
		gc()

		# OK, we know exactly which reads/alignments we will use
		mismatchs <- subset.data.frame( mismatchs, Align %in% isPrimary)

		# most types have the SeqID as is and the gene as a tags...
		seqIDs <- refID2seqID( refIDs, refData=refData)

		# but the splices have it embedded in the seqID
		if ( alignPhase == "Splice") {
			cat( " parseSplice..")
			terms <- strsplit( seqIDs, split="::")
			seqIDs <- sapply( terms, function(x) x[1])
			geneIDs <- sapply( terms, function(x) x[2])
		}


		# table of seqIDs
		if ( regexpr( "S", what) > 0) {
			cat( " seqs..")
			thisS <- evaluateSeqStats( seqIDs[isPrimary])
			ansS <- mergeSeqStats( rbind( ansS, thisS))
		}

		# gene-level stats...
		if ( regexpr( "G", what) > 0) {
			cat( " genes..")
			thisG <- evaluateGeneStats( seqIDs[isPrimary], geneIDs[isPrimary], verbose=verbose)
			ansG <- mergeGeneStats( rbind( ansG, thisG))
		}
	
		# MAR stats...
		if ( regexpr( "M", what) > 0) {
			cat( " MARs..")
			thisM <- evaluateMARstats( seqIDs, geneIDs, readFac, verbose=verbose)
			ansM <- mergeMARstats( ansM, thisM)
		}
	
		# ALIGN stats...
		if ( regexpr( "A", what) > 0) {
			cat( " AlignScores..")
			thisA <- evaluateALIGNstats( scores[isPrimary])
			ansA <- mergeALIGNstats( ansA, thisA)
		}
	
		# base flipping
		if ( regexpr( "B", what) > 0) {
			cat( " bases..")
			thisB <- evaluateBaseFlips( mismatchs, readLength=ReadLength, nAlign=nNow)
			ansB <- mergeBaseStats( ansB, thisB)
		}

		# indels
		if ( regexpr( "I|D", what) > 0) {
			cat( " indels..")
			thisI <- evaluateIndels( mismatchs, readLength=ReadLength, nAlign=nNow)
			ansI <- mergeIndelStats( ansI, thisI)
		}

		#after each buffer, show what we have...
		cat( " save..")
		alignStatsData <- list( "sampleID"= sampleID, "filename"=fileset, "nReads"=nReads, 
					"nAligns"=nAligns, "ReadLength"=ReadLength, "seqStatistics"=ansS, 
					"geneStatistics"=ansG, "marStatistics"=ansM, "baseStatistics"=ansB, 
					"indelStatistics"=ansI, "alignStatistics"=ansA, "banner"=banner, 
					"alignPhase"=alignPhase, "totalReads"=totalReads)
		save( alignStatsData, file=statsFile)

		if (plot) {
			cat( " plot.")
			plotAlignStats( alignStatsData, asPNG=TRUE, plotPath=statsPath, pause=pause, verbose=verbose)
		}

		rm( alignStatsData, mismatchs, seqIDs, refIDs, geneIDs)
		gc()
	} 

	cat( "\nDone.   Final Plots..\n")
	load( statsFile)
	plotAlignStats( alignStatsData, asPNG=TRUE, plotPath=statsPath, pause=pause, verbose=verbose)

	# close that graphics window
	dev.off()

	bamClose( reader)
	
	return()
}


evaluateIndels <- function( mismatchs, readLength=36, nAlign=max(mismatchs$Align)) {

	# we have a data frame with 'everything' we need...

	nIns <- nDel <- 0
	insPositionCnts <- delPositionCnts <- rep( 0, times=readLength)
	insPatternCnts <- delPatternCnts <- vector()
	nZeroInserts <- nZeroDeletes <- nAlign

	myTableFun <- base::table

	# look at only the insertions and deletions
	mis <- subset( mismatchs, Type == "I")
	if ( nrow(mis) > 0) {
		alignFac <- factor( mis$Align)
		nIns <- tapply( 1:nrow(mis), alignFac, length)
		posCnts <- myTableFun( mis$Position)
		posLocs <- as.integer( names( posCnts))
		where <- match( posLocs, 1:readLength, nomatch=0)
		insPositionCnts[ where] <- posCnts[ where > 0] 
		insPatternCnts <- myTableFun( mis$ReadBase)
		# we have to calculate how many had 'zero' events
		nZeroInserts <- nAlign - nlevels(alignFac)
		nIns <- c( nIns, rep( 0, nZeroInserts))
	}

	mis <- subset( mismatchs, Type == "D")
	if ( nrow(mis) > 0) {
		alignFac <- factor( mis$Align)
		nDel <- tapply( 1:nrow(mis), alignFac, length)
		posCnts <- myTableFun( mis$Position)
		posLocs <- as.integer( names( posCnts))
		where <- match( posLocs, 1:readLength, nomatch=0)
		delPositionCnts[ where] <- posCnts[ where > 0] 
		delPatternCnts <- myTableFun( mis$RefBase)
		# we have to calculate how many had 'zero' events
		nZeroDeletes <- nAlign - nlevels(alignFac)
		nDel <- c( nDel, rep( 0, nZeroDeletes))
	}

	avgNins <- mean.default( nIns)
	avgNdel <- mean.default( nDel)
	insCounts <- base::sort( myTableFun( nIns), decreasing=TRUE)
	delCounts <- base::sort( myTableFun( nDel), decreasing=TRUE)

	out <- list( "AvgN_Insert"=avgNins, "InsertFrequency"=insCounts, 
			"InsertPositions"=insPositionCnts, "InsertPatterns"=insPatternCnts,
			 "AvgN_Delete"=avgNdel, "DeleteFrequency"=delCounts, 
			"DeletePositions"=delPositionCnts, "DeletePatterns"=delPatternCnts)
	return( out)
}


mergeIndelStats <- function( stats1, stats2) {

	# given two baseStats results, re-combine into one.
	if ( is.null( stats1)) return( stats2)

	# positions are how many at each spot, (the lengths may be unevern
	len1 <- length( stats1$InsertPositions)
	len2 <- length( stats2$InsertPositions)
	lenNew <- max( c( len1, len2))
	newPos <- rep( 0, times=lenNew)
	if (len1 > 0) newPos[ 1:len1] <- newPos[ 1:len1] + stats1$InsertPositions
	if (len2 > 0) newPos[ 1:len2] <- newPos[ 1:len2] + stats2$InsertPositions
	newInsertPos <- newPos
	len1 <- length( stats1$DeletePositions)
	len2 <- length( stats2$DeletePositions)
	lenNew <- max( c( len1, len2))
	newPos <- rep( 0, times=lenNew)
	if (len1 > 0) newPos[ 1:len1] <- newPos[ 1:len1] + stats1$DeletePositions
	if (len2 > 0) newPos[ 1:len2] <- newPos[ 1:len2] + stats2$DeletePositions
	newDeletePos <- newPos
	
	# the patterns of bases
	newInsertPat <- mergeTables( stats1$InsertPatterns, stats2$InsertPatterns)
	newDeletePat <- mergeTables( stats1$DeletePatterns, stats2$DeletePatterns)

	# the frequencies are harder...
	allNames <- c( names( stats1$InsertFrequency), names( stats2$InsertFrequency))
	allNs <- c( stats1$InsertFrequency, stats2$InsertFrequency)
	f <- factor( allNames)
	newNs <- tapply( allNs, f, sum)
	newInsertCnts <- newNs <- base::sort( newNs, decreasing=TRUE)
	newNames <- names( newNs)
	sum <- sum( newNs)
	total <- sum(as.integer( newNames) * newNs)
	newInsertFraction <- total / sum
	allNames <- c( names( stats1$DeleteFrequency), names( stats2$DeleteFrequency))
	allNs <- c( stats1$DeleteFrequency, stats2$DeleteFrequency)
	f <- factor( allNames)
	newNs <- tapply( allNs, f, sum)
	newDeleteCnts <- newNs <- base::sort( newNs, decreasing=TRUE)
	newNames <- names( newNs)
	sum <- sum( newNs)
	total <- sum(as.integer( newNames) * newNs)
	newDeleteFraction <- total / sum

	out <- list( "AvgN_Insert"=newInsertFraction, "InsertFrequency"=newInsertCnts, 
			"InsertPositions"=newInsertPos, "InsertPatterns"=newInsertPat,
			 "AvgN_Delete"=newDeleteFraction, "DeleteFrequency"=newDeleteCnts, 
			"DeletePositions"=newDeletePos, "DeletePatterns"=newDeletePat)

	return( out)
}


evaluateBaseFlips <- function( mismatchs, readLength=36, nAlign=max(mismatchs$Align)) {

	# we have a data frame with 'everything' we need...

	nTotal <- nFlips <- 0
	positionCnts <- rep( 0, times=readLength)
	letterCnts <- rep(0, times=5)
	names( letterCnts) <- BASES <- c( "A","C","G","T","N")

	myTableFun <- base::table

	# look at only the base flips -- i.e. true mismatches
	mis <- subset( mismatchs, Type == "X")
	if ( nrow(mis) > 0) {
		alignFac <- factor( mis$Align)
		nFlips <- tapply( 1:nrow(mis), alignFac, length)

		posCnts <- myTableFun( mis$Position)
		posLocs <- as.integer( names( posCnts))
		where <- match( posLocs, 1:readLength, nomatch=0)
		positionCnts[ where] <- posCnts[ where > 0] 

		letCnts <- myTableFun( mis$ReadBase)
		letLocs <- names( letCnts)
		where <- match( letLocs, BASES, nomatch=0)
		letterCnts[ where] <- letCnts[ where > 0] 

		# we have to calculate how many had 'zero' events
		nZeroFlips <- nAlign - nlevels(alignFac)
		nFlips <- c( nFlips, rep( 0, nZeroFlips))
	}

	avgNflips <- mean.default( nFlips)
	flipCounts <- base::sort( myTableFun( nFlips), decreasing=TRUE)
	
	out <- list( "AvgN_Flips"=avgNflips, "Frequency"=flipCounts, "Positions"=positionCnts, 
			"Letters"=letterCnts)
	return( out)
}

mergeBaseStats <- function( stats1, stats2) {

	# given two baseStats results, re-combine into one.
	if ( is.null( stats1)) return( stats2)

	# positions are how many at each spot, (the lengths may be unevern
	len1 <- length( stats1$Positions)
	len2 <- length( stats2$Positions)
	lenNew <- max( c( len1, len2))
	newPos <- rep( 0, times=lenNew)
	newPos[ 1:len1] <- newPos[ 1:len1] + stats1$Positions
	newPos[ 1:len2] <- newPos[ 1:len2] + stats2$Positions
	
	# same for scores...and letters
	newLetters <- stats1$Letters + stats2$Letters

	# the frequencies are harder...
	allNames <- c( names( stats1$Frequency), names( stats2$Frequency))
	allNs <- c( stats1$Frequency, stats2$Frequency)
	f <- factor( allNames)
	newNs <- tapply( allNs, f, sum)
	newNs <- base::sort( newNs, decreasing=TRUE)
	newNames <- names( newNs)
	sum <- sum( newNs)
	total <- sum(as.integer( newNames) * newNs)
	newFraction <- total / sum

	out <- list( "AvgN_Flips"=newFraction, "Frequency"=newNs, "Positions"=newPos, 
			"Letters"=newLetters)

	return( out)
}


evaluateSeqStats <- function( seqids) {

	tb <- base::sort( base::table( seqids), decreasing=TRUE)
	nams <- names( tb)
	cnts <- as.vector( tb)
	seqidPcts <- 100 * cnts / length( seqids)
	seqidCumPcts <- cumsum( seqidPcts)
	outdf <- data.frame( nams, cnts, seqidPcts, seqidCumPcts, stringsAsFactors=FALSE)
	colnames( outdf) <- c( "SEQ_ID", "N_READS", "PCT_TOTAL", "CUMMULATIVE_PCT")
	rownames( outdf) <- 1:nrow(outdf)
	return( outdf)
}

mergeSeqStats <- function( mydf) {

	# given one dataframe with two sets of Seq stats, combine them
	f <- factor( mydf$SEQ_ID)
	cnts <- tapply( mydf$N_READS, f, sum)
	totalCounts <- sum(cnts)
	newdf <- data.frame( levels(f), cnts, stringsAsFactors=FALSE)
	colnames( newdf) <- c( "SEQ_ID", "N_READS")
	ord <- base::order( newdf$N_READS, decreasing=TRUE)
	newdf <- newdf[ ord, ]
	pct <- 100 * newdf$N_READS / totalCounts
	cpct <- cumsum( pct)
	newdf$PCT_TOTAL <- pct
	newdf$CUMMULATIVE_PCT <- cpct
	rownames( newdf) <- 1:nrow(newdf)
	return( newdf)
}


evaluateGeneStats <- function( seqids, geneids, verbose=TRUE) {

	tb <- base::sort( base::table( geneids), decreasing=TRUE)
	genenames <- names(tb)
	tbv <- as.vector( tb)
	pct <- 100 * tbv / length( geneids)
	cpct <- cumsum( pct)
	seqptr <- base::match( genenames, geneids)
	seqnames <- seqids[seqptr]
	gProd <- gene2ProductAllSpecies( genenames)
	species <- getSpeciesFromSeqID( seqnames, verbose=verbose)
	species[ is.na(species)] <- "unknownSpecies"

	outdf <- data.frame( genenames, seqnames, species, gProd, tbv, pct, cpct, 
				stringsAsFactors=FALSE)
	colnames( outdf) <- c( "GENE_ID", "SEQ_ID", "SPECIES_ID", "PRODUCT", "N_READS", 
				"PCT_TOTAL", "CUMMULATIVE_PCT")
	rownames( outdf) <- 1:nrow(outdf)
				
	return( outdf)
}


mergeGeneStats <- function( mydf) {

	# given one dataframe with two sets of gene stats, combine them
	f <- factor( mydf$GENE_ID)
	cnts <- tapply( mydf$N_READS, f, sum)
	totalCounts <- sum(cnts)
	ptr <- base::match( levels(f), mydf$GENE_ID)
	seqnames <- mydf$SEQ_ID[ ptr]
	species <- mydf$SPECIES_ID[ ptr]
	gProd <- mydf$PRODUCT[ ptr]

	newdf <- data.frame( levels(f), seqnames, species, gProd, cnts, stringsAsFactors=FALSE)
	colnames( newdf) <- c( "GENE_ID", "SEQ_ID", "SPECIES_ID", "PRODUCT", "N_READS")
	ord <- base::order( newdf$N_READS, decreasing=TRUE)
	newdf <- newdf[ ord, ]
	pct <- 100 * newdf$N_READS / totalCounts
	cpct <- cumsum( pct)
	newdf$PCT_TOTAL <- pct
	newdf$CUMMULATIVE_PCT <- cpct
	rownames( newdf) <- 1:nrow(newdf)
	return( newdf)
}


evaluateMARstats <- function( seqids, geneids, alignFac, verbose=TRUE) {

	# using the read factoring, let's know something about species/non-gene/gene breakdown
	# of the Multi Hit Reads
	species <- getSpeciesFromSeqID( seqids, verbose=verbose)
	species[ is.na(species)] <- "unknownSpecies"
	N <- nlevels(alignFac)
	spec <- gtype <- cnt <- vector( length=N)
	Nnow <- 0
	mySort <- base::sort
	myPaste <- base::paste

	tapply( 1:length(geneids), INDEX=alignFac, FUN=function(x) {
			if ( length(x) < 2) return()
			Nnow <<- Nnow + 1
			specSet <- mySort( unique.default( species[x]))
			spec[Nnow] <<- if ( length(specSet) > 1) myPaste( specSet, collapse="+") else specSet[1]
			geneSet <- unique.default( geneids[x])
			nng <- length( grep( "(ng)", geneSet, fixed=T))
			ng <- length( geneSet) - nng
			gtype[Nnow] <<- if ( ng == 1 && nng == 0) "Single Gene"
					else if ( ng > 1 && nng == 0) "Multiple Genes"
					else if ( nng > 0 && ng == 0) "Intergenic Only"
					else "Genes + Intergenic"
			cnt[Nnow] <<- length(x)
			return()
		})

	length(spec) <- length(gtype) <- length(cnt) <- Nnow
	if ( Nnow > 0) {
		geneTbl <- base::table( myPaste( spec, gtype, sep="::"))
		cntTbl <- base::table( cnt)
		return( list( "Genes"=geneTbl, "Counts"=cntTbl))
	}
	return( NULL)
}


mergeMARstats <- function( stats1, stats2) {

	# given two marStats results, re-combine into one.
	if ( is.null( stats1)) return( stats2)

	geneTbl <- mergeTables( stats1$Genes, stats2$Genes)
	cntTbl <- mergeTables( stats1$Counts, stats2$Counts)
	return( list( "Genes"=geneTbl, "Counts"=cntTbl))
}


evaluateALIGNstats <- function( scores) {

	cntTbl <- base::table( as.integer( scores))
	return( cntTbl)
}


mergeALIGNstats <- function( stats1, stats2) {

	# given two alignStats results, re-combine into one.
	if ( is.null( stats1)) return( stats2)

	cntTbl <- mergeTables( stats1, stats2)
	return( cntTbl)
}


plotAlignStats <- function( alignStatsData=NULL, statsFile=NULL, asPNG=FALSE, plotPath=".", pause=0,
				verbose=TRUE) {

	# wrapper to plot the various stats from an .bam file
	if ( !(capabilities()[ "png" ])) asPNG <- FALSE

	if ( ! is.null( statsFile)) {
		who <- load( file=statsFile)
		alignStatsData <- get( who[1])
	}
	if ( !is.null( alignStatsData$seqStatistics) && nrow(alignStatsData$seqStatistics) > 0) {
	   plotSeqStatData( alignStatsData, asPNG=asPNG, plotPath=plotPath, pause=pause, verbose=verbose)
	}
	if ( !is.null( alignStatsData$geneStatistics) && nrow(alignStatsData$geneStatistics) > 0) {
	   plotGeneStatData( alignStatsData, asPNG=asPNG, plotPath=plotPath, pause=pause)
	}
	if ( ! is.null( alignStatsData$marStatistics)) {
	   plotMARstatData( alignStatsData, asPNG=asPNG, plotPath=plotPath, pause=pause)
	}
	if ( ! is.null( alignStatsData$baseStatistics)) {
		plotBaseStatData( alignStatsData, asPNG=asPNG, plotPath=plotPath, pause=pause)
	}
	if ( ! is.null( alignStatsData$indelStatistics)) {
		plotIndelStatData( alignStatsData, asPNG=asPNG, plotPath=plotPath, pause=pause)
	}
	if ( ! is.null( alignStatsData$alignStatistics)) {
		plotAlignScoreData( alignStatsData, asPNG=asPNG, plotPath=plotPath, pause=pause)
	}
}


plotSeqStatData <- function( alignStatsData, asPNG=FALSE, plotPath=".", pause=0, verbose=TRUE) {

	# plot a few types...
	saveMAI <- par("mai")
	sampleID <- alignStatsData$sampleID
	banner <- alignStatsData$banner
	alignPhase <- alignStatsData$alignPhase
	nReads <- alignStatsData$nReads
	totalReads <- alignStatsData$totalReads

	seqStats <- alignStatsData$seqStatistics
	freqtbl <- seqStats
	species <- getSpeciesFromSeqID( freqtbl$SEQ_ID, verbose=verbose)
	species[ is.na(species)] <- "unknownSpecies"
	colorAns <- colorBySpecies( species)
	colset <- colorAns$colors
	legendColors <- colorAns$legend.colors
	legendNames <- names( legendColors)

	subText <- NULL
	if ( alignPhase == "RiboClear") {
		yLabel <- "Percent of Ribo Cleared Reads"
		if ( ! is.null( totalReads)) subText <-  paste( "Ribo Clearing removed ", 
			as.percent( nReads, big.value=totalReads), " of all reads.")
	} else if (alignPhase == "Splice") {
		yLabel <- "Percent of Splice Aligned Reads"
		if ( ! is.null( totalReads)) subText <-  paste( "Splice Alignment assigned ", 
			as.percent( nReads, big.value=totalReads), " of all reads.")
	} else {
		yLabel <- "Percent of Genomic Aligned Reads"
		if ( ! is.null( totalReads)) subText <-  paste( "Genomic Alignment assigned ", 
			as.percent( nReads, big.value=totalReads), " of all reads.")
	}

	# barplot of chromosome percentages
	par("mai"=c( 2.6, 0.92, 0.92, 0.3))
	N <- 20
	if ( nrow( freqtbl) < N) N <- nrow(freqtbl)
	x <- freqtbl$PCT_TOTAL[1:N]
	names(x) <- freqtbl$SEQ_ID[1:N]
	mainText <- paste( sampleID, "    ", banner, "\n Top", N, "Chromosome Hits     ", "\nN_Reads:",
			prettyNum( as.integer( alignStatsData$nReads), big.mark=","))
	barplot( x, las=3, col=colset[1:N], main=mainText, ylab=yLabel, sub=NULL, font.lab=2, font.axis=2)
	legend( "topright", legendNames, fill=legendColors, cex=1.05, bg="white")
	if ( asPNG) {
		pngfile <- paste( sampleID, alignPhase, "Seq.BarPlot.png", sep=".")
		dev.print( png, file=file.path( plotPath, pngfile), width=800, height=600, bg='white')
	} 
	if (pause > 0) Sys.sleep(pause)


	# pie of chromosomes
	par("mai"=c( 1.02, 0.52, 0.92, 0.5))
	mylabels <- paste( freqtbl$SEQ_ID[1:N], "  Pct= ", formatC( freqtbl$PCT_TOTAL[1:N], digits=2, 
			format="f")," %", sep="")
	pie( freqtbl$PCT_TOTAL[1:N], labels=mylabels, col=colset[1:N], radius=0.9, clockwise=FALSE, 
		main=mainText, cex.sub=1.2, font.sub=2, sub=subText)
	legend( "topleft", legendNames, fill=legendColors, cex=1.05, bg="white")
	if ( asPNG) {
		pngfile <- paste( sampleID, alignPhase, "Seq.Pie.png", sep=".")
		dev.print( png, file=file.path( plotPath, pngfile), width=800, height=600, bg='white')
	} 
	if (pause > 0) Sys.sleep(pause)

	# pie of organisms
	f <- factor( species)
	pcts <- tapply( freqtbl$PCT_TOTAL, INDEX=f, FUN=sum)
	lbls <- paste( levels(f), "  Pct= ", formatC( pcts, digits=2, format="f"),"%", sep="")
	colset2 <- legendColors
	ord <- base::order( pcts, decreasing=TRUE)
	mainText <- paste( sampleID, "      ", banner, "\nSpecies Hits     ", "\nN_Reads:",
			prettyNum( as.integer( alignStatsData$nReads), big.mark=","))
	pie( pcts[ord], labels=lbls[ord], radius=0.9, col=colset2[ord], clockwise=FALSE, main=mainText, 
			cex.sub=1.2, font.sub=2, sub=subText)
	legend( "topleft", legendNames, fill=legendColors, cex=1.05, bg="white")
	if ( asPNG) {
		pngfile <- paste( sampleID, alignPhase, "Species.Pie.png", sep=".")
		dev.print( png, file=file.path( plotPath, pngfile), width=800, height=600, bg='white')
	} 
	if (pause > 0) Sys.sleep(pause)

	par( "mai"=saveMAI)
	return()
}


plotGeneStatData <- function( alignStatsData, asPNG=FALSE, plotPath=".", pause=0) {

	# plot a few types...
	saveMAI <- par("mai")
	sampleID <- alignStatsData$sampleID
	banner <- alignStatsData$banner
	alignPhase <- alignStatsData$alignPhase
	nReads <- alignStatsData$nReads
	totalReads <- alignStatsData$totalReads

	# gather the detail and maps that we need...
	geneStats <- alignStatsData$geneStatistics
	freqtbl <- geneStats

	#grp <- gene2ProductAllSpecies( freqtbl$GENE_ID)
	grp <- freqtbl$PRODUCT
	grp[ grep( "(ng)", freqtbl$GENE_ID, fixed=TRUE) ] <- "intergenic"
	grp[ grp == ""] <- "no_product"
	groupFactor <- factor( grp)
	
	species <- freqtbl$SPECIES_ID
	colorAns <- colorBySpecies( species, intergenic=(grp == "intergenic"))
	colset <- colorAns$colors
	legendColors <- colorAns$legend.colors
	legendNames <- names( legendColors)

	# barplot of gene percentages
	par("mai"=c( 2.6, 0.92, 0.92, 0.3))
	N <- 20
	if ( nrow( freqtbl) < N) N <- nrow(freqtbl)
	x <- freqtbl$PCT_TOTAL[1:N]
	names(x) <- shortGeneName( freqtbl$GENE_ID[1:N])

	subText <- ""
	if ( !is.null(totalReads)) {
		subText <- paste( "Top", N, "Genes account for ", 
			as.percent( sum( freqtbl$N_READS[1:N]), big.value=totalReads), "of all reads.")
	}
	if ( alignPhase == "RiboClear") {
		yLabel <- "Percent of Ribo Cleared Reads"
		if ( ! is.null( totalReads)) subText <- paste( subText, "\n", "Ribo Clearing removed ", 
			as.percent( nReads, big.value=totalReads), " of all reads.")
	} else if (alignPhase == "Splice") {
		yLabel <- "Percent of Splice Aligned Reads"
		if ( ! is.null( totalReads)) subText <- paste( subText, "\n", "Splice Alignment assigned ", 
			as.percent( nReads, big.value=totalReads), " of all reads.")
	} else {
		yLabel <- "Percent of All Aligned Reads"
		if ( ! is.null( totalReads)) subText <- paste( subText, "\n", "Genomic Alignment assigned ", 
			as.percent( nReads, big.value=totalReads), " of all reads.")
	}

	mainText <- paste( sampleID, "    ", banner, "\n Top", N, "Gene Hits     ", "\nN_Reads:",
			prettyNum( as.integer( alignStatsData$nReads), big.mark=","))
	barplot( x, las=3, col=colset[1:N], main=mainText, ylab=yLabel, sub=NULL, font.lab=2, font.axis=2)
	legend( "topright", legendNames, fill=legendColors, cex=1.05, bg="white")
	if ( asPNG) {
		pngfile <- paste( sampleID, alignPhase, "Gene.BarPlot.png", sep=".")
		dev.print( png, file=file.path( plotPath, pngfile), width=800, height=600, bg='white')
	} 
	if (pause > 0) Sys.sleep(pause)

	# pie of genes
	par("mai"=c( 1.02, 0.52, 0.92, 0.5))
	mylabels <- paste( shortGeneName( freqtbl$GENE_ID[1:N]), "  Pct= ", 
				formatC( freqtbl$PCT_TOTAL[1:N], digits=2, format="f")," %", sep="")
	pie( freqtbl$PCT_TOTAL[1:N], labels=mylabels, col=colset[1:N], radius=0.8, clockwise=FALSE, main=mainText,
			cex.sub=1.2, font.sub=2, sub=subText, cex=0.9)
	legend( "topleft", legendNames, fill=legendColors, cex=1.05, bg="white")
	#if ( !is.null( subText)) text( 0.5, 0.1, subText, cex=1.2)
	if ( asPNG) {
		pngfile <- paste( sampleID, alignPhase, "Gene.Pie.png", sep=".")
		dev.print( png, file=file.path( plotPath, pngfile), width=800, height=600, bg='white')
	} 
	if (pause > 0) Sys.sleep(pause)


	# pie of products
	N <- 20
	mainText <- paste( sampleID, "    ", banner, "\n Top", N, " 'Gene Product' Hits     ", "\nN_Reads:",
			prettyNum( as.integer( alignStatsData$nReads), big.mark=","))
	pcts <- tapply( freqtbl$PCT_TOTAL, INDEX=groupFactor, FUN=sum)
	lbls <- paste( levels(groupFactor), "  Pct= ", formatC( pcts, digits=2, format="f"),"%", sep="")
	ord <- base::order( pcts, decreasing=TRUE)
	if ( nlevels( groupFactor) < N) N <- nlevels(groupFactor)
	pie( pcts[ord[1:N]], labels=lbls[ord[1:N]], radius=0.8, clockwise=FALSE, main=mainText,
			cex.sub=1.2, font.sub=2, sub=subText, cex=0.8)
	if ( asPNG) {
		pngfile <- paste( sampleID, alignPhase, "GeneProduct.Pie.png", sep=".")
		dev.print( png, file=file.path( plotPath, pngfile), width=800, height=600, bg='white')
	} 
	if (pause > 0) Sys.sleep(pause)

	par( "mai"=saveMAI)
	return()
}


plotBaseStatData <- function( alignStatsData, asPNG=FALSE, plotPath=".", pause=0) {

	# plot a few types...
	saveMAI <- par("mai")
	sampleID <- alignStatsData$sampleID
	banner <- alignStatsData$banner
	alignPhase <- alignStatsData$alignPhase
	nReads <- alignStatsData$nReads
	totalReads <- alignStatsData$totalReads

	baseStats <- alignStatsData$baseStatistics

	hasScores <- FALSE

	# number of bases being flipped
	par("mai"=c( 1.02, 0.92, 0.92, 0.5))
	avg <- baseStats$AvgN_Flips
	flipCnts <- baseStats$Frequency
	total <- sum( flipCnts, na.rm=TRUE)
	flipPcts <- 100 * flipCnts / total
	ord <- base::order( as.numeric( names( flipPcts)))
	flipPcts <- flipPcts[ ord]

	#if ( length( flipPcts) > 8) length(flipPcts) <- 8
	# find how many are 'worth' plotting...
	bigX <- max( which( flipPcts >= 0.1))
	length( flipPcts) <- bigX
	if (bigX < 5) bigX <- 5
	xLim <- ceiling( bigX * 1.3)

	yLim <- c( 0, max(flipPcts) * 1.2)
	
	yLabel <- "Percent of All Aligned Reads"
	if ( alignPhase == "RiboClear") {
		yLabel <- "Percent of Ribo Cleared Reads"
	}
	if ( alignPhase == "Splice") {
		yLabel <- "Percent of Splice Reads"
	}

	mainText <- paste( sampleID, "    ", banner, "\nBase Mismatch Counts", "\nN_Reads: ", 
			prettyNum( as.integer( alignStatsData$nReads), big.mark=","))

	barplot( flipPcts, main=mainText, xlab="Number of Base Mismatchs per Aligned Read", 
		ylab=yLabel, col="slateblue", xlim=c(0,xLim), ylim=yLim, 
		width=0.75, cex.lab=1.2, cex.names=1.2, cex.axis=1.2)
	legend( "top", legend=paste( "Average Mismatchs per Read: ", formatC( avg, digits=2, format="f")), 
		bg="white", cex=1.05)

	# the letters...
	lets <- baseStats$Letters
	letPcts <- 100 * lets / sum(lets, na.rm=TRUE)
	# scale to fit on existing the plot
	letPctsPlot <- letPcts * ( max(flipPcts)/100)
	cumPctsPlot <- cumsum( letPctsPlot)
	colset <- c( "red","blue","orange","green","brown")
	# draw it at the right edge
	xbl <- bigX * 1.02
	xbr <- bigX * 1.07
	for( j in 5:1) rect( xbl,0,xbr,cumPctsPlot[j],col=colset[j])
	txtLets <- base::paste( names(lets), " ", round(letPcts), "%",sep="")
	legend( "bottomright", legend=rev(txtLets), fill=rev(colset), cex=1.10)
	text( (bigX*1.14), max(flipPcts)*0.9, label="Base Letter", pos=4, cex=1.1)
	if ( asPNG) {
		pngfile <- paste( sampleID, alignPhase, "Base.Mismatchs.png", sep=".")
		dev.print( png, file=file.path( plotPath, pngfile), width=800, height=600, bg='white')
	} 
	if (pause > 0) Sys.sleep(pause)


	# where they are in the read
	posCnts <- baseStats$Positions
	names( posCnts) <- 1:length( posCnts)
	ttlCnts <- sum( posCnts, na.rm=TRUE)
	if ( ttlCnts > 0) {
		posPcts <- 100 * posCnts / sum( posCnts, na.rm=TRUE)
		yLim <- c( 0, max(posPcts) * 1.1)
	} else {
		posPcts <- posCnts
		yLim <- c( 0, 100)
	}
	# let's try using totoal reads instead
	ttlReads <- alignStatsData$nReads
	posPcts <- 100 * posCnts / ttlReads
	yLim <- c( 0, max(posPcts) * 1.5)

	mainText <- paste( sampleID, "    ", banner, "\nBase Mismatch Positions", "\nN_Reads: ", 
			prettyNum( as.integer( alignStatsData$nReads), big.mark=","))
	barplot( posPcts, main=mainText, ylim=yLim,
		xlab="Base Mismatch Position within Aligned Read", ylab="Percent of All Reads",
		col="steelblue", cex.lab=1.2, cex.names=1.2, cex.axis=1.2, border=(length(posPcts)>200))
	if ( asPNG) {
		pngfile <- paste( sampleID, alignPhase, "Base.Positions.png", sep=".")
		dev.print( png, file=file.path( plotPath, pngfile), width=800, height=600, bg='white')
	} 
	if (pause > 0) Sys.sleep(pause)

	return()
}


plotIndelStatData <- function( alignStatsData, asPNG=FALSE, plotPath=".", pause=0) {

	# plot a few types...
	saveMAI <- par("mai")
	sampleID <- alignStatsData$sampleID
	banner <- alignStatsData$banner
	alignPhase <- alignStatsData$alignPhase
	nReads <- alignStatsData$nReads
	totalReads <- alignStatsData$totalReads

	indelStats <- alignStatsData$indelStatistics

	# number of insertions and deletions per read
	par("mai"=c( 1.02, 0.92, 0.92, 0.2))
	avgI <- indelStats$AvgN_Insert
	avgD <- indelStats$AvgN_Delete
	insCnts <- indelStats$InsertFrequency
	delCnts <- indelStats$DeleteFrequency
	totalI <- sum( insCnts, na.rm=TRUE)
	totalD <- sum( delCnts, na.rm=TRUE)

	insPcts <- 100 * insCnts / totalI
	delPcts <- 100 * delCnts / totalD
	ord <- base::order( as.numeric( names( insPcts)))
	insPcts <- insPcts[ ord]
	ord <- base::order( as.numeric( names( delPcts)))
	delPcts <- delPcts[ ord]

	# find how many are 'worth' plotting...
	bigX <- max( which( insPcts >= 0.1), which( delPcts >= 0.1))
	if (bigX < 5) bigX <- 5
	length(insPcts) <- length(delPcts) <- bigX
	insPcts[ is.na(insPcts)] <- 0
	delPcts[ is.na(delPcts)] <- 0
	xLim <- ceiling( bigX * 1.13)

	yLim <- c( 0, max(insPcts, delPcts, na.rm=T) * 1.16)
	
	yLabel <- "Percent of Genomic Aligned Reads"
	if ( alignPhase == "RiboClear") {
		yLabel <- "Percent of Ribo Cleared Reads"
	}
	if ( alignPhase == "Splice") {
		yLabel <- "Percent of Splice Reads"
	}

	mainText <- paste( sampleID, "    ", banner, "\nInsertion/Deletion Counts", "\nN_Reads: ", 
			prettyNum( as.integer( alignStatsData$nReads), big.mark=","))

	m <- matrix( 0, nrow=2, ncol=bigX)
	m[ 1, ] <- insPcts
	m[ 2, ] <- delPcts
	colnames(m) <- 0:(bigX-1)

	barplot( m, main=mainText, xlab="Number of Indels per Aligned Read", 
		ylab=yLabel, beside=TRUE, col=c("red","steelblue"), xlim=c(0,xLim), ylim=yLim, 
		width=0.75, cex.lab=1.2, cex.names=1.2, cex.axis=1.2)
	legend( "top", legend=paste( "Average", c("Insertions","Deletions  "), "per Read: ", 
			formatC( c(avgI,avgD), digits=3, format="f")), fill=c("red","steelblue"), 
			bg="white", cex=1.05)

	if ( asPNG) {
		pngfile <- paste( sampleID, alignPhase, "Indel.Counts.png", sep=".")
		dev.print( png, file=file.path( plotPath, pngfile), width=800, height=600, bg='white')
	} 
	if (pause > 0) Sys.sleep(pause)


	# where they are in the read
	insPosCnts <- indelStats$InsertPositions
	names( insPosCnts) <- 1:length( insPosCnts)
	delPosCnts <- indelStats$DeletePositions
	names( delPosCnts) <- 1:length( delPosCnts)
	ttlReads <- alignStatsData$nReads
	insPosPcts <- 100 * insPosCnts / ttlReads
	delPosPcts <- 100 * delPosCnts / ttlReads
	yLim <- c( 0, max(insPosPcts,delPosPcts, na.rm=T) * 1.4)

	mainText <- paste( sampleID, "    ", banner, "\nInsertion/Deletion Positions", "\nN_Reads: ", 
			prettyNum( as.integer( alignStatsData$nReads), big.mark=","))

	m <- matrix( 0, nrow=2, ncol=length(insPosPcts))
	m[ 1, ] <- insPosPcts
	m[ 2, ] <- delPosPcts
	colnames(m) <- 1:length(insPosPcts)

	barplot( m, main=mainText, ylim=yLim,
		xlab="Indel Position within Aligned Read", ylab="Percent of All Reads",
		beside=TRUE, col=c("red","steelblue"), cex.lab=1.2, cex.names=1.2, cex.axis=1.2,
		border=(ncol(m)>125))
	legend( "topleft", c("Insertions","Deletions"), fill=c("red","steelblue"), bg="white",cex=1.05)

	if ( asPNG) {
		pngfile <- paste( sampleID, alignPhase, "Indel.Positions.png", sep=".")
		dev.print( png, file=file.path( plotPath, pngfile), width=800, height=600, bg='white')
	} 
	if (pause > 0) Sys.sleep(pause)

	# pie of indel bases
	insPatt <- indelStats$InsertPatterns
	delPatt <- indelStats$DeletePatterns
	if (length(insPatt)) names(insPatt) <- paste( "Ins:", names(insPatt), sep=" ")
	if (length(delPatt)) names(delPatt) <- paste( "Del:", names(delPatt), sep=" ")
	combo <- c( insPatt, delPatt)
	comboCol <- c( rep( "red", times=length(insPatt)), rep( "steelblue", times=length(delPatt)))
	ord <- order( combo, decreasing=T)
	combo <- combo[ord]
	comboCol <- comboCol[ord]
	ttlReads <- alignStatsData$nReads
	combo <- 100 * combo / ttlReads
	
	N <- 25
	if ( N > length( combo)) N <- length(combo)
	mainText <- paste( sampleID, "    ", banner, "\n Top", N, " 'Insertion/Deletion' Base Strings     ", "\nN_Reads:",
					prettyNum( as.integer(ttlReads), format="d", big.mark=","))
	if ( alignPhase == "RiboClear") {
		subText <-  paste( "Ribo Cleared Reads contain ", 
					prettyNum( as.integer(sum(insPatt)), format="d", big.mark=","), "Insertions and ",
					prettyNum( as.integer(sum(delPatt)), format="d", big.mark=","), "Deletions")
	} else if (alignPhase == "Splice") {
		subText <-  paste( "Splice Reads contain ", 
					prettyNum( as.integer(sum(insPatt)), format="d", big.mark=","), "Insertions and ",
					prettyNum( as.integer(sum(delPatt)), format="d", big.mark=","), "Deletions")
	} else {
		subText <-  paste( "Genomic Reads contain ", 
					prettyNum( as.integer(sum(insPatt)), format="d", big.mark=","), "Insertions and ",
					prettyNum( as.integer(sum(delPatt)), format="d", big.mark=","), "Deletions")
	}

	if ( N > 0) {
		pcts <- combo[ 1:N]
		comboCol <- comboCol[1:N]
		lbls <- paste( names(pcts), "  Pct= ", formatC( pcts, digits=3, format="f"),"%", sep="")
		pie( pcts, labels=lbls, radius=0.8, clockwise=FALSE, main=mainText, col=comboCol,
				cex.sub=1.2, font.sub=2, sub=subText, cex=0.85)
	}
	legend( "topleft", legend=c("Insert","Delete"), fill=c("red","steelblue"), cex=1.05, bg="white")
	if ( asPNG) {
		pngfile <- paste( sampleID, alignPhase, "Indel.BasePatterns.png", sep=".")
		dev.print( png, file=file.path( plotPath, pngfile), width=800, height=600, bg='white')
	} 
	par( "mai"=saveMAI)
	return()
}


plotMARstatData <- function( alignStatsData, asPNG=FALSE, plotPath=".", pause=0) {

	saveMAI <- par("mai")
	sampleID <- alignStatsData$sampleID
	banner <- alignStatsData$banner
	alignPhase <- alignStatsData$alignPhase
	nReads <- alignStatsData$nReads
	totalReads <- alignStatsData$totalReads

	# gather the detail from the names and counts
	marStats <- alignStatsData$marStatistics$Genes
	marNames <- names( marStats)
	marTerms <- strsplit( marNames, split="::")
	marSpecies <- sapply( marTerms, function(x) x[1])
	marGenes <- sapply( marTerms, function(x) x[2])
	marPcts <- marStats * 100 / sum( marStats)

	# color by species, uniqueness, etc
	speciesForColoring <- sub( "\\+.+", "", marSpecies)
	colorAns <- colorBySpecies( speciesForColoring, intergenic=(marGenes == "Intergenic Only"))
	marcolors <- colorAns$colors
	legendColors <- colorAns$legend.colors
	legendColors <- c( legendColors, "Mixed"="mediumpurple1")
	legendNames <- names(legendColors)

	# ones with mixed species get a different color
	isMixedSpecies <- grep( "+", marSpecies, fixed=T)
	marcolors[ isMixedSpecies] <- "mediumpurple1"

	# ones with a mix of intergenics get a bit darker
	dark1Colors <- adjustColorSet( marcolors, -0.2)
	dark2Colors <- adjustColorSet( marcolors, -0.4)
	isMulti <- grep( "Multiple", marGenes, fixed=T)
	marcolors[ isMulti] <- dark1Colors[ isMulti]
	isMixed <- grep( "+ Intergenic", marGenes, fixed=T)
	marcolors[ isMixed] <- dark2Colors[ isMixed]

	marCounts <- alignStatsData$marStatistics$Counts
	totalCounts <- sum( marCounts)
	subtotals <- as.integer( names( marCounts)) * marCounts
	avgCnt <- sum( subtotals) / totalCounts
	subText <- paste( "Average Aligments per MAR =  ", formatC( avgCnt, format="f", digits=2))

	# pie of genes
	par("mai"=c( 1.02, 0.52, 0.92, 0.5))
	marPcts <- as.numeric(marPcts)
	mylabels <- paste( marSpecies, "  ", marGenes, "  ", 
				formatC( marPcts, digits=2, format="f"),"%", sep="")
	mylabels[ marPcts < 0.25] <- ""
	ord <- order( marPcts, decreasing=T)
	mainText <- paste( sampleID, "    ", banner, "\n'Multiple Aligned Reads' ", "\nN_MARs:",
			prettyNum( as.integer( sum(marStats)), big.mark=","))

	pie( marStats[ord], labels=mylabels[ord], col=marcolors[ord], radius=0.8, clockwise=FALSE, main=mainText,
			cex=0.95, sub=subText, font.sub=2, cex.sub=1.2)
	legend( "topright", legendNames, fill=legendColors, cex=1.05, bg="white")

	if ( asPNG) {
		pngfile <- paste( sampleID, alignPhase, "MAR.Pie.png", sep=".")
		dev.print( png, file=file.path( plotPath, pngfile), width=800, height=600, bg='white')
	} 
	if (pause > 0) Sys.sleep(pause)

	par( "mai"=saveMAI)
	return()
}


plotAlignScoreData <- function( alignStatsData, asPNG=FALSE, plotPath=".", pause=0) {

	saveMAI <- par("mai")
	sampleID <- alignStatsData$sampleID
	banner <- alignStatsData$banner
	alignPhase <- alignStatsData$alignPhase
	nReads <- alignStatsData$nReads
	totalReads <- alignStatsData$totalReads

	alignStats <- alignStatsData$alignStatistics

	# number of scores of each AS value
	par("mai"=c( 1.02, 0.92, 0.92, 0.2))
	cnts <- alignStats
	scores <-names( cnts)
	pcts <- 100 * cnts / sum(cnts)

	# put these into a bargraph with all intervening pts
	scoreRange <- range( as.integer( scores))
	N <- diff( scoreRange) + 1
	if ( is.infinite(N)) return()
	ans <- rep( 0, times=N)
	names(ans) <- scoreRange[1]:scoreRange[2]
	where <- match( scores, names(ans))
	ans[ where] <- pcts

	# drop any that are very negative...
	keeps <- which( as.numeric( names(ans)) >= -40)
	ans <- ans[ keeps]

	yLim <- c( 0, max(pcts, na.rm=T) * 1.16)
	
	yLabel <- "Percent of Genomic Aligned Reads"
	if ( alignPhase == "RiboClear") yLabel <- "Percent of Ribo Cleared Reads"
	if ( alignPhase == "Splice") yLabel <- "Percent of Spliced Reads"

	mainText <- paste( sampleID, "    ", banner, "\nDistribution of Alignment Scores", "\nN_Reads: ", 
			prettyNum( as.integer( alignStatsData$nReads), big.mark=","))

	barplot( ans, main=mainText, xlab="Distribution of Alignment Scores",
		ylab=yLabel, col="lightgreen", ylim=yLim, width=0.75, cex.lab=1.2, cex.names=1.2, cex.axis=1.2)

	if ( asPNG) {
		pngfile <- paste( sampleID, alignPhase, "AlignScores.png", sep=".")
		dev.print( png, file=file.path( plotPath, pngfile), width=800, height=600, bg='white')
	} 
	if (pause > 0) Sys.sleep(pause)
	
	return()
}
