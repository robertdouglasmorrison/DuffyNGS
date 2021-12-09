# pipe.IsoformExpression.R -- compare splice junction detection depths

`pipe.IsoformExpression` <- function( sampleIDset, geneIDset=NULL, results.path=NULL, speciesID=getCurrentSpecies(), 
				annotationFile="Annotation.txt", optionsFile="Options.txt", visualize=F,
				verbose=TRUE) {

	# get options that we need
	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	splice.path <- file.path( results.path, "splicing")
	if ( speciesID != getCurrentSpecies()) setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()

	# make sure we have the splice summary files
	spliceFiles <- file.path( splice.path, paste( sampleIDset, "splice", prefix, "Summary.txt", sep="."))
	missing <- ! file.exists( spliceFiles)
	if ( any( missing)) {
		cat( "\nMissing Splice Summary files for some samples: ", sampleIDset[missing])
		sampleIDset <- sampleIDset[ ! missing]
		spliceFiles <- spliceFiles[ ! missing]
	}
	NS <- length( sampleIDset)
	if ( ! length( sampleIDset)) return(NULL)

	# get the set of genes to investigate
	# note that in the world of splice info, all GeneIDs are the short names, so adjust all
	exonMap <- getCurrentExonMap()
	exonMap$GENE_ID <- shortGeneName( exonMap$GENE_ID, keep=1)
	allGenes <- sort( unique( exonMap$GENE_ID))
	if ( is.null( geneIDset)) {
		geneIDset <- allGenes
	} else {
		geneIDset <- intersect( shortGeneName( geneIDset, keep=1), allGenes)
	}
	NG <- length( geneIDset)
	if ( NG < 1) {
		cat( "\nNo Genes selected")
		return(NULL)
	}


	`getSpliceCountsOneGene` <- function( gid, spliceDF) {

		emap <- subset.data.frame( exonMap, GENE_ID == gid)
		nExons <- nrow(emap)
		if ( nExons < 2) return(NULL)
		strand <- emap$STRAND[1]

		spliceDF <- subset( spliceDF, GENE_ID == gid)
		if ( ! nrow(spliceDF)) return(NULL)

		# build up the names of all the splice junctions that you should see
		if (strand == "+") {
			eLeftNum <- 1:(nExons-1)
			eRightNum <- 2:nExons
		} else {
			eRightNum <- 1:(nExons-1)
			eLeftNum <- 2:nExons
		}
		stdNames <- paste( "std:E", eLeftNum, "_E", eRightNum, sep="")
		skip1Names <- skip2Names <- character(0)
		if (nExons > 2) {
			if (strand == "+") {
				eLeftNum <- 1:(nExons-2)
				eRightNum <- 3:nExons
			} else {
				eRightNum <- 1:(nExons-2)
				eLeftNum <- 3:nExons
			}
			skip1Names <- paste( "alt:E", eLeftNum, "_E", eRightNum, sep="")
		}
		if (nExons > 3) {
			if (strand == "+") {
				eLeftNum <- 1:(nExons-3)
				eRightNum <- 4:nExons
			} else {
				eRightNum <- 1:(nExons-3)
				eLeftNum <- 4:nExons
			}
			skip2Names <- paste( "alt:E", eLeftNum, "_E", eRightNum, sep="")
		}
		altNames <- c( skip1Names, skip2Names)
		allNames <- c( stdNames, altNames)
		if ( ! length(allNames)) return(NULL)

		# OK, now access the counts of all the splice junctions that you should see
		stdCounts <- spliceDF$N_READS[ match( stdNames, spliceDF$SPLICE_ID)]
		stdCounts[ is.na( stdCounts)] <- 0
		names(stdCounts) <- stdNames
		NSTD <- length(stdCounts)
		altCounts <- spliceDF$N_READS[ match( altNames, spliceDF$SPLICE_ID)]
		altCounts[ is.na( altCounts)] <- 0
		names(altCounts) <- altNames
		NALT <- length(altCounts)

		# for each alt splice, evaluate its depth relative to the left exon
		stdExonNames <- sub( "(^std:)(E[0-9]+)(_E[0-9]+)", "\\2", names(stdCounts))
		altExonNames <- sub( "(^alt:)(E[0-9]+)(_E[0-9]+)", "\\2", names(altCounts))
		myStdPtr <- match( altExonNames, stdExonNames)
		myStdCounts <- stdCounts[ myStdPtr]
		# prevent divide by zero
		myStdCounts[ is.na(myStdCounts)] <- 0
		myStdCountsForRatio <- myStdCounts
		myStdCountsForRatio[ myStdCountsForRatio == 0] <- 1
		altRatios <- round( altCounts / myStdCountsForRatio, digits=3)
		names(altRatios) <- names(altCounts)

		# join these into one small data.frame
		out <- data.frame( "SPLICE_ID"=allNames, "STD_COUNT"=c( stdCounts, rep.int(0,NALT)), 
				"ALT_COUNT"=c( rep.int(0,NSTD), altCounts), "ALT_RATIO"=c( rep.int(0,NSTD), altRatios), 
				stringsAsFactors=F)
		rownames(out) <- 1:nrow(out)
		return( out)
	}	


	`getSpliceCountsOneSample` <- function( sid, geneIDset) {

		spliceFile <- spliceFiles[ match(sid, sampleIDset)]
		spliceDF <- read.delim( spliceFile, as.is=T)

		out <- data.frame()
		for ( g in geneIDset) {
			ans <- getSpliceCountsOneGene( g, spliceDF)
			if ( is.null( ans)) next

			# make a small data.frame of these calls
			sml <- data.frame( "GENE_ID"=g, ans, stringsAsFactors=F)
			out <- rbind( out, sml)
		}
		if ( nrow(out)) rownames(out) <- 1:nrow(out)
		return( out)
	}


	out <- vector( mode="list")
	for ( i in 1:NS) {

		sid <- sampleIDset[i]
		smlDF <- getSpliceCountsOneSample( sid, geneIDset)
		out[[i]] <- smlDF
		if (verbose) cat( "\n", i, sid, sum( smlDF$ALT_RATIO >= 1))
	}

	# if doing a single sample, return the result as it is
	if ( NS == 1) {
		tmp <- out[[i]]
		prod <- gene2Product( tmp$GENE_ID)
		out <- data.frame( "GENE_ID"=tmp$GENE_ID, "PRODUCT"=prod, tmp[,2:ncol(tmp)], stringsAsFactors=F)

		# allow showing the splice counts visually
		if ( visualize && length(sampleIDset) == 1 && length(geneIDset) == 1) {
			plotSpliceAlignmentCounts( sampleIDset[1], geneIDset[1], spliceCounts=out, plotPileups=T)
		}
		return( out)
	}

	# if we are given a list of samples, turn these into a matrix of results
	# create the union of all IDs we may see
	allGenes <- allSplices <- vector()
	for ( i in 1:NS) {
		smlDF <- out[[i]]
		if ( is.null(smlDF)) next
		allGenes <- c( allGenes, smlDF$GENE_ID)
		allSplices <- c( allSplices, smlDF$SPLICE_ID)
	}
	allSpliceKey <- sort( unique( paste( allGenes, allSplices, sep="::")))
	NSPL <- length( allSpliceKey)
	if ( ! NSPL) return(NULL)

	keyGene <- sub( "::.+", "", allSpliceKey)
	keySplice <- sub( ".+::", "", allSpliceKey)
	keyProd <- gene2Product( keyGene)

	countM <- ratioM <- stdCountM <- matrix( 0, nrow=NSPL, ncol=NS)
	colnames(countM) <- paste( "ALT_COUNT", sampleIDset, sep="_")
	colnames(ratioM) <- paste( "ALT_RATIO", sampleIDset, sep="_")
	colnames(stdCountM) <- paste( "STD_COUNT", sampleIDset, sep="_")
	rownames(countM) <- rownames(ratioM) <- rownames(stdCountM) <- 1:NSPL
	for ( i in 1:NS) {
		smlDF <- out[[i]]
		myKey <- paste( smlDF$GENE_ID, smlDF$SPLICE_ID, sep="::")
		wh <- match( myKey, allSpliceKey)
		stdCountM[ wh, i] <- smlDF$STD_COUNT
		countM[ wh, i] <- smlDF$ALT_COUNT
		ratioM[ wh, i] <- smlDF$ALT_RATIO
	}

	out <- data.frame( "GENE_ID"=keyGene, "PRODUCT"=keyProd, "SPLICE_ID"=keySplice, 
			stdCountM, "MEAN_STD_COUNT"=round(apply(stdCountM,1,mean,na.rm=T),digits=1), 
			countM, "MEAN_ALT_COUNT"=round(apply(countM,1,mean,na.rm=T),digits=1), 
			ratioM, "MEAN_ALT_RATIO"=round(apply(ratioM,1,mean,na.rm=T),digits=3), 
			stringsAsFactors=F)
	rownames(out) <- 1:nrow(out)

	return( out)
}


`plotSpliceAlignmentCounts` <- function( sampleID, geneID, spliceCounts=NULL, plotPileups=TRUE, 
					splice.col=1, pileup.col=2, ...) {

	# visualize counts of standard and alternate splices
	if (plotPileups) pipe.PlotGene( sampleID, geneID, pileup.col=pileup.col, addStrands=T, legend.cex=0, ...)

	# either calculate the splice counts or use what was given
	if ( is.null( spliceCounts)) spliceCounts <- pipe.IsoformExpression( sampleID, geneID, verbose=F)

	# show the counts, using the splice locations on the chromosome
	emap <- getCurrentExonMap()
	keep <- which( shortGeneName( emap$GENE_ID,keep=1) == shortGeneName(geneID,keep=1))
	emap <- emap[ keep, ]
	emap$MID_PT <- (emap$POSITION + emap$END) / 2
	
	# the splice info uses names like "E1_E2", so extract that to know where each splice belongs on the plot
	spliceExonLeft <- as.numeric( sub( "(^.+:E)([0-9]+)(_E.+)", "\\2", spliceCounts$SPLICE_ID))
	spliceExonRight <- as.numeric( sub( "(^.+_E)([0-9]+$)", "\\2", spliceCounts$SPLICE_ID))

	# use those to assign a location for each splice
	spliceCounts$LeftPt <- emap$MID_PT[ spliceExonLeft]
	spliceCounts$RightPt <- emap$MID_PT[ spliceExonRight]

	# show the standard splices first, when it has counts
	MIN_COUNT <- 1
	if ( (BigCount <- max( c( spliceCounts$STD_COUNT, spliceCounts$ALT_COUNT), na.rm=T)) > 100) MIN_COUNT <- BigCount * 0.01
	isStd <- which( abs( spliceExonLeft - spliceExonRight) == 1)
	hasCnts <- which( spliceCounts$STD_COUNT >= MIN_COUNT)
	toShow <- intersect( isStd, hasCnts)
	if ( length( toShow)) {
		for ( k in toShow) lines( c(spliceCounts$LeftPt[k], spliceCounts$RightPt[k]), rep.int( spliceCounts$STD_COUNT[k],2),
						col=splice.col, lwd=2, lty=1)
	}
	isAlt2 <- which( abs( spliceExonLeft - spliceExonRight) == 2)
	hasCnts <- which( spliceCounts$ALT_COUNT >= MIN_COUNT)
	toShow <- intersect( isAlt2, hasCnts)
	if ( length( toShow)) {
		for ( k in toShow) lines( c(spliceCounts$LeftPt[k], spliceCounts$RightPt[k]), rep.int( spliceCounts$ALT_COUNT[k],2),
						col=splice.col, lwd=3, lty=2)
	}
	isAlt3 <- which( abs( spliceExonLeft - spliceExonRight) == 3)
	toShow <- intersect( isAlt3, hasCnts)
	if ( length( toShow)) {
		for ( k in toShow) lines( c(spliceCounts$LeftPt[k], spliceCounts$RightPt[k]), rep.int( spliceCounts$ALT_COUNT[k],2),
						col=splice.col, lwd=4, lty=3)
	}
	legend( 'topright', c( "Standard", "Skip 1 Exon", "Skip 2 Exons"), lwd=c(2,3,4), lty=c(1,2,3), col=splice.col, 
		bg='white', title="Splice Junction Counts")
	dev.flush()
}
