# calcWIGmetrics.R

# calculate some metrics about how the WIG reads compare to the annotation

# changing from bins to the wiggle tracks compressed objects


`calcWIGmetrics` <- function( WIG, asDataFrame=FALSE, logFile=NULL ) {

	setCurrentSpecies( WIG$Info$Species)
	seqMap <- getCurrentSeqMap()
	geneMap <- getCurrentGeneMap()
	exonMap <- getCurrentExonMap()


	if ( ! is.null( logFile)) {
		sink( file=logFile, split=TRUE)
	}

	# count up how many reads hit the right and wrong strands, as per the genome annotation
	Nright <- Nwrong <- NinterG <- 0
	Nright_combo <- Nwrong_combo <- NinterG_combo <- 0
	ngenes <- nngs <- 0
	gOut <- readsOut <- correctPct <- vector()
	curSeqID <- ""

	nGenes1deep <- nGenes10deep <- nGenes100deep <- 0
	readLength <- WIG$Info$ReadLength
	
	for( ig in 1:nrow( geneMap)) {
		if ( ("REAL_G" %in% colnames( geneMap)) && geneMap$REAL_G[ig] == FALSE) {
			nngs <- nngs + 1
			isReal <- FALSE
		} else {
			ngenes <- ngenes + 1
			isReal <- TRUE
		}

		# get that chromosome's set of wiggle bins
		thisSeqID <- geneMap$SEQ_ID[ig]
		if ( curSeqID != thisSeqID) {
			wiggleChunk <- WIG_getWigglesOneSeq( WIG, thisSeqID)
			curSeqID <- thisSeqID
		}

		# get the chunks for this gene
		start <- geneMap$POSITION[ig]
		stop <- geneMap$END[ig]

		# get the read counts for this region
		counts <- WIG_getReadCountsInRegion( wiggleChunk, start, stop, readLength=readLength)
		nPlusMulti <- counts$Plus
		nPlusUnique <- counts$PlusUnique
		nMinusMulti <- counts$Minus
		nMinusUnique <- counts$MinusUnique

		if ( ig %% 100 == 0) cat( "\r", ig, geneMap$GENE_ID[ig], "\t", 
				nPlusMulti+nMinusMulti+nPlusUnique+nMinusUnique, "   ")

		# inter-genic regions may get special treatment
		if ( is.na( geneMap$STRAND[ig])) {
			NinterG <- NinterG + nPlusUnique + nMinusUnique
			NinterG_combo <- NinterG_combo + nPlusMulti + nMinusMulti
			# no other intergenic stuff
			next

		} else if ( geneMap$STRAND[ig] %in% c("","+")) {
			Nright <- Nright + nPlusUnique
			Nright_combo <- Nright_combo + nPlusMulti
			goodCnts <- nPlusUnique
			Nwrong <- Nwrong + nMinusUnique
			Nwrong_combo <- Nwrong_combo + nMinusMulti
			badCnts <- nMinusUnique
		} else {
			Nright <- Nright + nMinusUnique
			Nright_combo <- Nright_combo + nMinusMulti
			goodCnts <- nMinusUnique
			Nwrong <- Nwrong + nPlusUnique
			Nwrong_combo <- Nwrong_combo + nPlusMulti
			badCnts <- nPlusUnique
		}

		# other gene based metrics...
		if ( isReal) {
			nGenes1deep <- nGenes1deep + sum( goodCnts >= 1)
			nGenes10deep <- nGenes10deep + sum( goodCnts >= 10)
			nGenes100deep <- nGenes100deep + sum( goodCnts >= 100)
		}

		# build up some coverage numbers for each gene
		gOut[ ngenes] <- geneMap$GENE_ID[ ig]
		readsOut[ ngenes] <- goodCnts
		correctPct[ ngenes] <- ( goodCnts - badCnts) / ( goodCnts + badCnts) 
	}


	cat("\n\nStrand Correctness of 'In Gene' reads: \n")
	Nright <- round( Nright)
	Nwrong <- round( Nwrong)
	total <- Nright + Nwrong
	options( digits=4)
	cat( "\nUnique Reads hitting expected strand:    ", Nright, "\t Pct: ", as.percent( Nright / total))
	cat( "\nUnique Reads hitting wrong strand:       ", Nwrong, "\t Pct: ", as.percent( Nwrong / total))
	Nright_combo <- round( Nright_combo)
	Nwrong_combo <- round( Nwrong_combo)
	total_combo <- Nright_combo + Nwrong_combo
	cat( "\nAll Reads hitting expected strand:       ", Nright_combo, "\t Pct: ", 
			as.percent( Nright_combo / total_combo))
	cat( "\nAll Reads hitting wrong strand:          ", Nwrong_combo, "\t Pct: ", 
			as.percent( Nwrong_combo / total_combo))

	cat( "\n\nGene Specificity of reads: \n")
	NinterG <- round( NinterG)
	NinterG_combo <- round( NinterG_combo)
	totalAllReads <- total + NinterG
	cat( "\nUnique Reads hitting annotated genes:    ", total, "\t Pct: ", 
			as.percent( total / totalAllReads))
	cat( "\nUnique Reads hitting intergenic areas:   ", NinterG, "\t Pct: ", 
			as.percent( NinterG / totalAllReads))
	totalAllReads_combo <- total_combo + NinterG_combo
	cat( "\nAll Reads hitting annotated genes:       ", total_combo, "\t Pct: ", 
			as.percent( total_combo / totalAllReads_combo))
	cat( "\nAll Reads hitting intergenic areas:      ", NinterG_combo, "\t Pct: ", 
			as.percent( NinterG_combo / totalAllReads_combo))

	cat( "\n\nGene Sensitivity of reads: \n")
	cat( "\nGenes with at least    1 read:           ", nGenes1deep, "\t Pct: ", 
			as.percent( nGenes1deep / ngenes))
	cat( "\nGenes with at least   10 reads:          ", nGenes10deep, "\t Pct: ", 
			as.percent( nGenes10deep / ngenes))
	cat( "\nGenes with at least  100 reads:          ", nGenes100deep, "\t Pct: ", 
			as.percent( nGenes100deep / ngenes))

	cat( "\n")
	if ( ! is.null( logFile)) {
		sink( file=NULL)
		cat( "\n\nWrote Bin Metrics summary file: ", logFile, "\n")
	}

	# prep the output data.frame
	if ( asDataFrame) {
		out <- data.frame( gOut, readsOut, correctPct, stringsAsFactors=FALSE)
		colnames( out) <- c( "GENE_ID", "RAW_READS_U", "STRAND_CC")
		ord <- base::order( out$RAW_READS_U, decreasing=TRUE)
		out <- out[ ord, ]
		rownames( out) <- 1:nrow( out)
		return(out)
	}
}

