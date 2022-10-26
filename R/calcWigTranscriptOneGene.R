# calcTranscriptOneGene.R

# calculate the Read count and RPKM for one gene

`calcWigTranscriptOneGene` <- function( wiggleChunk, gene, gmapPtr=NULL, geneMap=getCurrentGeneMap(), 
				totalReads=list(  "Unique"=1000000, "Multi"=1000000), 
				readLength=32, useBothStrands=FALSE, exonsOnly=FALSE) {


	# define a roubust error result
	errorOut <- list( "rawReads"=0, "sigma"=0, "rpkm"=0, "strandness"=0, 
		"rawReads.Multi"=0, "sigma.Multi"=0, 
		"rpkm.Multi"=0, "strandness.Multi"=0, 
		"nBases"=100) 

	# calculate all the transcription terms for one gene, for speed we should know what gene already
	if ( is.null( gmapPtr)) {
		gmapPtr <- base::match( gene, geneMap$GENE_ID, nomatch=0)[1]
		if ( gmapPtr == 0) {
			cat( "\ncalcWigTranscriptOneGene:  not a found gene: ", gene, 
				"\tSpecies: ", getCurrentSpecies(),"n")
			return(errorOut)
		}
	}
	ig <- gmapPtr
	start <- geneMap$POSITION[ig]
	stop <- geneMap$END[ig]
	nBases <- max( geneMap$N_EXON_BASES[ig], 100)

	# allow a graceful exit if we have no wiggle data...
	if ( is.null( wiggleChunk)) {
		errorOut <- list( "rawReads"=0, "sigma"=0, "rpkm"=0, "strandness"=0, 
			"rawReads.Multi"=0, "sigma.Multi"=0, 
			"rpkm.Multi"=0, "strandness.Multi"=0, 
			"nBases"=nBases) 
		return( errorOut)
	}

	# get the read counts for this region
	if ( exonsOnly) {
		exonMap <- getCurrentExonMap()
		eptrs <- which( exonMap$GENE_ID == gene)
		if ( length(eptrs) < 1) exonsOnly <- FALSE
	}
	if ( exonsOnly) {
		nPlusMulti <- nPlusUnique <- nMinusMulti <- nMinusUnique <- nBases <- 0
		for ( j in eptrs) {
			counts <- WIG_getReadCountsInRegion( wiggleChunk, exonMap$POSITION[j], 
						exonMap$END[j], readLength=readLength)
			nPlusMulti <- nPlusMulti + counts$Plus
			nPlusUnique <- nPlusUnique + counts$PlusUnique
			nMinusMulti <- nMinusMulti + counts$Minus
			nMinusUnique <- nMinusUnique + counts$MinusUnique
			nBases <- nBases + (exonMap$END[j] - exonMap$POSITION[j] + 1)
		}
	} else {
		counts <- WIG_getReadCountsInRegion( wiggleChunk, start, stop, readLength=readLength)
		nPlusMulti <- counts$Plus
		nPlusUnique <- counts$PlusUnique
		nMinusMulti <- counts$Minus
		nMinusUnique <- counts$MinusUnique
	}
	if ( is.null( nPlusUnique)) nPlusUnique <- nPlusMulti
	if ( is.null( nMinusUnique)) nMinusUnique <- nMinusMulti

	# inter-genic regions may get special treatment
	interG <- (("REAL_G" %in% colnames( geneMap)) && (geneMap$REAL_G[ig] == FALSE))

	# only count the wiggles from the annotated strand, 
	# interGenic region use both strands
	if ( is.na( geneMap$STRAND[ig])) {
		tallyGoodUnique <- nPlusUnique + nMinusUnique
		tallyGoodMulti <- nPlusMulti + nMinusMulti
		tallyBadUnique <- tallyBadMulti <- 0
	} else if ( geneMap$STRAND[ig] %in% c("","+")) {
		tallyGoodUnique <- nPlusUnique
		tallyGoodMulti <- nPlusMulti
		tallyBadUnique <- nMinusUnique
		tallyBadMulti <- nMinusMulti
	} else {
		tallyGoodUnique <- nMinusUnique
		tallyGoodMulti <- nMinusMulti
		tallyBadUnique <- nPlusUnique
		tallyBadMulti <- nPlusMulti
	}
	if ( useBothStrands) {
		tallyGoodUnique <- tallyGoodUnique + tallyBadUnique
		tallyGoodMulti <- tallyGoodMulti + tallyBadMulti
	}

	# OK, we are ready to build the answer for this gene
	totalCnt <- tallyGoodUnique
	strandCC <- calcStrandCC( tallyGoodUnique, tallyBadUnique)
	sigma <- rosettaSigma( tallyGoodUnique)
	rpkm <- calculateRPKM( tallyGoodUnique, nBases, totalReads$Unique)

	# multis too
	combo_totalCnt <- totalCnt
	combo_sigma <- sigma
	combo_rpkm <- rpkm
	combo_strandCC <- strandCC
	if ( totalReads$Multi > 0) {
		combo_totalCnt <- tallyGoodMulti
		combo_badCnt <- tallyBadMulti
		combo_strandCC <- calcStrandCC( combo_totalCnt, combo_badCnt)
		combo_sigma <- rosettaSigma( combo_totalCnt)
		combo_rpkm <- calculateRPKM( combo_totalCnt, nBases, totalReads$Multi)
	}

	# make sure valid values if need be
	if ( is.null(rpkm) || is.nan(rpkm) || is.na(rpkm)) rpkm <- 0
	if ( is.null(combo_rpkm) || is.nan(combo_rpkm) || is.na(combo_rpkm)) combo_rpkm <- 0

	out <- list( "rawReads"=totalCnt, "sigma"=sigma, "rpkm"=rpkm, "strandness"=strandCC, 
			"rawReads.Multi"=combo_totalCnt, "sigma.Multi"=combo_sigma, 
			"rpkm.Multi"=combo_rpkm, "strandness.Multi"=combo_strandCC, 
			"nBases"=nBases) 
	return( out)
}

