# pipe.FindDeletions.R --  inspect the WIGGLE depth to find chromosomal deletions

`pipe.FindDeletions` <- function( sampleID, speciesID=getCurrentSpecies(), annotationFile="Annotation.txt", 
				optionsFile="Options.txt", results.path=NULL, min.depth=1, 
				min.pct.deleted=10) {
				
	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	if ( getCurrentSpecies() != speciesID) setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()

	wigfile <- file.path( results.path, "wig", paste( sampleID, prefix, "WIG.rda", sep="."))

	# load that WIG data
	load( wigfile)

	# see what type it is
	dataType <- wiggles$Info$DataType
	if ( is.null( dataType)) dataType <- "RNA-seq"
	if ( dataType %in% c( "RNA-seq", "DNA-seq")) {
		trackNames <- c("Plus", "Minus")
		NTracks <- length( trackNames)
	}
	if ( dataType %in% c( "ChIP-seq")) {
		trackNames <- c("Plus", "Minus")
		NTracks <- length( trackNames)
	}
	if ( ! exists( "NTracks")) stop( "Unknown WIG data type")


	# we will visit every real gene
	geneMap <- subset( getCurrentGeneMap(), REAL_G == TRUE)
	gids <- geneMap$GENE_ID
	sids <- geneMap$SEQ_ID
	curSeqID <- skipSeqID <- ""
	N <- length(gids)

	outAvgDep <- outPctDel <- vector( length=N)

	for ( i in 1:N) {
		# grab the wiggle data once for each chromosome
		thisS <- sids[i]
		if ( thisS != curSeqID) {
			# get the wanted chromosome
			# if we already tried and failed to load a chromosome, don't repeat...
			if ( thisS == skipSeqID) next
			wiggleChunk <- WIG_getWigglesOneSeq( wiggles, thisS)
			if ( is.null( wiggleChunk)) {
				cat( "\nFailed to load Wiggle data for seqID: ", thisS)
				skipSeqID <- thisS
				next
			}
			curSeqID <- thisS
		}

		# gather the set of all wiggle base depths in this range
		start <- geneMap$POSITION[i]
		stop <- geneMap$END[i]
		allBases <- start : stop
		NBases <- length(allBases)
		baseM <- matrix( 0, nrow=NBases, ncol=NTracks)
		rownames(baseM) <- allBases
		colnames(baseM) <- trackNames
		for ( j in 1:NTracks) {
			baseDepthTbl <- wiggleChunk[[ trackNames[j]]]
			if ( ! nrow( baseDepthTbl)) next
	
			# get the subset we want
			baseDepthTbl <- baseDepthTableSubset( baseDepthTbl, start, stop)
			if ( ! nrow( baseDepthTbl)) next
		
			# turn it into the 1-D vector of depths
			baseDepthVec <- baseDepthTableToVector( baseDepthTbl)
	
			# load those locations
			where <- match( rownames(baseM), names(baseDepthVec), nomatch=0)
			if ( any( where > 0)) baseM[ where > 0, j] <- baseDepthVec[ where]
		}

		# sum the strands and see how many bases are too low
		depth <- apply( baseM, 1, sum)
		nLow <- sum( depth < min.depth)
		pctLow <- nLow * 100 / NBases
		avgDepth <- sum(depth) / NBases
		outAvgDep[i] <- avgDepth
		outPctDel[i] <- pctLow
		if ( i %% 100 == 0) cat( "\r", i, gids[i], pctLow, avgDepth)
	}

	out <- data.frame( "GENE_ID"=gids, "PRODUCT"=geneMap$PRODUCT, "SEQ_ID"=sids, 
				"PCT_DELETED"=round( outPctDel), "AVG_DEPTH"=round(outAvgDep), 
				stringsAsFactors=F)
	ord <- order( -(out$PCT_DELETED), out$SEQ_ID, out$GENE_ID, decreasing=F)
	out <- out[ ord, ]

	# only keep/return those with some deletion...
	keep <- which( out$PCT_DELETED >= min.pct.deleted)
	out <- out[ keep, ]
	if ( nrow(out)) rownames(out) <- 1:nrow(out)
	out
}

