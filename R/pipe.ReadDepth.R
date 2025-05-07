# pipe.ReadDepth.R --  inspect the WIGGLE depth to find depth metrics per gene

`pipe.ReadDepth` <- function( sampleID, geneIDset=NULL, speciesID=getCurrentSpecies(), 
				annotationFile="Annotation.txt", optionsFile="Options.txt", results.path=NULL, 
				strandSpecific=FALSE, exonsOnly=FALSE) {
				
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

	wigReadLen <- wiggles$Info$ReadLength

	# we will visit every real gene
	geneMap <- subset( getCurrentGeneMap(), REAL_G == TRUE)
	if ( ! is.null( geneIDset)) {
		geneMap <- subset( geneMap, GENE_ID %in% geneIDset)
	}
	if (exonsOnly) {
		exonMap <- getCurrentExonMap()
		exonMap <- subset( exonMap, GENE_ID %in% geneMap$GENE_ID)
	}
	gids <- geneMap$GENE_ID
	sids <- geneMap$SEQ_ID
	curSeqID <- skipSeqID <- ""
	N <- length(gids)

	outAvgDep <- outReads <- vector( length=N)

	for ( i in 1:N) {
		# grab the wiggle data once for each chromosome
		thisG <- gids[i]
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
			cat( "\n  ",thisS, "\n")
		}

		# gather the set of all wiggle base depths in this range
		start <- geneMap$POSITION[i]
		stop <- geneMap$END[i]
		myStrand <- geneMap$STRAND[i]
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

		# all bases, or by exon?
		if ( exonsOnly) {
			emap <- subset.data.frame( exonMap, GENE_ID == thisG)
			if ( nrow(emap)) {
				#exonBases <- vector()
				#for (k in 1:nrow(emap)) exonBases <- c( exonBases, (emap$POSITION[k] : emap$END[k]))
				exonBases <- unlist( mapply( FUN=`:`, emap$POSITION, emap$END, SIMPLIFY=F))
				keep <- which( allBases %in% exonBases)
				baseM <- baseM[ keep, ]
				NBases <- nrow(baseM)
			}
		}

		# sum the strands and see how many bases are too low
		if ( ! strandSpecific) {
			depth <- apply( baseM, 1, sum)
		} else if ( myStrand == "-") {
			depth <- baseM[ , "Minus"]
		} else {
			depth <- baseM[ , "Plus"]
		} 
		
		avgDepth <- sum(depth) / NBases
		nReads <- sum(depth) / wigReadLen
		outAvgDep[i] <- avgDepth
		outReads[i] <- nReads
		if ( i %% 100 == 0) cat( "\r", i, gids[i], round(nReads), avgDepth)
	}

	out <- data.frame( "GENE_ID"=gids, "PRODUCT"=geneMap$PRODUCT, "SEQ_ID"=sids, 
				"N_READS"=round( outReads), "AVG_DEPTH"=round(outAvgDep, digits=2), 
				stringsAsFactors=F)
	ord <- order( -(out$AVG_DEPTH), out$SEQ_ID, out$GENE_ID, decreasing=F)
	out <- out[ ord, ]
	rownames(out) <- 1:nrow(out)

	# assign a rough copy number to every gene
	out$COPY_NUMBER <- round( out$AVG_DEPTH / median(out$AVG_DEPTH), digits=2)

	return( out)
}

