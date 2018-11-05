# pipe.ExtractWiggleDepths.R --  extract the base depths for all 4 wiggle tracks in a region

`pipe.ExtractWiggleDepths` <- function( sampleID, geneID=NULL, seqID=NULL, start=NULL, stop=NULL, 
				speciesID=getCurrentSpecies(), annotationFile="Annotation.txt", 
				optionsFile="Options.txt", results.path=NULL, noDepth=c("zero","drop")) {
				
	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	if ( getCurrentSpecies() != speciesID) setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()

	wigfile <- file.path( results.path, "wig", paste( sampleID, prefix, "WIG.rda", sep="."))

	# control what to do when no reads cover a region
	noDepth <- match.arg( noDepth)

	# let's be a bit more flexible with what we pass in
	# Could be a genomic range {seqID, start, stop}, a gene {geneID}, or any generic sequence...
	isRange <- isGene <- isGeneric <- FALSE
	geneMap <- getCurrentGeneMap()
	seqMap <- getCurrentSeqMap()

	# specify by either gene or range
	if ( is.null( geneID)) {
		if ( is.null(seqID)) stop( "One of 'seqID','geneID' must be not NULL")
		if ( seqID %in% seqMap$SEQ_ID) {
			isRange <- TRUE
			seqMap <- subset( seqMap, SEQ_ID == seqID)
			if ( is.null( start)) start <- 1
			if ( is.null( stop)) stop <- seqMap$LENGTH[1]
			gmap <- subset.data.frame( geneMap, SEQ_ID == seqID & stop >= POSITION & start <= END)
			geneID <- gmap$GENE_ID[1]
			geneID <- shortGeneName( geneID, keep=1)
		} else {
			isGeneric <- TRUE
			if ( any( c( is.null(start), is.null(stop)))) stop( "Generic consensus sequence needs explicit 'start,stop'")
		}
	} else {
		isGene <- TRUE
		gmap <- subset.data.frame( geneMap, GENE_ID == geneID)
		if ( ! nrow(gmap)) {
			cat( "\nGiven GeneID not found in current species: ", geneID, "\n")
			return(NULL)
		}
		if (is.null(start)) start <- gmap$POSITION[1]
		if (is.null(stop)) stop <- gmap$END[1]
		seqID <- gmap$SEQ_ID[1]
	}

	# load that WIG data
	load( wigfile)

	# see what type it is
	dataType <- wiggles$Info$DataType
	if ( is.null( dataType)) dataType <- "RNA-seq"
	if ( dataType %in% c( "RNA-seq", "DNA-seq")) {
		trackNames <- c("Plus", "Minus", "PlusUnique", "MinusUnique")
		NTracks <- length( trackNames)
	}
	if ( dataType %in% c( "ChIP-seq")) {
		trackNames <- c("Plus", "Minus")
		NTracks <- length( trackNames)
	}
	if ( ! exists( "NTracks")) stop( "Unknown WIG data type")

	# get the wanted chromosome
	wiggleChunk <- WIG_getWigglesOneSeq( wiggles, seqID)
	if ( is.null( wiggleChunk)) {
		cat( "\nFailed to load Wiggle data for seqID: ", seqID)
		return(NULL)
	}

	# gather the set of all wiggle base depths in this range
	allBases <- start : stop
	NBases <- length(allBases)
	baseM <- matrix( 0, nrow=NBases, ncol=NTracks)
	rownames(baseM) <- allBases
	colnames(baseM) <- trackNames

	for ( i in 1:NTracks) {
		baseDepthTbl <- wiggleChunk[[ trackNames[i]]]
		if ( ! nrow( baseDepthTbl)) next

		# get the subset we want
		baseDepthTbl <- baseDepthTableSubset( baseDepthTbl, start, stop)
		if ( ! nrow( baseDepthTbl)) next
	
		# turn it into the 1-D vector of depths
		baseDepthVec <- baseDepthTableToVector( baseDepthTbl)

		# load those locations
		where <- match( rownames(baseM), names(baseDepthVec), nomatch=0)
		if ( any( where > 0)) baseM[ where > 0, i] <- baseDepthVec[ where]
	}

	if ( noDepth == "drop") {
		maxDepth <- apply( baseM, 1, max)
		toDrop <- which( maxDepth == 0)
		if ( length( toDrop)) baseM <- baseM[ -toDrop, ]
		if( ! nrow(baseM)) {
			out <- matrix( NA, nrow=0, ncol=6)
			colnames(out) <- c( "SEQ_ID", "POSITION", trackNames)
			return( as.data.frame(out))
		}
	}

	out <- data.frame( "SEQ_ID"=seqID, "POSITION"=rownames(baseM), baseM, stringsAsFactors=F)
	rownames(out) <- 1:nrow(out)
	out
}
