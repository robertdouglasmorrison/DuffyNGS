# pipe.MARsurvey.R

pipe.MARsurvey <- function( sampleID, geneID=NULL, seqID=NULL, start=NULL, stop=NULL, speciesID=getCurrentSpecies(),
				annotationFile="Annotation.txt", optionsFile="Options.txt", 
				results.path=NULL, maxAlignments=NULL, mode=c("genomic","riboClear"),
				short.gene.names=TRUE) {
				
	# make sure we can see all species we aligned against
	setCurrentTarget( optionsFile=optionsFile)
	speciesSet <- getCurrentTargetSpecies()
	setCurrentSpecies( speciesID)

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}

	# which BAM file depends on the mode
	mode <- match.arg( mode)
	if ( mode == "genomic") {
		bamfile <- file.path( results.path, "align", paste( sampleID, "genomic.bam", sep="."))
	} else {
		bamfile <- file.path( results.path, "riboClear", paste( sampleID, "ribo.converted.bam", sep="."))
	}
	if ( ! file.exists( bamfile)) {
		cat( "\nUnsorted BAM file not found:  ", bamfile)
		cat( "\nMAR alignments only detectable from pre-sorted BAM file...")
		return(NULL)
	}

	# specify by either gene or range
	geneMap <- getCurrentGeneMap()
	seqMap <- getCurrentSeqMap()
	if ( is.null( geneID)) {
		if ( is.null(seqID)) stop( "One of 'seqID','geneID' must be not NULL")
		if ( seqID %in% seqMap$SEQ_ID) {
			gmap <- subset.data.frame( geneMap, SEQ_ID == seqID & stop >= POSITION & start <= END)
			geneID <- gmap$GENE_ID[1]
			geneID <- shortGeneName( geneID, keep=1)
		}
	} else {
		gmap <- subset.data.frame( geneMap, GENE_ID == geneID)
		if ( nrow(gmap) == 0) gmap <- subset.data.frame( geneMap, NAME == geneID)
		if ( nrow(gmap) == 0) {
			cat( "\nNo Gene with that ID:  ", geneID)
			return(NULL)
		}
		if (is.null(start)) start <- gmap$POSITION[1]
		if (is.null(stop)) stop <- gmap$END[1]
		seqID <- gmap$SEQ_ID[1]
	}

	# how to interpret depends on whether we started with PairedEnd data or not
	pairedEnd <- getAnnotationTrue( annotationFile, key=sampleID, columnArg="PairedEnd", notfound=FALSE)

	chunkSize <- as.numeric( getOptionValue( optT, "readBufferSize", notfound="1000000", verbose=T))

	reader <- bamReader( bamfile)
	on.exit( bamClose( reader))
	refData <- getRefData( reader)
	out <- data.frame()
	nreads <- 0
	cat( "\nScanning BAM file for MAR alignments..\n")

	WHICH <- base::which
	PASTE <- base::PASTE

	repeat {
		if ( !is.null( maxAlignments) && nreads >= maxAlignments) break
		chunk <- getNextChunk( reader, n=chunkSize)
		nreads <- nreads + size(chunk)
		if ( size(chunk) < 1) break
		thisDF <- as.data.frame(chunk)
		refid <- thisDF$refid
		pos <- thisDF$position
		sid <- refID2seqID( refid, reader=reader, refData=refData)
	
		hitTarget <- WHICH( sid == seqID & pos >= start & pos <= stop)
		if ( length(hitTarget) < 1) {
			cat(".")
			next
		}
	
		# we got at least one read in this region
		myRIDs <- thisDF$name[ hitTarget]
		# keep all the alignments that have these readIDs
		keep <- WHICH( thisDF$name %in% myRIDs)
		smlDF <- thisDF[ keep, ]
		smlSID <- sid[ keep]
		out <- rbind( out, cbind( "SEQ_ID"=smlSID, smlDF[,c(1:5,8,9)], stringsAsFactors=FALSE))
		cat( "\rAlignments:", formatC(nreads,format="d",big.mark=","), "  Hits:", nrow(out))
	}
	cat( "  Done scanning.\n")

	# augment these alignments with annotation facts
	ans <- fastSP2GP( seqid=out$SEQ_ID, seqbase=out$position)
	out$GENE_ID <- ans$GENE_ID
	out$PRODUCT <- gene2ProductAllSpecies(ans$GENE_ID)
	isblank <- which( out$PRODUCT == "")
	out$PRODUCT[isblank] <- out$GENE_ID[isblank]
	if ( short.gene.names) out$GENE_ID <- shortGeneName( out$GENE_ID, keep=1)
	out <- out[ order( out$SEQ_ID, out$position), ]
	rownames(out) <- 1:nrow(out)

	# visit each MAR, and assess
	cat( "\nAssessing MARs..\n")
	idFac <- factor( out$name)
	seqText <- productText <- geneText <- nAlign <- vector()
	nout <- 0
	tapply( 1:nrow(out), idFac, function(x) { 
		nout <<- nout + 1
		nAlign[nout] <<- length(x)
		geneSet <- out$GENE_ID[x]
		if (pairedEnd) {
			geneSet <- unique.default( geneSet)
			nAlign[nout] <<- length(geneSet)
		}
		geneText[nout] <<- PASTE( geneSet, collapse="; ")
		productText[nout] <<- PASTE( unique.default(out$PRODUCT[x]), collapse="; ")
		seqText[nout] <<- PASTE( unique.default(out$SEQ_ID[x]), collapse="; ")
		return(NULL)
	})
	cat( "  Done.\n")

	# summarize
	gTextFac <- factor( geneText)
	nGtext <- nlevels(gTextFac)
	geneCnt <- tapply( geneText, gTextFac, length)
	geneTbl <- tapply( geneText, gTextFac, function(x) x[1])
	productTbl <- tapply( productText, gTextFac, function(x) x[1])
	seqTbl <- tapply( seqText, gTextFac, function(x) x[1])
	alignPerTbl <- tapply( nAlign, gTextFac, function(x) round( sum(x)/length(x), digits=2))
	genePcts <- round( geneCnt * 100 / sum(geneCnt), digits=2)

	out <- data.frame( "N_Reads"=as.vector(geneCnt), "Pct_Reads"=as.vector(genePcts), 
			"AlignsPerRead"= alignPerTbl, "SeqList"=as.vector(seqTbl),
			"GeneList"=as.vector(geneTbl), "ProductList"=as.vector(productTbl), 
			stringsAsFactors=F)
	ord <- order( out$N_Reads, decreasing=T)
	out <- out[ ord, ]
	rownames(out) <- 1:nrow(out)
	
	if ( pairedEnd) {
		isMAR <- (out$AlignsPerRead > 1.05)
		out <- data.frame( "Type"=ifelse( isMAR, "MAR", "UAR"), out, stringsAsFactors=F)
		rownames(out) <- 1:nrow(out)
	}

	return( out)
}
