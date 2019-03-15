# pipe.BAMpileup.R -- get the mpileup data from a BAM file

`pipe.BAMpileup` <- function( sampleID, geneID=NULL, seqID=NULL, start=NULL, stop=NULL, 
				annotationFile="Annotation.txt", optionsFile="Options.txt", 
				max.depth=10000, results.path=NULL, summarize.calls=TRUE, 
				verbose=FALSE) {
				
	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	genomicFastaFile <- getOptionValue( optT, "genomicFastaFile", verbose=F)

	bamfile <- file.path( results.path, "align", paste( sampleID, "genomic.bam", sep="."))

	# let's be a bit more flexible with what we pass in
	# Could be a genomic range {seqID, start, stop}, a gene {geneID}, or any generic sequence...
	geneMap <- getCurrentGeneMap()
	exonMap <- getCurrentExonMap()
	seqMap <- getCurrentSeqMap()

	# specify by either gene or range
	if ( is.null( geneID)) {
		if ( is.null(seqID)) stop( "One of 'seqID','geneID' must be not NULL")
		if ( seqID %in% seqMap$SEQ_ID) {
			gmap <- subset.data.frame( geneMap, SEQ_ID == seqID & stop >= POSITION & start <= END)
			emap <- subset.data.frame( exonMap, GENE_ID %in% gmap$GENE_ID)
			geneID <- gmap$GENE_ID[1]
			geneID <- shortGeneName( geneID, keep=1)
		} else {
			stop( "SeqID not in current species")
		}
	} else {
		gmap <- subset.data.frame( geneMap, GENE_ID == geneID)
		emap <- subset.data.frame( exonMap, GENE_ID == geneID)
		if ( nrow(gmap) != 1) stop( paste( "GeneID not in current species: ", geneID))
		if (is.null(start)) start <- gmap$POSITION[1]
		if (is.null(stop)) stop <- gmap$END[1]
		seqID <- gmap$SEQ_ID[1]
	}

	# load that portion of the PILEUPS
	curMPU <- BAM.mpileup( bamfile, seqID=seqID, fastaFile=genomicFastaFile, start=start, stop=stop, 
			max.depth=max.depth, summarize.calls=summarize.calls, verbose=verbose)

	return( curMPU)
}
