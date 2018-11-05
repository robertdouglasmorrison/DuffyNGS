# genericTranscriptome.R -- make files that look like NGS transcriptomes from generic expression data


`pipe.GenericTranscriptomes` <- function( m, annotationFile="Annotation.txt", optionsFile="Options.txt", 
				speciesID=NULL, results.path=NULL, readsPerSample=1000000, verbose=TRUE) {

	annT <- readAnnotationTable( annotationFile)
	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	if ( ! file.exists( results.path)) dir.create( results.path, recursive=TRUE, showWarnings=FALSE)
	pathOut <- file.path( results.path, "transcript")
	if ( ! file.exists( pathOut)) dir.create( pathOut, showWarnings=FALSE)

	if ( is.null(speciesID)) speciesID <- getCurrentSpecies()
	speciesPrefix <- getCurrentSpeciesFilePrefix()

	# sample IDs are the column names
	sampleIDs <- colnames(m)
	NS <- ncol(m)
	NG <- nrow(m)

	gids <- rownames(m)
	prods <- gene2Product( gids)

	# try to find/estimate a number of exonic bases
	gmap <- getCurrentGeneMap()
	nexonb <- rep.int( 1000, NG)
	wh <- match( gids, gmap$GENE_ID, nomatch=0)
	wh2 <- match( gids, gmap$NAME, nomatch=0)
	use2 <- which( wh == 0 & wh2 > 0)
	wh[use2] <- wh2[use2]
	nexonb[ wh > 0] <- gmap$N_EXON_BASES[wh]

	# also we can give 'real' IDs to the genes if we want?...
	# not yet done...

	fileOutTrans <- paste( sampleIDs, speciesPrefix, "Transcript.txt", sep=".")
	fileOutTrans <- file.path( pathOut, fileOutTrans)

	for ( i in 1:NS) {

		# grab each column, treat the value given as our RPKM, and generate a proxy for N_READS
		v <- m[ ,i]
		v[ is.na(v)] <- 0

		rpkm <- v
		reads <- rpkm * readsPerSample / sum(rpkm)
		sigma <- rosettaSigma( reads)
		sml <- data.frame( "GENE_ID"=gids, "PRODUCT"=prods,
				"RPKM_M"=rpkm, "READS_M"=reads, "SIGMA_M"=sigma, "STRAND_M"=1,
				"RPKM_U"=rpkm, "READS_U"=reads, "SIGMA_U"=sigma, "STRAND_U"=1,
				"N_EXON_BASES"=nexonb, stringsAsFactors=F)
		ord <- order( rpkm, decreasing=T)
		sml <- sml[ ord, ]

		write.table( sml, fileOutTrans[i], sep="\t", quote=F, row.names=F)
		cat( "\r", i, sampleIDs[i])
	}
	return( sampleIDs)
}

