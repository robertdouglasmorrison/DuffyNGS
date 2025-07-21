# calcWigTranscriptome.R

# turn a Wiggles object into a transcriptome.   Converting from WB to 'wiggles' data strucure...

`calcWigTranscriptome` <- function( WIG, geneMap=NULL, useBothStrands=FALSE, keepIntergenics=FALSE, 
			exonsOnly=FALSE, fileout=NA, minReadsPerSpecies=1000, verbose=!interactive() ) {

	setCurrentSpecies( WIG$Info$Species)

	if ( is.null( geneMap)) {
		geneMap <- getCurrentGeneMap()
	} else {
		sameSpecies <- all( unique.default( geneMap$SEQ_ID) %in% getCurrentSeqMap()$SEQ_ID)
		if ( ! sameSpecies) stop( "calcTranscriptome:  given 'geneMap' does not go with current Species")
	}

	cat( "\n\nCalculating Transcriptome:");
	cat( "\nSpecies:         \t", getCurrentSpecies());
	cat( "\nAlignment files: \n  ", paste( WIG$Info$FileName, collapse="\n  "), "\n", sep="")

	# allocate vectors for result
	geneSet <- gProdSet <- vector( length=nrow( geneMap))
	totCntSet <- sigmaSet <- rpkmSet <- strandSet <- vector( length=nrow( geneMap))
	totCntSetM <- sigmaSetM <- rpkmSetM <- strandSetM <- vector( length=nrow( geneMap))
	nBaseSet <- vector( length=nrow( geneMap))

	# pre evaluate what we can
	hasNonGenes <- ("REAL_G" %in% colnames( geneMap))
	totalWIGreads <- WIG_getTotalReadsForRPKM( WIG, minReadsPerSpecies=minReadsPerSpecies)
	# catch the very rare case of no reads
	if ( totalWIGreads$Unique < 1) totalWIGreads$Unique <- 1
	if ( totalWIGreads$Multi < 1) totalWIGreads$Multi <- 1
	readLength <- WIG$Info$ReadLength
	if ( is.na(readLength) || readLength < 32) readLength <- 32

	# visit every gene
	ngenes <- nTranscribed <- 0
	curSeq <- ""

	for( ig in 1:nrow( geneMap)) {

		# bypass this one ?
		if ( !keepIntergenics && hasNonGenes && ( geneMap$REAL_G[ig] == FALSE)) next

		gene <- geneMap$GENE_ID[ig]
		seqID <- geneMap$SEQ_ID[ig]
		if ( seqID != curSeq) {
			wiggleChunk <- WIG_getWigglesOneSeq( WIG, seqID)
			curSeq <- seqID
		}

		ans <- calcWigTranscriptOneGene( wiggleChunk, gene, gmapPtr=ig, geneMap=geneMap, 
				totalReads=totalWIGreads, readLength=readLength,
				useBothStrands=useBothStrands, exonsOnly=exonsOnly)
		if ( is.null(ans)) next

		# load the vectors
		ngenes <- ngenes + 1
		geneSet[ ngenes] <- gene
		gProdSet[ ngenes] <- geneMap$PRODUCT[ig]
		rpkmSet[ ngenes] <- ans$rpkm
		totCntSet[ ngenes] <- ans$rawReads
		sigmaSet[ ngenes] <- ans$sigma
		strandSet[ ngenes] <- ans$strandness
		rpkmSetM[ ngenes] <- ans$rpkm.Multi
		totCntSetM[ ngenes] <- ans$rawReads.Multi
		sigmaSetM[ ngenes] <- ans$sigma.Multi
		strandSetM[ ngenes] <- ans$strandness.Multi
		nBaseSet[ ngenes] <- ans$nBases

		if (ngenes %% 1000 == 0) cat( "\n", ngenes, "\t", gene,"\t", ans$rawReads)

		if ( ans$rpkm.Multi > 0) nTranscribed <- nTranscribed + 1
	}

	# now trim to the true size
	length(geneSet) <- length(gProdSet) <- length(nBaseSet) <- ngenes
	length(rpkmSetM) <- length(totCntSetM) <- length(sigmaSetM) <- length(strandSetM) <- ngenes
	length(rpkmSet) <- length(totCntSet) <- length(sigmaSet) <- length(strandSet) <- ngenes

	# we can now do the TPM metric of gene expression
	tpmM <- tpm( totCntSetM, nBaseSet, minReadsPerSpecies=minReadsPerSpecies)
	tpmU <- tpm( totCntSet, nBaseSet, minReadsPerSpecies=minReadsPerSpecies)

	# calc gene rankings
	rnkM <- rank( tpmM, ties.method="min")
	rnkU <- rank( tpmU, ties.method="min")

	# trim to sensible digits resolution
	rpkmSetM <- round( rpkmSetM, digits=4)
	rpkmSet <- round( rpkmSet, digits=4)
	tpmM <- round( tpmM, digits=4)
	tpmU <- round( tpmU, digits=4)
	totCntSetM <- round( totCntSetM, digits=2)
	totCntSet <- round( totCntSet, digits=2)
	sigmaSetM <- round( sigmaSetM, digits=2)
	sigmaSet <- round( sigmaSet, digits=2)
	strandSetM <- round( strandSetM, digits=3)
	strandSet <- round( strandSet, digits=3)

	out <- data.frame( geneSet, gProdSet, rpkmSetM, totCntSetM, tpmM, rnkM, sigmaSetM, strandSetM, 
			rpkmSet, totCntSet, tpmU, rnkU, sigmaSet, strandSet, nBaseSet, stringsAsFactors=FALSE)
	colnames( out) <- c("GENE_ID", "PRODUCT", 
			"RPKM_M", "READS_M", "TPM_M", "RANK_M", "SIGMA_M", "STRAND_M", 
			"RPKM_U", "READS_U", "TPM_U", "RANK_U", "SIGMA_U", "STRAND_U", 
			"N_EXON_BASES")
			
	# now sort into intensity order
	ord <- base::order( out$TPM_M, decreasing=TRUE)
	out <- out[ ord, ]
	rownames( out) <- 1:nrow(out)
	
	if ( ! is.na( fileout)) {

		if ( getCurrentSpecies() %in% MAMMAL_SPECIES) out <- addHumanIDterms( out)
		#if ( getCurrentSpecies() %in% ORIGID_PARASITE_SPECIES) out <- addOrigIDterms( out)
		if ( getCurrentSpecies() %in% BACTERIA_SPECIES) out <- addNameIDterms( out)

		write.table( out, file=fileout, sep="\t", quote=FALSE, row.names=FALSE)
		cat( "\nWrote Transcriptome file:  \t", fileout)
	}
	
	if (verbose) {
		cat( "\nN_Gene regions processed:         \t", nrow(out))
		cat( "\nN_Genes regions with RPKM > 0:     \t", nTranscribed, "\n")
	}

	return( out)
}

