# pipe.ExtractRiboClearedReadCounts.R -- mine the RiboClear BAM file to get cleared read counts

`pipe.ExtractRiboClearedReadCounts` <- function( sampleIDset, optionsFile="Options.txt", 
				speciesID=getCurrentSpecies(), results.path=NULL, verbose=TRUE) {

	NS <- length( sampleIDset)
	# set up for one species
	if (speciesID != getCurrentSpecies()) setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()
	sMap <- getCurrentSeqMap()
	rrnaMap <- subset( getCurrentRrnaMap(), SEQ_ID %in% sMap$SEQ_ID)

	# get the genes we will count, and set up storage
	rrnaGenes <- sort( rrnaMap$GENE_ID)
	NG <- length(rrnaGenes)
	countM <- matrix( 0, nrow=NS, ncol=NG)
	rownames(countM) <- sampleIDset
	colnames(countM) <- rrnaGenes

	for ( i in 1:NS) {
	
		# get all the alignments from the BAM file
		s <- sampleIDset[i]
		tmp <- pipe.GatherGeneAlignments( s, genes=rrnaGenes, optionsFile=optionsFile,
		 		results.path=results.path, stages="riboClear", mode="all",
		 		asFASTQ=F, verbose=FALSE)
		 		
		# we want read counts, not alignments, so only keep one instance of each
		keep <- which( ! duplicated( tmp$name))
		tmp <- tmp[ keep, ]
		
		# gather those gene counts
		cntsTbl <- table( factor( tmp$geneid, levels=rrnaGenes))
		countM[ i, ] <- as.numeric( cntsTbl)
		if (verbose) cat( "\n", s, "  \tRibo Reads: ", sum(cntsTbl))
	}

	out <- data.frame( "SampleID"=sampleIDset, countM, stringsAsFactors=FALSE)
	rownames(out) <- 1:nrow(out)
	return( out)
}
