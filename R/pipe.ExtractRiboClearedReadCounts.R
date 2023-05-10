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
	countM <- matrix( 0, nrow=NG, ncol=NS)
	colnames(countM) <- sampleIDset
	rownames(countM) <- rrnaGenes

	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	
	for ( i in 1:NS) {
	
		# get all the alignments from the BAM file
		s <- sampleIDset[i]
		
		# the ribo conversion step seems to misplaced some reads, such that 'position' may
		# be off.  Use a different barebones method for now
		#tmp <- pipe.GatherGeneAlignments( s, genes=rrnaGenes, optionsFile=optionsFile,
		# 		results.path=results.path, stages="riboClear", mode="all",
		# 		asFASTQ=F, verbose=FALSE) 
		# we want read counts, not alignments, so only keep one instance of each
		#keep <- which( ! duplicated( tmp$name))
		#tmp <- tmp[ keep, ]
		# gather those gene counts
		#cntsTbl <- table( factor( tmp$geneid, levels=rrnaGenes))
		#countM[ , i] <- as.numeric( cntsTbl)
		
		# read the BAM file directly
		bamfile <- file.path( results.path, "riboClear", paste( s, ".ribo.converted.bam", sep=""))
		bamfile <- BAM.verifySorted( bamfile, verbose=verbose, threads=2)
		br <- bamReader( bamfile)
		
		# grab gene ID tag elements in chunks, keeping only one read per MAR
		geneTags <- vector()
		nRead <- 0
		while (TRUE) { 
			ch <- getNextChunk( br, n=1000000); 
			if (size(ch) < 1) break; 
			nRead <- nRead + size(ch);
			tags <- getTag(ch, "GI"); 
			# primary alignments only
			tags <- tags[ ! secondaryAlign(ch)]; 
			# only keep those in this species
			tags <- tags[ tags %in% rrnaGenes]
			geneTags <- c( geneTags, tags); 
			if (verbose) cat( "\r", s, nRead, "\tRibo Reads: ", length(geneTags))
		}
		# done with the BAM file
		bamClose(br)
		cntsTbl <- table( factor( geneTags, levels=rrnaGenes))
		countM[ , i] <- as.numeric( cntsTbl)
		if (verbose) cat( "\n")
	}

	out <- data.frame( "GENE_ID"=rrnaGenes, "PRODUCT"=gene2Product(rrnaGenes), countM, 
			stringsAsFactors=FALSE)
	rownames(out) <- 1:nrow(out)
	return( out)
}
