# buildRiboClearingFiles.R

`buildRiboClearingFiles` <- function( speciesID, genomicFastaFile, outPath=".", mapFile="riboMap.txt",
				tailSize=24, verbose=TRUE) {
	
	# load that species
	setCurrentSpecies( speciesID=speciesID)
	outfilePrefix <- getCurrentSpeciesFilePrefix()
	outfileFasta <- file.path( outPath, base::paste( outfilePrefix, "riboClear.fasta", sep="."))
	outfileMap <- file.path( outPath, mapFile)
	
	# remove any old files
	file.delete( outfileFasta)
	file.delete( outfileMap)

	ans <- buildRiboClearData( speciesID, genomicFastaFile, tailSize)

	# the new FASTA is made from descriptors and the DNA fragments
	cat( "\nWriting riboClearing files...")

	# write the FASTA
	writeFasta( as.Fasta( ans$desc, toupper(ans$seq)), file=outfileFasta)
	if ( verbose) cat( "\nWrote riboClear FASTA file:  ", outfileFasta, "\nN_riboClear genes:  ", length( ans$desc), "\n")

	# write the riboClear Map
	write.table( ans$riboclearMap, file=outfileMap, sep="\t", quote=FALSE, row.names=FALSE)
	if ( verbose) cat( "\nWrote riboClear Map file:   ", outfileMap, "\n")

	out <- list( "FastaFile"=outfileFasta, "MapFile"=outfileMap)
	return( out)
}


`buildRiboClearData` <- function( speciesID, genomicFastaFile, tailSize=24) {

	rrnaMap <- getCurrentRrnaMap()
	exonMap <- getCurrentExonMap()

	allocSize <- nrow(rrnaMap) * 2
	desc <- gid <- seqID <- seq <- grp <- seqBeg <- seqEnd <- faBeg <- faEnd <- vector( mode="character", 
			length=allocSize)

	curSeqID <- ""
	curSeqDNA <- ""
	nout <- 0

	# visit every ribo RNA entry
	cat( "\nbuilding riboClear data...")

	if ( ! ("CLEAR" %in% colnames( rrnaMap))) {
		cat( "\nNo 'CLEAR' column found in riboMap... flagging all ribo genes for clearing.")
		rrnaMap$CLEAR <- rep( TRUE, times=nrow(rrnaMap))
	}

	for ( ig in 1:nrow(rrnaMap)) {

		# only use entries that are flagged for clearing
		if ( ! as.logical( rrnaMap$CLEAR[ig])) next

		# watch for change of chromosome
		gene <- rrnaMap$GENE_ID[ ig]
		seqid <- rrnaMap$SEQ_ID[ ig]
		if ( seqid != curSeqID) {
			curSeqID <- seqid
			curSeqDNA <- getFastaSeqFromFilePath( filePath=genomicFastaFile, seqID=seqid)
			emap <- subset( exonMap, SEQ_ID == seqid)
		}

		if ( (ig %% 10 == 0) || (nout %% 10 == 0)) cat( "\n{N_genes, Gene, N_clears} = ", ig, gene, nout)

		# grab the DNA for this one gene, plus the tails...
		gBegin <- rrnaMap$POSITION[ig]
		gEnd <- rrnaMap$END[ig]
		gGroup <- rrnaMap$GROUP[ig]
		geneDNA <- as.character( base::substr( curSeqDNA, (gBegin-tailSize), (gEnd+tailSize)))

		# add it to the set
		nout <- nout + 1
		gid[nout] <- gene
		seqID[nout] <- seqid
		desc[nout] <- base::paste( seqid, gene, gGroup, sep="::")
		seq[nout] <- geneDNA
		grp[nout] <- gGroup
		seqBeg[nout] <- gBegin
		seqEnd[nout] <- gEnd
		faBeg[nout] <- tailSize + 1
		faEnd[nout] <- tailSize + 1 + (gEnd - gBegin)

		# if this gene has exon splices, build that construct too
		exPtrs <- which( emap$GENE_ID == gene)
		if ( length( exPtrs) < 2) next

		geneDNA <- ""
		gGroup <- paste( gGroup, "splice", sep="_")
		for ( j in exPtrs) {
			smallDNA <- as.character( base::substr( curSeqDNA, emap$POSITION[j], emap$END[j]))
			# when the exon is long enough, 'N' out the interior so we end up with just the splice boundaries
			firstN <- tailSize + 5
			lastN <- nchar( smallDNA) - tailSize - 5
			if ( (nN <- lastN - firstN + 1) > 0) {
				filler <- paste( rep.int( "N", nN), collapse="")
				substr( smallDNA, firstN, lastN) <- filler
			}
			geneDNA <- paste( geneDNA, smallDNA, sep="")
		}
		nout <- nout + 1
		gid[nout] <- gene
		seqID[nout] <- seqid
		desc[nout] <- base::paste( seqid, gene, gGroup, sep="::")
		seq[nout] <- geneDNA
		grp[nout] <- gGroup
		seqBeg[nout] <- emap$POSITION[ exPtrs[1]]
		seqEnd[nout] <- emap$END[ exPtrs[ length(exPtrs)]]
		faBeg[nout] <- 1
		faEnd[nout] <- nchar( geneDNA)
	}
	cat( "\nOrganizing...")

	# trim to actual length
	length( gid) <- nout
	length( seqID) <- nout
	length( desc) <- nout
	length( seq) <- nout
	length( grp) <- nout
	length( seqBeg) <- nout
	length( seqEnd) <- nout
	length( faBeg) <- nout
	length( faEnd) <- nout

	outDF <- data.frame( gid, grp, seqID, seqBeg, seqEnd, faBeg, faEnd, stringsAsFactors=FALSE)
	colnames( outDF) <- RIBOCLEAR_COLUMNS
	rownames( outDF) <- 1:nrow(outDF)

	return( list( "desc"=desc, "seq"=seq, "riboclearMap"=outDF))
}

