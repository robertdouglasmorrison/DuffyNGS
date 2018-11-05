# buildSpliceJunctionFiles.R

`buildSpliceJunctionFiles` <- function( speciesID, genomicFastaFile, outPath=".", 
					spliceMapPrefix="spliceMap", fragmentSize=60, 
					maxExonSkip=2, intronSplicesToo=TRUE, verbose=TRUE) {
	
	require( Biostrings)

	# load that species
	setCurrentSpecies( speciesID=speciesID)
	outfilePrefix <- getCurrentSpeciesFilePrefix()
	outfileFasta <- file.path( outPath, paste( outfilePrefix, spliceMapPrefix, "fasta", sep="."))
	outfileMap <- file.path( outPath, paste( outfilePrefix, spliceMapPrefix, "txt", sep="."))

	# remove any old file copies
	file.delete( outfileFasta)
	file.delete( outfileMap)

	# build a FASTA file and junctionMap of all the possible exon splicing combinations. 
	ans <- buildSpliceJunctions( genomicFastaFile, fragmentSize, maxExonSkip, intronSplicesToo)

	# the new FASTA is made from descriptors and the DNA fragments
	cat( "\nWriting splice files...")

	# write the FASTA
	writeLongFasta( ans$desc, ans$seq, file=outfileFasta)
	if ( verbose) cat( "\nWrote splice FASTA file:  ", outfileFasta, "\nN_Splices:  ", 
			length( ans$desc), "\n")

	# write the splice junction Map
	write.table( ans$spliceMap, file=outfileMap, sep="\t", quote=FALSE, row.names=FALSE)
	if ( verbose) cat( "\nWrote spliceMap file:     ", outfileMap, "\n")

	out <- list( "FastaFile"=outfileFasta, "MapFile"=outfileMap)
	return( out)
}


`buildSpliceJunctions` <- function( genomicFastaFile, fragmentSize=60, maxExonSkip=2, intronSplicesToo=FALSE) {

	seqMap <- getCurrentSeqMap()
	
	# pre-fetch one seq of genomic DNA to load that memory prior to the "multicore.lapply'
	curSeqDNA <- getFastaSeqFromFilePath( filePath=genomicFastaFile, seqID=seqMap$SEQ_ID[1])
	rm( curSeqDNA)

	# do each chromosome in parallel
	if ( nrow(seqMap) > 1) {
		ans <- try( multicore.lapply( seqMap$SEQ_ID, FUN=buildSpliceJunctionsOneSeqID, 
				genomicFastaFile=genomicFastaFile, fragmentSize=fragmentSize,
				maxExonSkip=maxExonSkip, intronSplicesToo=intronSplicesToo))
		names(ans) <- seqMap$SEQ_ID
		descOut <- seqOut <- vector()
		mapOut <- data.frame()
		cat( "\nCombining...\n")
		for ( i in 1:length(ans)) {
	
			if ( is.null( ans[[i]])) next
			cat( "\r", i, names(ans)[i])
			descOut <- c( descOut, ans[[i]]$desc)
			seqOut <- c( seqOut, ans[[i]]$seq)
			mapOut <- rbind( mapOut, ans[[i]]$spliceMap)
		}
		return( list( "desc"=descOut, "seq"=seqOut, "spliceMap"=mapOut))
	} else {
		ans <- buildSpliceJunctionsOneSeqID( seqMap$SEQ_ID[1], genomicFastaFile=genomicFastaFile, 
				fragmentSize=fragmentSize, maxExonSkip=maxExonSkip, intronSplicesToo=intronSplicesToo)
		return( ans)
	}
}


`buildSpliceJunctionsOneSeqID` <- function( SeqID, genomicFastaFile, fragmentSize=60, maxExonSkip=2, 
					intronSplicesToo=FALSE) {

	geneMap <- subset.data.frame( getCurrentGeneMap(), SEQ_ID == SeqID, 
					select=c( GENE_ID, POSITION, END, SEQ_ID, REAL_G))
	if ( nrow( geneMap) < 1) return(NULL)
	exonMap <- subset.data.frame( getCurrentExonMap(), SEQ_ID == SeqID, 
					select=c( GENE_ID, POSITION, END, STRAND))
	if ( nrow( exonMap) < 2) return(NULL)

	minFragSize <- 8

	cat( "\nMaking Splices:  \t", SeqID)
	# allocate enough space for all the columns, trim it later
	allocChunk <- ceiling( maxExonSkip * nrow(exonMap) / 1000) * 1000
	if ( allocChunk < 1000) allocChunk <- 1000
	allocSize <- allocChunk
	cat( "\nInitial splice memory allocation: ", allocSize)
	desc <- gid <- splid <- seqID <- seq <- exSize <- exBeg <- exEnd <- vector( mode="character", 
			length=allocSize)

	curSeqID <- ""
	curSeqDNA <- ""
	veryLastBase <- 0
	nout <- 0

	# visit every gene
	hasNonGenes <- ("REAL_G" %in% colnames( geneMap))
	for ( ig in 1:nrow(geneMap)) {

		# only real genes
		if ( hasNonGenes && geneMap$REAL_G[ig] == FALSE) next

		# watch for change of chromosome
		gene <- geneMap$GENE_ID[ ig]
		seqid <- geneMap$SEQ_ID[ ig]
		if ( seqid != curSeqID) {
			curSeqID <- seqid
			curSeqDNA <- getFastaSeqFromFilePath( filePath=genomicFastaFile, seqID=seqid)
			veryLastBase <- nchar( curSeqDNA)
		}

		# less than 2 exons: nothing to do
		exons <- which( exonMap$GENE_ID == gene)
		NX <- length( exons)
		if ( NX < 2) next

		if ( (ig > 0) && (nout > 0) && ((ig %% 1000 == 0) || (nout %% 1000 == 0))) 
			cat( "\n{N_genes, gene, N_Splices} = ", ig, gene, nout)

		if ( nout > (allocSize * 0.95)) {
			allocSize <- allocSize + allocChunk
			cat( "\nExtending splice memory allocation to: ", allocSize, "\n")
			length( desc) <- length( gid) <- length( splid) <- length( seqID) <- allocSize
			length( seq) <- length( exSize) <- length( exBeg) <- length( exEnd) <- allocSize
		}

		# grab the DNA for this one gene
		gBegin <- geneMap$POSITION[ig]
		gEnd <- geneMap$END[ig]
		if ( gEnd > veryLastBase) gEnd <- veryLastBase
		geneDNA <- strsplit( as.character( base::substr( curSeqDNA, gBegin, gEnd)), split="")[[1]]
		goffset <- gBegin - 1


		# force the exons to be in genomic order
		exonPOSITION <- exonMap$POSITION[ exons]
		ord <- order( exonPOSITION)
		exons <- exons[ ord]

		# make the exon fragments that we can
		exonPOSITION <- exonMap$POSITION[ exons]
		exonEND <- exonMap$END[ exons]
		exonLen <- exonEND - exonPOSITION + 1
		exonName <- if ( exonMap$STRAND[ exons[1]] == "+") 1:NX else NX:1
		exonFragEnd <- exonPOSITION + base::pmin( fragmentSize, exonLen) - 1
		exonFragSize <- exonFragEnd - exonPOSITION + 1
		exonFragSeq <- vector( length=NX)
		for (j in 1:NX) exonFragSeq[j] <- base::paste( geneDNA[ (exonPOSITION[j]-goffset) : (exonFragEnd[j]-goffset)], 
							collapse="")

		# make explicit introns too
		intronPOSITION <- exonEND[1:(NX-1)] + 1
		intronEND <- exonPOSITION[2:NX] - 1
		NIN <- NX - 1
		intronLen <- intronEND - intronPOSITION + 1
		intronName <- if ( exonMap$STRAND[ exons[1]] == "+") 1:NIN else NIN:1
		intronFragEnd <- intronPOSITION + base::pmin( fragmentSize, intronLen) - 1
		intronFragSize <- intronFragEnd - intronPOSITION + 1
		intronFragSeq <- vector( length=NIN)
		for (j in 1:NIN) intronFragSeq[j] <- base::paste( geneDNA[ (intronPOSITION[j]-goffset) : (intronFragEnd[j]-goffset)], 
							collapse="")

		# visit all in-order pairs of exons, till too many are skipped
		for ( ix in 1:(NX-1)) {
			# determine the fragments
		    end1 <- exonEND[ix]
		    from1 <- exonEND[ix] - min( c(fragmentSize, exonLen[ix])) + 1
		    siz1 <- end1 - from1 + 1
		    if ( siz1 < minFragSize) next
		    frag1 <- base::paste( geneDNA[ (from1-goffset) : (end1-goffset)], collapse="")

		    # loop over down stream exons
		    for ( ix2 in (ix+1):NX) {
			if ( (ix2 - ix - 1) > maxExonSkip) break
			from2 <- exonPOSITION[ix2]
			end2 <- exonFragEnd[ix2]
			siz2 <- exonFragSize[ix2]
		        if ( siz2 < minFragSize) next
			frag2 <- exonFragSeq[ix2]

			# add it to the set
			nout <- nout + 1
			type <- if ( ix2-ix == 1) "std" else "alt"
			spliceTxt <- base::paste( type, ":E", exonName[ix], "_E", exonName[ix2], sep="")
			splid[nout] <- spliceTxt
			gid[nout] <- gene
			seqID[nout] <- seqid
			desc[nout] <- base::paste( seqid, gene, spliceTxt, sep="::")
			seq[nout] <- base::paste( frag1, frag2, sep="")
			exBeg[nout] <- base::paste( from1, from2, sep=",")
			exEnd[nout] <- base::paste( end1, end2, sep=",")
			exSize[nout] <- base::paste( siz1, siz2, sep=",")
		    }
		
		    if ( ix >= NIN) next
		    # loop over down stream introns
		    if ( ! intronSplicesToo) next

		    for ( ix2 in (ix+1):NIN) {
			if ( (ix2 - ix - 1) > maxExonSkip) break
			from2 <- intronPOSITION[ix2]
			end2 <- intronFragEnd[ix2]
			siz2 <- intronFragSize[ix2]
		        if ( siz2 < minFragSize) next
			frag2 <- intronFragSeq[ix2]

			# add it to the set
			nout <- nout + 1
			type <- "alt"
			spliceTxt <- base::paste( type, ":E", exonName[ix], "_I", intronName[ix2], sep="")
			splid[nout] <- spliceTxt
			gid[nout] <- gene
			seqID[nout] <- seqid
			desc[nout] <- base::paste( seqid, gene, spliceTxt, sep="::")
			seq[nout] <- base::paste( frag1, frag2, sep="")
			exBeg[nout] <- base::paste( from1, from2, sep=",")
			exEnd[nout] <- base::paste( end1, end2, sep=",")
			exSize[nout] <- base::paste( siz1, siz2, sep=",")
		    }
		}
		
		# now visit all in-order pairs from introns!, till too many are skipped
		if ( NX < 3) next
		if ( ! intronSplicesToo) next
		
		for ( ix in 1:(NIN-1)) {
			# determine the fragments
		    end1 <- intronEND[ix]
		    from1 <- intronEND[ix] - min( c(fragmentSize, intronLen[ix])) + 1
		    siz1 <- end1 - from1 + 1
		    if ( siz1 < minFragSize) next
		    frag1 <- base::paste( geneDNA[ (from1-goffset) : (end1-goffset)], collapse="")

		    # loop over down stream exons
		    for ( ix2 in (ix+2):NX) {
			if ( (ix2 - ix - 1) > maxExonSkip) break
			from2 <- exonPOSITION[ix2]
			end2 <- exonFragEnd[ix2]
			siz2 <- exonFragSize[ix2]
		        if ( siz2 < minFragSize) next
			frag2 <- exonFragSeq[ix2]

			# add it to the set
			nout <- nout + 1
			type <- "alt"
			spliceTxt <- base::paste( type, ":I", intronName[ix], "_E", exonName[ix2], sep="")
			splid[nout] <- spliceTxt
			gid[nout] <- gene
			seqID[nout] <- seqid
			desc[nout] <- base::paste( seqid, gene, spliceTxt, sep="::")
			seq[nout] <- base::paste( frag1, frag2, sep="")
			exBeg[nout] <- base::paste( from1, from2, sep=",")
			exEnd[nout] <- base::paste( end1, end2, sep=",")
			exSize[nout] <- base::paste( siz1, siz2, sep=",")
		    }

		    # loop over down stream introns
		    for ( ix2 in (ix+1):NIN) {
			if ( (ix2 - ix - 1) > maxExonSkip) break
			from2 <- intronPOSITION[ix2]
			end2 <- intronFragEnd[ix2]
			siz2 <- intronFragSize[ix2]
		        if ( siz2 < minFragSize) next
			frag2 <- intronFragSeq[ix2]

			# add it to the set
			nout <- nout + 1
			type <- "alt"
			spliceTxt <- base::paste( type, ":I", intronName[ix], "_I", intronName[ix2], sep="")
			splid[nout] <- spliceTxt
			gid[nout] <- gene
			seqID[nout] <- seqid
			desc[nout] <- base::paste( seqid, gene, spliceTxt, sep="::")
			seq[nout] <- base::paste( frag1, frag2, sep="")
			exBeg[nout] <- base::paste( from1, from2, sep=",")
			exEnd[nout] <- base::paste( end1, end2, sep=",")
			exSize[nout] <- base::paste( siz1, siz2, sep=",")
		    }
		}
		
	}

	# there is a tiny chance that we have none...
	if ( nout < 1) return(NULL)


	# trim to actual length
	cat( "\nOrganizing...")
	length( gid) <- nout
	length( splid) <- nout
	length( seqID) <- nout
	length( desc) <- nout
	length( seq) <- nout
	length( exSize) <- nout
	length( exBeg) <- nout
	length( exEnd) <- nout

	outDF <- data.frame( gid, splid, seqID, exSize, exBeg, exEnd, stringsAsFactors=FALSE)
	colnames( outDF) <- SPLICEMAP_COLUMNS
	rownames( outDF) <- 1:nrow(outDF)

	return( list( "desc"=desc, "seq"=seq, "spliceMap"=outDF))
}

