# buildRiboClearingFiles.R

`buildRiboClearingFiles` <- function( speciesID, genomicFastaFile, outPath=".", mapFile="riboMap.txt",
				tailSize=24, as.contigs=TRUE, max.contig.gap=1000, verbose=TRUE) {
	
	# load that species
	setCurrentSpecies( speciesID=speciesID)
	outfilePrefix <- getCurrentSpeciesFilePrefix()
	outfileFasta <- file.path( outPath, base::paste( outfilePrefix, "riboClear.fasta", sep="."))
	outfileMap <- file.path( outPath, mapFile)
	
	# remove any old files
	file.delete( outfileFasta)
	file.delete( outfileMap)

	ans <- buildRiboClearData( speciesID, genomicFastaFile, tailSize, as.contigs=as.contigs, 
				max.contig.gap=max.contig.gap)

	# the new FASTA is made from descriptors and the DNA fragments
	cat( "\nWriting riboClearing files...")

	# write the FASTA
	writeFasta( as.Fasta( ans$desc, toupper(ans$seq)), file=outfileFasta)
	if ( verbose) cat( "\nWrote riboClear FASTA file:  ", outfileFasta, "\nN_riboClear contigs:  ", length( ans$desc), "\n")

	# write the riboClear Map
	write.table( ans$riboclearMap, file=outfileMap, sep="\t", quote=FALSE, row.names=FALSE)
	if ( verbose) cat( "\nWrote riboClear Map file:   ", outfileMap, "\nN_riboClear genes:  ", nrow(ans$riboclearMap), "\n")

	out <- list( "FastaFile"=outfileFasta, "MapFile"=outfileMap)
	return( out)
}


`buildRiboClearData` <- function( speciesID, genomicFastaFile, tailSize=24, as.contigs=TRUE, max.contig.gap=1000) {

	rrnaMap <- getCurrentRrnaMap()
	exonMap <- getCurrentExonMap()
	geneMap <- getCurrentExonMap()

	# in case we make contigs, get the Rrna map ordered
	ord <- order( rrnaMap$SEQ_ID, rrnaMap$POSITION)
	rrnaMap <- rrnaMap[ ord, ]

	allocSize <- nrow(rrnaMap) * 2
	desc <- gid <- seqID <- seq <- grp <- seqBeg <- seqEnd <- faBeg <- faEnd <- cid <- vector( mode="character", 
			length=allocSize)

	curSeqID <- ""
	curSeqDNA <- ""
	nout <- 0

	# also track what we need to make contigs on the fly
	contigSet <- vector()
	nContigs <- 0

	# visit every ribo RNA entry
	cat( "\nbuilding riboClear data...")

	if ( ! ("CLEAR" %in% colnames( rrnaMap))) {
		cat( "\nNo 'CLEAR' column found in riboMap... flagging all ribo genes for clearing.")
		rrnaMap$CLEAR <- rep( TRUE, times=nrow(rrnaMap))
	}


	# local function to make a contig from 2+ genes
	contigify <- function() {
			if ( length(contigSet) >= 2) {
				cBegin <- as.numeric( seqBeg[ contigSet[1]])
				cEnd <- as.numeric( seqEnd[ contigSet[length(contigSet)]])
				contigDNA <- as.character( base::substr( curSeqDNA, (cBegin-tailSize), (cEnd+tailSize)))
				nContigs <<- nContigs + 1
				myShortCID <- paste("Contig", nContigs, sep=".")
				myLongCID <- paste( curSeqID, myShortCID, sep="::")
				# give all these genes this ContigID
				cid[ contigSet] <<- myShortCID
				# blank out all but one of the desc/seq set, and give the contig to the first one
				desc[ contigSet] <<- ""
				seq[ contigSet] <<- ""
				desc[ contigSet[1]] <<- myLongCID
				seq[ contigSet[1]] <<- contigDNA
				# and adjust all there FASTA location data
				tmpSeqBeg <- as.numeric( seqBeg[contigSet])
				tmpSeqEnd <- as.numeric( seqEnd[contigSet])
				tmpFaBeg <- as.numeric( faBeg[contigSet])
				tmpFaEnd <- as.numeric( faEnd[contigSet])
				deltaSeqBeg <- tmpSeqBeg - tmpSeqBeg[1]
				deltaSeqEnd <- tmpSeqBeg - tmpSeqBeg[1]
				tmpFaBeg <- tmpFaBeg + deltaSeqBeg
				tmpFaEnd <- tmpFaEnd + deltaSeqEnd
				faBeg[ contigSet] <<- as.character( tmpFaBeg)
				faEnd[ contigSet] <<- as.character( tmpFaEnd)
			}
			# that's all.  Always clean up for next contig set...
			contigSet <<- vector()
			return()
	}


	# ready to visit every gene in the ribo clear map
	# and be ready to make contigs any time we transition
	for ( ig in 1:nrow(rrnaMap)) {

		# only use entries that are flagged for clearing
		if ( ! as.logical( rrnaMap$CLEAR[ig])) {
			contigify()
			next
		}

		# watch for change of chromosome
		gene <- rrnaMap$GENE_ID[ ig]
		seqid <- rrnaMap$SEQ_ID[ ig]
		if ( seqid != curSeqID) {
			contigify()
			curSeqID <- seqid
			curSeqDNA <- getFastaSeqFromFilePath( filePath=genomicFastaFile, seqID=seqid)
			emap <- subset( exonMap, SEQ_ID == seqid)
		}

		if ( (ig %% 10 == 0) || (nout %% 10 == 0)) cat( "\n{N_genes, Gene, N_clears} = ", ig, gene, nout)

		# grab the DNA for this one gene, plus the tails...
		gBegin <- rrnaMap$POSITION[ig]
		gEnd <- rrnaMap$END[ig]
		gGroup <- rrnaMap$GROUP[ig]
		# since the gene group will be in the FASTA descriptor, make sure it will always parse correctly
		gGroup <- gsub( " ", ".", gGroup)
		geneDNA <- as.character( base::substr( curSeqDNA, (gBegin-tailSize), (gEnd+tailSize)))

		# add it to the set
		nout <- nout + 1
		gid[nout] <- gene
		seqID[nout] <- seqid
		desc[nout] <- base::paste( seqid, gene, gGroup, sep="::")
		seq[nout] <- geneDNA
		grp[nout] <- gGroup
		seqBeg[nout] <- as.character(gBegin)
		seqEnd[nout] <- as.character(gEnd)
		faBeg[nout] <- as.character( tailSize + 1)
		faEnd[nout] <- as.character( tailSize + 1 + (gEnd - gBegin))

		# either start or maybe add to a contig
		if ( ! length(contigSet)) {
			contigSet[1] <- nout
		} else {
			# we have 2 criteria:   first is close enough together
			prevEnd <- as.numeric( seqEnd[ nout-1])
			thisGap <- abs( gBegin - prevEnd)
			closeEnough <- ( thisGap <= max.contig.gap)
			# second is there are no other non-cleared genes in between
			otherGenes <- FALSE
			whGmap1 <- match( gid[nout-1], geneMap$GENE_ID)
			whGmap2 <- match( gid[nout], geneMap$GENE_ID)
			if ( ! any( is.na( c( whGmap1, whGmap2)))) {
				if ( (nOthers <- (whGmap2 - whGmap1)) > 1) {
					nReal <- sum( geneMap$REAL_G[(whGmap1+1):(whGmap2-1)])
					otherGenes <- (nReal > 0)
				}
			}
			if ( closeEnough && !otherGenes) {
				contigSet <- c( contigSet, nout)
			} else {
				contigify()
				contigSet[1] <- nout
			}
		}

		# if this gene has exon splices, build that construct too
		exPtrs <- which( emap$GENE_ID == gene)
		if ( length( exPtrs) < 2) next

		# but never let the multi-exon construct go into a contig
		contigify()

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
		seqBeg[nout] <- as.character( emap$POSITION[ exPtrs[1]])
		seqEnd[nout] <- as.character( emap$END[ exPtrs[ length(exPtrs)]])
		faBeg[nout] <- as.character( 1)
		faEnd[nout] <- as.character( nchar( geneDNA))
	}
	# all done, maybe make the last contig?
	contigify()

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
	length( cid) <- nout

	# because of possibly making contigs, some of the desc/seq entries are now unused
	keepFA <- which( desc != "")

	outDF <- data.frame( gid, grp, seqID, seqBeg, seqEnd, faBeg, faEnd, cid, stringsAsFactors=FALSE)
	colnames( outDF) <- c( RIBOCLEAR_COLUMNS, "CONTIG_ID")
	rownames( outDF) <- 1:nrow(outDF)

	return( list( "desc"=desc[keepFA], "seq"=seq[keepFA], "riboclearMap"=outDF))
}

