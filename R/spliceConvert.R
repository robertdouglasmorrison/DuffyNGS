# spliceConvert.R

# translate the alignments that came from a Splice Junction Index, back into their genomic coordinate system

`spliceConvert` <-
function( filein, fileout=sub( "bam$", "converted.bam", filein), genomefile="genomic.bam", 
		spliceMapPath=bowtie2Par("IndexPath"), spliceMapPrefix="spliceMap", 
		rawReadCount=NULL, readSense="sense", readBufferSize=1000000, sampleID="", verbose=TRUE) {

	timestart <- proc.time()

	# read one existing BAM
	if ( ! file.exists( filein)) {
		cat( "\nSplice BAM file not found: ", filein)
		return(NULL)
	}
	conSplice <- bamReader( filein)
	spliceRefData <- getRefData( conSplice)

	# we also need a bam file from the genome that we are converting back too
	if ( ! file.exists( genomefile)) {
		cat( "\nGenomic BAM file not found: ", genomefile)
		return(NULL)
	}
	conGenome <- bamReader( genomefile)
	headerGenome <- getHeader( conGenome)
	genomeRefData <- getRefData( conGenome)

	# open a new BAM file to hold the fragments after resolving the splices back to the genome
	conOut <- bamWriter( headerGenome, filename=fileout)

	ans <- calcAlignSummary( mode="setup", filename=filein, rawReadCount=rawReadCount, alignPhase="Splicing")

	# set up to do a buffer at a time
	nReads <- nReadsOut <- nUniqueReads <- nMultiReads <- nBad <- 0;
	hasMore <- TRUE

	# using local function...
	my_paste <- base::paste
	my_substr <- base::substr

	# grab a buffer
	repeat {
		if ( ! hasMore) break
		cat( "\nReadBAM..")
		chunk <- getNextChunk( conSplice, n=readBufferSize, alignedOnly=TRUE)
		nNow <- size(chunk)
		if ( nNow < 2) break
		if ( nNow < readBufferSize) hasMore <- FALSE
		nReads <- nReads + nNow
		cat( "  N_Splice: ", prettyNum( as.integer( nReads), big.mark=","))

		# some things can be done as a chunk
		#myRefID <- refID( chunk)
		#spliceStart <- position( chunk)
		#spliceSeq <- alignSeq( chunk)
		#spliceQual <- alignQual( chunk)
		#myRIDs <- readID( chunk)
		#myMateRefID <- mateRefID( chunk)
		#mateSpliceStart <- matePosition( chunk)

		# do it as a data frame to interate over the file less times
		chunkDF <- as.data.frame( chunk)
		myRefID <- chunkDF$refid
		spliceStart <- chunkDF$position
		spliceSeq <- chunkDF$seq
		spliceQual <- chunkDF$qual
		myRIDs <- chunkDF$name
		myMateRefID <- chunkDF$mrefid
		mateSpliceStart <- chunkDF$mposition
		myWts <- weight.align( chunk)

		# names for the 2 new splice pieces
		newReadID1 <- my_paste( myRIDs, "::splice1", sep="")
		newReadID2 <- my_paste( myRIDs, "::splice2", sep="")

		# the SeqID field is in spliceMap notation  'seqID::geneID::spliceID'
		spliceIndexID <- refID2seqID( myRefID, refData=spliceRefData)
		spliceIDterms <- strsplit( spliceIndexID, split="::")
		allSeqIDs <- sapply( spliceIDterms, `[`, 1)
		allGeneIDs <- sapply( spliceIDterms, `[`, 2)
		allSpliceIDs <- sapply( spliceIDterms, `[`, 3)

		# gather up the splice junction mapping data for these splices
		cat( "  junctionLookup..")
		spliceMap <- spliceJunctionLookup( spliceIndexID, spliceMapPath, spliceMapPrefix);

		# map from the true SeqIDs back to the genomic BAM refIDs
		mySIDs <- spliceMap$SEQ_ID
		myNewRefIDs <- seqID2refID( mySIDs, refData=genomeRefData)

		# some things get parsed out of the splice map details
		sizeTerms <- strsplit( spliceMap$SIZES, split=",")
		startTerms <- strsplit( spliceMap$POSITIONS, split=",")
		endTerms <- strsplit( spliceMap$ENDS, split=",")
		allSizes <- lapply( sizeTerms, as.integer)
		allStarts <- lapply( startTerms, as.integer)
		allEnds <- lapply( endTerms, as.integer)
		allLengths <- nchar( spliceSeq)

		# we may also have/need the mate pair info
		isPaired <- paired(chunk)
		DO_PAIRS <- length( (whoPaired <- which( isPaired)))
		if ( DO_PAIRS) {
			spliceIndexID2 <- refID2seqID( myMateRefID[whoPaired], refData=spliceRefData)
			spliceMap2 <- spliceJunctionLookup( spliceIndexID2, spliceMapPath, spliceMapPrefix);
			#spliceID2terms <- strsplit( spliceIndexID2, split="::")
			mySID2s <- spliceMap2$SEQ_ID
			myNewRefID2s <- seqID2refID( mySID2s, refData=genomeRefData)
			# stuff these converted refIDs back where they go
			myMateRefID[ whoPaired] <- myNewRefID2s
		}

		# OK, ready to visit each splice
		cat( "  converting..")
		newChunk <- bamRange()
		newGID <- newSPID <- wtStrs <- rep( "", times=nNow*2)
		nOutNow <- 0


		for ( i in 1:nNow) {

			# for speed, try using the .Call layer directly to avoid method lookup/load
			tmpAns <- .Call( "bam_range_get_next_align", chunk@range, FALSE, FALSE, PACKAGE="DuffyNGS")
			thisAlign <- new( "bamAlign", tmpAns)

			# invalid splices get ignored
			if ( spliceMap$GENE_ID[i] == "") next
			thisSeq <- allSeqIDs[i]
			thisGene <- allGeneIDs[i]
			thisSplice <- allSpliceIDs[i]
			if ( is.na(thisSeq) || is.na(thisGene) || is.na(thisSplice)) {
				warning( paste( "Invalid spliceID: ", spliceTable$SPLICE_INDEX_ID[i]))
				nBad <- nBad + 1
				next
			}

			# Create the splice in terms of chromosome coordinates
			sizes <- allSizes[[i]]
			starts <- allStarts[[i]]
			ends <- allEnds[[i]]
			if ( length(starts) != 2) {
				nBad <- nBad + 1
				next
			}

			#get the splice-centric coordinates of the read
			nbases <- allLengths[i]
			read.start = starts[1] + spliceStart[i] - 1;
			read.end = read.start + nbases - 1;
			read.frag1end <- ends[1]
			read.frag2start <- starts[2]
			loca.start <- 1
			loca.end <- nbases
			loca.frag1end <- read.frag1end - read.start + 1
			loca.frag2start <- loca.frag1end + 1
		
			# prep all the new fields we will write
			newPos1 <- read.start
			newPos2 <- read.frag2start
			newAlignSeq1 <- my_substr( spliceSeq[i], loca.start, loca.frag1end);
			newAlignQual1 <- my_substr( spliceQual[i], loca.start, loca.frag1end);
			newAlignSeq2 <- my_substr( spliceSeq[i], loca.frag2start, loca.end);
			newAlignQual2 <- my_substr( spliceQual[i], loca.frag2start, loca.end);

			# all the parts are ready if not a mate pair
			# do separate calls to save time when we can
			if ( !DO_PAIRS) {
				# write these 2 new alignments, calling C directly for speed
				tmpAns <- .Call( "bam_align_modify", thisAlign@align, as.integer(myNewRefIDs[i]), 
						as.integer(newPos1), NULL, NULL,
						newReadID1[i], newAlignSeq1, newAlignQual1, PACKAGE="DuffyNGS")
				new1 <- new( "bamAlign", tmpAns)
				.Call( "bam_range_push_back", newChunk@range, new1@align, PACKAGE="DuffyNGS")
	
				tmpAns <- .Call( "bam_align_modify", thisAlign@align, as.integer(myNewRefIDs[i]), 
						as.integer(newPos2), NULL, NULL,
						newReadID2[i], newAlignSeq2, newAlignQual2, PACKAGE="DuffyNGS")
				new2 <- new( "bamAlign", tmpAns)
				.Call( "bam_range_push_back", newChunk@range, new2@align, PACKAGE="DuffyNGS")
			} else {
				# may have mate pair to worry about
				myMatePos1 <- myMatePos2 <- mateSpliceStart[i]
				if ( myMateRefID[i] >= 0) {
					diffPos <- mateSpliceStart[i] - spliceStart[i]
					myMatePos1 <- newPos1 + diffPos
					myMatePos2 <- newPos2 + diffPos
				}

				# write these 2 new alignments, calling C directly for speed
				tmpAns <- .Call( "bam_align_modify", thisAlign@align, as.integer(myNewRefIDs[i]), 
						as.integer(newPos1), as.integer(myMateRefID[i]), as.integer(myMatePos1),
						newReadID1[i], newAlignSeq1, newAlignQual1, PACKAGE="DuffyNGS")
				new1 <- new( "bamAlign", tmpAns)
				.Call( "bam_range_push_back", newChunk@range, new1@align, PACKAGE="DuffyNGS")
	
				tmpAns <- .Call( "bam_align_modify", thisAlign@align, as.integer(myNewRefIDs[i]), 
						as.integer(newPos2), as.integer(myMateRefID[i]), as.integer(myMatePos2),
						newReadID2[i], newAlignSeq2, newAlignQual2, PACKAGE="DuffyNGS")
				new2 <- new( "bamAlign", tmpAns)
				.Call( "bam_range_push_back", newChunk@range, new2@align, PACKAGE="DuffyNGS")
			}
			nReadsOut <- nReadsOut + 2

			newGID[nOutNow + (1:2)] <- thisGene
			newSPID[nOutNow + (1:2)] <- thisSplice
			wtStrs[nOutNow + (1:2)] <- "1"
			if ( myWts[i] < 1) wtStrs[nOutNow + (1:2)] <- formatC( myWts[i], format="f", digits=2)
			nOutNow <- nOutNow + 2
		}

		# with all the new alignments in this new range, we can slam in some new tags
		length(newGID) <- length(newSPID) <- nOutNow
		setTag( newChunk, BAMTAG_GENEID, newGID)
		setTag( newChunk, BAMTAG_SPLICEID, newSPID)
		setTag( newChunk, BAMTAG_READWEIGHT, wtStrs)

		# now we can accumulate summary facts too
		who <- seq.int( 1, nOutNow, 2)
		mySpliceIDs <- paste( shortGeneName( newGID[who], keep=1), newSPID[who], sep="::")
		mySpeciesIDs <- getSpeciesFromSeqID( mySIDs)
		ans <- calcAlignSummary( mode="addData", chunk=chunk, geneIDs=mySpliceIDs, 
						speciesIDs=mySpeciesIDs)

		cat( "  writeBAM..")
		bamSave( conOut, newChunk)
		rm( newChunk)
		gc()
	} # end of each buffer...

	bamClose( conOut)
	bamClose( conSplice)
	bamClose( conGenome)

	# tally results
	ans <- calcAlignSummary( mode="report")
	cat( "\n", ans$textSummary, "\n")

	return( c( ans, list( "Alignments"=nReadsOut, "SplicedReads"=nReads, "NonSpliceAlignments"=nBad)))
}
