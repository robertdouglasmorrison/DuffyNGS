# riboConvert.R

# translate the alignments that came from a Ribo Clearing Index, back into their genomic coordinate system

`riboConvert` <- function( filein, fileout=filein, sampleID="", genomefile="genomic.bam",
			rawReadCount=NULL, riboMapFile="HsPf.riboClearMap.txt", 
			readBufferSize=1000000, verbose=TRUE) {

	# catch and adjust for same filename
	if ( fileout == filein) {
		tmpFileOut <- paste( filein, ".tmpBufferFile", sep="")
	} else {
		tmpFileOut <- fileout
	}
	if ( ! file.exists( filein)) {
		cat( "\nRiboClear BAM file not found: ", filein)
		return(NULL)
	}
	conRibo <- bamReader( filein)
	riboRefData <- getRefData( conRibo)

	# we also need a bam file from the genome that we are converting back too
	if ( ! file.exists( genomefile)) {
		cat( "\nGenomic BAM file not found: ", genomefile)
		return(NULL)
	}
	conGenome <- bamReader( genomefile)
	headerGenome <- getHeader( conGenome)
	genomeRefData <- getRefData( conGenome)
	
	# open a new BAM file to hold the reads after resolving back to the genome
	conOut <- bamWriter( headerGenome, filename=tmpFileOut)

	# gather up the Ribo conversion Map
	if ( ! file.exists( riboMapFile)) {
		cat( "\n\nCan not find riboMap file:  ", riboMapFile,"    skipping ribo comversion")
		return()
	}
	riboMap <- read.delim( riboMapFile, as.is=TRUE);
	if ( ! ("SEQ_POSITION" %in% colnames( riboMap))) riboMap$SEQ_POSITION <- riboMap$START

	ans <- calcAlignSummary( mode="setup", filename=filein, rawReadCount=rawReadCount, alignPhase="RiboClear")

	# read in the alignment file one buffer at a time
	hasMore <- TRUE
	nReads <- 0
	repeat {
		if ( ! hasMore) break
		if (verbose) cat( "\nReadBAM..")
		chunk <- getNextChunk( conRibo, n=readBufferSize, alignedOnly=TRUE)
		nNow <- size(chunk)
		if ( nNow < 1) break
		if ( nNow < readBufferSize) hasMore <- FALSE
		nReads <- nReads + nNow
		if (verbose) cat( "  N_RiboClear: ", prettyNum( as.integer( nReads), big.mark=","))

		# the SeqID field is in riboMap notation
		riboIndexID <- refID2seqID( refID( chunk), refData=riboRefData)
		#readIDs <- readID( chunk)

		# take apart the riboClear seqID,  they're in the form <seqid>::<geneID>::<group>
		if (verbose) cat( "  riboClearLookup..")

		# these are highly redundant, so factor them and lookup the smaller set
		ridFac <- factor( riboIndexID)
		ridPtrs <- tapply( 1:nNow, ridFac, FUN=NULL)
		rids <- levels( ridFac)
		terms <- strsplit( rids, split="::", fixed=TRUE)
		ridSID <- sapply( terms, function(x) return( x[1]))
		ridGID <- sapply( terms, function(x) return( x[2]))
		rm( terms)
		ridSPID <- getSpeciesFromSeqID( ridSID)

		# hand these unique parts back to all the duplicate locations
		# map these seqIDs back to genomic RefIDs
		newSeqIDs <- seqID2refID( ridSID, refData=genomeRefData)[ ridPtrs]
		geneIDs <- ridGID[ ridPtrs]
		speciesIDs <- ridSPID[ ridPtrs]

		# map from riboIndex units back to genomic position
		riboPosition <- position( chunk)
		whereInMap <- base::match( geneIDs, riboMap$GENE_ID, nomatch=0)
		seqOffset <- rep( NA, times=length(riboPosition))
		seqOffset[ whereInMap > 0] <- riboMap$SEQ_POSITION[ whereInMap]
		newPosition <- seqOffset + riboPosition - 1
		rm( seqOffset, whereInMap, riboPosition, ridSPID, ridGID, ridSID, ridPtrs, ridFac)
		gc()
	
		# ready to modify the refID & position
		cat( "  converting..")
		newChunk <- modifyAlign( chunk, refID=newSeqIDs, pos=newPosition)
		setTag( newChunk, BAMTAG_GENEID, geneIDs)

		cat( "  WriteBAM..")
		bamSave( conOut, newChunk)
		rm( newChunk)
		gc()

		calcAlignSummary( mode="addData", chunk=chunk, geneIDs=geneIDs, speciesIDs=speciesIDs)
		rm( chunk, geneIDs, speciesIDs, newPosition)
		gc()
	} # end of each buffer...

	bamClose( conOut)
	bamClose( conRibo)
	bamClose( conGenome)

	if ( filein == fileout) file.rename( tmpFileOut, fileout)

	# update the line count...
	if ( nReads > 0) quickFileLineCountRecord( fileout, sampleID, lineCount= nReads)

	# tally results
	ans <- calcAlignSummary( mode="report")
	cat( "\n", ans$textSummary, "\n")

	return( c( ans, "Alignments"=nReads ))
}
