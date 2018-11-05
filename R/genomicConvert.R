# genomicConvert.R

# translate the alignments that came from a Ribo Clearing Index, back into their genomic coordinate system

`genomicConvert` <- function( filein, fileout=filein, sampleID="", rawReadCount=NULL, 
			readBufferSize=1000000, verbose=TRUE) {

	# catch and adjust for same filename
	if ( fileout == filein) {
		tmpFileOut <- paste( filein, ".tmpBufferFile", sep="")
	} else {
		tmpFileOut <- fileout
	}
	if ( ! file.exists( filein)) {
		cat( "\nGenomic BAM file not found: ", filein)
		return(NULL)
	}
	conGenomic <- bamReader( filein)
	genomicRefData <- getRefData( conGenomic)
	headerGenome <- getHeader( conGenomic)
	
	# open a new BAM file to hold the reads after conversion
	conOut <- bamWriter( headerGenome, filename=tmpFileOut)

	ans <- calcAlignSummary( mode="setup", filename=filein, rawReadCount=rawReadCount, alignPhase="Genomic")

	# read in the alignment file one buffer at a time
	hasMore <- TRUE
	nReads <- 0
	repeat {
		if ( ! hasMore) break
		if (verbose) cat( "\nReadBAM..")
		chunk <- getNextChunk( conGenomic, n=readBufferSize, alignedOnly=TRUE)
		nNow <- size(chunk)
		if ( nNow < 1) break
		if ( nNow < readBufferSize) hasMore <- FALSE
		nReads <- nReads + nNow
		if (verbose) cat( "  N_Genomic: ", prettyNum( as.integer( nReads), big.mark=","))

		# convert the refID to a SeqID to figure the gene locations
		cat( "  geneIDs..")
		seqIDs <- refID2seqID( refID( chunk), refData=genomicRefData)
		positions <- position( chunk)
		lens <- nchar( readSeq( chunk))
		middles <- positions + round( lens/2)

		ans <- fastSP2GP( seqIDs, middles)
		geneIDs <- ans$GENE_ID
		speciesIDs <- ans$SPECIES

		# add a tag that has the GeneID terms
		setTag( chunk, BAMTAG_GENEID, geneIDs)

		# figure the unique vs multi read weights
		cat( "  readWts..")
		wts <- weight.align( chunk)
		wtStrs <- rep.int( "1", nNow)
		conv <- which( wts < 1)
		wtStrs[ conv] <- formatC( wts[conv], format="f", digits=2)
		setTag( chunk, BAMTAG_READWEIGHT, wtStrs)

		ans <- calcAlignSummary( mode="addData", chunk=chunk, geneIDs=geneIDs, speciesIDs=speciesIDs)

		cat( "  WriteBAM..")
		bamSave( conOut, chunk)
		gc()

	} # end of each buffer...

	bamClose( conOut)
	bamClose( conGenomic)

	if ( filein == fileout) file.rename( tmpFileOut, fileout)

	# update the line count...
	if (nReads > 0) quickFileLineCountRecord( fileout, sampleID, lineCount=nReads)

	# tally results
	ans <- calcAlignSummary( mode="report")
	cat( "\n", ans$textSummary, "\n")

	return( c( ans, "Alignments"=nReads))
}
