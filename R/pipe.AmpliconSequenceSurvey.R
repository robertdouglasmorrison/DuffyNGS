# pipe.AmpliconSequenceSurvey.R -- find true AA sequences from mate pair amplicon FASTQ data

require( DuffyNGS)

`pipe.AmpliconSequenceSurvey` <- function( sampleID, referenceAA=NULL, annotationFile="Annotation.txt",
				optionsFile="Options.txt", results.path=NULL, 
				min.count=10, verbose=TRUE) {

	require( Biostrings)
	require( pwalign)
	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\n\nStarting 'Amplicon Sequence Survey' for sample:     ", sampleID, "\n\n")
	}
	startT <- proc.time()

	# setup to get at the fastq and alignment files we need...
	annT <- readAnnotationTable( annotationFile)
	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=FALSE)
	}
	fastq.path <- getOptionValue( optT, "fastqData.path", notfound=".")
	if ( ! ( sampleID %in% annT$SampleID)) stop( paste( "SampleID not in Annotation file: ", sampleID))
	pairedEnd <- getAnnotationTrue( annT, sampleID, "PairedEnd", notfound=FALSE)
	if ( ! pairedEnd) stop( "Amplicon Sequence pipe requires mate pair FASTQ data")
	
	# make .fastq readers for both raw files
	fastqFile <- annT$Filename[ match( sampleID, annT$SampleID)]
	fastqFileSet <- strsplit( fastqFile, split=", *")[[1]]
	fastqFile1 <- file.path( fastq.path, fastqFileSet[1])
	fastqFile2 <- file.path( fastq.path, fastqFileSet[2])
	fastqFile1 <- allowCompressedFileName( fastqFile1)
	rawFQ1 <- fastqReader()
	rawFQ1$initialize( fastqFile1)
	fastqFile2 <- allowCompressedFileName( fastqFile2)
	rawFQ2 <- fastqReader()
	rawFQ2$initialize( fastqFile2)

	# set up for doing fast overlap queries
	DNA_MATRIX <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
	tieBreakMode <- if ( is.null( referenceAA)) "evalue" else "reference"
	
	# done with all setup...
	# get ready to read thru the files
	nTotal <- nBadNN <- nMismatch <- nGoodPair <- 0
	aaSeqs <- vector()
	nAAseq <- 0
	if (verbose) cat( "\nEvaluating Amplicon Read Pairs for perfect overlap..\n")
	repeat {

		# get the next pair of reads
		raw1 <- rawFQ1$read(1)
		raw2 <- rawFQ2$read(1)
		if ( is.null( raw1) || is.null( raw2)) break
		nTotal <- nTotal + 1
		if ( nTotal %% 1000 == 0 && verbose) cat( "\nN_Pairs:", nTotal, " \tN_Good:", nAAseq, " \tN_BadN:", nBadNN, " \tN_Mismatch:", nMismatch)

		# for speed, turn them into DNASTRING objects
		dna1 <- DNAString( raw1$seq)
		dna2 <- DNAString( raw2$seq)
		if ( any( grepl( "N", c(dna1,dna2), fixed=T))) {
			nBadNN <- nBadNN + 1
			next
		}
		
		# we expect them to be a perfect overlapping mate pair, with no discrepancies where they overlap
		rc2 <- reverseComplement( dna2)
		pa <- pwalign::pairwiseAlignment( dna1, rc2, type="overlap", substitutionMatrix=DNA_MATRIX)
		frag1 <- as.character( pattern(pa))
		frag2 <- as.character( subject(pa))
		if ( frag1 != frag2) {
			nMismatch <- nMismatch + 1
			next
		}
		
		# we get a good overlap pair, so build the composite:  the non-overlap leader from #1 with all of #2
		overlapStart <- start( pattern(pa))
		out1 <- 
		bigdna <- paste( substr( raw1$seq, 1, overlapStart-1), as.character(rc2), sep="")
		
		# find the best coding frame
		thisAA <- DNAtoBestPeptide( bigdna, tieBreakMode=tieBreakMode, reference=referenceAA)
		nAAseq <- nAAseq + 1
		aaSeqs[nAAseq] <- thisAA
	}

	if (verbose) cat( "\nDone.")
	cat( "\nN_RawReadPairs: ", nTotal, " \tN_GoodPairs:  ", nAAseq, " \tN_BadNN:  ", nBadNN, " \tN_MisMatch:  ", nMismatch, "\n")
	rawFQ1$finalize()
	rawFQ2$finalize()
	
	# last item to do is to tabulate and decide a real vs noise cutoff
	aaTbl <- table( aaSeqs)
	drops <- which( as.numeric(aaTbl) < min.count)
	if ( length( drops)) {
		if (verbose) cat( "\nRemoving ", length(drops), " low count sequences")
		aaTbl <- aaTbl[ -drops]
	}
	aaTbl <- sort( aaTbl, decreasing=T)
	
	out <- data.frame( "AA.Sequence"=names(aaTbl), "Count"=as.numeric(aaTbl), stringsAsFactors=F)
	out$Frequency <- round( out$Count * 100 / sum( out$Count), digits=2)
	out$Pct.Raw.FASTQ <- round( out$Count * 100 / nTotal, digits=3)
	rownames(out) <- 1:nrow(out)
	return( out)
}
