# pipe.AmpliconSequenceSurvey.R -- find true AA sequences from mate pair amplicon FASTQ data

require( DuffyNGS)

`pipe.AmpliconSequenceSurvey` <- function( sampleID, reference=NULL, annotationFile="Annotation.txt",
				optionsFile="Options.txt", results.path=NULL, min.count=5, min.confidence=75,
				chunkSize=1000, 
				primers=NULL, verbose=TRUE) {

	require( Biostrings)
	require( pwalign)

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nStarting 'Amplicon Sequence Survey' for sample:     ", sampleID, "\n")
	}
	startT <- proc.time()

	# setup to inspect the fastq files we need...
	annT <- readAnnotationTable( annotationFile)
	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=FALSE)
	}
	fastq.path <- getOptionValue( optT, "fastqData.path", notfound=".", verbose=F)
	if ( ! ( sampleID %in% annT$SampleID)) stop( paste( "SampleID not in Annotation file: ", sampleID))
	pairedEnd <- getAnnotationTrue( annT, sampleID, "PairedEnd", notfound=FALSE, verbose=F)
	if ( ! pairedEnd) stop( "Amplicon Sequence pipe requires mate pair FASTQ data")

	out.path <- file.path( results.path, "AmpliconSequences")
	if ( ! file.exists( out.path)) dir.create( out.path, recursive=T)
	
	# the tool can be given reference proteins to match against, as either sequences or GeneIDs
	referenceAA <- reference
	if ( ! is.null( reference)) {
		gmap <- subset( getCurrentGeneMap(), REAL_G == TRUE)
		isGene <- which( reference %in% gmap$GENE_ID)
		if ( length( isGene)) {
			genomeFasta <- getOptionValue( optT, "genomicFastaFile", notfound="", verbose=FALSE)
			for ( ig in isGene) {
				tmpFA <- gene2Fasta( reference[ig], genomeFasta, mode="aa", verbose=F)
				referenceAA[ig] <- tmpFA$seq
				names(referenceAA)[ig] <- tmpFA$desc
			}
		}
		if ( is.null( names(referenceAA))) cat( "\nWarning:  'reference' proteins should have a names attribute")
	}

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
	if ( is.null( reference)) {
		tieBreakMode <- "evalue"
		referenceAAstring <- NULL
	} else {
		tieBreakMode <- "reference"
		# for finding the best reading frame, concatenate multiple proteins into one string
		referenceAAstring <- paste( referenceAA, collapse="XXX")
	}
	
	# done with all setup...
	# get ready to read thru the files
	nTotal <- nBadNN <- nMismatch <- nGoodPair <- 0
	joinedSeq <- vector()
	nOut <- 0
	if (verbose) cat( "\n\nEvaluating Amplicon Read Pairs for perfect overlap..\n")
	repeat {

		# get the next pair of reads
		raw1 <- rawFQ1$read(1)
		raw2 <- rawFQ2$read(1)
		if ( is.null( raw1) || is.null( raw2)) break
		nTotal <- nTotal + 1
		if ( nTotal %% chunkSize == 0 && verbose) cat( "\nN_Pairs:", nTotal, " \tN_Good:", nOut, " \tN_BadN:", nBadNN, " \tN_Mismatch:", nMismatch)

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
		bigdna <- paste( substr( raw1$seq, 1, overlapStart-1), as.character(rc2), sep="")
		nOut <- nOut + 1
		joinedSeq[nOut] <- bigdna
	}

	if (verbose) cat( "\nDone.")
	cat( "\nN_RawPairs:", nTotal, " \tN_GoodPairs:", nOut, " \tN_BadN:", nBadNN, " \tN_Mismatch:", nMismatch, "\n")
	rawFQ1$finalize()
	rawFQ2$finalize()
	
	# now tabulate and decide a real vs noise cutoff
	dnaTbl <- table( joinedSeq)
	drops <- which( as.numeric(dnaTbl) < min.count)
	if ( length( drops)) {
		if (verbose) cat( "\nRemoving ", length(drops), " low count sequences having less than ", min.count, "copies.")
		dnaTbl <- dnaTbl[ -drops]
	}
	dnaTbl <- sort( dnaTbl, decreasing=T)
	# turn these unique DNA strings into AA, after finding/trimming any unwanted primer/adaptor tails
	dnaSeqs <- names(dnaTbl)
	dnaCounts <- as.numeric(dnaTbl)

	# if we were given forward/reverse primer sequences, find them now and use them to trim edges
	if ( ! is.null( primers)) {
		cat( "\nTrimmer primer edges..\n")
		for ( k in 1:length(dnaTbl)) {
			for ( i in 1:length(primers)) {
				primerSeq <- primers[i]
				minScore <- nchar(primerSeq) * 0.9
				pa <- pwalign::pairwiseAlignment( primerSeq, dnaSeqs[k], type="global-local", substitutionMatrix=DNA_MATRIX)
				if (score(pa) > minScore) {
					keepStart <- start( subject(pa))
					if ( keepStart > 1 && keepStart < nchar(dnaSeqs[k])/2) {
						trimmedSeq <- substr( dnaSeqs[k], 1, keepStart-1)
						dnaSeqs[k] <- substring( dnaSeqs[k], keepStart)
						if (verbose) cat( "\nTrimmed from front: ", k, nchar(trimmedSeq), trimmedSeq)
					}
				}
				revcompSeq <- myReverseComplement( primerSeq)
				pa <- pwalign::pairwiseAlignment( revcompSeq, dnaSeqs[k], type="global-local", substitutionMatrix=DNA_MATRIX)
				if (score(pa) > minScore) {
					seqStart <- start( subject(pa))
					seqEnd <- seqStart + width( subject(pa)) - 1
					if ( seqEnd < nchar(dnaSeqs[k]) && seqEnd > nchar(dnaSeqs[k])/2) {
						trimmedSeq <- substring( dnaSeqs[k], seqEnd+1)
						dnaSeqs[k] <- substr( dnaSeqs[k], 1, seqEnd)
						if (verbose) cat( "\nTrimmed from end: ", k, nchar(trimmedSeq), trimmedSeq)
					}
				}
			}
		}
	}

	# now we can find the best coding frame for each DNA
	aaSeqs <- vector()
	for ( k in 1:length(dnaSeqs)) {
		thisAA <- DNAtoBestPeptide( dnaSeqs[k], tieBreakMode=tieBreakMode, reference=referenceAAstring)
		aaSeqs[k] <- thisAA
	}
	
	# because we may have done base trimming, and because of translation to AA and reading frame, 
	# we may now have duplicate AA strings.  Combine them and recount.
	aaFac <- factor( aaSeqs)
	newCounts <- tapply( dnaCounts, aaFac, sum)
	ord <- order( newCounts, decreasing=T)
	aaSeqsFinal <- levels(aaFac)[ord]
	aaCounts <- newCounts[ord]
	
	out <- data.frame( "SampleID"=sampleID, "AA.Sequence"=aaSeqsFinal, "Length"=nchar(aaSeqsFinal), "Count"=aaCounts, stringsAsFactors=F)
	rownames(out) <- 1:nrow(out)

	# lastly, append extra info we need/want about the gene and where the sequence hit
	out$Frequency <- round( out$Count * 100 / sum( out$Count), digits=2)
	out$Pct.Raw.FASTQ <- round( out$Count * 100 / nTotal, digits=3)
	# calculate how similar/dissimilar the sequence is to the reference
	out$GENE_ID <- NA
	out$PRODUCT <- NA
	out$Start <- NA
	out$Stop <- NA
	out$Ref.Sequence <- NA
	out$EditDistance <- NA
	out$Confidence <- NA
	out$Mutations <- ""
	if ( ! is.null(reference)) {
		data( BLOSUM62)
		for ( i in 1:nrow(out)) {
			thisSeq <- out$AA.Sequence[i]
			thisLen <- nchar(thisSeq)
			pa <- pwalign::pairwiseAlignment( referenceAA, thisSeq, type="local-global", substitutionMatrix=BLOSUM62)
			best <- if (length(referenceAA) > 1) which.max( score(pa)) else 1
			out$Ref.Sequence[i] <- refFrag <- as.character( pattern(pa)[best])
			out$Start[i] <- start( pattern(pa)[best])
			out$Stop[i] <- width( pattern(pa)[best]) + out$Start[i] - 1
			editDist <- adist( thisSeq, refFrag)
			out$GENE_ID[i] <- myGID <- names(referenceAA)[best]
			out$PRODUCT[i] <- gene2Product( myGID)
			out$EditDistance[i] <- editDist
			out$Confidence[i] <- round( (thisLen - editDist) * 100 / thisLen, digits=1)
			if ( editDist > 0) {
				thisAAvec <- strsplit( thisSeq, split="")[[1]]
				refAAvec <- strsplit( refFrag, split="")[[1]]
				nCheck <- min( length(thisAAvec), length(refAAvec))
				isMUT <- which( thisAAvec[1:nCheck] != refAAvec[1:nCheck])
				refOffset <- start( pattern(pa)[best]) - 1
				mutationStr <- paste( refAAvec[isMUT], (isMUT+refOffset), thisAAvec[isMUT], sep="", collapse="; ")
				out$Mutations[i] <- mutationStr
			}
		}
	}
	
	# append the raw data to support what we found
	where <- match( aaSeqsFinal, aaSeqs)
	out$Trimmed.Joined.DNA <- dnaSeqs[where]
	
	# lastly, drop those with terrible confidence, as remnants of junk amplification
	drops <- which( out$Confidence < min.confidence)
	if ( length(drops)) {
		if (verbose) cat( "\nRemoving ", length(drops), " low quality sequences having less than ", min.confidence, "confidence.")
		out <- out[ -drops, ]
		if (nrow(out)) rownames(out) <- 1:nrow(out)
	}

	# all done.  Write it out, and return it
	outfile <- file.path( out.path, paste( sampleID, "AmpliconSequences.csv", sep="."))
	write.csv( out, outfile, row.names=F)
	return( out)
}
