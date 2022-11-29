# pipe.TargetSearch.R -- mine No Hits (or raw fASTQ) for evidence of arbitrary Target

`pipe.TargetSearch` <- function( sampleIDset, targetName="Virus.Genes", fastaName=targetName, max.mismatch.per.100=5,
				max.low.complexity=80, input.mode=c("nohits","fastq"), 
				annotationFile="Annotation.txt", optionsFile="Options.txt", 
				k.bowtie=1, chunk.size=10000, min.reads.RPKM=10000, verbose=FALSE) {

	# get options that we need
	optT <- readOptionsTable( optionsFile)
	results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=verbose)
	fastq.path <- getOptionValue( optT, "fastqData.path", verbose=verbose)
	bowtie.path <- getOptionValue( optT, "bowtie2Index.path", notfound=".", verbose=verbose)

	# make sure we see and can load the target data
	fastaFile <- file.path( bowtie.path, "FastaFiles", paste( fastaName, "fasta", sep="."))
	if ( ! file.exists( fastaFile)) stop( paste( "Required Target FASTA file not found: ", fastaFile))
	targetFA <- loadFasta( fastaFile, short.desc=FALSE, verbose=verbose)
	descTerms <- strsplit( targetFA$desc, split=" | ", fixed=T)
	targetIDs <- sapply( descTerms, `[`, 1)
	targetNames <- sapply( descTerms, `[`, 2)
	targetSizes <- sapply( targetFA$seq, nchar)
	genomeSize <- sum( targetSizes, na.rm=T)
	if ( any( is.na( targetNames))) cat( "\nWarning:  some Fasta headers missing ' | ' product term delimitors")

	# same with the target Bowtie index
	bowtieFile <- file.path( bowtie.path, paste( targetName, "1.bt2", sep="."))
	if ( ! file.exists( bowtieFile)) stop( paste( "Required", targetName, "Bowtie2 target index file not found: ", bowtieFile))

	# make the place/file for our results
	bam.path <- file.path( results.path, "TargetSearch", targetName)
	if ( ! file.exists( bam.path)) dir.create( bam.path, recursive=T)

	# since we don't expect much hits (if at all), allow doing a set of samples into one result
	bigDF <- data.frame()
	nRejectMismatch <- nRejectLowComplex <- nFound <- 0

	# detecting multiple alignments is based on assumptions about adjacent records in the BAM file.
	# make sure we can see those
	if ( k.bowtie > (chunk.size/1000)) stop( "Increase 'chunk.size' for searching BAM file. Must be bigger than 'k.size' * 1000")

	for (sid in sampleIDset) {

		cat( "\n\nStarting '", targetName, "' Search for Sample:    ", sid, 
			"\n  Max base mismatches per 100bp:  ", max.mismatch.per.100, 
			"\n  Max low complexity base %:      ", max.low.complexity, 
			"\n  Max alignments per read 'K':    ", k.bowtie, 
			"\n  Min reads for RPKM calculation: ", min.reads.RPKM, sep="")
	
		# get the file of noHit reads, as target reads should still be in there
		input.mode <- match.arg( input.mode)
		if ( input.mode == "nohits") {
			infile <- paste( sid, "noHits", "fastq", sep=".") 
			infile <- file.path( results.path, "fastq", infile)
			infile <- allowCompressedFileName( infile)
			if ( ! file.exists( infile)) {
				notgenomicfile <- paste( sid, "not.genomic", "fastq", sep=".") 
				notgenomicfile <- file.path( results.path, "fastq", notgenomicfile)
				notgenomicfile <- allowCompressedFileName( notgenomicfile)
				cat( "\n'No Hits' File not found:  ", infile, "\nTrying: ", notgenomicfile)
				infile <- notgenomicfile
			}
			nNoHits <- getFileLineCount( infile, sid) / 4
			cat( "\n  Searching the 'No Hits' bin of unaligned reads")
		} else {
			# use the raw FASTQ files instead of the NoHits
			rawFastq <- getRawFastqFileNames( sid, annotationFile, optionsFile, verbose=FALSE)
			infile <- rawFastq$files
			nNoHits <- getFileLineCount( infile[1], sid) / 4
			cat( "\n  Searching the FASTQ file(s) of all raw reads")
		}
		if ( ! all( file.exists( infile))) {
			cat( "\nFASTQ file not found:  ", infile)
			next
		}

		# call Bowtie
		cat( "\n  Aligning..")
		bamFile <- file.path( bam.path, paste( sid, targetName, "Hits.bam", sep="."))
		ans1 <- fastqToBAM( inputFastqFile=infile, outputFile=bamFile, k=k.bowtie, sampleID=sid, optionsFile=optionsFile,
					alignIndex=targetName, index.path=bowtie.path, keepUnaligned=FALSE, 
					verbose=verbose, label=sid)
		nReadsIn <- ans1$RawReads

		# now see what got hit
		cat( "  Scan for high quality hits..")
		reader <- bamReader( bamFile)
		refData <- getRefData( reader)
		targetHits <- targetWtsM <- targetWtsU <- vector()
		nSeen <- 0
		repeat {
			chunk <- getNextChunk( reader, n=chunk.size, alignedOnly=T)
			nReads <- size( chunk)
			if ( nReads < 1) break
			nSeen <- nSeen + nReads
			refIDs <- refID( chunk)
			seqID <- refID2seqID( refIDs, reader=reader, refData=refData)
			# look at how many mismatches per read
			misMatchCnt <- as.numeric( getTag( chunk, "XM"))
			readLen <- nchar( readSeq( chunk))
			badMatch <- which( (misMatchCnt * 100 / readLen) > max.mismatch.per.100)
			# also look for low complexity
			bases <- strsplit( readSeq( chunk), split="")
			pctLow <- sapply( bases, function(x) {
					myPct <- table(x) / length(x)
					return( max(myPct))
				})
			badComplex <- which( (pctLow * 100) > max.low.complexity)
			goodHits <- setdiff( 1:nReads, union( badMatch, badComplex))
			targetHits <- c( targetHits, seqID[goodHits])
			# also keep the weights if K more than one
			if ( k.bowtie == 1) {
				targetWts <- c( targetWts, rep.int(1,length(goodHits)))
			} else {
				# we need to calc the wts
				readIDs <- readID( chunk)
				readWtsM <- readWtsU <- rep.int( 1, nReads)
				tapply( 1:nReads, factor(readIDs), function(x) { 
					readWtsM[x] <<- (1/length(x)); 
					readWtsU[x] <<- if ( length(x) > 1) 0 else 1
					return(NULL)
				})
				targetWtsM <- c( targetWtsM, readWtsM[goodHits])
				targetWtsU <- c( targetWtsU, readWtsU[goodHits])
			}

			nFound <- nFound + length(goodHits)
			nRejectMismatch <- nRejectMismatch + length(badMatch)
			nRejectLowComplex <- nRejectLowComplex + length(badComplex)
		}
		bamClose( reader)
	
		if ( ! length(targetHits)) {
			cat( "\nNo high quality", targetName, "hits detected..")
			next
		}

		# now tablify
		if ( k.bowtie > 1) {
			hitTblM <- tapply( targetWtsM, factor( targetHits), sum, na.rm=T)
			hitTblU <- tapply( targetWtsU, factor( targetHits), sum, na.rm=T)
		} else {
			hitTblM <- table( targetHits)
			hitTblU <- rep.int( 0, length(hitTblM))
		}
		hitCntsM <- round( as.numeric( hitTblM))
		hitCntsU <- round( as.numeric( hitTblU))
		pctHitsM <- round( hitCntsM * 100 / sum(hitCntsM), digits=2)
		pctHitsU <- round( hitCntsU * 100 / sum(hitCntsU), digits=2)
		where <- match( names(hitTblM), targetIDs)
		products <- targetNames[where]
		# make a RPKM and TPM call too, both use gene lengths
		hitLens <- targetSizes[where]
		hitLens[ is.na( hitLens)] <- 100
		totalReadsM <- sum( hitCntsM, na.rm=T)
		totalReadsU <- sum( hitCntsU, na.rm=T)
		totalReadsForRPKMm <- max( totalReadsM, min.reads.RPKM)
		totalReadsForRPKMu <- max( totalReadsU, min.reads.RPKM)
		hitsRPKM_M <- round( calculateRPKM( hitCntsM, exonBases=hitLens, totalReads=totalReadsForRPKMm), digits=2)
		hitsRPKM_U <- round( calculateRPKM( hitCntsU, exonBases=hitLens, totalReads=totalReadsForRPKMu), digits=2)
		hitsTPM_M <- round( tpm( readCount=hitCntsM, geneLen=hitLens, readLen=100, minReadsPerSpecies=min.reads.RPKM), digits=2)
		hitsTPM_U <- round( tpm( readCount=hitCntsU, geneLen=hitLens, readLen=100, minReadsPerSpecies=min.reads.RPKM), digits=2)

		smlDF <- data.frame( "SampleID"=sid, "TargetID"=names(hitTblM), "Product"=products, 
				"READS_M"=hitCntsM, "PCT_M"=pctHitsM, "RPKM_M"=hitsRPKM_M, "TPM_M"=hitsTPM_M, 
				"READS_U"=hitCntsU, "PCT_U"=pctHitsU, "RPKM_U"=hitsRPKM_U, "TPM_U"=hitsTPM_U,
				stringsAsFactors=F)
		ord <- order( hitCntsU, hitCntsM, decreasing=T)
		smlDF <- smlDF[ ord, ]
		rownames(smlDF) <- 1:nrow(smlDF)
		outfile <- file.path( bam.path, paste( sid, targetName, "Summary.csv", sep="."))
		write.table( smlDF, outfile, sep=",", quote=T, row.names=F)

		cat( "\nDone.  N_Aligned: ", nSeen, "  N_Good_Hits: ", sum(hitTblM >= 1), "\n")
		print( head( smlDF, 3))

		# accumulate the results
		bigDF <- rbind( bigDF, smlDF)
	}
	if ( ! nrow(bigDF)) return( bigDF)

	# final ordering is...
	ord <- order( bigDF$READS_U, bigDF$READS_M, decreasing=T)
	bigDF <- bigDF[ ord, ]
	rownames(bigDF) <- 1:nrow(bigDF)

	# report the overall picture
	cat( "\n\nSummary:")
	cat( "\n  N Good Target Reads Found:      ", nFound)
	cat( "\n  N Rejected for Mismatches:     ", nRejectMismatch)
	cat( "\n  N Rejected for Low Complexity: ", nRejectLowComplex, "\n")

	return( bigDF)
}


`pipe.TargetProteinPileups` <- function( sampleID, geneID, targetName="Virus.Genes", fastaName=targetName, 
					optionsFile="Options.txt", results.path=NULL, max.depth=60, txt.cex=0.21, 
					forceSetup=FALSE, plotOnly=FALSE, max.drawnPerSite=3,
					trim5.aligns=0, trim3.aligns=0, draw.box=FALSE, chunkSize.pileup=50000, 
					startFromReference=FALSE, showFrameShiftPeptides=TRUE, verbose=TRUE) {

	# a version of the CPP tool, specialized for our arbitrary "TargetSearch" workflow
	
	require(Biostrings)

	# make sure we were given valid parameters..
	if ( length(sampleID) > 1) {
		sampleID <- sampleID[1]
		cat( "\nWarning:  only one sampleID allowed")
	}

	optT <- readOptionsTable( optionsFile)
	if ( is.null(results.path)) results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	bowtie.path <- getOptionValue( optT, "bowtie2Index.path", notfound=".", verbose=verbose)

	# make sure we see and can load the target data
	targetFastaFile <- file.path( bowtie.path, "FastaFiles", paste( fastaName, "fasta", sep="."))
	if ( ! file.exists( targetFastaFile)) stop( paste( "Required Target FASTA file not found: ", targetFastaFile))
	targetFA <- loadFasta( targetFastaFile, short.desc=FALSE, verbose=verbose)
	descTerms <- strsplit( targetFA$desc, split=" | ", fixed=T)
	targetIDs <- sapply( descTerms, `[`, 1)
	targetNames <- sapply( descTerms, `[`, 2)
	targetSizes <- sapply( targetFA$seq, nchar)
	target.path <- file.path( results.path, "TargetSearch", targetName)
	if ( ! file.exists( target.path)) stop( paste( "Target Search results folder not found: ", target.path))
	
	# make sure the wanted gene is in that target
	whereTarget <- match( geneID, targetIDs, nomatch=0)
	if (whereTarget == 0) {
		whereTarget <- grep( paste( "^", geneID, sep=""), targetIDs)[1]
		if ( is.na( whereTarget)) stop( paste( "Given GeneID not in target: ", geneID))
		geneID <- targetIDs[ whereTarget]
	}
	geneLengthDNA <- targetSizes[whereTarget]
	
	# these gene names may be MHC/HLA, with non-standard characters.  Make a filename safe version..
	geneFileID <- file.cleanSpecialCharactersFromFileName( geneID)

	# step 1: make sure all needed peptide and pileup files are in place
	geneReadsFile <- file.path( target.path, paste( sampleID, geneFileID, "fastq.gz", sep="."))
	genePeptidesFile <- file.path( target.path, paste( sampleID, geneFileID, "RawReadPeptides.txt", sep="."))
	if ( forceSetup || ! file.exists( geneReadsFile)) {
		cat( "\nExtracting Gene alignments..")
		nAligns <- pipe.GatherTargetAlignments( sampleID, geneID, targetName=targetName, 
						geneLength=geneLengthDNA, asFASTQ=T, fastq.keyword=geneFileID)
		if (nAligns < 1) {
			cat( "\nZero reads aligned to gene..")
			return( 0)
		}
		# new file, so make sure we remake its decendants
		file.delete( genePeptidesFile)
	}
	if ( forceSetup || ! file.exists( genePeptidesFile)) {
		cat( "\nConverting Gene alignments to peptides..")
		nPeptides <- fastqToPeptides( geneReadsFile, genePeptidesFile, chunk=100000, lowComplexityFilter=FALSE,
				trim5=trim5.aligns, trim3=trim3.aligns, clipAtStop=FALSE)
		if (nPeptides < 1) {
			cat( "\nGene alignments gave zero peptides..")
			return( 0)
		}
	}

	# step 2:  run the tool that aligns peptides to a construct
	bamFile <- paste( sampleID, targetName, "Hits.bam", sep=".")
	bamFile <- file.path( target.path, bamFile)
	consensusAAfile <- file.path( target.path, paste( sampleID, geneFileID, "ConsensusAA.fasta", sep="."))	
	consensusDNAfile <- file.path( target.path, paste( sampleID, geneFileID, "ConsensusDNA.fasta", sep="."))
	consensusBASEfile <- file.path( target.path, paste( sampleID, geneFileID, "ConsensusBaseCalls.txt", sep="."))
	if ( forceSetup || ! all( file.exists( c( consensusBASEfile, consensusAAfile, consensusDNAfile)))) {
		cat( "\nTurning alignment pileups into consensus base calls..")
		ans <- pipe.ConsensusBaseCalls( sampleID, geneID=NULL, seqID=geneID, start=1, stop=geneLengthDNA, 
				genomicFastaFile=targetFastaFile, bamfile=bamFile, as.cDNA=TRUE, noReadCalls="blank")
		# since the target is not part of a standard gene map, do the AA step manually
		ans <- consensusBaseCallsToCDNA( ans, geneID=NULL, start=1, stop=geneLengthDNA, 
					best.frame=TRUE, verbose=verbose)
		callsTable <- ans$callsTable
		myAAvec <- callsTable$AA
		myAAvec[ is.na(myAAvec)] <- ""
		myDNAvec <- callsTable$DNA
		myDNAvec[ is.na(myDNAvec)] <- ""
		myAA <- paste( myAAvec, collapse="")
		myDNA <- paste( myDNAvec, collapse="")
		myDesc <- paste( sampleID, geneID, sep="_")
		writeFasta( as.Fasta( myDesc, myAA), consensusAAfile, line.width=100)
		writeFasta( as.Fasta( myDesc, myDNA), consensusDNAfile, line.width=100)
		write.table( callsTable, consensusBASEfile, sep="\t", quote=F, row.names=FALSE)
		nAA <- nchar( myAA)
	} else {
		fa <- loadFasta( consensusAAfile, verbose=F)
		nAA <- nchar( fa$seq)
	}

	# better luck drawing if we close the current graphics window first
	if ( dev.cur() > 1) dev.off()
	x11Width <- max( round( nAA/550 * 25), 10)
	x11Height <- max( round( max.depth/10, digits=1), 6)
	pdfWidth <- max( round( nAA/550 * 20), 10)
	pdfHeight <- max( round( max.depth/12, digits=1), 5)
	checkX11( bg='white', width=x11Width, height=x11Height)

	consensusAns <- proteinConstructPeptidePileups( sampleID, geneName=geneID, constructFile=consensusAAfile, 
				peptide.path=target.path, txt.cex=txt.cex, maxNoHits=0, 
				max.depth=max.depth, max.drawnPerSite=max.drawnPerSite, draw.box=draw.box,
				chunkSize=chunkSize.pileup, showFrameShiftPeptides=showFrameShiftPeptides)

	pdfFile <- file.path( target.path, paste( sampleID, geneFileID, "PeptidePileups.pdf", sep="."))
	dev.print( pdf, pdfFile, width=pdfWidth, height=pdfHeight)

	consensusSaveFile <- file.path( target.path, paste( sampleID, geneFileID, "ConsensusAnswer.rda", sep="."))
	save( consensusAns, file=consensusSaveFile)
	
	# step 3:  summarize how the consensus differs from what it used as its construct
	differences <- proteinConstructPileupSummary( consensusSaveFile, sampleID=sampleID, geneName=geneID, 
				constructName=paste(sampleID,geneFileID,sep="."), txt.cex=txt.cex, doPlot=FALSE)
	return( invisible( differences))
}


`pipe.GatherTargetAlignments` <- function( sampleID, genes, targetName="Virus.Genes", optionsFile="Options.txt", 
				results.path=NULL, tail.width=0, mode=c("all","best.one"), geneLength=10000, 
				asFASTQ=FALSE, fastq.keyword="Genes", verbose=TRUE) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	NG <- length( genes)

	# decide if we want all reads or just best one from each location
	mode <- match.arg( mode)

	target.path <- file.path( results.path, "TargetSearch", targetName)
	bamFile <- paste( sampleID, targetName, "Hits.bam", sep=".")
	bamFile <- file.path( target.path, bamFile)

	outrefid <- outpos <- outncig <- outcig <- outflag <- outseq <- outqual <- vector()
	outname <- outrev <- outsize <- outgid <- vector()
	nout <- 0

	# make sure we have that BAM file sorted and indexed
	cat( "\nFile: ", basename( bamFile))
	bamf <- BAM.verifySorted( bamFile, verbose=F, threads=2)
	if ( is.null( bamf)) {
		cat( "\n  Problem: BAM file not found: ", bamFile)
		return(NULL)
	}
	bamidx <- paste( bamf, "bai", sep=".")
	reader <- bamReader( bamf, indexname=bamidx)
	refData <- getRefData( reader)

	# visit every gene we were given
	for ( ig in 1:NG) {
		thisGene <- genes[ig]
		# extract the chunk of reads for this gene's loci
		refid <- seqID2refID( thisGene, refData=refData)
		chunk <- bamRange( reader, coords=c(refid, 1, geneLength))
		if ( size(chunk) < 1) next
	
		smallDF <- as.data.frame( chunk)
		if (asFASTQ) {
			smallDF$seq <- readSeq( chunk)
			smallDF$qual <- readQual( chunk)
		}
		smallDF$geneid <- thisGene
		 
		# if we want the raw reads, don't keep MARs
		if (asFASTQ) {
			dups <- which( duplicated( smallDF$name))
			if ( length(dups) > 0) {
				smallDF <- smallDF[ -dups, ]
			}
		}

		if ( nrow(smallDF) && mode == "best.one") {
			# for tools like denovo assembly, limit the read depth by keeping just one copy per
			# location.  Use count if possible, else random
			posFac <- factor( smallDF$position)
			keep <- tapply( 1:nrow(smallDF), posFac, function(x) {
					if ( length(x) == 1) return(x)
					seqTbl <- sort( table( smallDF$seq[x]), decreasing=T)
					bestCnt <- seqTbl[1]
					isBest <- which( seqTbl == bestCnt)
					if ( length(isBest) > 1) isBest <- sample( isBest, size=1)
					# return the first row that has this chosen sequence
					bestHit <- match( names(seqTbl)[isBest], smallDF$seq[x])
					return( x[bestHit])
				})
			smallDF <- smallDF[ keep, ]
		}

		if ( verbose) cat( "\n", thisGene, "\tN_Alignments: ", nrow(smallDF))

		if ( nrow(smallDF) < 1) next

		now <- (nout + 1) : (nout + nrow(smallDF))
		outrefid[now] <- smallDF$refid
		outpos[now] <- smallDF$position
		outncig[now] <- smallDF$nCigar
		outcig[now] <- smallDF$cigar
		outflag[now] <- smallDF$flag
		outseq[now] <- smallDF$seq
		outqual[now] <- smallDF$qual
		outname[now] <- smallDF$name
		outrev[now] <- smallDF$revstrand
		outsize[now] <- smallDF$insertsize
		outgid[now] <- smallDF$geneid
		nout <- max( now)
	}
	bamClose( reader)
	
	if ( verbose) cat( "\nTotal Alignments: ", nout)

	# put into chromosomal order
	cat( "\nSorting alignments..")
	ord <- order( outrefid, outpos)
	outrefid <- outrefid[ ord]
	outpos <- outpos[ ord]
	outncig <- outncig[ ord]
	outcig <- outcig[ ord]
	outflag <- outflag[ ord]
	outseq <- outseq[ ord]
	outqual <- outqual[ ord]
	outname <- outname[ ord]
	outrev <- outrev[ ord]
	outsize <- outsize[ ord]
	outgid <- outgid[ ord]
	out <- data.frame( "refid"=outrefid, "position"=outpos, "nCigar"=outncig, "cigar"=outcig,
				"flag"=outflag, "seq"=outseq, "qual"=outqual, "name"=outname,
				"revstrand"=outrev, "insertsize"=outsize, "geneid"=outgid, 
				stringsAsFactors=FALSE)
	if (nrow(out)) rownames(out) <- 1:nrow(out)
	cat( "  Done.")

	if ( asFASTQ) {
		if (verbose) cat( "\nConverting Alignments back to FASTQ..")
		outfile <- paste( sampleID, fastq.keyword, "fastq.gz", sep=".")
		outfile <- file.path( target.path, outfile)
		fqDF <- data.frame( "READ_ID"=out$name, "READ_SEQ"=out$seq, "SCORE"=out$qual,
				stringsAsFactors=FALSE)

		# there may be duplicate readIDs, that mapped to more than one location in the genome
		# don't let them be written out more than once...
		dups <- which( duplicated( fqDF$READ_ID))
		if ( length(dups) > 0) {
			if (verbose) cat( "\nDropping redundant MAR alignments from FASTQ: ", length(dups))
			fqDF <- fqDF[ -dups, ]
		}
		writeFastq( fqDF, outfile, compress=T)
		return( nrow(fqDF))
	} else {
		return( out)
	}
}



`pipe.TargetProteinExtraction` <- function( sampleID, geneID, targetName="Virus.Genes", optionsFile="Options.txt",
					results.path=NULL, max.proteins=NULL, min.minor.pct=6, 
					min.mutation.pct=1.0, min.minor.read.count=2, verbose=TRUE ) {
						
	# min.minor.pct - at any one amino acid, how frequent must a minor AA call be to not be treated as just noise
	# min.mutation.pct -  for the full length protein, what fraction of AA must be mutated to call it a separate protein call

	# make sure we were given valid parameters..
	if ( length(sampleID) > 1) {
		sampleID <- sampleID[1]
		cat( "\nWarning:  only one sampleID allowed")
	}
	if ( is.null(results.path)) results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	target.path <- file.path( results.path, "TargetSearch", targetName)
	if ( ! file.exists( target.path)) {
		cat( "\nNo 'Target Search' results found for target: ", targetName)
		return(NULL)
	}

	# build the filename of the expected results, and get the consensus protein call details 
	geneName <- geneID
	geneFileName <- file.cleanSpecialCharactersFromFileName( geneName)

	summaryFile <- file.path( target.path, paste( sampleID, geneFileName, "ConsensusProteinSummary.txt", sep="."))						
	if ( ! file.exists( summaryFile)) {
		cat( "\nProtein summary file not found: ", summaryFile)
		return(NULL)
	}
	protDF <- read.delim( summaryFile, as.is=T)
	NAA <- nrow( protDF)
	
	# get the consensus protein(s) themselves, and extract the initial guess about proportions from their names
	topProteins <- cpp.ExtractTopProteins( summaryFile, min.minor.pct=min.minor.pct, min.mutation.pct=min.mutation.pct,
						drop.gaps=FALSE)
	NPROT <- length( topProteins)
	protNames <- names( topProteins)

	# allow the user to put a upper limit on how many proteins we consider
	if ( ! is.null( max.proteins)) {
		topProteins <- topProteins[ 1:min(NPROT,max.proteins)]
		NPROT <- length( topProteins)
		protNames <- names( topProteins)
	}
	protProportions <- rep.int( NA, NPROT)
	protProportions[1] <- 100
	if (NPROT > 1) {
		for ( i in 2:NPROT) protProportions[i] <- as.numeric( sub( "(^MinorVariant.+Pct=)([0-9]+)(%$)", "\\2", protNames[i]))
		protProportions[1] <- 100 - sum( protProportions[2:NPROT], na.rm=T)
	}
	if (verbose) {
		cat( "\nInitial Proportions: \n")
		print( paste( protNames, protProportions, sep=" = "))
	}
	# as defined, all these constructs should be the same length.  Make sure
	if ( any( nchar(topProteins) != NAA)) {
		stop( "\nSize error:  consensus protein lengths not all equal..")
	}
	# turn all the protein calls into one AA matrix
	protV <- strsplit( topProteins, split="")
	aaM <- matrix( unlist(protV), nrow=NAA, ncol=NPROT)
	colnames(aaM) <- protNames
	
	# prep all the protein calls and all the summary percentages & counts
	protV <- strsplit( topProteins, split="")
	pcts <- protDF$Percentages
	cnts <- protDF$Counts
	pctTerms <- strsplit( pcts, split="; ")
	cntTerms <- strsplit( cnts, split="; ")
	pctAAs <- lapply( pctTerms, function(x) sub( ":[0-9]+", "", x))
	pctVals <- lapply( pctTerms, function(x) as.numeric( sub( "[-*A-Z]:", "", x)))
	cntAAs <- lapply( cntTerms, function(x) sub( ":[0-9]+", "", x))
	cntVals <- lapply( cntTerms, function(x) as.numeric( sub( "[-*A-Z]:", "", x)))
	
	# we can tell who is conserved or not by how many AA were called
	nAAterms <- sapply( pctAAs, length)
	useForScoring <- which( nAAterms > 1)
	
	# we will ignore very small read counts as being not seen enough to be considered a true minor variant
	for ( k in useForScoring) {
		myCnts <- cntVals[[k]]
		tooLow <- which( myCnts < min.minor.read.count)
		if ( length(tooLow)) {
			# remove these from both percents and counts
			cntVals[[k]] <- cntVals[[k]][ -tooLow]
			cntAAs[[k]] <- cntAAs[[k]][ -tooLow]
			pctVals[[k]] <- pctVals[[k]][ -tooLow]
			pctAAs[[k]] <- pctAAs[[k]][ -tooLow]
			# re-evaluate the percentages after this trim of low count minor variants
			pctVals[[k]] <- round( cntVals[[k]] * 100 / sum( cntVals[[k]], na.rm=T))
		}
	}
	
	# now we need to reassess who is conserved or not by how many AA were called
	nAAterms <- sapply( pctAAs, length)
	useForScoring <- which( nAAterms > 1)
	if (verbose) cat( "\nN AA with minor variants: ", length(useForScoring))
	
	# use one universal set of all possible AA for doing our compares and distance scoring
	ALL.AA <- sort( unique( c( unlist( pctAAs), unlist(cntAAs))))
	
	# prebuild the static distributions of the called AAs
	aaPctTables <- vector( mode="list", length=NAA)
	for ( k in useForScoring) {
		myAA <- rep( pctAAs[[k]], times=pctVals[[k]])
		myTbl <- table( factor(myAA, levels=ALL.AA))
		aaPctTables[[k]] <- myTbl
	}
	
	
	# local function that scores how well the consensus proteins match/explain the call details
	# puts some results into global storage for speed / later use...
	aaDist <- rep.int( 0, NAA)
	aaM.out <- aaM
	
	consensus.AA.Distance <- function( aaM, protPcts) {
	
		# given one matrix of all the protein AA calls for every location, and the percentage 
		# that each protein contributes to the total, see how this current consensus and percentage
		# fits to the called data from the CPP tool
		protPcts <- as.integer( protPcts)

		# visit each AA that has any minor forms, and see how distant is current from what was called
		for ( k in useForScoring) {
			myAA <- rep( aaM[ k, ], times=protPcts)
			myTbl <- table( factor(myAA, levels=ALL.AA))
			calledTbl <- aaPctTables[[k]]
			myDistTbl <- myTbl - calledTbl
			myDist <- sum( abs( myDistTbl), na.rm=T) 
			# store these where the caller can use them later
			aaDist[k] <<- myDist 
		}
		return( sum(aaDist) / nrow(aaM))
	}
	
	
	optimize.Proportions <- function( aaM, protPcts, verbose=T) {
	
		# given one matrix of all the protein AA calls for every location, and the percentage 
		# that each protein contributes to the total, optimize those percentages to give the best
		# fit to the called data
		
		# ready to do the optimization.   Evaluate at the initial estimate
		startingDist <- consensus.AA.Distance( aaM, protProportions)
		if (verbose) cat( "\nInitial Distance: ", startingDist)

		bestDist <- startingDist
		bestProportion <- protProportions
		iterCount <- 0
		while ( TRUE ) {
			# step over all alternative pairs of proteins
			thisProportion <- lastProportion <- bestProportion
			iterCount <- iterCount + 1
			bestI <- bestJ <- 0
			for (i in 1:NPROT) {
				for ( j in 1:NPROT) {
					if ( i == j) next
					# tweak the proportions between all 2 choices of protein
					thisProportion <- lastProportion
					thisProportion[i] <- thisProportion[i] + 1
					thisProportion[j] <- thisProportion[j] - 1
					# never let a proportion go to zero
					if ( any( thisProportion < 1)) next
					# rescore the distance
					thisDist <- consensus.AA.Distance( aaM, thisProportion)
					if ( thisDist < bestDist) {
						bestDist <- thisDist
						bestProportion <- thisProportion
						bestI <- i
						bestJ <- j
					}
				}
			}
			# we did all possible steps, was any better?
			if ( bestI == 0 || bestJ == 0) break
			if (verbose) cat( "\nIter: ", iterCount, "  Dist: ", bestDist, "  Proportions: ", bestProportion)
		}
		return( list( "Distance"=bestDist, "Proportions"=bestProportion))
	}


	optimize.ResidueCalls <- function( aaM, protPcts, verbose=T) {
	
		# given one matrix of all the protein AA calls for every location, and the percentage 
		# that each protein contributes to the total, optimize the AA residues between the proteins
		# to give the best fit to the called data
		
		# the current distance for each minor site already sits in global storage...

		aaM.out <- aaM
		dist.out <- aaDist
		
		# visit each AA that has any minor forms, and see how distant is affected by permuting the residues
		for ( k in useForScoring) {
		
			# start with the set of AA the proteins were on function entry
			myAAorig <- aaM[ k, ]
			aaSet <- sort( unique( myAAorig))
			nChoice <- length( aaSet)
			nNeed <- ncol(aaM)
			
			# step over all possible combinations of AA
			bestDist <- origDist <- aaDist[k]
			setPtr <- rep.int( 1, nNeed)
			bestAAtry <- NULL
			while (TRUE) {
				# use the set pointers to make this possible choice of AA calls
				myAAtry <- aaSet[ setPtr]
				myAA <- rep( myAAtry, times=protPcts)
				myTbl <- table( factor(myAA, levels=ALL.AA))
				calledTbl <- aaPctTables[[k]]
				myDistTbl <- myTbl - calledTbl
				thisDist <- sum( abs( myDistTbl), na.rm=T) 
				if ( thisDist < bestDist) {
					bestDist <- thisDist
					bestAAtry <- myAAtry
				}
				# now iterate by incrementing the pointers to eventually try all possible combinations
				setPtr[1] <- setPtr[1] + 1
				for (ich in 1:(nNeed-1)) {
					if (setPtr[ich] > nChoice) {
						setPtr[ich] <- 1
						setPtr[ich+1] <- setPtr[ich+1] + 1
					}
				}
				# we are done when the very last pointer is too big
				if ( setPtr[nNeed] > nChoice) break
			}
			# we did all possible steps, was any better?
			if ( is.null(bestAAtry)) next
			if (verbose && any( myAAorig != bestAAtry)) cat( "\nAA Swap: ", k, "  OldDist:", origDist, "  OldAA:", myAAorig, 
									"|  NewDist:", bestDist, "  NewAA:", bestAAtry)
			
			# update the proteins
			aaM.out[ k, ] <- bestAAtry
			dist.out[k] <- bestDist
		}
		
		newDist <- sum( dist.out) / nrow(aaM)
		return( list( "Distance"=newDist, "AA_Matrix"=aaM.out))
	}
	
	
	topDistSites <- function( aaM, aaDist, cutoffDist, nExtraAA=2) {
	
		# gather the top residual distance sites as extra info
		keep <- which( aaDist > cutoffDist)
		if ( ! length(keep)) return(NULL)
		badDist <- aaDist[ keep]
		badMotif <- aaM[ keep, 1]
		for ( i in 1:length(keep)) {
			k <- keep[i]
			badMotif[i] <- paste( aaM[(max(1,k-nExtraAA)):(min(nrow(aaM),k+nExtraAA)), 1], collapse="")
		}
		sml <- data.frame( "AA.Position"=keep, "AA.Motif"=badMotif, "Residual.Distance"=badDist, stringsAsFactors=F)
		ord <- order( sml$Residual.Distance, decreasing=T)
		sml <- sml[ ord, ]
		rownames(sml) <- 1:nrow(sml)
		return(sml)
	}

	
	# ready to do the work...
	# Step 1:  optimize the starting proportions, given what the extractor guessed
	myProportions <- protProportions
	ans1 <- optimize.Proportions( aaM, myProportions, verbose=verbose)
	myProportions <- ans1$Proportions
	myDist <- ans1$Distance
	# note that the optimization may have decreased the number of proteins!...
	NPROT <- length( myProportions)
	
	if ( NPROT > 1) {
		# Step 2:  visit each minor mutation site, and see if swapping residues improves the distance scoring
		ans2 <- optimize.ResidueCalls( aaM, myProportions, verbose=verbose)
		aaM.out <- ans2$AA_Matrix
	
		# Step 3:  re-optimize the  proportions, given the improved residue calls
		ans3 <- optimize.Proportions( aaM.out, myProportions, verbose=verbose)
		myProportions <- ans3$Proportions
		myDist <- ans3$Distance
	}
	
	# make the final answer
	finalProportions <- myProportions
	finalDist <- myDist
	finalAAseqs <- apply( aaM.out, MARGIN=2, FUN=function(x) paste( x, collapse=""))
	finalNames <- rep.int( paste( sampleID, geneName, sep="_"), NPROT)
	finalNames[1] <- paste( finalNames[1], "_Consensus", sep="")
	if (NPROT > 1) {
		finalNames[2:NPROT] <- paste( finalNames[2:NPROT], "_Variant", 2:NPROT, sep="")
	}
	finalSuffix <- paste( round( finalProportions), "%", sep="")
	finalNames <- paste( finalNames, finalSuffix, sep="_")
	
	if (verbose) cat( "\nFinal Extraction: \n", finalNames, "\nFinal Distance: ", finalDist)

	outSeqs <- as.Fasta( finalNames, finalAAseqs)
	finalFile <- file.path( target.path, paste( sampleID, geneFileName, "FinalExtractedAA.fasta", sep="."))
	writeFasta( outSeqs, finalFile, line=100)

	alnFile <- file.path( target.path, paste( sampleID, geneFileName, "FinalExtractedAA.aln", sep="."))
	if ( NPROT > 1) {
		aln <- mafft( finalFile, alnFile)
		writeALN( aln, alnFile, line=100)	
	} else {
		# make sure any earlier result gets removed
		file.delete( alnFile)
	}
	
	ans4 <- topDistSites( aaM.out, aaDist, finalDist*2)
	if ( ! is.null(ans4)) {
		residSiteFile <- file.path( target.path, paste( sampleID, geneFileName, "FinalExtracted.ResidualSites.csv", sep="."))
		write.table( ans4, residSiteFile, sep=",", quote=T, row.names=F)	
	}
	
	out <- list( "AA.Fasta"=outSeqs, "Residual"=finalDist, "Problem.Sites"=ans4)
	return( invisible(out))
}
	
