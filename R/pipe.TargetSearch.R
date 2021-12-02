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
	bam.path <- file.path( results.path, "TargetSearch", paste( targetName, "BAMS", sep="."))
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
		ord <- order( hitCntsM, decreasing=T)
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
	ord <- order( bigDF$READS_M, decreasing=T)
	bigDF <- bigDF[ ord, ]
	rownames(bigDF) <- 1:nrow(bigDF)

	# report the overall picture
	cat( "\n\nSummary:")
	cat( "\n  N Good Target Reads Found:      ", nFound)
	cat( "\n  N Rejected for Mismatches:     ", nRejectMismatch)
	cat( "\n  N Rejected for Low Complexity: ", nRejectLowComplex, "\n")

	return( bigDF)
}

