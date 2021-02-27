# pipe.VirusSearch.R -- mine No Hits for evidence of Viruses

`pipe.VirusSearch` <- function( sampleIDset, virusTargetName="Virus.Genes", max.mismatch.per.100=5,
				max.low.complexity=80, annotationFile="Annotation.txt", 
				optionsFile="Options.txt", verbose=FALSE) {

	# get options that we need
	optT <- readOptionsTable( optionsFile)
	results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=verbose)
	fastq.path <- getOptionValue( optT, "fastqData.path", verbose=verbose)
	bowtie.path <- getOptionValue( optT, "bowtie2Index.path", notfound=".", verbose=verbose)

	# make sure we see and can load the virus data
	viralFastaFile <- file.path( bowtie.path, "FastaFiles", paste( virusTargetName, "fasta", sep="."))
	if ( ! file.exists( viralFastaFile)) stop( paste( "Required Virus FASTA file not found: ", viralFastaFile))
	viralFA <- loadFasta( viralFastaFile, short.desc=FALSE, verbose=verbose)
	descTerms <- strsplit( viralFA$desc, split=" | ", fixed=T)
	viralIDs <- sapply( descTerms, `[`, 1)
	viralNames <- sapply( descTerms, `[`, 2)
	if ( any( is.na( viralNames))) cat( "\nWarning:  some Fasta headers missing ' | ' product term delimitors")

	# same with the viral Bowtie index
	bowtieFile <- file.path( bowtie.path, paste( virusTargetName, "1.bt2", sep="."))
	if ( ! file.exists( bowtieFile)) stop( paste( "Required", virusTargetName, "Bowtie2 target index file not found: ", bowtieFile))

	# make the place/file for our results
	bam.path <- file.path( results.path, paste( virusTargetName, "BAMS", sep="."))
	if ( ! file.exists( bam.path)) dir.create( bam.path, recursive=T)

	# since we don't expect much hits (if at all), allow doing a set of samples into one result
	bigDF <- data.frame()
	nRejectMismatch <- nRejectLowComplex <- nFound <- 0

	for (sid in sampleIDset) {

		cat( "\n\nStarting '", virusTargetName, "' Search for Sample:    ", sid, 
			"\n  Max base mismatches per 100bp:  ", max.mismatch.per.100, 
			"\n  Max low complexity base %:      ", max.low.complexity, sep="")
	
		# get the file of noHit reads, as viral reads should still be in there
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
		if ( ! file.exists( infile)) {
			cat( "\nFASTQ file not found:  ", infile)
			next
		}

		# call Bowtie
		cat( "\n  Aligning..")
		bamFile <- file.path( bam.path, paste( sid, virusTargetName, "Hits.bam", sep="."))
		ans1 <- fastqToBAM( inputFastqFile=infile, outputFile=bamFile, k=1, sampleID=sid, optionsFile=optionsFile,
					alignIndex=virusTargetName, index.path=bowtie.path, keepUnaligned=FALSE, 
					verbose=verbose, label=sid)
		nReadsIn <- ans1$RawReads

		# now see what got hit
		cat( "  Scan for high quality hits..")
		reader <- bamReader( bamFile)
		refData <- getRefData( reader)
		viralHits <- vector()
		nSeen <- 0
		repeat {
			chunk <- getNextChunk( reader, n=10000, alignedOnly=T)
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
			viralHits <- c( viralHits, seqID[goodHits])
			nFound <- nFound + length(goodHits)
			nRejectMismatch <- nRejectMismatch + length(badMatch)
			nRejectLowComplex <- nRejectLowComplex + length(badComplex)
		}
		bamClose( reader)
	
		if ( ! length(viralHits)) {
			cat( "\nNo high quality", virusTargetName, "hits detected..")
			next
		}

		# now tablify
		hitTbl <- table( viralHits)
		hitCnts <- as.numeric( hitTbl)
		pctHits <- round( hitCnts * 100 / sum(hitCnts), digits=2)
		perMilNoHits <- round( hitCnts / (nReadsIn/1000000), digits=4)
		where <- match( names(hitTbl), viralIDs)
		products <- viralNames[where]

		smlDF <- data.frame( "SampleID"=sid, "VirusID"=names(hitTbl), "Product"=products, 
				"Percent.Sample"=pctHits, "N_Hits"=hitCnts, "PerMillionNoHits"=perMilNoHits, 
				stringsAsFactors=F)
		ord <- order( hitCnts, decreasing=T)
		smlDF <- smlDF[ ord, ]
		rownames(smlDF) <- 1:nrow(smlDF)
		colnames(smlDF)[2] <- paste( virusTargetName, "ID", sep=".")
		outfile <- file.path( bam.path, paste( sid, virusTargetName, "Summary.csv", sep="."))
		write.table( smlDF, outfile, sep=",", quote=T, row.names=F)

		cat( "\nDone.  N_Aligned: ", nSeen, "  N_Good_Hits: ", sum(hitTbl), "\n")
		print( head( smlDF, 3))

		# accumulate the results
		bigDF <- rbind( bigDF, smlDF)
	}
	if ( ! nrow(bigDF)) return( bigDF)

	# final ordering is...
	ord <- order( bigDF$N_Hits, decreasing=T)
	bigDF <- bigDF[ ord, ]
	rownames(bigDF) <- 1:nrow(bigDF)

	# report the overall picture
	cat( "\n\nSummary:")
	cat( "\n  N Good Virus Reads Found:      ", nFound)
	cat( "\n  N Rejected for Mismatches:     ", nRejectMismatch)
	cat( "\n  N Rejected for Low Complexity: ", nRejectLowComplex, "\n")

	return( bigDF)
}

