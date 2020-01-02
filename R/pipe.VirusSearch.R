# pipe.VirusSearch.R -- mine No Hits for evidence of Viruses

`pipe.VirusSearch` <- function( sampleIDset, annotationFile="Annotation.txt",
		optionsFile="Options.txt", virusTargetName="Virus.Genes", verbose=FALSE) {

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

	# same with the viral Bowtie index
	bowtieFile <- file.path( bowtie.path, paste( virusTargetName, "1.bt2", sep="."))
	if ( ! file.exists( bowtieFile)) stop( paste( "Required Virus Bowtie2 target index file not found: ", bowtieFile))

	# make the place/file for our results
	bam.path <- file.path( results.path, "Viral.BAMS")
	if ( ! file.exists( bam.path)) dir.create( bam.path, recursive=T)

	# since we don't expect much hits (if at all), allow doing a set of samples into one result
	bigDF <- data.frame()
	for (sid in sampleIDset) {

		cat( "\n\nStarting 'Virus Search' for Sample:     ", sid, "\n")
	
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
		bamFile <- file.path( bam.path, paste( sid, "ViralHits.bam", sep="."))
		ans1 <- fastqToBAM( inputFastqFile=infile, outputFile=bamFile, k=30, sampleID=sid, optionsFile=optionsFile,
					alignIndex=virusTargetName, index.path=bowtie.path, keepUnaligned=FALSE, 
					verbose=verbose, label=sid)
		nReadsIn <- ans1$RawReads

		# now see what got hit
		reader <- bamReader( bamFile)
		refData <- getRefData( reader)
		viralHits <- vector()
		repeat {
			chunk <- getNextChunk( reader, n=10000, alignedOnly=T)
			nReads <- size( chunk)
			if ( nReads < 1) break
			refIDs <- refID( chunk)
			seqID <- refID2seqID( refIDs, reader=reader, refData=refData)
			viralHits <- c( viralHits, seqID)
		}
		bamClose( reader)
	
		if ( ! length(viralHits)) {
			cat( "\nNo Viral gene hits detected..")
			next
		}

		# now tablify
		hitTbl <- table( viralHits)
		hitCnts <- as.numeric( hitTbl)
		pctHits <- round( hitCnts * 100 / sum(hitCnts), digits=2)
		perMilNoHits <- round( hitCnts / (nReadsIn/1000000), digits=4)
		where <- match( names(hitTbl), viralIDs)
		products <- viralNames[where]

		smlDF <- data.frame( "SampleID"=sid, "VirusID"=names(hitTbl), "Product"=products, "N_Hits"=hitCnts,
				"PercentHits"=pctHits, "PerMillionNoHits"=perMilNoHits, stringsAsFactors=F)
		ord <- order( hitCnts, decreasing=T)
		smlDF <- smlDF[ ord, ]
		rownames(smlDF) <- 1:nrow(smlDF)
		outfile <- file.path( bam.path, paste( sid, "ViralHits.Summary.csv", sep="."))
		write.table( smlDF, outfile, sep=",", quote=T, row.names=F)

		cat( "\nDone.  N_Hits: ", sum(hitTbl), "\n")
		print( head( smlDF, 3))

		# accumulate the results
		bigDF <- rbind( bigDF, smlDF)
	}

	# final ordering is...
	ord <- order( bigDF$N_Hits, decreasing=T)
	bigDF <- bigDF[ ord, ]
	rownames(bigDF) <- 1:nrow(bigDF)

	return( bigDF)
}

