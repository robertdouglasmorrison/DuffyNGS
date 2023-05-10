# pipe.GatherGeneAlignments.R -- collect up the reads that align to some genes, and 
#				optionally repackage in their original FASTQ format

`pipe.GatherGeneAlignments` <- function( sampleID, genes, optionsFile="Options.txt", 
				results.path=NULL, tail.width=0, 
				stages=c("genomic", "splice"), mode=c("all","best.one"), 
				asFASTQ=FALSE, fastq.keyword="Genes", verbose=TRUE) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}

	NG <- length( genes)
	gmap <- getCurrentGeneMap()
	smap <- getCurrentSeqMap()
	where <- match( genes, gmap$GENE_ID, nomatch=0)
	if ( any( where == 0)) {
		cat( "\nSome genes not found in current species: ", genes[ where == 0])
		where <- where[ where > 0]
		genes <- genes[ where]
		NG <- length( genes)
	}
	gptrs <- where

	# determine the set of BAM files to visit
	Stages <- c( "riboClear", "genomic", "splice")
	if ( ! all( stages %in% Stages)) {
		cat( "\nAllowed pipeline stages: ", Stages)
		stop()
	}

	# decide if we want all reads or just best one from each location
	mode <- match.arg( mode)

	bamFiles <- vector()
	if ( "riboClear" %in% stages) {
		bamFiles <- c( bamFiles, file.path( results.path, "riboClear", 
				paste( sampleID, ".ribo.converted.bam", sep="")))
	}
	if ( "genomic" %in% stages) {
		bamFiles <- c( bamFiles, file.path( results.path, "align", 
				paste( sampleID, ".genomic.bam", sep="")))
	}
	if ( "splice" %in% stages) {
		bamFiles <- c( bamFiles, file.path( results.path, "splicing", 
				paste( sampleID, ".splice.converted.bam", sep="")))
	}

	outrefid <- outpos <- outncig <- outcig <- outflag <- outseq <- outqual <- vector()
	outname <- outrev <- outsize <- outgid <- outstage <- vector()
	nout <- 0

	for ( f in bamFiles) {

		thisStage <- "genomic"
		if ( regexpr( "ribo", f) > 0) thisStage <- "riboClear"
		if ( regexpr( "splice", f) > 0) thisStage <- "splice"

		# make sure we have that BAM file sorted and indexed
		if (verbose) cat( "\nFile: ", basename(f))
		bamf <- BAM.verifySorted( f, verbose=F, threads=2)
		if ( is.null( bamf)) next
		bamidx <- paste( bamf, "bai", sep=".")

		reader <- bamReader( bamf, indexname=bamidx)
		refData <- getRefData( reader)

		# visit every gene we were given
		for ( ig in 1:NG) {
			sml <- gmap[ gptrs[ ig], ]

			# extract the chunk of reads for this gene's loci
			refid <- seqID2refID( sml$SEQ_ID, refData=refData)
			start <- sml$POSITION - tail.width
			end <- sml$END + tail.width
			# quietly catch case of that chromosome not in the data
			if (refid < 0) next

			# force clip to this chromosome's limits
			start <- max( start, 1)
			end <- min( end, smap$LENGTH[ match( sml$SEQ_ID, smap$SEQ_ID)])
			chunk <- bamRange( reader, coords=c(refid, start, end))
			if ( size(chunk) < 1) next
	
			smallDF <- as.data.frame( chunk)
			if (asFASTQ) {
				smallDF$seq <- readSeq( chunk)
				smallDF$qual <- readQual( chunk)
			}
			smallDF$geneid <- sml$GENE_ID
			smallDF$stage <- thisStage
			 
			# splices have the reads broken by the splice junction, and a modified readID
			# we need to rebuild the originals
			if ( thisStage == "splice") {
				smallDF <- rejoinSplicedReads( smallDF)
			}

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

			if ( verbose) cat( "\n", sml$GENE_ID, "\tN_Alignments: ", nrow(smallDF))

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
			outstage[now] <- smallDF$stage
			nout <- max( now)
		}
		bamClose( reader)
	}
	if ( verbose) cat( "\nTotal Alignments: ", nout)

	# put into chromosomal order
	if (verbose) cat( "\nSorting alignments..")
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
	outstage <- outstage[ ord]
	out <- data.frame( "refid"=outrefid, "position"=outpos, "nCigar"=outncig, "cigar"=outcig,
				"flag"=outflag, "seq"=outseq, "qual"=outqual, "name"=outname,
				"revstrand"=outrev, "insertsize"=outsize, "geneid"=outgid, 
				"stage"=outstage, stringsAsFactors=FALSE)
	if (nrow(out)) rownames(out) <- 1:nrow(out)
	cat( "  Done.")

	if ( asFASTQ) {
		if (verbose) cat( "\nConverting Alignments back to FASTQ..")
		outfile <- paste( sampleID, fastq.keyword, "fastq.gz", sep=".")
		outfile <- file.path( results.path, "fastq", outfile)
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


`pipe.GatherRegionAlignments` <- function( sampleID, seqids, starts, stops,
				optionsFile="Options.txt", 
				results.path=NULL, stages=c("genomic", "splice"), 
				asFASTQ=FALSE, fastq.keyword="Region", verbose=TRUE) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}

	# determine the set of BAM files to visit
	Stages <- c( "riboClear", "genomic", "splice")
	if ( ! all( stages %in% Stages)) {
		cat( "\nAllowed pipeline stages: ", Stages)
		stop()
	}

	bamFiles <- vector()
	if ( "riboClear" %in% stages) {
		bamFiles <- c( bamFiles, file.path( results.path, "riboClear", 
				paste( sampleID, ".ribo.converted.bam", sep="")))
	}
	if ( "genomic" %in% stages) {
		bamFiles <- c( bamFiles, file.path( results.path, "align", 
				paste( sampleID, ".genomic.bam", sep="")))
	}
	if ( "splice" %in% stages) {
		bamFiles <- c( bamFiles, file.path( results.path, "splicing", 
				paste( sampleID, ".splice.converted.bam", sep="")))
	}

	# can have more than one region...
	nRegions <- length( starts)
	if (length(stops) != nRegions) stop( "'starts' and 'stops' must be of same length")
	if (length(seqids) < nRegions) seqids <- rep( seqids, length.out=nRegions)

	out <- data.frame()

	for ( f in bamFiles) {

		thisStage <- "genomic"
		if ( regexpr( "ribo", f) > 0) thisStage <- "riboClear"
		if ( regexpr( "splice", f) > 0) thisStage <- "splice"

		# make sure we have that BAM file sorted and indexed
		bamf <- BAM.verifySorted( f, verbose=FALSE, threads=2)
		if ( is.null( bamf)) next
		bamidx <- paste( bamf, "bai", sep=".")

		reader <- bamReader( bamf, indexname=bamidx)
		refData <- getRefData( reader)

		# visit this region
		for ( iregion in 1:nRegions) {
		 	seqid <- seqids[iregion]
		 	start <- starts[iregion]
		 	stop <- stops[iregion]

			# extract the chunk of reads for this gene's loci
			refid <- seqID2refID( seqid, refData=refData)
			chunk <- bamRange( reader, coords=c(refid, start, stop))
			if ( size(chunk) < 1) next
		
			smallDF <- as.data.frame( chunk)
			if ( asFASTQ) {
				smallDF$seq <- readSeq( chunk)
				smallDF$qual <- readQual( chunk)
			}
			smallDF$stage <- thisStage
	
			# splices have the reads broken by the splice junction, and a modified readID
			# we need to rebuild the originals
			if ( thisStage == "splice") {
				smallDF <- rejoinSplicedReads( smallDF)
			}


			# if we want the raw reads, don't keep MARs
			if (asFASTQ) {
				dups <- which( duplicated( smallDF$name))
				if ( length(dups) > 0) {
					smallDF <- smallDF[ -dups, ]
				}
			}

			if ( verbose) cat( "\n", basename(f), "\nSeqID, Start, Stop: ", seqid, start, stop, 
					"\tN_Alignments: ", nrow(smallDF))

			if ( nrow(smallDF) < 1) next

			out <- rbind( out, smallDF)
		}
		bamClose( reader)
	}
	if ( verbose) cat( "\nTotal Aignments: ", nrow(out))

	# put into chromosomal order
	ord <- order( out$seq, out$position)
	out <- out[ ord, ]
	rownames(out) <- 1:nrow(out)

	if ( asFASTQ) {
		if (verbose) cat( "\nConverting Alignments back to FASTQ..")
		outfile <- paste( sampleID, fastq.keyword, "fastq.gz", sep=".")
		outfile <- file.path( results.path, "fastq", outfile)
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


`rejoinSplicedReads` <- function( tbl) {

	# given a data frame of alignments from a splice BAM file, put the halve back together
	if ( ! all( c( "position", "seq", "name", "qual") %in% colnames(tbl))) stop( "Not given a splice BAM alignment data frame")

	posIn <- tbl$position
	nameIn <- tbl$name
	seqIn <- tbl$seq
	qualIn <- tbl$qual

	# all the first halves say 'splice1'
	isFront <- grep( "::splice1", nameIn, fixed=T)
	isBack <- grep( "::splice2", nameIn, fixed=T)

	# grab all the partial items we will need
	nameOut <- sub( "::splice[12]", "", nameIn)
	nameFront <- nameOut[ isFront]
	nameBack <- nameOut[ isBack]

	# to be a usable read, we need to see both halves of the same readID
	frontHitsBack <- match( nameFront, nameBack, nomatch=0)
	keepers <- isFront[ frontHitsBack > 0]
	keepBack <- isBack[ frontHitsBack]

	# grab that subset of the given table as the result, then update the bits that need it
	out <- tbl[ keepers, ]
	out$name <- nameOut[ keepers]
	out$seq <- paste( seqIn[keepers], seqIn[keepBack], sep="")
	out$qual <- paste( qualIn[keepers], qualIn[keepBack], sep="")

	# all done, just these pairs that got resolved go back
	out
}
