# BAM.fileTools.R -- wrappers to SAMTOOLS calls, for various BAM file manipulations


`BAM.sort` <- function( file, what=c("position", "readID"), index=TRUE, memory=2147000000, 
			threads=2, verbose=TRUE) {

	if ( ! file.exists( file)) {
		if (verbose) cat( "\nFile not found:  ", file, "\nFailed to find existing BAM file.")
		return(NULL)
	}

	bamOutFile <- sub( "bam$", "sorted.bam", file)
	prefix <- sub( ".bam", "", basename(file))

	what <- match.arg( what)
	sortOpt <- if ( what == "readID") " -n " else ""
	memValue <- as.integer(memory)
	if ( is.na( memValue)) memValue <- 500000000
	memOpt <- if ( is.null( memory)) "" else paste( " -m", as.integer(memValue))

	samtools <- Sys.which( "samtools")
	if ( samtools == "") stop( "Executable not found on search path:  'samtools'")
	cmdline <- paste( samtools, " sort ", sortOpt, memOpt, " -o", bamOutFile, " -T", prefix, " -@", threads, file)
	if (verbose) cat( "\nSorting BAM file (threads=",threads, "):   ", basename(file), sep="")
	system( cmdline)
	if (verbose) cat( "\nDone.")

	if ( index) BAM.index( bamOutFile, verbose=verbose)
	return( bamOutFile)
}


`BAM.index` <- function( file, verbose=TRUE) {

	if ( ! file.exists( file)) {
		if (verbose) cat( "\nFile not found:  ", file, "\nFailed to find existing BAM file.")
		return(NULL)
	}
	idxfile <- sub( "bam$", "bam.bai", file)
	samtools <- Sys.which( "samtools")
	if ( samtools == "") stop( "Executable not found on search path:  'samtools'")
	if (verbose) cat( "\nIndexing BAM file:  ", basename(file))
	cmdline <- paste( samtools, " index ", file)
	system( cmdline)
	if (verbose) cat( "\nDone.")
	return( idxfile)
}


`BAM.verifySorted` <- function( file, index=TRUE, threads=2, verbose=TRUE) {

	# make sure that all the names end up with the sort prefix
	file <- sub( "sorted.bam$", "bam", file)
	bamfile <- sub( "bam$", "sorted.bam", file)

	# we can allow the sorted file to exist without the original!!
	if ( ! file.exists( file)) {
		if ( ! file.exists( bamfile)) {
			if (verbose) cat( "\nFile not found:  ", file, "\nFailed to find existing BAM file.")
			return(NULL)
		}
	}

	# see if the 'sorted' one is there yet
	if ( ! file.exists( bamfile)) {
		ans <- BAM.sort( file, index=index, verbose=verbose, threads=threads)
	}

	# see if the 'index' is there too
	if (index) {
		idxfile <- sub( "bam$", "bam.bai", bamfile)
		if ( ! file.exists( idxfile)) ans <- BAM.index( bamfile, verbose=verbose)
	}
	return( bamfile)
}


`BAM.merge` <- function( files, newfile, index=TRUE, verbose=TRUE) {

	allFound <- TRUE
	for ( file in files) {
		if ( ! file.exists( file)) {
			if (verbose) cat( "\nFile not found:  ", file, "\nFailed to find existing BAM file.")
			allFound <- FALSE
		}
		if ( ! allFound) return(NULL)
	}

	if ( regexpr( "sorted.bam$", newfile) < 1) {
		newfile <- paste( newfile, "sorted.bam", sep=".")
	}

	samtools <- Sys.which( "samtools")
	if ( samtools == "") stop( "Executable not found on search path:  'samtools'")
	cmdline <- paste( samtools, " merge -f -c ", newfile, " ", paste( files, collapse=" "))
	if (verbose) cat( "\nMerging BAM files:  ", basename(files))
	system( cmdline)
	if (verbose) cat( "\nDone.  Created new merged file: ", basename(newfile), "\n")

	if ( index) BAM.index( newfile, verbose=verbose)
	return( newfile)
}


`BAM.indexFASTA` <- function( file, verbose=TRUE) {

	if ( ! file.exists( file)) {
		if (verbose) cat( "\nFile not found:  ", file, "\nFailed to index FASTA file.")
		return(NULL)
	}
	newfile <- paste( file, "fai", sep=".")
	if ( file.exists( newfile)) return(newfile)

	samtools <- Sys.which( "samtools")
	if ( samtools == "") stop( "Executable not found on search path:  'samtools'")
	cmdline <- paste( samtools, " faidx ", file)
	if (verbose) cat( "\nIndexing FASTA file:  ", basename(file))
	system( cmdline)
	if (verbose) cat( "\nDone.")
	return( newfile)
}


`BAM.mpileup` <- function( files, seqID, fastaFile, start=NULL, stop=NULL, min.depth=1, max.depth=10000,
			min.gap.fraction=0.25, mpileupArgs="", summarize.calls=FALSE, verbose=TRUE) {

	if ( ! file.exists( fastaFile)) {
		cat( "\nGenomic Fasta File not found:  ", fastaFile, "\nUnable to call SAMTOOLS MPILEUP..\n")
		return(NULL)
	}

	N <- length( files)
	allFound <- TRUE
	# SAMTOOLS expects these to be '.sorted.bam' files that already have '.sorted.bam.bai' files...  Check!
	for ( i in 1:N) {
		thisFile <- BAM.verifySorted( files[i], index=TRUE)
		if ( is.null(thisFile)) {
			allFound <- FALSE
			next
		} else {
			files[i] <- thisFile
			ff <- paste( thisFile, "bai", sep=".")
			if ( ! file.exists( ff)) {
				cat( "\nBAM index file not found: ", ff)
				allFound <- FALSE
			}
		}
	}
	if ( ! allFound) {
		cat( "\nCall to 'SAMTOOLS MPILEUP' failed.\n")
		return(NULL)
	}

	fileArg <- files
	if ( N > 1) fileArg <- paste( files, collapse="  ")
	tmpFile <- tempfile()

	region <- seqID
	if ( ! is.null( start)) region <- paste( region, as.integer(start), sep=":")
	if ( ! is.null( stop)) region <- paste( region, as.integer(stop), sep="-")

	samtools <- Sys.which( "samtools")
	if ( samtools == "") stop( "Executable not found on search path:  'samtools'")
	cmdline <- paste( samtools, " mpileup -A -B -r ", region, " -f ", fastaFile, " -m ", min.depth, 
			" -d ", max.depth, " -F ", min.gap.fraction, " -L ", max.depth, 
			" ", mpileupArgs, " ", fileArg, " > ", tmpFile)
	if ( !verbose) cmdline <- paste( cmdline, "  2> /dev/null")

	if (verbose) cat( "\nGenerating Pileups for ", region, " of:  ", basename(files), "\n")
	system( cmdline)
	# force the reference bases to be characters
	ans <- try( read.delim( tmpFile, header=FALSE, comment.char="", quote="", as.is=T, 
				colClasses=c("V3"="character")), silent=TRUE)
	if ( class( ans) == "try-error") return( data.frame())
	file.remove( tmpFile)

	if ( N == 1) {
		colnames( ans)[1:6] <- c( "SEQ_ID", "POSITION", "REF_BASE", "DEPTH", "CALL_BASE", "CALL_SCORE")
	} else {
		# extract sampleIDs from the file names...
		sampleIDs <- sub( "(^.+)(\\.{1}?)(.+)", "\\1", basename(files))
		colnames( ans)[1:3] <- c( "SEQ_ID", "POSITION", "REF_BASE")
		colnames( ans)[4:ncol(ans)] <- paste( rep( c( "DEPTH", "CALL_BASE", "CALL_SCORE"), times=N),
							rep( sampleIDs, each=3), sep="_")
	}

	# we may want the base call summarized
	if ( summarize.calls) {
		if (verbose) cat( "\nSummarizing Pileups into final Base Calls..")
		calledBaseColumns <- grep( "CALL_BASE", colnames(ans))
		for ( baseColumn in calledBaseColumns) {

			rawBases <- ans[[ baseColumn]]
			callAns <- MPU.callBases( rawBases, ans$REF_BASE)
			ans[[ baseColumn]] <- callAns$call
			ans[[ baseColumn + 1]] <- MPU.callTableToString( callAns$depth.table)
			colnames(ans)[ baseColumn + 1] <- sub( "CALL_SCORE", "BASE_TABLE", colnames(ans)[baseColumn+1])
		}
	}

	if (verbose) cat( "\nDone.\nN_Lines: ", nrow(ans))
	return( ans)
}


`BAM.variantCalls` <- function( files, seqID, fastaFile, start=NULL, stop=NULL, 
				prob.variant=0.5, min.depth=1, max.depth=10000, min.gap.fraction=0.25,
				mpileupArgs="", vcfArgs="", ploidy=1, geneMap=getCurrentGeneMap(), 
				snpCallMode=c("all","multiallelic","consensus"), verbose=TRUE) {

	N <- length( files)
	fileArg <- files
	if ( N > 1) fileArg <- paste( files, collapse="  ")

	tmpFile <- tempfile()

	region <- seqID
	if ( ! is.null( start)) region <- paste( region, as.integer(start), sep=":")
	if ( ! is.null( stop)) region <- paste( region, as.integer(stop), sep="-")

	samtools <- Sys.which( "samtools")
	if ( samtools == "") stop( "Executable not found on search path:  'samtools'")
	bcftools <- Sys.which( "bcftools")
	if ( samtools == "") stop( "Executable not found on search path:  'bcftools'")

	# as of SAMTOOLS ~1.6, the ploidy is not numeric again...
	ploidyArg <- if ( ploidy == 1) " --ploidy 1 " else ""

	snpCallMode <- match.arg( snpCallMode)
	callModes <- c( " -m ", " -c ")
	if (snpCallMode == "multiallelic") callModes <- " -m "
	if (snpCallMode == "consensus") callModes <- " -c "

	out <- data.frame()
	for (thisCallMode in callModes) {

		# as of SAMTOOLS ~1.6 and up, the MPILEUP call for doing variant calling has changed!
		#cmdline <- paste( samtools, " mpileup -A -B -v -u -t DP -r ", region, " -f ", fastaFile, 
		#		" -d ", max.depth, " -m ", min.depth, " -F", min.gap.fraction,
		#		" -L", max.depth, mpileupArgs, "  ", fileArg, 
		cmdline <- paste( bcftools, " mpileup -A -B -r ", region, " -f ", fastaFile, 
				" -d ", max.depth, " -m ", min.depth, " -F", min.gap.fraction,
				" -L", max.depth, mpileupArgs, "  ", fileArg, 
				" | ", bcftools, " call -v ", thisCallMode, " -p ", prob.variant, ploidyArg,
				" -P 0 ", vcfArgs, " -O v -o ", tmpFile)
		if (verbose) {
			cat( "\nGenerating Variant Calls for ", region, " of:  ", basename(files), "\n")
			cat( "Command Line:   ", cmdline, "\n")
		}
		system( cmdline)
		ans <- try( read.delim( tmpFile, header=FALSE, comment.char="#", as.is=T), silent=TRUE)
		file.remove( tmpFile)

		if ( class( ans) == "try-error") {
			next
		} else {
	
			colnames( ans)[1:9] <- c( "SEQ_ID", "POSITION", "GENE_ID", "REF_BASE", "ALT_BASE", "QUAL", "FILTER", "INFO", "FORMAT")
	
			# tiny but non-zero chance that either of the 2 base call columns have nothing but 'T'
			# and R table reader may treat them as locical 'TRUE'
			if ( ncol(ans) < 1000) {
				ans$REF_BASE <- as.character( ans$REF_BASE)
				ans$REF_BASE[ ans$REF_BASE == "TRUE"] <- "T"
				ans$ALT_BASE <- as.character( ans$ALT_BASE)
				ans$ALT_BASE[ ans$ALT_BASE == "TRUE"] <- "T"
			}

			# there could be 'N's in the reference, that are of no use to us...
			isN <- which( ans$REF_BASE == "N")
			if ( length(isN) > 0) ans <- ans[ -isN, ]
	
			# extract sampleIDs from the file names...
			sampleIDs <- sub( "(^.+)(\\.{1}?)(.+)", "\\1", basename(files))
			colnames( ans)[10:(10+length(sampleIDs)-1)] <- sampleIDs
	
			# let's push the GeneID in...
			geneMap <- dropAntiSenseGenes( geneMap)
			gmap <- subset( geneMap, SEQ_ID == seqID)
			ptrs <- findInterval( ans$POSITION, gmap$POSITION)
			# tiny chance that no gene at very front edge of this map
			ptrs[ ptrs < 1] <- 1
			ans$GENE_ID <- gmap$GENE_ID[ ptrs]
		}

		# either combine two sets or not...
		if ( ! nrow(out)) {
			out <- ans
		} else {
			cat( "\nMerging 2 SNP call methods..")
			keyOut <- paste( out$SEQ_ID, out$POSITION, sep="::")
			keyAns <- paste( ans$SEQ_ID, ans$POSITION, sep="::")
			both <- sort( union( keyOut, keyAns))
			whOut <- match( both, keyOut, nomatch=0)
			whAns <- match( both, keyAns, nomatch=0)
			# don't take the ones that are lower scoring in the other
			inBoth <- which( whOut > 0 & whAns > 0)
			tieScoreOut <- out$QUAL[ whOut[ inBoth]]
			tieScoreAns <- ans$QUAL[ whAns[ inBoth]]
			whOut[ inBoth[ tieScoreAns > tieScoreOut]] <- 0
			whAns[ inBoth[ tieScoreOut >= tieScoreAns]] <- 0
			# OK, grab those 2 chunks
			partOut <- out[ whOut, ]
			partAns <- ans[ whAns, ]
			out <- rbind( partOut, partAns)
			ord <- order( out$SEQ_ID, out$POSITION)
			out <- out[ ord, ]
			rownames(out) <- 1:nrow(out)
		}
	}
	return( out)
}


cleanupBAMfiles <- function( path="results") {

	# delete any large unneeded BAM files left behind by the DuffyNGS pipeline
	filePatterns <- c( "genomic.bam$", "ribo.bam$", "ribo.converted.bam$", "splice.bam$", 
			"splice.converted.bam$")
	filePaths <- c( "align", "riboClear", "riboClear", "splicing", "splicing")
	sortedMustExist <- c( TRUE, FALSE, FALSE, FALSE, TRUE)

	totalBytes <- 0
	nfiles <- 0

	for (i in 1:length( filePatterns)) {
		patt <- filePatterns[i]
		thisPath <- filePaths[i]
		checkIfSorted <- sortedMustExist[i]
		files <- dir( file.path( path,thisPath), pattern=patt, recursive=T, full.name=T)
		if ( length( files) < 1) next
		for ( f in files) {
			if (checkIfSorted) {
				sortedBAM <- sub( "bam$", "sorted.bam", f)
				if ( ! file.exists(sortedBAM)) next
			}
			bytes <- file.info( f)$size
			nfiles <- nfiles + 1
			cat( "\r", nfiles, "  Size: ", bytes, "  File: ", f)
			totalBytes <- totalBytes + bytes
			file.delete( f)
		}
		cat( "\n")
	}
	cat( "\nDeleted_Files: ", nfiles, "\tDeleted Bytes: ", prettyNum( totalBytes, big.mark=","), "\n")
}


`BAM.fieldTable` <- function( file, field, chunkSize=100000, maxReads=NULL, verbose=TRUE) {

	if ( ! file.exists( file)) {
		if (verbose) cat( "\nFile not found:  ", file, "\nFailed to find existing BAM file.")
		return(NULL)
	}
	samtools <- Sys.which( "samtools")
	if ( samtools == "") stop( "Executable not found on search path:  'samtools'")

	# the field term can be numeric, or a tag name
	isNUMBER <- (length( grep( "^[0-9]+$", as.character(field))) == 1)
	if (isNUMBER) {
		field <- as.integer(field)
	} else {
		field <- paste( "^", field, sep="")
	}

	# make a pipe that readsfrom the BAM file
	cmdString <- paste( samtools, "view", file, sep=" ")
	conIn <- pipe( cmdString, open="r+t")

	# read the file in chunks, and accumulate thet field
	valueTable <- vector()
	nRead <- 0
	repeat {
		txt <- readLines( conIn, n=chunkSize)
		if ( ! length(txt)) break
		nRead <- nRead + length(txt)
		terms <- strsplit( txt, split="\t")

		if (isNUMBER) {
			values <- sapply( terms, `[`, field)
		} else {
			# find that field by patterm match
			values <- sapply( terms, function(x) {
					wh <- grep( field, x)
					if (length(wh))  x[wh[1]] else NA
				})
			# remove the field from the value we found
			values <- sub( field, "", values)
		}
		values <- values[ !is.na(values)]
		smallTable <- table( values)
		valueTable <- mergeTables( valueTable, smallTable)
		if ( verbose) cat( "\rN_Lines: ", nRead, "  N_Terms: ", length( valueTable))
		# if we read less than expected, the file is done
		if ( length(txt) < chunkSize) break
		# or quit early by read count
		if ( !is.null(maxReads) && nRead >= maxReads) break
	}

	close( conIn)

	nams <- names(valueTable)
	cnts <- as.numeric(valueTable)
	pcts <- round( cnts * 100 / sum(cnts), digits=4)
	out <- data.frame( "Value"=nams, "Count"=cnts, "Percentage"=pcts, stringsAsFactors=F)
	ord <- order( out$Count, decreasing=T)
	out <- out[ ord, ]
	if ( nrow(out)) rownames(out) <- 1:nrow(out)
	out
}


`SAM2BAM` <- function( file, verbose=TRUE, delete.SAM=TRUE) {

	if ( ! file.exists( file)) {
		if (verbose) cat( "\nFile not found:  ", file, "\nFailed to find existing SAM file.")
		return(NULL)
	}
	outfile <- sub( "sam$", "bam", file)
	if ( outfile == file) {
		if (verbose) cat( "\nFile not .sam:  ", file, "\nFailed to find SAM suffix.")
		return(NULL)
	}
	samtools <- Sys.which( "samtools")
	if ( samtools == "") stop( "Executable not found on search path:  'samtools'")
	if (verbose) cat( "\nCompressing SAM file:  ", basename(file))
	cmdline <- paste( samtools, " view -b -o ", outfile, " ", file)
	system( cmdline)
	if (verbose) cat( "\nDone.")

	if (delete.SAM) {
		# verify the new file exists first...
		if ( file.exists( outfile)) file.delete( file)
	}

	return( outfile)
}
