# pipe.Kmerize.R

# turn raw FASTQ reads into a 1-D table of Kmer counts


`pipe.Kmerize` <- function( sampleID=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt", 
				kmer.size=33, doCutadapt=TRUE, cutadaptProgram=Sys.which("cutadapt"), forceMatePairs=NULL, 
				buffer.size=1000000, maxReads=NULL, min.count=2, verbose=TRUE) {

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nStarting  'Kmerize Pipeline' on Sample:     ", sampleID,
			"\nStart Date/Time:   \t", date(), 
			"\nKmer Size:         \t", kmer.size, "\n")
	}

	results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)

	# file(s) to process comes from annotation and options...
	rawFastq <- getRawFastqFileNames( sampleID, annotationFile, optionsFile, verbose=FALSE)
	inputFastqFiles <- rawFastq$files
	asMatePairs <- rawFastq$asMatePairs

	# if we are doing the CutAdapt pass, check for those files, and call it if we need.
	if (doCutadapt) {
		# see if the files we need are already there
		cat( "\n")
		if( grepl( "\\.fastq", inputFastqFiles[1])) expectedFiles <- sub( "\\.fastq", ".trimmed.fastq", inputFastqFiles)
		if( grepl( "\\.fq", inputFastqFiles[1])) expectedFiles <- sub( "\\.fq", ".trimmed.fq", inputFastqFiles)
		if ( all( file.exists( expectedFiles))) {
			cat( "\nUsing existing 'cutadapt' results..")
			filesToDo <- expectedFiles
		} else {
			# is the data paired end or not
			if (asMatePairs || (!is.null(forceMatePairs) && forceMatePairs)) {
				if (verbose) cat( "\nRunning 'CutAdapt' on mate pair files..")
				filesToDo <- cutadapt( file1=inputFastqFiles[1], file2=inputFastqFiles[2], 
                        				cutadaptProgram=cutadaptProgram, min.length=kmer.size)
			} else {
				filesToDo <- vector()
				for ( i in 1:length(inputFastqFiles)) {
					if (verbose) cat( "\nRunning 'CutAdapt' on single file: ", basename(inputFastqFiles[i]))
					filesToDo[i] <- cutadapt( file1=inputFastqFiles[i], file2=NULL, cutadaptProgram=cutadaptProgram, 
								min.length=kmer.size)
				}
			}
		}
	} else { 
		filesToDo <- inputFastqFiles
	}

	# set up to write out the results as they grow
	kmer.path <- file.path( results.path, "Kmers")
	if ( ! file.exists( kmer.path)) dir.create( kmer.path, recursive=T)
	outfile <- file.path( kmer.path, paste( sampleID, "Kmer", kmer.size, "Table.rda", sep="."))
	# removed any previous results
	file.delete( outfile)

	startT <- proc.time()

	nReadsIn <- nDistinctKmers <- nTotalKmers <- 0
	nFiles <- length( filesToDo)
	kmer.size <- as.integer( kmer.size)

	# let's try refactoring how we do this, to minimize memory usage
	# use one GLOBAL table, and just extend as needed
	bigKmerTable <<- vector()

	cat( "\nN_Files to Kmerize: ", nFiles, "\n")
	for ( i in 1:nFiles) {
		ans <- kmerizeOneFastqFile( filesToDo[i], kmer.size=kmer.size, buffer.size=buffer.size, 
					sampleID=sampleID, kmer.path=kmer.path, 
					maxReads=maxReads, min.count=min.count)
		nReadsNow <- ans$nReadsIn

		cat( "  save..")
		save( bigKmerTable, file=outfile)

		nReadsIn <- nReadsIn + nReadsNow
		nDistinctKmers <- length( bigKmerTable)
		nTotalKmers <- sum( bigKmerTable)
		cat( "\nFile: ", i, basename(filesToDo[i]), "  N_Distinct: ", nDistinctKmers, "  N_Total: ", nTotalKmers, "\n")
		if ( ! is.null(maxReads) && nReadsIn >= maxReads) {
			cat( "\nReached 'maxReads' count. Stopping early..")
			break
		}
	}

	# now do the cleanup...
	rm( bigKmerTable, ans, inherits=T)
	gc()

	cat( verboseOutputDivider)
	cat( "\n\nFinished 'Kmerize Pipeline' on Sample:     ", sampleID, "\n")
	cat( "\nTiming Stats: \n")
	myTime <- elapsedProcTime( startT, proc.time(), N=nTotalKmers, what="Kmer")
	print( myTime)

	out <- list( "nReadsIn"=nReadsIn, "nDistinctKmers"=nDistinctKmers, "nTotalKemrs"=nTotalKmers)
	return( out)
}


`pipe.KmerTable` <- function( sampleIDset, annotationFile="Annotation.txt", optionsFile="Options.txt", 
				kmer.size=33, min.count=5, min.samples=2) {

	results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	kmer.path <- file.path( results.path, "Kmers")

	fileSet <- file.path( kmer.path, paste( sampleIDset, "Kmer", kmer.size, "Table.rda", sep="."))
	found <- file.exists( fileSet)
	if ( any( !found)) {
		cat( "\nSome Kmer files not found: ", fileSet[!found])
		fileSet <- fileSet[ found]
		sampleIDset <- sampleIDset[ found]
	}
	NS <- length(fileSet)
	if ( NS < 1) return(NULL)

	# we need to build up one giant vector of Kmers
	# try to be memory efficient
	expectNkmers <- NS * 10000000
	allKmers <- vector( mode="character", length=expectNkmers)
	nKmers <- 0

	cat( "\nFind union of all Sample Kmers..\n")
	for ( f in fileSet) {
		cat( "  load: ", basename(f))
		load(f)
		myKmers <- names(bigKmerTable)
		nNow <- length( myKmers)
		cat( "  N:", nNow)
		if ( ! nKmers) {
			allKmers[1:nNow] <- myKmers
			nKmers <- nNow
		} else { 
			cat( "  lookup..")
			where <- match( myKmers, allKmers[1:nKmers], nomatch=0)
			# only add the new ones
			myKmers <- myKmers[ where == 0]
			nNow <- length( myKmers)
			if ( nNow) {
				cat( "  extend..")
				now <- (nKmers+1) : (nKmers+nNow)
				allKmers[now] <- myKmers
				nKmers <- nKmers + nNow
			}
			rm( where)
		}
		cat( "  N_Kmer:", nKmers)
	}
	NK <- nKmers
	# trim the memory allocation if we need
	if ( expectNkmers > NK) length(allKmers) <- NK
	rm( myKmers, bigKmerTable)
	gc()

	# now fill the table
	cat( "\nReloading to fill table..\n")
	kmerTbl <- matrix( 0, nrow=NK, ncol=NS)
	colnames(kmerTbl) <- sampleIDset
	rownames(kmerTbl) <- allKmers

	for ( i in 1:NS) {
		f <- fileSet[i]
		cat( "  load: ", basename(f))
		load(f)
		myKmers <- names(bigKmerTable)
		myCounts <- as.vector( bigKmerTable)
		cat( "  lookup..")
		wh <- match( myKmers, allKmers)
		v <- rep.int(0,NK)
		cat( "  store..")
		v[ wh] <- myCounts
		kmerTbl[ , i] <- v
		rm( bigKmerTable, myKmers, myCounts, wh)
		gc()
	}
	cat( "\nDone loading.\n")
	
	cat( "\nChecking for low coverage Kmers to drop: \n  At least", min.count, "Kmers in at least", min.samples, "samples..")
	nGood <- apply( kmerTbl, 1, function(x) sum( x >= min.count))
	drops <- which( nGood < min.samples)
	if ( length(drops)) {
		cat( "  Removing", length(drops), "Kmers rows..")
		kmerTbl <- kmerTbl[ -drops, ]
		cat( "  N_Kmer: ", nrow(kmerTbl))
	}
	
	return( kmerTbl)
}


`pipe.KmerCompare` <- function( kmerTbl, sampleIDset, groupSet, normalize=c("LKPTM"), n.remove=100) {

	# just just the samples asked for
	where <- match( sampleIDset, colnames(kmerTbl))
	if ( any( is.na(where))) stop( "Some SampleIDs not in Kmer table")
	useTbl <- kmerTbl[ , where]
	NC <- ncol(useTbl)
	NR <- nrow(useTbl)
	
	# get our group factors
	grpFac <- factor( groupSet)
	grpLvls <- levels( grpFac)
	nGrps <- nlevels( grpFac)
	if ( nGrps != 2) stop( "Expected exactly 2 Sample Groups")
	cat( "\nBreakdown by Group:\n")
	print( table( groupSet))
	
	# drop the super high count Kmers
	if (n.remove > 0) {
		useTbl <- remove.HighCountKmers( useTbl, n.remove=n.remove)
		NR <- nrow(useTbl)
	}

	# convert everything to normalized unit
	normalize <- match.arg( normalize)
	if ( normalize == "LKPTM") {
		cat( "\nConverting to LKPTM (Log2 Kmers Per Ten Million)..")
		useTbl <- as.LKPTM( useTbl)
	}
	
	# we are now ready to do tests per Kmer
	cat( "\nRunning Linear Model on ", NR, " Kmers..\n")
	fold <- pval <- vector( mode="numeric", length=NR)
	avg <- matrix( 0, 2, NR)
	rownames(avg) <- paste( "Avg", grpLvls, sep="_")

	for ( i in 1:NR) {
		v <- useTbl[ i, ]
		ans <- t.test( v ~ grpFac)
		pval[i] <- ans$p.value
		avg[ ,i] <- ans$estimate
		fold[i] <- log2( (avg[2,i]+1) / (avg[1,i]+1))
		if ( i %% 10000 == 0) cat( "\r", i, rownames(useTbl)[i], fold[i], pval[i])
	}
	
	out <- data.frame( "Kmer"=rownames(useTbl), t(avg), "Log2.Fold"=fold, "P.Value"=pval,
			stringsAsFactors=F)
	ord <- diffExpressRankOrder( out$Log2.Fold, out$P.Value)
	out <- out[ ord, ]
	rownames(out) <- 1:NR
	
	return(out)
}


`alignKmersToGenome` <- function( kmers, optionsFile="Options.txt", quiet=TRUE, verbose=!quiet) {

	if ( is.data.frame(kmers) && ("Kmer" %in% colnames(kmers))) {
		kmers <- kmers$Kmer
	}

	# map some Kmers onto the reference genome
	N <- length( kmers)
	
	# make a temp FASTA file from the Kmer sequences
	kmerIDs <- paste( "Kmer", 1:N, sep="_")
	kmerFastaFile <- "Temp.Kmers.BowtieInput.fasta"
	writeFasta( as.Fasta( kmerIDs, kmers), kmerFastaFile)
	kmerBamFile <- "Temp.Kmers.BowtieOutput.bam"
	
	# build a call to Bowtie2
	ans <- kmerCallBowtie( kmerFastaFile, kmerBamFile, optionsFile=optionsFile, quiet=quiet, verbose=verbose)
	
	# read up the BAM file contents to extract what we need
	ans2 <- kmerReadBam( kmerBamFile, verbose=verbose)
	
	# force the results to match the input order
	wh <- match( kmers, ans2$Kmer)
	out <- ans2[ wh, ]
	return( out)
}


`mapKmersToProteins` <- function( kmerAlignments, verbose=TRUE) {

	# takes the output from 'alignKmersToGenome', and looks up where those sites 
	# land on proteins
	if ( ! all( c( "SEQ_ID", "POSITION", "GENE_ID") %in% colnames( kmerAlignments))) {
		cat( "\nWarning:  expected a data.frame with SEQ_ID, POSITION, GENE_ID columns.")
		cat( "\nPerhaps run 'alignKmersToGenome()' first...")
		return( NULL)
	}

	N <- nrow(kmerAlignments)
	
	# for Kmers that hit genes, let's map to AA location and guess the protein fragment
	aaPos <- rep.int( NA, N)
	aaFrag <- rep.int( "", N)
	geneHits <- setdiff( 1:N, grep( "(ng)", kmerAlignments$GENE_ID, fixed=T)) 
	if (verbose) cat( "\nMapping", length(geneHits), "Kmer Gene hits to AA location and fragment seqs..\n")
	for (j in geneHits) {
		if ( any( is.na( c( kmerAlignments$SEQ_ID[j], kmerAlignments$POSITION[j])))) next
		smlAns <- convertGenomicDNApositionToAAposition( kmerAlignments$SEQ_ID[j], kmerAlignments$POSITION[j])
		aaPos[j] <- smlAns$AA_POSITION
		if ( is.na( aaPos[j])) next
		if ( is.na( kmerAlignments$STRAND[j])) next
		readFrame <- if ( kmerAlignments$STRAND[j] == "+") 1:3 else if (kmerAlignments$STRAND[j] == "-") 4:6 else 1:6
		aaFrag[j] <- DNAtoBestPeptide( kmerAlignments$Kmer[j], clip=F, readingFrame=readFrame)
		if (verbose && j %% 1000 == 0) cat( "\r", j, kmerAlignments$GENE_ID[j], aaPos[j], aaFrag[j])
	}
	
	out <- kmerAlignments
	out$AA_POSITION <- aaPos
	out$PROTEIN_FRAGMENT <- aaFrag
	return( out)
}


`kmerizeOneFastqFile` <- function( filein, kmer.size=33, buffer.size=1000000, 
				sampleID="SampleID", kmer.path=".", maxReads=NULL, min.count=2) {

	# get ready to read the fastq file in chuncks...
	filein <- allowCompressedFileName( filein)
	if ( ! file.exists( filein)) stop( paste("Can't find input file: ", filein))
	conIn <- openCompressedFile( filein, open="r")
	on.exit( close( conIn))

	# 4 lines per read
	chunkSize <- buffer.size
	readChunkSize <- chunkSize * 4
	nread <- 0
	hasMore <- TRUE
	firstPass <- TRUE
	#bigTable <- vector()

	repeat {
		if ( ! hasMore) break
		if ( ! is.null(maxReads) && nread >= maxReads) break
		cat( "\nBuffer Read..")
		chunk <- readLines( conIn, n=readChunkSize)
		if ( length( chunk) < 1) break
		if ( length( chunk) < readChunkSize) hasMore <- FALSE

		# in a FASTQ file, it takes 4 lines for each read;  1 is ID, 2nd is the seq, 4th is scores
		# all we want is the raw reads
		seqTxt <- chunk[ seq.int( 2, length(chunk), by=4)]
		nread <- nread + length( seqTxt)
		cat( "  N_reads: ", prettyNum( as.integer(nread), big.mark=","))
		rm( chunk)

		timeStat <- system.time( smlTable <- kmerizeOneChunk( seqTxt, kmer.size=kmer.size))
		cat( "  time(usr,sys)=",timeStat[3],"(",timeStat[1],",",timeStat[2],"}", sep="")
		rm( seqTxt)

		# do the merge in place in GLOBAL storage
		if ( ! length(bigKmerTable)) {
			bigKmerTable <<- smlTable
		} else {
			cat( "  store..")
			where <- match( names(smlTable), names(bigKmerTable), nomatch=0)
			yesNo <- (where > 0)
			hitsTo <- where[ yesNo]
			cntsWas <- as.vector( bigKmerTable[ hitsTo])
			cntsNew <- smlTable[ yesNo]
			bigKmerTable[ hitsTo] <<- cntsWas + cntsNew
			cat( "  extend..")
			newFrom <- which( ! yesNo)
			lenNew <- length( newFrom)
			lenWas <- length( bigKmerTable)
			length(bigKmerTable) <<- lenWas + lenNew
			newTo <- (lenWas+1) : (lenWas+lenNew)
			bigKmerTable[ newTo] <<- as.vector( smlTable[ newFrom])
			names(bigKmerTable)[ newTo] <<- names( smlTable)[ newFrom]
		}
		cat( "  N_Kmer:", length( bigKmerTable))
		gc()
	}
	rm( smlTable)

	# any K-mers with too few observations are not to be kept
	cat( "\nDrop if count < ", min.count)
	tooFew <- which( bigKmerTable < min.count)
	if ( length( tooFew)) {
		cat( "  Found:", length(tooFew))
		tmpTable <- bigKmerTable[ -tooFew]
		rm( bigKmerTable, inherits=T)
		bigKmerTable <<- tmpTable
		rm( tmpTable, tooFew)
		cat( "  N_Kmer:", length( bigKmerTable))
	}
	gc()
	
	# also do any RevComp resolving now too
	cat( "\nSearch for RevComp pairs to join..")
	kmerRevComp.GlobalTable( sampleID=sampleID, kmer.path=kmer.path, kmer.size=kmer.size)

	#return( list( "nReadsIn"=nread, "Kmer.Table"=bigTable))
	return( list( "nReadsIn"=nread))
}


`kmerizeOneChunk` <- function( seqs, kmer.size=33) {

	# make Kmers of every raw read in this chunk
	# set up to get all N-mers for all the sequences in this chunk
	LAPPLY <- base::lapply
	SUBSTR <- base::substr

	sizeM1 <- kmer.size - 1

	# high chance of having duplicates, so only do each one once
	seqTbl <- table.nosort( seqs)
	seqCounts <- as.vector( seqTbl)
	seqs <- names( seqTbl)
	nSeq <- length( seqs)
	cat( "  N_Unique:", nSeq)
	seqLens <- base::nchar( seqs)
	rm( seqTbl)
	gc()

	# pre-guess a size of output, based on size of input
	expectNout <- length(seqs) * (max(seqLens) - sizeM1)
	kmersOut <- vector( mode="character", length=expectNout)
	nOut <- 0

	# we can save a bit of time by doing all reads of the same length at one time
	lenFac <- factor( seqLens)
	cat( "  ReadLen=")
	tapply( 1:nSeq, lenFac, function(x) {
			
		# given all the pointers to seqs of the same length
		lenNow <- seqLens[ x[1]]
		cat( ".", lenNow, sep="")
		veryLastStart <- lenNow - sizeM1
		if ( veryLastStart < 1) return(NULL)
		fromSet <- 1:veryLastStart
		toSet <- fromSet + sizeM1
		nSubstrs <- length( fromSet)
		#for ( j in x) {
		LAPPLY( x, function(j) {
			thisCount <- seqCounts[j]
			hasN <- grepl( "N", seqs[j], fixed=T)
			kmer <- SUBSTR( rep.int( seqs[j], nSubstrs), fromSet, toSet)
			# we only want to keep valid DNA chunks, so discard any 'N's
			if (hasN) {
				isNNN <- grep( "N", kmer, fixed=T)
				if ( length(isNNN)) kmer <- kmer[ -isNNN]
				if ( ! length(kmer)) return(NULL)  #next
			}
			# account for the multiplicity of this read
			if ( thisCount > 1) kner <- rep.int( kmer, thisCount)
			n <- length(kmer)
			now <- (nOut+1) : (nOut+n)
			if ( nOut+n > length(kmersOut)) {
				cat( "  MemRealloc..")
				length(kmersOut) <<- (expectNout <<- expectNout + 1000000)
			}
			kmersOut[now] <<- kmer
			nOut <<- nOut + n
			return(NULL)
		})
		return(NULL)
	})
	# lastly trim size if we over estimated
	if ( expectNout > nOut) length(kmersOut) <- nOut
	gc()

	cat( "  Tablulate..")
	return( table.nosort( kmersOut))
}


`kmerRevComp.file` <- function( kmerFile, sampleID="SampleID", kmer.path=".", kmer.size=33) {

	# the Kmers may be from both strands, and we only want to keep one form of each
	# so merge those where both are present
	
	# allow being given a saved file
	load( kmerFile)
	kmerTable <- bigKmerTable
	isFile <- TRUE
	
	kmers <- names( kmerTable)
	cnts <- as.numeric( kmerTable)
	N <- length(kmerTable)
	cat( "\nN_Kmers in: ", N)

	# see what the RevComp of every Kmer is
	rcKmer <- findKmerRevComp( kmers, sampleID=sampleID, kmer.path=kmer.path, kmer.size=kmer.size)
	
	# then see if those Rev Comps are already in the table
	cat( "\nLocate RevComp pairs..")
	where <- match( rcKmer, kmers, nomatch=0)

	# we will make a new 1-D table that has just the first form of each Kmer
	# we can combine those that have both forms
	cat( "  join..")
	out <- rep.int( 0, N)
	hasRC <- which( where > 0)
	if ( length(hasRC)) {
		myLoc <- (1:N)[hasRC]
		rcLoc <- where[hasRC]
		firstLoc <- pmin( myLoc, rcLoc)
		out[firstLoc] <- cnts[myLoc] + cnts[rcLoc]
		names(out)[firstLoc] <- kmers[firstLoc]
	}

	# then those without their RevComp just stay as is
	noRC <- which( where == 0)
	if ( length(noRC)) {
		out[noRC] <- cnts[noRC]
		names(out)[noRC] <- kmers[noRC]
	}

	# lastly thown away the empty slots
	drops <- which( out == 0)
	if ( length( drops)) out <- out[ -drops]
	Nout <- length(out)
	cat( "\nN_Kmers out: ", Nout)
	
	cat( "  Pct Reduction: ", round( (N-Nout) * 100 / N), "%")

	# done
	bigKmerTable <- out
	cat( "  re-save Kmers file..")
	save( bigKmerTable, file=kmerFile)
	return( kmerFile)
}


`kmerRevComp.GlobalTable` <- function( sampleID="SampleID", kmer.path=".", kmer.size=33) {

	# the Kmers may be from both strands, and we only want to keep one form of each
	# so merge those where both are present
	
	# faster memory efficient to operate on one global vector
	kmers <- names( bigKmerTable)
	cnts <- as.vector( bigKmerTable)
	N <- length(bigKmerTable)
	cat( "\nN_Kmers in: ", N)

	rm( bigKmerTable, inherits=T)
	gc()

	# see what the RevComp of every Kmer is
	rcKmer <- findKmerRevComp( kmers, sampleID=sampleID, kmer.path=kmer.path, kmer.size=kmer.size)
	gc()
	
	# then see if those Rev Comps are already in the table
	cat( "\nLocate RevComp pairs..")
	where <- match( rcKmer, kmers, nomatch=0)

	# we will make a new 1-D table that has just the first form of each Kmer
	# we can combine those that have both forms
	cat( "  join..")
	out <- rep.int( 0, N)
	hasRC <- which( where > 0)
	if ( length(hasRC)) {
		myLoc <- (1:N)[hasRC]
		rcLoc <- where[hasRC]
		firstLoc <- pmin( myLoc, rcLoc)
		out[firstLoc] <- cnts[myLoc] + cnts[rcLoc]
		names(out)[firstLoc] <- kmers[firstLoc]
		rm( myLoc, rcLoc)
	}

	# then those without their RevComp just stay as is
	noRC <- which( where == 0)
	if ( length(noRC)) {
		out[noRC] <- cnts[noRC]
		names(out)[noRC] <- kmers[noRC]
	}

	# lastly thown away the empty slots
	drops <- which( out == 0)
	if ( length( drops)) out <- out[ -drops]
	Nout <- length(out)
	cat( "\nN_Kmers out: ", Nout)
	
	cat( "  Pct Reduction: ", round( (N-Nout) * 100 / N), "%")

	# done
	bigKmerTable <<- out

	rm( out, kmers, cnts, hasRC, noRC, drops)
	gc()
	return( length( bigKmerTable))
}


findKmerRevComp <- function( kmers, sampleID=sampleID, kmer.path=".", kmer.size=33) {

	# set up to allow 'faster' lookup of previously calculated Rev Comps
	N <- length( kmers)
	out <- rep.int( "", N)
	
	# allow 'fast' lookup via a file of previously done kmers
	# we are switching to separate files for every sample, to try to minimize file
	# access collisions.  You can read fron any sample's file, but only write to
	# your own.
	myKmerXrefFile <- file.path( kmer.path, paste( sampleID, kmer.size, "Kmer.RevComp.Xref.rda", sep="."))
	allKmerXrefFiles <- dir( kmer.path, pattern="Kmer.RevComp.Xref.rda$", full=T)
	
	if ( length(allKmerXrefFiles)) {
	    cat( "  searching", length(allKmerXrefFiles), "RevComp Xref files:")
	    for ( f in allKmerXrefFiles) {
		xref <- data.frame()
		status <- tryCatch( load( f), error=function(e) { xref <<- data.frame()})
		if ( ! nrow( xref)) next
	
		# we got some from one Xref file, see if we get any hits
		# but only check on the ones we still need on this pass around
		toTest <- which( out == "")
		cat( "  lookup..")
		where1 <- match( kmers[toTest], xref$Kmer, nomatch=0)
		out[ toTest[ where1 > 0]] <- xref$RevComp[ where1]
		where2 <- match( kmers[toTest], xref$RevComp, nomatch=0)
		out[ toTest[ where2 > 0]] <- xref$Kmer[ where2]
		cat( "  found=", sum(where1 > 0 | where2 > 0))
		rm( xref, toTest, where1, where2)
	    }
	}
	
	# calculate the rest
	toCalc <- which( out == "")
	if ( length( toCalc)) {
		cat( "  calc RevComps.. N=", length(toCalc), " (", round( length(toCalc)*100/N, digits=1),"%)", sep="")
		rcKmer <- vapply( kmers[ toCalc], FUN=myReverseComplement, FUN.VALUE="ACGT", USE.NAMES=F)
		out[ toCalc] <- rcKmer
	
		# store these for the future
		smlXref <- data.frame( "Kmer"=kmers[toCalc], "RevComp"=rcKmer, stringsAsFactors=FALSE)
		
		# but reload this sample's Xref first
		xref <- data.frame()
		if ( file.exists( myKmerXrefFile)) load( myKmerXrefFile)
		cat( "  save updated RevComp Xref file..")
		xref <- rbind( xref, smlXref)
		save( xref, file=myKmerXrefFile)
		rm( xref, rcKmer, smlXref, toCalc)
	}
	return( out)
}


remove.HighCountKmers <- function( kmerTbl, n.remove=100) {

	# there are some Kmers that are pure noise, like poly-A and ATAT.
	# allow removal of the super-high count ones, to give a more mornalized distribution
	# between all samples
	whoHigh <- vector()
	n.remove <- min( n.remove, round(nrow(kmerTbl)*0.01))

	for ( i in 1:ncol(kmerTbl)) {
		ord <- order( kmerTbl[,i], decreasing=T)
		whoHigh <- c( whoHigh, ord[1:n.remove])
	}
	whoHigh <- sort( unique( whoHigh))
	nDrop <- length(whoHigh)

	mDrop <- kmerTbl[ whoHigh, ]
	kmerCountDrop <- sum( mDrop)
	kmerCountIn <- sum( kmerTbl)

	cat( "\nHigh Count Removal flagged", nDrop, "Kmers")
	cat( "\n  Being ", round( kmerCountDrop * 100 / kmerCountIn, digits=2), "% of all Kmer counts")
	cat( "\nTop Culprits:\n")
	n.show <- min( nDrop, 20)
	print( head( mDrop, n.show))

	out <- kmerTbl[ -whoHigh, ]
	out
}


as.LKPTM <- function( kmerTbl) {

	kSums <- apply( kmerTbl, 2, sum)
	NC <- ncol(kmerTbl)
	lkpbTbl <- kmerTbl
	for ( i in 1:NC) {
		myCnts <- kmerTbl[ ,i]
		myKPB <- myCnts * 1e7 / kSums[i]
		myLKPTM <- round( log2( myKPB + 1), digits=4)
		lkpbTbl[ , i] <- myLKPTM
	}
	return( lkpbTbl)
}


kmerCallBowtie <- function( kmerFastaFile, kmerBamFile, optionsFile="Options.txt", quiet=T, verbose=F) {

	
	optT <- readOptionsTable( optionsFile)
	# force the look at options file to initialize bowtie.
	bowtie2Par.defaults( optionsFile, verbose=FALSE)
	## build the Bowtie alignment program's command line
	progName <- getOptionValue( optT, "bowtie2Program", verbose=FALSE)
	cmdLine <- progName
	
	# next is "input options" using FASTA instead of default FASTQ
	cmdLine <- paste( cmdLine, " -f")
	# we can know if a Kmer is unique or not by asking for up to 2 alignments
	cmdLine <- paste( cmdLine, " -k 2")
	# next is explicit alignment policy
	cmdLine <- paste( cmdLine, " --very-sensitive")
	if ( quiet) cmdLine <- paste( cmdLine, " --quiet")
	# when not catching the unaligned, allow explicit throw-away of unaligned
	cmdLine <- paste( cmdLine, " --no-unal")
	# next is threads
	nCores <- as.integer( getOptionValue( optT, "nCores", notfound="4", verbose=FALSE))
	cmdLine <- paste( cmdLine, " --threads", nCores)
	
	# let's use a slightly more lax scoring threshold, to get a few more aligned sites for variant Kmers
	cmdLine <- paste( cmdLine, " --score-min L,-0.7,-0.7 ")

	# the index
	index.path <- getOptionValue( optT, "bowtie2Index.path", notfound=".", verbose=FALSE)
	alignIndex <- getOptionValue( optT, "GenomicIndex", verbose=F)
	myIndexFile <- file.path( index.path, alignIndex)
	# test to see that a file is really there...
	if ( ! file.exists(  paste( myIndexFile, ".1.bt2", sep=""))) {
	    # new long indexex are possible
	    if ( ! file.exists(  paste( myIndexFile, ".1.bt2l", sep=""))) {
		stop( paste( "Bowtie2 index Path and/or File not found. \n", "Filename as given:  ", myIndexFile))
	    }
	}
	cmdLine <- paste( cmdLine, " -x", myIndexFile)
	# the input file
	cmdLine <- paste( cmdLine, " -U", kmerFastaFile)
	
	# lastly, the output destination
	outputTerm <- paste ( " | samtools view -bS -o", kmerBamFile, " - ")
	cmdLine <- paste( cmdLine, outputTerm, sep="  ")
	
	#OK, now ready to cal Bowtie
	return( callBowtie2( cmdLine, verbose=verbose))
}


kmerReadBam <- function( kmerBamFile, chunkSize=100000, verbose=T) {


	# set up to read the BAM file in chunks
	con <- bamReader( kmerBamFile)
	refData <- getRefData( con)
	
	kmerOut <- sidOut <- posOut <- gidOut <- strandOut <- uniqOut <- vector()
	nOut <- 0
	
	# read in the alignment file one buffer at a time
	hasMore <- TRUE
	nReads <- 0
	repeat {
		if ( ! hasMore) break
		if (verbose) cat( "\rReadBAM..")
		chunk <- getNextChunk( con, n=chunkSize, alignedOnly=TRUE)
		nNow <- size(chunk)
		if ( nNow < 1) break
		if ( nNow < chunkSize) hasMore <- FALSE
		nReads <- nReads + nNow
		if (verbose) cat( "  N_Kmers: ", prettyNum( as.integer( nReads), big.mark=","))

		# convert the refID to a SeqID to figure the gene locations
		if (verbose) cat( "  geneIDs..")
		seqIDs <- refID2seqID( refID( chunk), refData=refData)
		positions <- position( chunk)
		kmerSeq <- readSeq( chunk)
		lens <- nchar( kmerSeq)
		middles <- positions + round( lens/2)
		ans <- fastSP2GP( seqIDs, middles)
		geneIDs <- ans$GENE_ID
		
		# also get the strand hit
		strand <- ! reverseStrand( chunk)
		
		# and a yes/no call about who had 2+ alignments
		isDup <- duplicated( kmerSeq)
		if (any( isDup)) {
			# we only keep one copy of each
			dupSeqs <- kmerSeq[isDup]
			keep <- which( ! isDup)
			kmerSeq <- kmerSeq[keep]
			seqIDs <- seqIDs[keep]
			middles <- middles[keep]
			geneIDs <- geneIDs[keep]
			strand <- strand[keep]
			nNow <- length(keep)
			uniq <- rep.int(TRUE,nNow)
			notUniq <- which( kmerSeq %in% dupSeqs)
			uniq[notUniq] <- FALSE
		}
		
		# save what we need to send back
		now <- (nOut+1) : (nOut+nNow)
		kmerOut[now] <- kmerSeq
		sidOut[now] <- seqIDs
		posOut[now] <- middles
		gidOut[now] <- geneIDs
		strandOut[now] <- ifelse( strand, "+", "-")
		uniqOut[now] <- uniq
		nOut <- nOut + nNow

	} # end of each buffer...
	if (verbose) cat("\n")
	bamClose( con)
	
	# package up what we send back
	out <- data.frame( "Kmer"=kmerOut, "UNIQUE"=uniqOut, "SEQ_ID"=sidOut, "POSITION"=posOut,  
					"GENE_ID"=gidOut, "STRAND"=strandOut, stringsAsFactors=F)
	if ( nOut) {
		# because of the way we worked in buffers, there is a tiny chance of duplicates.  Check explicitly once more.
		dups <- duplicated( out$Kmer)
		if (any( dups)) {
			isDup <- which( dups)
			dupSeqs <- out$Kmer[ isDup]
			out <- out[ -isDup, ]
			where <- match( unique.default(dupSeqs), out$Kmer)
			out$UNIQUE[ where] <- FALSE
		}
		rownames(out) <- 1:nrow(out)
	}
	return( out)
}
