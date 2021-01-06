# pipe.Kmerize.R

# turn raw FASTQ reads into a 1-D table of Kmer counts

# using Biostrings 'DNAStrings' to improve speed/capacity.  THere is is strict limit on
# XVector hash table size, that we can never exceed
MAX_KMERS <- 250000000


# turn raw FASTQ into Kmers for one sample

`pipe.Kmerize` <- function( sampleID=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt", 
				kmer.size=33, doCutadapt=TRUE, cutadaptProgram=Sys.which("cutadapt"), forceMatePairs=NULL, 
				buffer.size=1000000, maxReads=NULL, min.count=2, verbose=TRUE) {

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nStarting  'Kmerize Pipeline' on Sample:     ", sampleID,
			"\nStart Date/Time:   \t", date(), 
			"\nKmer Size:         \t", kmer.size, "\n")
	}
	require(Biostrings)

	# file(s) to process comes from annotation and options...
	results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	rawFastq <- getRawFastqFileNames( sampleID, annotationFile, optionsFile, verbose=FALSE)
	inputFastqFiles <- rawFastq$files
	asMatePairs <- rawFastq$asMatePairs

	# if we are doing the Cutadapt pass, check for those files, and call it if we need.
	if (doCutadapt) {
		filesToDo <- checkOrCallCutadapt( inputFastqFiles, asMatePairs=asMatePairs, forceMatePairs=forceMatePairs,
						cutadaptProgram=cutadaptProgram, kmer.size=kmer.size, verbose=verbose)
	} else { 
		filesToDo <- inputFastqFiles
	}
	nFiles <- length( filesToDo)

	# set up to write out the results as they grow
	kmer.path <- file.path( results.path, "Kmers")
	if ( ! file.exists( kmer.path)) dir.create( kmer.path, recursive=T)
	outfile <- file.path( kmer.path, paste( sampleID, "Kmer", kmer.size, "Table.rda", sep="."))
	# removed any previous results
	file.delete( outfile)

	startT <- proc.time()
	nReadsIn <- nDistinctKmers <- nTotalKmers <- 0
	kmer.size <- as.integer( kmer.size)

	# let's try refactoring how we do this, to minimize memory usage, by switching to use Biostrings DNAString package
	require( Biostrings)

	# use one GLOBAL list of DNA strings and counts
	bigKmerStrings <<- bigKmerCounts <<- vector( mode="list")
	gc()

	cat( "\nN_Files to Kmerize: ", nFiles, "\n")
	for ( i in 1:nFiles) {
		# do one file, and save what we got
		ans <- kmerizeOneFastqFile( filesToDo[i], kmer.size=kmer.size, buffer.size=buffer.size, 
					sampleID=sampleID, kmer.path=kmer.path, maxReads=maxReads, min.count=min.count)
		nReadsNow <- ans$nReadsIn
		nReadsIn <- nReadsIn + nReadsNow
		saveKmers( outfile)

		# report those matrics
		nDistinctKmers <- length( bigKmerCounts[[1]])
		nTotalKmers <- sum( bigKmerCounts[[1]])
		cat( "\nSample: ", sampleID, "  File: ", i, basename(filesToDo[i]), "  N_Distinct: ", nDistinctKmers, 
				"  N_Total: ", nTotalKmers, "\n")
		if ( ! is.null(maxReads) && nReadsIn >= maxReads) {
			cat( "\nReached 'maxReads' count. Stopping early..")
			break
		}
	}

	# now do the cleanup...
	if ( exists("bigKmerStrings", envir=.GlobalEnv)) rm( bigKmerStrings, envir=.GlobalEnv)
	if ( exists("bigKmerCounts", envir=.GlobalEnv)) rm( bigKmerCounts, envir=.GlobalEnv)
	gc()

	cat( verboseOutputDivider)
	cat( "\n\nFinished 'Kmerize Pipeline' on Sample:     ", sampleID, "\n")
	cat( "\nTiming Stats: \n")
	myTime <- elapsedProcTime( startT, proc.time(), N=nTotalKmers, what="Kmer")
	print( myTime)

	# it is possible that the sample ends up with zero kmers...
	if ( nTotalKmers < 100) {
		cat( "\nSample produced too few Kmers..  Deleting Kmer file.")
		file.delete( outfile)
	}

	out <- list( "nReadsIn"=nReadsIn, "nDistinctKmers"=nDistinctKmers, "nTotalKemrs"=nTotalKmers)
	return( out)
}


# Join all Kmers from many samples into a matrix of counts

`pipe.KmerTable` <- function( sampleIDset, annotationFile="Annotation.txt", optionsFile="Options.txt", 
				kmer.size=33, min.count=5, min.samples=2) {

	results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	kmer.path <- file.path( results.path, "Kmers")
	require(Biostrings)

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
	# try to be memory efficient, using DNAStrings instead of characters
	allKmers <- DNAStringSet()
	allCounts <- integer(0)
	nKmers <- 0

	cat( "\nFind union of all Sample Kmers..\n")
	for ( i in 1:NS) {
		f <- fileSet[i]
		loadSavedKmers(f)
		myKmers <- bigKmerStrings[[1]]
		myCounts <- bigKmerCounts[[1]]
		nNow <- length( myKmers)
		cat( "  N_In:", nNow)

		# if not enough Kmers, just discard it
		if ( sum(myCounts) < 100) {
			cat( "  Too Few, Ignore!")
			next
		}

		# during the scan phase, only Kmers that meet the minimum count
		# would get kept later, so no need to keep them at this point
		cat( "  scanDrops..")
		lowCnts <- which( myCounts < min.count)
		if ( length( lowCnts)) {
			cat( "  drop:", length(lowCnts))
			myKmers <- myKmers[ -lowCnts]
			myCounts <- myCounts[ -lowCnts]
		}
		if ( ! length( myKmers)) next

		# if the very first batch, just stow them
		if ( nKmers == 0) {
			allKmers <- myKmers
			allCounts <- myCounts
			nKmers <- nNow
			next
		}

		cat( "  lookup..")
		where <- match( myKmers, allKmers, nomatch=0)
		# those we do find get counts increased
		allCounts[where] <- allCounts[where] + myCounts[where > 0]
		# only add the new ones
		myKmers <- myKmers[ where == 0]
		myCounts <- myCounts[ where == 0]
		nNow <- length( myKmers)
		if ( nNow) {
			cat( "  N_New:", nNow)
			allKmers <- c( allKmers, myKmers)
			allCounts <- c( allCounts, myCounts)
			nKmers <- nKmers + nNow
		}
		cat( "  N_Kmer:", nKmers)

		# watch for too many Kmers, so we don't break XStrings
		min.count.now <- min.count 
		while ( nKmers > MAX_KMERS) {
			min.count.now <- min.count.now + 1
			tooFew <- which( allCounts < min.count.now)
			if ( length( tooFew)) {
				cat( "\n  Exceeded MAX_KMERS: drop low count Kmers < ", min.count.now)
				# if not enough flagged, just go around again right now
				if ( (nKmers - length(tooFew)) > MAX_KMERS) next
				allKmers <- allKmers[ -tooFew]
				allCounts <- allCounts[ -tooFew]
			}
			cat( " Now:", length(allCounts))
			nKmers <- length(allCounts)
		}
		if ( exists( "tooFew")) {
			rm( tooFew)
			gc()
		}

		rm( myKmers, myCounts, where)
		if ( i %% 4 == 0) gc()
	}
	NK <- nKmers
	# do some cleanup...
	if ( exists("bigKmerStrings", envir=.GlobalEnv)) rm( bigKmerStrings, envir=.GlobalEnv)
	if ( exists("bigKmerCounts", envir=.GlobalEnv)) rm( bigKmerCounts, envir=.GlobalEnv)
	gc()

	# now fill the table
	cat( "\nReloading to fill table..\n")
	kmerTbl <- matrix( 0, nrow=NK, ncol=NS)
	colnames(kmerTbl) <- sampleIDset
	rownames(kmerTbl) <- 1:NK

	for ( i in 1:NS) {
		f <- fileSet[i]
		loadSavedKmers(f)
		myKmers <- bigKmerStrings[[1]]
		myCounts <- bigKmerCounts[[1]]
		cat( "  lookup..")
		wh <- match( myKmers, allKmers, nomatch=0)
		# these may have been removed due to MAX_KMERS
		v <- rep.int(0,NK)
		cat( "  store..")
		v[ wh] <- myCounts[ wh > 0]
		kmerTbl[ , i] <- v
		# do some cleanup...
		rm( myKmers, myCounts, wh)
		if ( exists("bigKmerStrings", envir=.GlobalEnv)) rm( bigKmerStrings, envir=.GlobalEnv)
		if ( exists("bigKmerCounts", envir=.GlobalEnv)) rm( bigKmerCounts, envir=.GlobalEnv)
		if ( i %% 4 == 0) gc()
	}
	cat( "\nDone loading.\n")
	
	cat( "\nChecking for low coverage Kmers to drop: \n  At least", min.count, "counts in at least", min.samples, "samples..")
	nGood <- apply( kmerTbl, 1, function(x) sum( x >= min.count))
	drops <- which( nGood < min.samples)
	if ( length(drops)) {
		cat( "  Removing", length(drops), "Kmers rows..")
		kmerTbl <- kmerTbl[ -drops, ]
		allKmers <- allKmers[ -drops]
		cat( "  N_Kmer: ", nrow(kmerTbl))
	}
	if ( exists("bigKmerStrings", envir=.GlobalEnv)) rm( bigKmerStrings, envir=.GlobalEnv)
	if ( exists("bigKmerCounts", envir=.GlobalEnv)) rm( bigKmerCounts, envir=.GlobalEnv)
	gc()
	
	# lastly, put the kmers as names on the rows
	rownames(kmerTbl) <- as.character( allKmers)
	return( kmerTbl)
}


# do differential Kmer analysis by sample groups

`pipe.KmerCompare` <- function( kmerTbl, sampleIDset, groupSet, normalize=c("LKPTM"), 
				kmer.size=33, min.count=NULL, min.samples=NULL, n.remove=100) {

	# just just the samples asked for
	where <- match( sampleIDset, colnames(kmerTbl))
	if ( any( is.na(where))) stop( "Some SampleIDs not in Kmer table")
	if ( length(where) != ncol(kmerTbl) || !all( where == (1:ncol(kmerTbl)))) {
		useTbl <- kmerTbl[ , where]
	} else {
		useTbl <- kmerTbl
	}
	NC <- ncol(useTbl)
	NR <- nrow(useTbl)
	
	# get our group factors
	grpFac <- factor( groupSet)
	grpLvls <- levels( grpFac)
	nGrps <- nlevels( grpFac)
	if ( nGrps != 2) stop( "Expected exactly 2 Sample Groups")
	grp1 <- which( as.numeric(grpFac) == 1)
	grp2 <- which( as.numeric(grpFac) == 2)
	cat( "\nBreakdown by Group:\n")
	print( table( groupSet))
	
	# allow a further trimming of Kmers for low detection, after it was already done during table creation
	if ( ! ( is.null(min.count) && is.null(min.samples))) {
		if ( is.null(min.count) || is.null(min.samples)) stop( "Error: must give both 'min.count' and 'min.samples' values")
		cat( "\nChecking for low coverage Kmers to drop: \n  At least", min.count, "counts in at least", min.samples, "samples..")
		nGood <- apply( useTbl, 1, function(x) sum( x >= min.count))
		drops <- which( nGood < min.samples)
		if ( length(drops)) {
			cat( "  Removing", length(drops), "Kmers rows..")
			useTbl <- useTbl[ -drops, ]
			NR <- nrow(useTbl)
			cat( "  N_Kmer: ", NR)
		}
		rm( nGood)
	}

	# drop the super high count Kmers
	if (n.remove > 0) {
		useTbl <- remove.HighCountKmers( useTbl, n.remove=n.remove)
		NR <- nrow(useTbl)
	}

	# convert everything to normalized unit
	normalize <- match.arg( normalize)
	if ( normalize == "LKPTM") {
		cat( "\nNormalizing to LKPTM (Log2 Kmers Per Ten Million)..")
		useTbl <- as.LKPTM( useTbl)
	}
	gc()
	
	# we are now ready to do tests per Kmer
	cat( "\nRunning Linear Model on ", NR, " Kmers..\n")
	fold <- pval <- vector( mode="numeric", length=NR)
	avg <- matrix( 0, 2, NR)
	rownames(avg) <- paste( "Avg", grpLvls, "LKPTM", sep="_")
	cnt <- matrix( 0, 2, NR)
	rownames(cnt) <- paste( "N", grpLvls, "Samples", sep="_")

	for ( i in 1:NR) {
		v <- useTbl[ i, ]
		ans <- t.test( v[grp1], v[grp2])
		pval[i] <- ans$p.value
		avg[ ,i] <- ans$estimate
		fold[i] <- log2( (avg[2,i]+1) / (avg[1,i]+1))
		cnt[1,i] <- sum( v[grp1] > 0)
		cnt[2,i] <- sum( v[grp2] > 0)
		if ( i %% 10000 == 0) cat( "\r", i, rownames(useTbl)[i], fold[i], pval[i])
	}
	cat( "\nDone.\n")
	
	out <- data.frame( "Kmer"=rownames(useTbl), t(cnt), t(avg), "Log2.Fold"=fold, "P.Value"=pval,
			row.names=seq_len(NR), stringsAsFactors=F)
	ord <- diffExpressRankOrder( out$Log2.Fold, out$P.Value)
	out <- out[ ord, ]
	rownames(out) <- 1:NR
	return(out)
}


# add all the details end user may want about the Kmers, typically after the 'Compare' Step

`pipe.KmerAnnotation` <- function( kmerTbl, optionsFile="Options.txt", speciesID=getCurrentSpecies(), 
				quiet=TRUE, verbose=!quiet) {

	neededColumns <- c( "Kmer")
	if ( ! all( neededColumns %in% colnames(kmerTbl))) {
		cat( "\nError: Kmer table is missing some needed columns: ", neededColumns)
		return( NULL)
	}

	# 3 main steps

	# STep 1:  Align to Genome
	cat( "\n\nStep 1:  Align Kmers to Genome via Bowtie2..\n")
	ansAlign <- alignKmersToGenome( kmerTbl$Kmer, optionsFile=optionsFile, speciesID=speciesID,
				quiet=quiet, verbose=verbose)

	# Step 2:  map those aligned locations onto proteins
	cat( "\n\nStep 2:  Map Kmer Genome hits to Protein Fragments..\n")
	ansProt <- mapKmersToProteins( ansAlign, optionsFile=optionsFile, verbose=verbose)

	# Step 3:  discern protein fragment difference from expected reference
	cat( "\n\nStep 3:  Compare Protein Fragments to Reference..\n")
	ansDiff <- mapProteinFragmentToReference( ansProt, optionsFile=optionsFile, verbose=verbose)

	# lastly, sew it all back together in a sensible order...
	cat( "\n\nStep 4:  Combine and organize all results..")
	out <- cbind( ansProt, "DIF_FROM_REF"=ansDiff, kmerTbl[ ,2:ncol(kmerTbl)], stringsAsFactors=F)
	cat( "\nDone.\n")
	return( out)
}


# align Kmers to the genome

`alignKmersToGenome` <- function( kmers, optionsFile="Options.txt", speciesID=getCurrentSpecies(), 
				quiet=TRUE, verbose=!quiet) {

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

	# to interpret the Bowtie results, make sure we are the corrent target and species
	# force knowing the species before we look at target
	mySpecies <- speciesID
	target <- getAndSetTarget( optionsFile, verbose=TRUE)
	setCurrentSpecies( mySpecies)
	
	# build a call to Bowtie2
	ans <- kmerCallBowtie( kmerFastaFile, kmerBamFile, optionsFile=optionsFile, quiet=quiet, verbose=verbose)
	
	# read up the BAM file contents to extract what we need
	ans2 <- kmerReadBam( kmerBamFile, verbose=verbose)
	
	# force the results to match the input order
	wh <- match( kmers, ans2$Kmer)
	out <- ans2[ wh, ]
	return( out)
}


# further map aligned Kmers to protein location & sequence

`mapKmersToProteins` <- function( kmerAlignments, optionsFile="Options.txt", verbose=TRUE) {


	# takes the output from 'alignKmersToGenome', and looks up where those sites 
	# land on proteins
	if ( ! all( c( "SEQ_ID", "POSITION", "GENE_ID") %in% colnames( kmerAlignments))) {
		cat( "\nWarning:  expected a data.frame with SEQ_ID, POSITION, GENE_ID columns.")
		cat( "\nPerhaps run 'alignKmersToGenome()' first...")
		return( NULL)
	}

	N <- nrow(kmerAlignments)
	
	genome.file <- getOptionValue( optionsFile, "genomicFastaFile", verbose=verbose)
	# pre-fetch one chromosome for speed and better diagnostics
        gdna <- getFastaSeqFromFilePath( filePath=genome.file, seqID=getCurrentSeqMap()$SEQ_ID[1],
					verbose=verbose)

	# for Kmers that hit genes, let's map to AA location and guess the protein fragment
	aaPos <- rep.int( NA, N)
	aaFrag <- rep.int( "", N)
	#refFrag <- rep.int( "", N)
	geneHits <- setdiff( 1:N, grep( "(ng)", kmerAlignments$GENE_ID, fixed=T)) 
	if (verbose) cat( "\nMapping", length(geneHits), "Kmer Gene hits to AA location and protein fragment sequence..\n")

	# visit them in gene order to make it faster
	geneFac <- factor( kmerAlignments$GENE_ID[ geneHits])
	require(Biostrings)
	data(BLOSUM62)

	nDone <- 0
	tapply( geneHits, geneFac, function(x) {
		
		# given all the rows that belong to one gene
		i <- x[1]
		mySeqID <- kmerAlignments$SEQ_ID[i]
		myGeneID <- kmerAlignments$GENE_ID[i]
		if ( is.na(mySeqID) || is.na(myGeneID)) return()
		nDone <<- nDone + length(x)

		myPos <- kmerAlignments$POSITION[x]
		if ( any( is.na(myPos))) {
			x <- x[ ! is.na(myPos)]
			myPos <- kmerAlignments$POSITION[x]
		}

		# pre fetch some maps
		curGmap <- subset.data.frame( getCurrentGeneMap(), GENE_ID == myGeneID)
		curCDSmap <- subset.data.frame( getCurrentCdsMap(), GENE_ID == myGeneID)

		smlAns <- convertGenomicDNApositionToAAposition( mySeqID, myPos, geneID=myGeneID,
								genemap=curGmap, cdsmap=curCDSmap)
		aaPos[x] <<- myAApos <- smlAns$AA_POSITION

		# fetch the reference protein
		refProtein <- gene2Fasta( myGeneID, genome.file, mode="aa")$seq[1]
		if ( !is.na(refProtein)) {
			aaFrag[x] <<- myFrags <- DNAtoBestPeptide( kmerAlignments$Kmer[x], clip=F, readingFrame=1:6, 
							tieBreakMode="reference", reference=refProtein,
							substitutionMatrix=BLOSUM62)
		} else {
			for (j in x) {
				readFrame <- if ( kmerAlignments$STRAND[j] == "+") 1:3 else if (kmerAlignments$STRAND[j] == "-") 4:6 else 1:6
				aaFrag[j] <<- DNAtoBestPeptide( kmerAlignments$Kmer[j], clip=F, readingFrame=readFrame)
			}
			myFrags <- aaFrag[x]
		}

		if (verbose) cat( "\r", nDone, myGeneID, length(x), aaFrag[i])
		return(NULL)
	})
	cat( "\nDone.\n")
	
	out <- kmerAlignments
	out$AA_POSITION <- aaPos
	out$KMER_FRAGMENT <- aaFrag
	#out$REF_FRAGMENT <- refFrag
	return( out)
}


# try to show how the protein fragment differs from reference protein

`mapProteinFragmentToReference` <- function( kmerAlignments, optionsFile="Options.txt", verbose=TRUE) {


	# takes the output from 'alignKmersToGenome', and looks up where those sites 
	# land on proteins
	if ( ! all( c( "SEQ_ID", "GENE_ID", "AA_POSITION", "KMER_FRAGMENT") %in% colnames( kmerAlignments))) {
		cat( "\nWarning:  expected a data.frame with SEQ_ID, AA_POSITION, GENE_ID, KMER_FRAGMENT columns.")
		cat( "\nPerhaps run 'alignKmersToGenome()' and/or 'mapKmersToProteins()' first...")
		return( NULL)
	}

	N <- nrow(kmerAlignments)
	
	genome.file <- getOptionValue( optionsFile, "genomicFastaFile", verbose=verbose)
	# pre-fetch one chromosome for speed and better diagnostics
        gdna <- getFastaSeqFromFilePath( filePath=genome.file, seqID=getCurrentSeqMap()$SEQ_ID[1],
					verbose=verbose)

	# for Kmers that hit genes, let's map to AA location and guess the protein fragment
	aaDiff <- rep.int( "N/A", N)
	geneHits <- setdiff( 1:N, grep( "(ng)", kmerAlignments$GENE_ID, fixed=T)) 
	if (verbose) cat( "\nMapping", length(geneHits), "Kmer Gene hits to AA location and protein fragment sequence..\n")

	# only try to map those with a valid fragment
	noFrag <- which( kmerAlignments$KMER_FRAGMENT == "")
	geneHits <- setdiff( geneHits, noFrag)

	# visit them in gene order to make it faster
	geneFac <- factor( kmerAlignments$GENE_ID[ geneHits])
	require(Biostrings)
	data(BLOSUM62)

	nDone <- 0
	tapply( geneHits, geneFac, function(x) {
		
		# given all the rows that belong to one gene
		i <- x[1]
		mySeqID <- kmerAlignments$SEQ_ID[i]
		myGeneID <- kmerAlignments$GENE_ID[i]
		if ( is.na(mySeqID) || is.na(myGeneID)) return()
		nDone <<- nDone + length(x)

		# pre fetch some maps
		curGmap <- subset( getCurrentGeneMap(), GENE_ID == myGeneID)
		curCDSmap <- subset( getCurrentCdsMap(), GENE_ID == myGeneID)

		# fetch the reference protein
		refProtein <- gene2Fasta( myGeneID, genome.file, mode="aa")$seq[1]
		if ( !is.na(refProtein)) {
			myFrags <- kmerAlignments$KMER_FRAGMENT[x]
			pa <- pairwiseAlignment( myFrags, refProtein, type="global-local", scoreOnly=F)
			refFrags <- as.character( subject(pa))

			# try to tell how the Kmer differs from the reference
			# put a dot notation to show just what's different
			refAAvec <- strsplit( refFrags, split="")
			myAAvec <- strsplit( myFrags, split="")
			lapply( 1:length(myFrags), function(i) {
				chUse <- 1 : min( length(refAAvec[[i]]), length(myAAvec[[i]]))
				same <- which( refAAvec[[i]][chUse] == myAAvec[[i]][chUse])
				if ( length(same)) {
					tmp <- myAAvec[[i]]
					tmp[same] <- "."
					myFrags[i] <<- PASTE(tmp,collapse="")
				}
				return(NULL)
			})
			aaDiff[x] <<- myFrags
		}

		if (verbose) cat( "\r", nDone, myGeneID, length(x), aaDiff[i])
		return(NULL)
	})
	cat( "\nDone.\n")
	return( aaDiff)
}


# Kmerize a single FASTQ file


# Kmerize a single FASTQ file

`kmerizeOneFastqFile` <- function( filein, kmer.size=33, buffer.size=1000000, 
				sampleID="SampleID", kmer.path=".", maxReads=NULL, min.count=2) {

	# get ready to read the fastq file in chuncks...
	filein <- allowCompressedFileName( filein)
	if ( ! file.exists( filein)) stop( paste("Can't find input file: ", filein))
	conIn <- openCompressedFile( filein, open="r")
	on.exit( close( conIn))
	require(Biostrings)

	# 4 lines per read
	chunkSize <- buffer.size
	readChunkSize <- chunkSize * 4
	nread <- 0
	hasMore <- TRUE
	firstPass <- TRUE

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
		cat( "  N_Reads: ", prettyNum( as.integer(nread), big.mark=","))
		rm( chunk)

		timeStat <- round( system.time( kmerizeOneChunk.as.DNAString( seqTxt, kmer.size=kmer.size)), digits=1)
		cat( "  Time(usr,sys)=",timeStat[3],"(",timeStat[1],",",timeStat[2],")", sep="")
		rm( seqTxt)

		# each call just extends the global lists, just clean up memory
		gc()
	}

	# with the file done, merge all the chunk results

	# any K-mers with too few observations are not to be kept
	mergeKmerChunks( min.count)

	# also do any RevComp resolving now too
	if ( bigKmerCounts[[1]] > 0) {
		cat( "\nSearch for RevComp pairs to join..")
		kmerRevComp.GlobalTable( sampleID=sampleID, kmer.path=kmer.path, kmer.size=kmer.size)
	}

	return( list( "nReadsIn"=nread))
}


mergeKmerChunks <- function( min.count) {

	# the FASTA file is now a big list of DNAStrings and Counts, in global storage
	# step one is merge.   Within each chunk we are unique, 
	# but find where we have dups in later chunks.  And Join and Extend as needed...
	bigStrings <- bigKmerStrings[[1]]
	bigCounts <- bigKmerCounts[[1]]

	if ( (N <- length(bigKmerStrings)) > 1) {
	    cat( "\nMerge Chunks: ")
	    for ( j in 2:N) {
		strs2 <- bigKmerStrings[[j]]
		cnts2 <- bigKmerCounts[[j]]
		cat( "  Chunk",j)
		where <- match( strs2, bigStrings, nomatch=0)
		hitsFrom <- which( where > 0)
		hitsTo <- where[ hitsFrom]
		if ( length( hitsFrom)) {
			cat( " Hit:",length(hitsFrom))
			bigCounts[hitsTo] <- bigCounts[hitsTo] + cnts2[hitsFrom]
		}
		missFrom <- which( where == 0)
		if ( length( missFrom)) {
			cat( " New:",length(missFrom))
			bigStrings <- c( bigStrings, strs2[missFrom])
			bigCounts <- c( bigCounts, cnts2[missFrom])
		}
		cat( " Now:", length(bigCounts))
		
		# any time we exceed the maximum number of Kmers, trim away low counts to prevent crashing XVector hash size
		min.count.now <- 1
		while ( length(bigCounts) > MAX_KMERS) {
			min.count.now <- min.count.now + 1
			tooFew <- which( bigCounts < min.count.now)
			if ( length( tooFew)) {
				cat( "\n  Exceeded MAX_KMERS: drop low count Kmers < ", min.count.now)
				bigStrings <- bigStrings[ -tooFew]
				bigCounts <- bigCounts[ -tooFew]
			}
			cat( " Now:", length(bigCounts))
		}
		if ( exists( "tooFew")) {
			rm( tooFew)
			gc()
		}

		# remove the old elements after we merged them
		bigKmerStrings[[j]] <<- vector()
		bigKmerCounts[[j]] <<- vector()
		if ( j %% 5 == 0) gc()
		
	    }
	    rm( strs2, cnts2, where, hitsFrom, hitsTo, missFrom)
	    length(bigKmerStrings) <<- length(bigKmerCounts) <<- 1
	    gc()
	}

	cat( "\nDrop low Kmer counts < ", min.count)
	tooFew <- which( bigCounts < min.count)
	if ( length( tooFew)) {
		cat( "  Found:", length(tooFew))
		bigStrings <- bigStrings[ -tooFew]
		bigCounts <- bigCounts[ -tooFew]
		cat( "  N_Kmer:", length( bigCounts))
	}

	# stuff these merged results back now
	bigKmerStrings[[1]] <<- bigStrings
	bigKmerCounts[[1]] <<- bigCounts

	rm( bigStrings, bigCounts)
	gc()
	return(NULL)
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


`kmerizeOneChunk.as.DNAString` <- function( seqs, kmer.size=33) {

	# make Kmers of every raw read in this chunk
	# set up to get all N-mers for all the sequences in this chunk
	sizeM1 <- kmer.size - 1
	nSeq <- length( seqs)
	seqLens <- base::nchar( seqs)

	# ignore any with N for building the Kmers
	hasN <- grep( "N", seqs, fixed=T)
	
	# initial where we will store the growing output
	kmerSets <- vector( mode="list", length=nSeq)
	nOut <- 0

	# we can save a bit of time by doing all reads of the same length at one time
	lenFac <- factor( seqLens)
	cat( "  Kmerize..")
	tapply( 1:nSeq, lenFac, function(x) {
			
		# given all the pointers to seqs of the same length
		lenNow <- seqLens[ x[1]]
		veryLastStart <- lenNow - sizeM1
		if ( veryLastStart < 1) return(NULL)
		fromSet <- 1:veryLastStart
		toSet <- fromSet + sizeM1
		nSubstrs <- length( fromSet)

		x <- setdiff( x, hasN)
		
		base::lapply( x, function(j) {
			kmers <- subseq( rep.int( seqs[j], nSubstrs), fromSet, toSet)
			#if ( j %in% hasN) kmers <- kmers[ -grep("N",kmers)]
			nOut <<- nOut + 1
			kmerSets[[nOut]] <<- kmers
			return(NULL)
		})
		return(NULL)
	})
	rm( seqLens, hasN, lenFac)

	# turn from list back to one giant DNAStringSet
	cat( " tabulate..")
	length(kmerSets) <- nOut
	kmersOut <- do.call( c, kmerSets)
	rm( kmerSets)
	
	# then do the counting
	cntsV <- table.nosort( kmersOut)
	rm( kmersOut)

	cat( " store..")
	nNow <- length(bigKmerStrings) + 1
	bigKmerCounts[[nNow]] <<- as.vector(cntsV)
	bigKmerStrings[[nNow]] <<- DNAStringSet( names(cntsV))
	cat( "  N_Kmer:", length(cntsV))
	return( NULL)
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
	require(Biostrings)
	
	# faster memory efficient to operate on one global vector
	# stored as a list of DNASTrings and Counts
	kmers <- bigKmerStrings[[1]]
	cnts <- bigKmerCounts[[1]]
	N <- length(kmers)
	cat( "\nN_Kmers in: ", N)

	# see what the RevComp of every Kmer is
	rcKmer <- DNAStringSet( findKmerRevComp( as.character(kmers), sampleID=sampleID, kmer.path=kmer.path, kmer.size=kmer.size))
	gc()
	
	# then see if those Rev Comps are already in the table
	cat( "\nLocate RevComp pairs..")
	where <- match( rcKmer, kmers, nomatch=0)
	rm( rcKmer)

	# we will make a new 1-D table that has just the first form of each Kmer
	# the Kmers are not ordered anymore, but this test of which to keep needs/assumes ordering
	# so generate it now
	ord <- order( kmers)
	rankOrder <- seq_len(N)
	rankOrder[ ord] <- rankOrder
	rm( ord)

	# we can combine those that have both forms
	cat( "  join..")
	cntOut <- rep.int( 0, N)
	hasRC <- which( where > 0)
	if ( length(hasRC)) {
		myLoc <- (1:N)[hasRC]
		rcLoc <- where[hasRC]
		myRank <- rankOrder[ myLoc]
		rcRank <- rankOrder[ rcLoc]
		#firstLoc <- pmin( myLoc, rcLoc)
		firstLoc <- ifelse( myRank < rcRank, myLoc, rcLoc)
		cntOut[firstLoc] <- cnts[myLoc] + cnts[rcLoc]
		rm( myLoc, rcLoc, firstLoc, myRank, rcRank)
	}
	rm( rankOrder)

	# then those without their RevComp just stay as is
	noRC <- which( where == 0)
	if ( length(noRC)) {
		cntOut[noRC] <- cnts[noRC]
	}
	rm(where)

	# lastly thown away the empty slots
	drops <- which( cntOut == 0)
	if ( length( drops)) {
		cntOut <- cntOut[ -drops]
		kmerOut <- kmers[ -drops]
	} else {
		kmerOut <- kmers
	}
	Nout <- length(cntOut)
	cat( "\nN_Kmers out: ", Nout)
	cat( "  Pct Reduction: ", round( (N-Nout) * 100 / N), "%")

	# done, stash these back in global storage
	bigKmerStrings[[1]] <<- kmerOut
	bigKmerCounts[[1]] <<- cntOut

	rm( cntOut, kmerOut, cnts, hasRC, noRC, drops)
	gc()
	return( Nout)
}


findKmerRevComp <- function( kmers, sampleID=sampleID, kmer.path=".", kmer.size=33) {

	# set up to allow 'faster' lookup of previously calculated Rev Comps
	N <- length( kmers)
	out <- rep.int( "", N)
	require(Biostrings)
	
	# allow 'fast' lookup via a file of previously done kmers
	# we are switching to separate files for every sample, to try to minimize file
	# access collisions.  You can read fron any sample's file, but only write to
	# your own.
	myKmerXrefFile <- file.path( kmer.path, paste( sampleID, kmer.size, "Kmer.RevComp.Xref.rda", sep="."))
	allKmerXrefFiles <- dir( kmer.path, pattern="Kmer.RevComp.Xref.rda$", full=T)
	
	if ( length(allKmerXrefFiles)) {
	    cat( "  searching", length(allKmerXrefFiles), "RevComp Xref files:")
	    nLook <- 0
	    for ( j in 1:length(allKmerXrefFiles)) {
	    	f <- allKmerXrefFiles[j]
		xref <- data.frame()
		status <- tryCatch( load( f), error=function(e) { xref <<- data.frame()})
		if ( ! nrow( xref)) next
	
		# we got some from one Xref file, see if we get any hits
		# but only check on the ones we still need on this pass around
		toTest <- which( out == "")
		cat( "  Lookup:", nLook <- nLook+1)
		where1 <- match( kmers[toTest], xref$Kmer, nomatch=0)
		out[ toTest[ where1 > 0]] <- xref$RevComp[ where1]
		where2 <- match( kmers[toTest], xref$RevComp, nomatch=0)
		out[ toTest[ where2 > 0]] <- xref$Kmer[ where2]
		cat( "  Found:", sum(where1 > 0 | where2 > 0))
		rm( xref, toTest, where1, where2)
		if ( j %% 5 == 0) gc()
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
	if ( ! nDrop) return( kmerTbl)

	mDrop <- kmerTbl[ whoHigh, ]
	kmerCountDrop <- sum( mDrop)
	kmerCountIn <- sum( kmerTbl)

	cat( "\nHigh Count Removal flagged", nDrop, "Kmers")
	cat( "\n  Being ", round( kmerCountDrop * 100 / kmerCountIn, digits=2), "% of all Kmer counts")
	cat( "\nTop Culprits:\n")
	n.show <- min( nDrop, 20)
	n.col.show <- min( ncol(mDrop), 10)
	# sort to show the biggest values
	ord <- order( apply( mDrop, 1, sum), decreasing=T)
	print( head( mDrop[ ord, 1:n.col.show], n.show))

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
	# next is explicit alignment policy.   Since these are short Kmers that may not match the genome well
	# we need to be quite lax to get the best chance of getting some hit for every kmer
	# cmdLine <- paste( cmdLine, " --very-sensitive")  # default is " -D 20 -R 3 -N 0 -L 20 -i S,1,0.50"
	cmdLine <- paste( cmdLine, " -D 30 -R 5 -N 0 -L 7 -i C,1,0")
	if ( quiet) cmdLine <- paste( cmdLine, " --quiet")
	# when not catching the unaligned, allow explicit throw-away of unaligned
	cmdLine <- paste( cmdLine, " --no-unal")
	# next is threads
	nCores <- as.integer( getOptionValue( optT, "nCores", notfound="4", verbose=FALSE))
	cmdLine <- paste( cmdLine, " --threads", nCores)
	
	# let's use a slightly more lax scoring threshold, to get a few more aligned sites for variant Kmers
	cmdLine <- paste( cmdLine, " --score-min L,-0.9,-0.9 ")

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
		if (verbose) cat( "  N_Alignments: ", prettyNum( as.integer( nReads), big.mark=","))

		# convert the refID to a SeqID to figure the gene locations
		if (verbose) cat( "  geneIDs..")
		seqIDs <- refID2seqID( refID( chunk), refData=refData)
		positions <- position( chunk)
		kmerSeq <- readSeq( chunk)
		lens <- nchar( kmerSeq)
		middles <- positions + floor( lens/2)
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
	if (verbose) cat("\nDone.\n")
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


`saveKmers` <- function(outfile, mode="first.time") {

	# push the global results to disk
	# turn the DNAString data to character
	require(Biostrings)
	cat( "  saving..")
	bigKmerCounts <- as.vector( bigKmerCounts[[1]])

	# the first time we make Kmers from larger DNAStrings, there seems to be excess memory use.
	# so convert them to character and then back to DNAStrings
	if (mode == "first.time") {
		bigKmerStrings <- as.character( bigKmerStrings[[1]])
		bigKmerStrings <- DNAStringSet( bigKmerStrings)
	}
	save( bigKmerStrings, bigKmerCounts, file=outfile)

	if ( exists("bigKmerStrings")) rm( bigKmerStrings, bigKmerCounts)
	gc()
}


`loadSavedKmers` <- function(file) {

	# load the one global results variables from disk
	# turn the character data back to DNAString 
	require(Biostrings)
	cat( "  loading:", basename(file))
	load( file)

	# our convention is to use lists, of DNAStrings and Counts
	tmpCounts <- as.integer( bigKmerCounts)
	bigKmerCounts <<- list( tmpCounts)

	# the Kmers may already be S4 clases, or they may be just character strings
	if ( is.character( bigKmerStrings)) {
		tmpStrings <- DNAStringSet( bigKmerStrings)
	} else {
		tmpStrings <- bigKmerStrings
	}
	bigKmerStrings <<- list( tmpStrings)

	# remove the local copies
	rm( bigKmerStrings, bigKmerCounts, tmpCounts, tmpStrings)
	gc()
	return( length( bigKmerCounts[[1]]))
}


`checkOrCallCutadapt` <- function( inputFastqFiles, kmer.size=33, cutadaptProgram=Sys.which("cutadapt"), 
				asMatePairs=FALSE, forceMatePairs=NULL, verbose=TRUE) {
		
	# see if the files we need are already there
	if (verbose) cat( "\n")
	expectedFiles <- inputFastqFiles
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
	return( filesToDo)
}

