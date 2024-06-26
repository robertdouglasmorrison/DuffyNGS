# pipe.Kmerize.R

# turn raw FASTQ reads into a 1-D table of Kmer counts

# using Biostrings 'DNAStrings' to improve speed/capacity.  THere is is strict limit on
# XVector hash table size, that we can never exceed
MAX_KMERS <- 250000000


# turn raw FASTQ into Kmers for one sample

`pipe.Kmerize` <- function( sampleID=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt", 
				kmer.size=33, doCutadapt=TRUE, cutadaptProgram=Sys.which("cutadapt"), forceMatePairs=NULL, 
				buffer.size=1000000, maxReads=NULL, min.count=2, kmer.debug=NULL, verbose=TRUE) {

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nStarting  'Kmerize Pipeline' on Sample:     ", sampleID,
			"\nStart Date/Time:   \t", date(), 
			"\nKmer Size:         \t", kmer.size, "\n")
	}
	require(Biostrings)
	if (version$major == "4" && as.numeric( version$minor) >= 4) require( pwalign)

	# file(s) to process comes from annotation and options...
	results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	rawFastq <- getRawFastqFileNames( sampleID, annotationFile, optionsFile, verbose=FALSE)
	inputFastqFiles <- rawFastq$files
	asMatePairs <- rawFastq$asMatePairs

	# allow trimming of raw reads before Kmerize (but afgter Cutadapt..)
	trim <- getReadTrimming( sampleID=sampleID, annotationFile=annotationFile, optionsFile=optionsFile, verbose=verbose)
	trim5 <- trim$trim5
	trim3 <- trim$trim3

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
	# use one GLOBAL list of DNA strings and counts
	bigKmerStrings <<- bigKmerCounts <<- vector( mode="list")
	gc()

	cat( "\nN_Files to Kmerize: ", nFiles, "\n")
	for ( i in 1:nFiles) {
		# do one file, and save what we got
		ans <- kmerizeOneFastqFile( filesToDo[i], kmer.size=kmer.size, buffer.size=buffer.size, 
					maxReads=maxReads, min.count=min.count, kmer.debug=kmer.debug,
					trim5=trim5, trim3=trim3)
		nReadsNow <- ans$nReadsIn
		nReadsIn <- nReadsIn + nReadsNow
		saveKmers( outfile)

		# report those matrics
		nDistinctKmers <- length( bigKmerCounts[[1]])
		nTotalKmers <- sum( bigKmerCounts[[1]])
		cat( "\nSample: ", sampleID, "  File: ", i, basename(filesToDo[i]), "  N_Distinct: ", nDistinctKmers, 
				"  N_Total: ", nTotalKmers, "\n")
		# treat the 'maxReads' as a "per file" max, so both Fwd & Rev get equal treatment...
		#if ( ! is.null(maxReads) && nReadsIn >= maxReads) {
		#	cat( "\nReached 'maxReads' count. Stopping early..")
		#	break
		#}
	}

	# now do the cleanup...
	if ( exists("bigKmerStrings", envir=.GlobalEnv)) rm( bigKmerStrings, envir=.GlobalEnv)
	if ( exists("bigKmerCounts", envir=.GlobalEnv)) rm( bigKmerCounts, envir=.GlobalEnv)
	gc()

	# one new last step:  Force all Kmers to be in RevComp alphabetical ordering, so when we later do table compares,
	# we can garauntee that we don't have both a Kmer and its own RevComp present in 2 different files
	kmerAlphaRevCompOrder( outfile)

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
	if (version$major == "4" && as.numeric( version$minor) >= 4) require( pwalign)

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
			nNow <- length( myKmers)
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
	
	if ( min.count > 1 || min.samples > 1) {
		cat( "\nChecking for low coverage Kmers to drop: \n  At least", min.count, "counts in at least", min.samples, "samples..")
		nGood <- apply( kmerTbl, 1, function(x) sum( x >= min.count))
		drops <- which( nGood < min.samples)
		if ( length(drops)) {
			cat( "  Removing", length(drops), "Kmers rows..")
			kmerTbl <- kmerTbl[ -drops, ]
			allKmers <- allKmers[ -drops]
			cat( "  N_Kmer: ", nrow(kmerTbl))
		}
	}

	# clean up
	if ( exists("bigKmerStrings", envir=.GlobalEnv)) rm( bigKmerStrings, envir=.GlobalEnv)
	if ( exists("bigKmerCounts", envir=.GlobalEnv)) rm( bigKmerCounts, envir=.GlobalEnv)
	gc()
	
	# lastly, put the kmers as names on the rows
	rownames(kmerTbl) <- as.character( allKmers)
	return( kmerTbl)
}


# do differential Kmer abundance analysis by sample groups, using normalized data

`pipe.KmerCompare` <- function( kmerTbl, sampleIDset, groupSet, levels=sort(unique(groupSet)), normalize=c("LKPTM"), 
				min.count=NULL, min.samples=NULL, n.remove=100, offset=0.1) {

	# use just the samples asked for, and make sure the columns are in sample & group order
	where <- match( sampleIDset, colnames(kmerTbl), nomatch=0)
	if ( any( where == 0)) stop( "Some SampleIDs not in Kmer table")
	if ( length(where) != ncol(kmerTbl) || !all( where == (1:ncol(kmerTbl)))) {
		useTbl <- kmerTbl[ , where]
	} else {
		useTbl <- kmerTbl
	}
	NC <- ncol(useTbl)
	NR <- nrow(useTbl)
	
	# get our group factors
	if ( ! is.factor( groupSet)) {
		grpFac <- factor( groupSet, levels=levels)
	} else {
		grpFac <- groupSet
		groupSet <- as.character(groupSet)
	}
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
		normTbl <- as.LKPTM( useTbl)
	} else {
		stop( "Unknown or Unsupported Kmer normalization method: ", normalize)
	}
	gc()
	
	# we are now ready to do tests per Kmer
	cat( "\nRunning Differential Detection on ", NR, " Kmers..\n")
	fold <- pval <- vector( mode="numeric", length=NR)
	avg <- avg.cnt <- matrix( 0, 2, NR)
	rownames(avg) <- paste( "Avg", grpLvls, "LKPTM", sep="_")
	rownames(avg.cnt) <- paste( "Avg", grpLvls, "Depth", sep="_")
	cnt <- matrix( 0, 2, NR)
	rownames(cnt) <- paste( "N", grpLvls, "Samples", sep="_")

	for ( i in 1:NR) {
		v <- normTbl[ i, ]
		vCnt <- useTbl[ i, ]
		ans <- sparse.t.test( v[grp1], v[grp2], min.random.value=0)
		pval[i] <- ans$p.value
		avg[ ,i] <- ans$estimate
		fold[i] <- log2( (avg[2,i]+offset) / (avg[1,i]+offset))
		cnt[1,i] <- sum( v[grp1] > 0)
		cnt[2,i] <- sum( v[grp2] > 0)
		avg.cnt[1,i] <- mean( vCnt[grp1])
		avg.cnt[2,i] <- mean( vCnt[grp2])
		if ( i %% 10000 == 0) cat( "\r", i, rownames(useTbl)[i], fold[i], pval[i])
	}
	cat( "\nDone.\n")
	rm( normTbl)
	
	out <- data.frame( "Kmer"=rownames(useTbl), t(cnt), round(t(avg.cnt)), round(t(avg),digits=3), 
			"Log2.Fold"=round(fold,digits=4), "P.Value"=pval,
			row.names=seq_len(NR), stringsAsFactors=F)
	#ord <- diffExpressRankOrder( out$Log2.Fold, out$P.Value)
	#out <- out[ ord, ]
	#rownames(out) <- 1:NR
	return(out)
}


# do differential Kmer abundance analysis using the EdgeR counts tool

`pipe.KmerEdgeR` <- function( kmerTbl, sampleIDset, groupSet, levels=sort(unique(groupSet)), 
				min.count=NULL, min.samples=NULL, n.remove=100, dispersion=0.05) {

	require( edgeR)

	# use just the samples asked for, and make sure the columns are in sample & group order
	where <- match( sampleIDset, colnames(kmerTbl), nomatch=0)
	if ( any( where == 0)) stop( "Some SampleIDs not in Kmer table")
	if ( length(where) != ncol(kmerTbl) || !all( where == (1:ncol(kmerTbl)))) {
		useTbl <- kmerTbl[ , where]
	} else {
		useTbl <- kmerTbl
	}
	NC <- ncol(useTbl)
	NR <- nrow(useTbl)
	
	# get our group factors
	if ( ! is.factor( groupSet)) {
		grpFac <- factor( groupSet, levels=levels)
	} else {
		grpFac <- groupSet
		groupSet <- as.character(groupSet)
	}
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
		cat( "  Removing high count Kmers..")
		useTbl <- remove.HighCountKmers( useTbl, n.remove=n.remove)
		NR <- nrow(useTbl)
	}
	gc()

	# we are now ready to do tests per Kmer
	cat( "\nRunning EdgeR Differential Detection on ", NR, " Kmers..\n")
	
	# step 1:  build the EdgeR object
	# the column flags for EdgeR are to end up as 1,2, where 1 is the one we want...
	# as in Group 1 divided by Group 2
	grpNumber <- rep( 2, times=length(groupSet))
	grpNumber[ groupSet == groupSet[2]] <- 1
	cat( "  Build DEGList..")
	ans <- DGEList( counts=useTbl, group=grpNumber)
	
	# step 2: filtering of low expression genes.  Not for Kmers
	
	# step 3: normalize the library sizes
	cat( "  Normalize..")
	ans <-  normLibSizes( ans)
	
	# step 4:  Dispersion -- EdgeR can take very long to calculate.  Allow a simple override
	if ( is.null(dispersion) || is.character(dispersion)) {
		cat( "  EstimateDispersion..")
		if ( is.null(dispersion)) {
			dispersion <- "auto"
		} else {
			dispersion <- dispersion[1]
		}
		ans <- estimateDisp( ans)
	} else {
		dispersion <- as.numeric(dispersion)
		cat( "  UseExplicitDispersion=", dispersion)
	}

	# step 5: do the DE assessment
	cat( "  ExactTest..")
	etAns <- exactTest( ans, pair=c(2,1), dispersion=dispersion)
	rm( ans); gc()
	
	cat( "  ExtractTags..")
	NR <- nrow(useTbl)
	edgeRout <- topTags( etAns, n=NR)$table
	# get what we want/need out
	fout <- edgeRout[ , 1]
	pout <- edgeRout[ , 3]
	qout <- edgeRout[ , 4]
	gnames <- rownames( edgeRout)
	rm( edgeRout);  gc()
	
	cat( "  BuildResults..")
	where <- match( rownames(useTbl), gnames, nomatch=0)
	fold <- pval <- qval <- rep.int(NA, NR)
	fold[ where > 0] <- fout[where]
	pval[ where > 0] <- pout[where]
	qval[ where > 0] <- qout[where]
	nPos <- avgCnt <- matrix( NA, nrow=NR, ncol=2)
	colnames(nPos) <- paste( "N.Samples", grpLvls, sep=".")
	colnames(avgCnt) <- paste( "Avg.Count", grpLvls, sep=".")
	for( i in 1:NR) {
		x <- useTbl[ i, ]
		nPos[ i, ] <- tapply( x, grpFac, function(y) sum( y > 0))
		avgCnt[ i, ] <- tapply( x, grpFac, mean)
	}
	cat( "\nDone.\n")
	
	out <- data.frame( "Kmer"=rownames(useTbl), nPos, round(avgCnt, digits=3), "Log2.Fold"=round(fold,digits=4), 
			"P.Value"=pval, "Q.Value"=round(qval,digits=5), row.names=seq_len(NR), stringsAsFactors=F)
	#ord <- diffExpressRankOrder( out$Log2.Fold, out$P.Value)
	#out <- out[ ord, ]
	#rownames(out) <- 1:NR
	return(out)
}


# add all the details end user may want about the Kmers, typically after the 'Compare' Step

`pipe.KmerAnnotate` <- function( kmerTbl, optionsFile="Options.txt", speciesID=getCurrentSpecies(), 
				quiet.bowtie=TRUE, verbose=TRUE) {

	neededColumns <- c( "Kmer")
	if ( ! all( neededColumns %in% colnames(kmerTbl))) {
		cat( "\nError: Kmer table is missing some needed columns: ", neededColumns)
		return( NULL)
	}

	# 3 main steps

	# Step 1:  Align to Genome
	cat( "\n\nStep 1:  Align Kmers to Genome via Bowtie2..\n")
	ansAlign <- alignKmersToGenome( kmerTbl$Kmer, optionsFile=optionsFile, speciesID=speciesID,
				quiet=quiet.bowtie, verbose=verbose)

	# Step 2:  map those aligned locations onto proteins
	cat( "\n\nStep 2:  Map Kmer Genome hits to Protein Fragments..\n")
	ansProt <- mapKmersToProteins( ansAlign, optionsFile=optionsFile, verbose=verbose)

	# Step 3:  discern protein fragment difference from expected reference
	cat( "\n\nStep 3:  Compare Protein Fragments to Reference..\n")
	ansDiff <- mapKmerProteinFragmentsToReference( ansProt, optionsFile=optionsFile, verbose=verbose)

	# lastly, sew it all back together in a sensible order...
	# but force it all to comply exactly with the order we got as input
	cat( "\n\nStep 4:  Combine and organize all results..")
	kmer.out <- kmerTbl$Kmer
	whereProt <- match( kmer.out, ansProt$Kmer)
	outProt <- ansProt[ whereProt, 2:ncol(ansProt)]
	outDiff <- ansDiff[ whereProt]
	out <- cbind( "Kmer"=kmer.out, outProt, "DIF_FROM_REF"=outDiff, kmerTbl[ ,2:ncol(kmerTbl)], stringsAsFactors=F)
	cat( "\nDone.\n")
	return( out)
}


# build DNA and AA contigs from Kmers unique to one comparison group

`pipe.KmerNovelContigs` <- function( kmerTbl, group, folder.keyword=group, min.count=5, min.kmers=100, 
				max.copies.out=5,
				velvet.path=dirname(Sys.which("velveth")),  velveth.args="", velvetg.args="", 
				min.contig.length=100, proteinsFastaFile=NULL, 
				min.score=100, verbose=FALSE) {

	wantedDepthColumn <- paste( "Avg", group, "Depth", sep="_")
	neededColumns <- c( "Kmer", wantedDepthColumn)
	if ( ! all( neededColumns %in% colnames(kmerTbl))) {
		cat( "\nError: Kmer table is missing some needed columns: ", neededColumns)
		return( NULL)
	}

	results.path <- paste( "Novel.Kmer.Contigs_", folder.keyword, sep="")
	if ( ! file.exists(results.path)) dir.create( results.path, recursive=T)
	
	# 3 main steps

	# Step 1:  Find the Kmers only seen by this group, and make a FASTA file for de novo
	cat( "\n\nStep 1:  Gather Novel Kmers..\n")
	otherDepthColumn <- setdiff( grep( "^Avg_.+_Depth$", colnames(kmerTbl), value=T), wantedDepthColumn)
	if ( length(otherDepthColumn) != 1) {
		cat( "\nFailed to find the 2 expected Kmer Depth columns..")
		return( NULL)
	}
	wantedCounts <- as.numeric( kmerTbl[[ wantedDepthColumn]])
	otherCounts <- as.numeric( kmerTbl[[ otherDepthColumn]])
	novelRows <- which( wantedCounts >= min.count & otherCounts == 0)
	if ( (NN <- length(novelRows)) < min.kmers) {
		cat( "\nToo few novel Kmers to bother making contigs: ", NN)
		return(NULL)
	}
	cat( "  N_Novel Kmers for group: ", group, "  = ", NN)
	# make multiple copies of each Kmer, using the depth, up to some max
	kmerOut <- kmerTbl$Kmer[novelRows]
	nEach <- kmerTbl[[wantedDepthColumn]][ novelRows]
	nEach[ nEach > max.copies.out] <- max.copies.out
	kmersOut <- rep( kmerOut, times=nEach)
	idsOut <- paste( "Kmer", rep(1:NN,times=nEach), "_copy", unlist(lapply(nEach,function(x) 1:x)), sep="")
	fa <- as.Fasta( idsOut, kmersOut)
	faFile <- file.path( results.path, paste( "Novel.Kmers_", group, ".fasta", sep=""))
	writeFasta( fa, faFile, line=100)
	
	# Step 2:  hand those Kmer to Velvet
	cat( "\n\nStep 2:  Do Vevet de novo assembly on Kmers..")
	kmer.size <- nchar( fa$seq[1])
	# when Velvet chps up the Kmers, it generates non-unique contigs, which violates the entire premise
	# so force it to use the same size
	#velvet.K <- round(kmer.size * 0.33) * 2 + 1	# make Velvet use an odd size about 2/3 of Kmer size
	velvet.K <- kmer.size - 2
	
	makeVelvetContigs( faFile, outpath=results.path,
		velvet.path=velvet.path, buildHash=TRUE, buildContigs=TRUE,
		kmerSize=velvet.K, minLength=min.contig.length, doPairedEnd=FALSE, minCoverage=NULL, 
		velveth.args=velveth.args, velvetg.args=velvetg.args, verbose=verbose)

	# Step 3:  turn those contigs into proteins
	cat( "\n\nStep 3:  Map de novo contigs to best Proteins..\n")
	contigs <- loadFasta( file.path( results.path, "contigs.fa"), short=FALSE)
	cat( "\nN_Velvet_Contigs: ", length( contigs$desc))
	if ( ! (NC <- length(contigs$desc))) {
		cat( "\n\nWarning:  Velvet returned no contigs..")
		return(0)
	}

	# convert to peptides
	peps <- DNAtoBestPeptide( contigs$seq, clipAtStop=FALSE)
	pepNames <- contigs$desc
	# drop those that are too short to be interesting
	min.aa.len <- min.contig.length / 3
	keep <- which( nchar(peps) >= min.aa.len)
	peps <- peps[keep]
	pepNames <- pepNames[keep]
	NC <- length(keep)
	faOut <- as.Fasta( paste( group, pepNames, sep="_"), peps)
	outfile <- paste( group, "_Novel.Kmer.Velvet.Peptides.fasta", sep="")
	outfile <- file.path( results.path, outfile)
	writeFasta( faOut, file=outfile, line.width=100) 
	cat( "\nWrote Kmer Novel Peptides .FASTA file: ", outfile)

	# perhaps also lastly map these novel peptides onto reference proteins
	if ( ! is.null( proteinsFastaFile)) {
		cat( "\nNow mapping peptides onto reference proteins..")
		nCoresWas <- multicore.currentCoreCount()
		if ( nCoresWas < 2) multicore.setup( 4)
		ans <- peptide2BestProtein( peps, proteinsFastaFile=proteinsFastaFile, 
				substitutionMatrix=BLOSUM62, tieBreakMode="all", details=T)
		bigAns <- data.frame()
		for (i in 1:length(ans)) { 
			if (is.null(ans[[i]])) next
			smlAns <- data.frame( "ContigID"=pepNames[i], ans[[i]], stringsAsFactors=F)
			bigAns <- rbind( bigAns, smlAns)
		}
		# put the biggest/best at the top
		ord <- order( bigAns$Score, decreasing=T)
		bigAns <- bigAns[ ord, ]
		keep <- which( bigAns$Score >= min.score)
		bigAns <- bigAns[ keep, ]
		rownames(bigAns) <- 1:nrow(bigAns)
		# write it out
		outfile <- paste( group, "_Novel.Kmer.Velvet.ProteinHits.csv", sep="")
		outfile <- file.path( results.path, outfile)
		write.csv( bigAns, outfile, row.names=F)
		cat( "\nWrote Kmer Novel Protein Hits .CSV file: ", outfile)
		cat( "\nN_High Scoring Protein Fragments: ", nrow(bigAns))
		if ( nCoresWas < 2) multicore.setup( 1)
	}
	cat( "\nDone.\n")
	return( NC)
}


# search all Kmers for hits to one given CDS protein sequence

`pipe.KmerSearchForProtein` <- function( kmerTbl, cds.seq, chunk.size=1000) {

	neededColumns <- c( "Kmer")
	if ( ! all( neededColumns %in% colnames(kmerTbl))) {
		cat( "\nError: Kmer table is missing some needed columns: ", neededColumns)
		return( NULL)
	}
	require( Biostrings)
	if (version$major == "4" && as.numeric( version$minor) >= 4) require( pwalign)
	
	# prep the cds protein sequence we will search for hits to
	cds.seq <- as.character( cds.seq[1])
	aa <- DNAtoAA( cds.seq, clipAtStop=F, readingFrame=1)
	cdsDNAstr <- DNAString( cds.seq)
	
	# prep all the Kmers
	kmers <- as.character( kmerTbl$Kmer)
	NKmer <- length(kmers)
	kmer.size <- nchar( kmers[1])
	kmer.hw.cds <- floor( kmer.size/2)
	kmer.hw.aa <- floor( kmer.size/6)
	
	# prep for doing a bunch of alignmetns
	dnaMatrix <- nucleotideSubstitutionMatrix()
	data(BLOSUM62)
	minScore <- round( kmer.size * 0.80)
	
	# we will visit in chunks, and collect Kmers that score well enough
	outKmer <- outAAfrag <- outAAdiff <- vector()
	outRow <- outCDSloc <- outAAloc <- outScore <- vector()
	nOut <- 0
	nDone <- 0
	cat( "\nSearching Kmers..\n")
	while (nDone < NKmer) {
		first <- (nDone+1)
		last <- min( NKmer, nDone+chunk.size)
		nextChunk <- first:last
		nDone <- last
		myKmers <- DNAStringSet( kmers[ nextChunk])
		myRcKmers <- reverseComplement( myKmers)
		
		# get the scores only for now
		scoreFwd <- pairwiseAlignment( myKmers, cdsDNAstr, type="global-local", scoreOnly=T, substitutionMatrix=dnaMatrix)
		scoreRev <- pairwiseAlignment( myRcKmers, cdsDNAstr, type="global-local", scoreOnly=T, substitutionMatrix=dnaMatrix)
		myScores <- pmax( scoreFwd, scoreRev)
		myHits <- which( myScores >= minScore)
		if ( ! (NHits <- length(myHits))) next
		
		# grab the right form of the high scoring Kmers
		bestKmers <- ifelse( scoreFwd[myHits] > scoreRev[myHits], as.character(myKmers[myHits]), as.character(myRcKmers[myHits]))
		bestScores <- myScores[myHits]
		
		# and see exactly where they land
		bestKmerStrs <- DNAStringSet( bestKmers)
		pa <- pairwiseAlignment( bestKmerStrs, cdsDNAstr, type="global-local", scoreOnly=F, substitutionMatrix=dnaMatrix)
		cdsStarts <- start( subject(pa))
		
		# calc all the facts we need
		aaStarts <- floor((cdsStarts-1)/3) + 1
		aaFrags <- sapply( 1:NHits, function(x) DNAtoBestPeptide( bestKmers[x], clipAtStop=F, tieBreak="reference", ref=aa))
		pa2 <- pairwiseAlignment( aaFrags, aa, type="global-local", scoreOnly=F, substitutionMatrix=BLOSUM62)
		aaCDS <- as.character( alignedSubject(pa2))
		aaChars <- strsplit( aaFrags, split="")
		cdsChars <- strsplit( aaCDS, split="")
		aaDiffs <- sapply( 1:NHits, function(x) {
				n <- min( length( ch1 <- aaChars[[x]]), length(ch2 <- cdsChars[[x]]))
				isSame <- which( ch1[1:n] == ch2[1:n])
				ch1[ isSame] <- "."
				return( paste( ch1, collapse=""))
		})
		
		# and finally save them up
		now <- (nOut+1) : (nOut+NHits)
		outKmer[now] <- bestKmers
		outScore[now] <- bestScores
		# convention is to report the center of the kmer, not the head
		outCDSloc[now] <- cdsStarts + kmer.hw.cds
		outAAloc[now] <- aaStarts + kmer.hw.aa
		outAAfrag[now] <- aaFrags
		outAAdiff[now] <- aaDiffs
		outRow[now] <- nextChunk[myHits]
		nOut <- nOut + NHits
		
		# go get anouther chunk
		cat( "\r(Kmers,Hits): ", nDone, nOut)
	}
	
	# pull it all together
	out1 <- kmerTbl[ outRow, ]
	out2 <- data.frame( "CDS_SCORE"=outScore, "CDS_POS"=outCDSloc, "AA_POS"=outAAloc, "AA_FRAGMENT"=outAAfrag, 
				"DIF_FROM_CDS"=outAAdiff, stringsAsFactors=FALSE)
	out <- cbind( out1, out2, stringsAsFactors=FALSE)
	return( out)
}


`plotKmerSearchForProteinDepth` <- function( kmerTbl, countColumns, col=rainbow(length(countColumns),end=0.7), 
						featureMap=NULL, label="", lwd=2, legend.cex=1, feature.cex=1,
						mode=c("depth","KPM"), min.score=NULL, forceYmax=NULL) {

	neededColumns <- c( "CDS_POS", "CDS_SCORE")
	if ( ! all( neededColumns %in% colnames(kmerTbl))) {
		cat( "\nError: Kmer table is missing some needed columns: ", neededColumns)
		return( NULL)
	}

	# extract the numeric counts we will show
	NR <- nrow( kmerTbl)
	cdsP <- as.numeric( kmerTbl$CDS_POS)
	cdsS <- as.numeric( kmerTbl$CDS_SCORE)
	cntsM <- as.matrix( kmerTbl[ , countColumns])
	NS <- ncol(cntsM)
	cdsOrd <- order( cdsP)
	if ( ! all( cdsOrd == 1:NR)) {
		cntsM <- cntsM[ ord, ]
		cdsP <- cdsP[ ord]
		cdsS <- cdsS[ ord]
	}

	# raw depth, or nornalized
	yLabel <- "Kmer Depth"
	mode <- match.arg( mode)
	if ( mode == "KPM") {
		for ( i in 1:NS) {
			v <- cntsM[ ,i]
			cntsM[ ,i] <- v * 1000000 / sum(v)
		}
		yLabel <- "Normalized Kmer Depth (KPM)"
	}

	# perhaps discard some low scoring Kmers?..
	if ( ! is.null( min.score)) {
		drops <- which( cdsS < min.score)
		if ( length(drops)) {
			cntsM <- cntsM[ -drops, ]
			cdsP <- cdsP[ -drops]
			NR <- nrow( cntsM)
		}
	}

	# there can be 2+ values at each CDS, so trim to just the deepest
	if ( any( duplicated( cdsP))) {
		cntsM <- apply( cntsM, MARGIN=2, function(x) tapply( x, factor( cdsP), sum))
		cdsP <- unique( cdsP)
		NR <- nrow(cntsM)
	}
	
	# set up to do the plotting
	xLimits <- range( c( 1, cdsP))
	bigX <- xLimits[2]
	yLimits <- range( c( 0, cntsM))
	if ( ! is.null( forceYmax)) yLimits[2] <- as.numeric( forceYmax)
	col <- rep( col, length.out=NS)

	if ( ! is.null( featureMap)) {
		featureColumns <- c( "GENE_ID", "POSITION", "END", "COLOR")
		if ( ! all( featureColumns %in% colnames(featureMap))) {
			cat( "\nWarning: Feature Map is missing some needed columns: ", featureColumns, "\nIgnoring..")
			featureMap <- NULL
		} else {
			yFeatureLow <- -yLimits[2] * 0.1
			yFeatureHigh <- yFeatureLow * 0.25
			yFeatureMid <- (yFeatureLow + yFeatureHigh) / 2
			yLimits[1] <- yFeatureLow * 1.05
		}
	}

	mainText <- paste( "Kmer Search: ", label)
	plot( 1, 1, type="n", main=mainText, xlab="Position on Target Sequence", ylab=yLabel,
			xlim=xLimits*c(1,1.1), ylim=yLimits*c(1,1.1))

	if ( ! is.null( featureMap)) {
		for ( j in 1:nrow(featureMap)) {
			rect( featureMap$POSITION[j], yFeatureLow, featureMap$END[j], yFeatureHigh, col=featureMap$COLOR[j], border=1)
			text( (featureMap$POSITION[j] + featureMap$END[j])/2, yFeatureMid, featureMap$GENE_ID[j], col=1, cex=feature.cex)
		}
	}

	for ( i in 1:NS) {
		vIn <- cntsM[ , i]
		vOut <- rep.int( 0, bigX)
		vOut[ cdsP] <- vIn
		lines( 1:bigX, vOut, lwd=lwd, lty=1, col=col[i])
	}

	legend( 'topright', colnames(cntsM), lwd=lwd, col=col, bg='white', cex=legend.cex)
	dev.flush()
	return(NULL)
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
	if (version$major == "4" && as.numeric( version$minor) >= 4) require( pwalign)
	data(BLOSUM62)

	geneMap <- getCurrentGeneMap()
	cdsMap <- getCurrentCdsMap()

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
		curGmap <- subset.data.frame( geneMap, GENE_ID == myGeneID)
		curCDSmap <- subset.data.frame( cdsMap, GENE_ID == myGeneID)

		smlAns <- convertGenomicDNApositionToAAposition( mySeqID, myPos, geneID=myGeneID,
								genemap=curGmap, cdsmap=curCDSmap)
		aaPos[x] <<- myAApos <- smlAns$AA_POSITION

		# fetch the reference protein
		refProtein <- gene2Fasta( myGeneID, genome.file, mode="aa")$seq[1]
		if ( !is.na(refProtein)) {
			aaFrag[x] <<- myFrags <- DNAtoBestPeptide( kmerAlignments$Kmer[x], clip=F, readingFrame=1:6, 
							tieBreakMode="reference", reference=refProtein,
							substitutionMatrix=BLOSUM62, breakAtStops=F)
		} else {
			for (j in x) {
				readFrame <- if ( kmerAlignments$STRAND[j] == "+") 1:3 else if (kmerAlignments$STRAND[j] == "-") 4:6 else 1:6
				aaFrag[j] <<- DNAtoBestPeptide( kmerAlignments$Kmer[j], clip=F, readingFrame=readFrame, breakAtStops=F)
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

`mapKmerProteinFragmentsToReference` <- function( kmerAlignments, optionsFile="Options.txt", verbose=TRUE) {


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
	if (version$major == "4" && as.numeric( version$minor) >= 4) require( pwalign)
	data(BLOSUM62)

	geneMap <- getCurrentGeneMap()
	cdsMap <- getCurrentCdsMap()

	LAPPLY <- base::lapply
	STRSPLIT <- base::strsplit
	WHICH <- base::which
	MIN <- base::min
	PASTE <- base::paste

	nDone <- 0
	tapply( geneHits, geneFac, function(x) {
		
		# given all the rows that belong to one gene
		i <- x[1]
		mySeqID <- kmerAlignments$SEQ_ID[i]
		myGeneID <- kmerAlignments$GENE_ID[i]
		if ( is.na(mySeqID) || is.na(myGeneID)) return()
		nDone <<- nDone + length(x)

		# pre fetch some maps
		curGmap <- subset.data.frame( geneMap, GENE_ID == myGeneID)
		curCDSmap <- subset.data.frame( cdsMap, GENE_ID == myGeneID)

		# fetch the reference protein
		refProtein <- gene2Fasta( myGeneID, genome.file, mode="aa")$seq[1]
		if ( !is.na(refProtein)) {
			myFrags <- kmerAlignments$KMER_FRAGMENT[x]
			pa <- pairwiseAlignment( myFrags, refProtein, type="global-local", scoreOnly=F, 
						substitutionMatrix=BLOSUM62)

			# the 'alignedSubject' function seems broke/slow.  Try to do by hand
			#refFrags <- as.character( alignedSubject(pa))
			refFrags <- as.character( subject(pa))
			# when the pattern was not all used, append gap to start of subject
			refStartOffsets <- start( pattern(pa))
			toFix <- which( refStartOffsets > 1)
			if ( NFIX <- length(toFix)) {
				nGaps <- refStartOffsets[toFix] - 1
				gapstring <- "------------------------------------------------------------------"
				refFrags[toFix] <- paste( substr( rep.int(gapstring,NFIX), 1, nGaps), refFrags[toFix], sep="")
			}

			# try to tell how the Kmer differs from the reference
			# put a dot notation to show just what's different
			refAAvec <- STRSPLIT( refFrags, split="")
			myAAvec <- STRSPLIT( myFrags, split="")
			LAPPLY( 1:length(myFrags), function(i) {
				chUse <- 1 : MIN( length(refAAvec[[i]]), length(myAAvec[[i]]))
				same <- WHICH( refAAvec[[i]][chUse] == myAAvec[[i]][chUse])
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

`kmerizeOneFastqFile` <- function( filein, kmer.size=33, buffer.size=1000000, 
				maxReads=NULL, min.count=2, kmer.debug=NULL, trim5=0, trim3=0) {

	# get ready to read the fastq file in chuncks...
	filein <- allowCompressedFileName( filein)
	if ( ! file.exists( filein)) stop( paste("Can't find input file: ", filein))
	conIn <- openCompressedFile( filein, open="r")
	on.exit( close( conIn))
	require(Biostrings)
	if (version$major == "4" && as.numeric( version$minor) >= 4) require( pwalign)

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

		# allow base trimming before we kmerize
		if ( trim5 > 0 || trim3 > 0) {
			cat( "  Trim..")
			firstBP <- trim5 + 1
			lastBP <- max( nchar( seqTxt)) - trim3
			seqTxt <- substr( seqTxt, firstBP, lastBP)
			drops <- which( nchar(seqTxt) < kmer.size)
			if ( length(drops)) seqTxt <- seqTxt[ -drops]
		}

		timeStat <- round( system.time( kmerizeOneChunk.as.DNAString( seqTxt, kmer.size=kmer.size)), digits=1)
		cat( "  Time(usr,sys)=",timeStat[3],"(",timeStat[1],",",timeStat[2],")", sep="")
		rm( seqTxt)

		# each call just extends the global lists, just clean up memory
		gc()
	}
	
	if ( ! is.null(kmer.debug)) {
		cat( "\n\nDebug Kmer Counts:  Before Merge:\n")
		debugBigKmerTable( kmer.debug)
	}

	# with the file done, merge all the chunk results
	# any K-mers with too few observations are not to be kept
	mergeKmerChunks( min.count)
	if ( ! is.null(kmer.debug)) {
		cat( "\n\nDebug Kmer Counts:  After Merge, before RevComp Combine:\n")
		debugBigKmerTable( kmer.debug)
	}

	# also do any RevComp resolving now too
	if ( sum(bigKmerCounts[[1]]) > 0) {
		cat( "\nSearch for RevComp pairs to join..")
		kmerRevComp.GlobalTable()
		if ( ! is.null(kmer.debug)) {
			cat( "\n\nDebug Kmer Counts:  After RevComp Combine:\n")
			debugBigKmerTable( kmer.debug)
		}
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


`kmerizeOneChunk.as.DNAString` <- function( seqs, kmer.size=33) {

	# make Kmers of every raw read in this chunk
	# set up to get all N-mers for all the sequences in this chunk
	sizeM1 <- kmer.size - 1
	nSeq <- length( seqs)
	seqLens <- base::nchar( seqs)

	# ignore any with N for building the Kmers
	# long reads may have more N's, must do this after chopping, not before
	#hasN <- grep( "N", seqs, fixed=T)
	
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

		# x <- setdiff( x, hasN)
		
		base::lapply( x, function(j) {
			kmers <- subseq( rep.int( seqs[j], nSubstrs), fromSet, toSet)
			#if ( j %in% hasN) kmers <- kmers[ -grep("N",kmers,fixed=T)]
			hasN <- grep( "N", kmers, fixed=T)
			if ( length(hasN)) kmers <- kmers[ -hasN]
			nOut <<- nOut + 1
			kmerSets[[nOut]] <<- kmers
			return(NULL)
		})
		return(NULL)
	})
	rm( seqLens, lenFac)

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


`kmerRevComp.GlobalTable` <- function() {

	# the Kmers may be from both strands, and we only want to keep one form of each
	# so merge those where both are present
	require(Biostrings)
	if (version$major == "4" && as.numeric( version$minor) >= 4) require( pwalign)
	
	# faster memory efficient to operate on one global vector
	# stored as a list of DNASTrings and Counts
	kmers <- bigKmerStrings[[1]]
	cnts <- bigKmerCounts[[1]]
	N <- length(kmers)
	cat( "\nN_Kmers in: ", N)

	# see what the RevComp of every Kmer is.  Now with DNAStrings, it is easy
	rcKmer <- reverseComplement(kmers)
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


`kmerAlphaRevCompOrder` <- function( kmerFile) {

	# the Kmers within any one file may be from either strand.  We have already forced reduction of 
	# if both a Kmer and its RevComp were present within a single file.
	#
	# BUT: before we can combine/compare between files, we need to make sure all files are using just
	# the same set of Kmers.  Solution:  keep only the alphabetically first for every Kmer
	
	# allow being given a saved file
	loadSavedKmers( kmerFile)
	myKmers <- bigKmerStrings[[1]]
	N <- length( myKmers)
	cat( "\nForcing Kmers into RevComp alphabetical order.  N_Kmers in: ", N)

	# see what the RevComp of every Kmer is.  Now with DNAStrings, it is easy
	rcKmer <- reverseComplement( myKmers)
	
	# we will keep the 'lower' one by sort order
	needRC <- which( rcKmer < myKmers)
	if (length(needRC)) {
		cat( "\nConverting to RevComp Kmer: ", length(needRC))
		myKmers[needRC] <- rcKmer[needRC]
		cat( "  re-save Kmers file..")
		bigKmerStrings[[1]] <<- myKmers
		saveKmers( kmerFile)
	} else {
		cat( "  All Kmers already in alphabetic Kmer order..")
	}
	rm( myKmers, rcKmer, needRC)

	# now do the global cleanup...
	if ( exists("bigKmerStrings", envir=.GlobalEnv)) rm( bigKmerStrings, envir=.GlobalEnv)
	if ( exists("bigKmerCounts", envir=.GlobalEnv)) rm( bigKmerCounts, envir=.GlobalEnv)
	if ( exists("xref", envir=.GlobalEnv)) rm( xref, envir=.GlobalEnv)
	gc()
	
	return( kmerFile)
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
		uniq <- rep.int(TRUE,nNow)
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


`saveKmers` <- function(outfile) {

	# push the global results to disk
	# turn the DNAString data to character
	require(Biostrings)
	cat( "  saving..")
	bigKmerCounts <- as.vector( bigKmerCounts[[1]])

	# the first time we make Kmers from larger DNAStrings, there seems to be excess memory use.
	# so convert them to character and then back to DNAStrings
	bigKmerStrings <- as.character( bigKmerStrings[[1]])
	bigKmerStrings <- DNAStringSet( bigKmerStrings)
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
	filesToDo <- expectedFiles
	if ( all( file.exists( expectedFiles))) {
		cat( "\nUsing existing 'cutadapt' results..")
	} else {
		# is the data paired end or not
		if (asMatePairs || (!is.null(forceMatePairs) && forceMatePairs && length(filesToDo) == 2)) {
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


`plotSignificantKmers` <- function( kmerTbl, speciesID=getCurrentSpecies(), nKmersPerGene=100,
					nGenesToLabel=10, gene.cex=0.8, gene.pos=3) {

	neededColumns <- c( "SEQ_ID", "POSITION", "GENE_ID", "Log2.Fold", "P.Value")
	if ( ! all( neededColumns %in% colnames(kmerTbl))) {
		cat( "\nError: Kmer table is missing some of these needed columns: ", neededColumns)
		return( NULL)
	}
	NK <- nrow(kmerTbl)

	if ( speciesID != getCurrentSpecies()) setCurrentSpecies( speciesID)
	gmap <- getCurrentGeneMap()

	# force the data to be numeric
	kmerTbl$P.Value <- as.numeric( kmerTbl$P.Value)
	kmerTbl$Log2.Fold <- as.numeric( kmerTbl$Log2.Fold)

	# only use the Kmers that did show some change
	#cat( "\nDropping Kmers with minimal difference.  cut.fold=", cut.fold, "  cut.pvalue=", cut.pvalue)
	#keepFC <- which( abs( kmerTbl$Log2.Fold) >= cut.fold)
	#keepPV <- which( kmerTbl$P.Value <= cut.pvalue)
	#keep <- sort( union( keepFC, keepPV))
	#cat( "\nDropping: ", NK-length(keep), "   Keeping: ", length(keep))
	#kmerTbl <- kmerTbl[ keep, ]
	#NK <- nrow(kmerTbl)

	# get the names of the two comparison groups. Take it from the 'N.Samples' columns
	groupColumns <- grep( "N.Samples", colnames(kmerTbl), value=T)
	groupColumns <- sub( "^N.Samples.", "", groupColumns)
	grp1Name <- groupColumns[1]
	grp2Name <- groupColumns[2]

	# we want to visit every gene in chromosomal order, so make a key that preserves that
	# use the position of the gene's start in the annotation, not the position of the Kmer
	cat( "\nFinding Kmers per Gene..")
	myKmerGenePos <- gmap$POSITION[ match( kmerTbl$GENE_ID, gmap$GENE_ID)]
	myKmerGenePos <- formatC( myKmerGenePos, format="d", width=12)
	kmerKey <- paste( kmerTbl$SEQ_ID, myKmerGenePos, shortGeneName(kmerTbl$GENE_ID,keep=1), sep="::")
	kmerFac <- factor( kmerKey)
	NG <- nlevels( kmerFac)
	cat( "   N_Genes=", NG, "\n")

	# for each gene, we will gather a few metrics that summarize how big the fold and strong the P-value
	# trap any zero P-values at something small but not zero
	smallP <- min( kmerTbl$P.Value[ kmerTbl$P.Value > 0])
	if (smallP < 1e-300) smallP <- 1e-300
	kmerTbl$P.Value[ kmerTbl$P.Value < smallP] <- smallP
	gUpFold <- gDownFold <- rep.int( 0, NG)
	gUpPval <- gDownPval <- rep.int( 1, NG)
	gSID <- gName <- rep.int( "", NG)
	nOut <- 0
	tapply( 1:NK, kmerFac, function(x) {
		# given all the rows from one gene
		myGene <- kmerTbl$GENE_ID[ x[1]]
		if ( is.na(myGene) || myGene == "") return(NULL)
		nOut <<- nOut + 1
		gName[nOut] <<- shortGeneName(myGene,keep=1)
		gSID[nOut] <<- kmerTbl$SEQ_ID[x[1]]

		myFC <- kmerTbl$Log2.Fold[x]
		myPV <- kmerTbl$P.Value[x]
		isUP <- which( myFC > 0)
		isDOWN <- which( myFC < 0)
		# use our traditional way of ranking by both Fold and Pvalue
		ord <- diffExpressRankOrder( myFC, myPV)
		# use the number of Kmers asked for, for each group
		nUp <- min( length(isUP), nKmersPerGene)
		nDown <- min( length(isDOWN), nKmersPerGene)
		isUP <- ord[ 1:nUp]
		isDOWN <- rev(ord)[ 1:nDown]
		if ( length(isUP)) {
			gUpFold[nOut] <<- mean( myFC[ isUP], na.rm=T)
			gUpPval[nOut] <<- logmean( myPV[ isUP], na.rm=T)
		}
		if ( length(isDOWN)) {
			gDownFold[nOut] <<- mean( myFC[ isDOWN], na.rm=T)
			gDownPval[nOut] <<- logmean( myPV[ isDOWN], na.rm=T)
		}
		if ( nOut %% 100 == 0) cat( "\rDebug: ", nOut, myGene, gUpFold[nOut], gUpPval[nOut], gDownFold[nOut], gDownPval[nOut])
		return(NULL)
	})
	length(gUpFold) <- length(gDownFold) <- nOut
	length(gUpPval) <- length(gDownPval) <- nOut
	length(gSID) <- length(gName) <- nOut
	NG <- nOut

	# now with all the values known, set up to plot and color each gene
	log10pvUp <- -log10( gUpPval)
	log10pvDown <- log10( gDownPval)
	#log10pvUp[ is.nan(log10pvUp) | is.infinite(log10pvUp)] <- NA
	#log10pvDown[ is.nan(log10pvDown) | is.infinite(log10pvDown)] <- NA
	bigFC <- max( gUpFold, abs(gDownFold), na.rm=T)
	bigPV <- max( log10pvUp, abs(log10pvDown), na.rm=T)

	colorRamp <- heatMapColors( 41, palette="red-white-blue", inflex=0.5, rampExponent=1.0)
	myColorUp <- colorRamp[ 21 + round( log10pvUp*20/bigPV)]
	myColorDown <- colorRamp[ 21 + round( log10pvDown*20/bigPV)]

	# add some space for legends
	bigPV <- bigPV * 1.15
	plot( 1,1, type="n", main="Manhattan plot:  Genes with Most Significant Kmers", 
		xlab="Genes and Intergenic Gaps, in chromosomal order",
		ylab=paste( "Log10 of EdgeR P-Values "), xlim=c(1,NG), ylim=c(-bigPV,bigPV))
	
	ord <- order( log10pvUp)
	points( (1:NG)[ord], log10pvUp[ord], pch=19, col=myColorUp[ord], cex=1)
	ord <- order( log10pvDown, decreasing=T)
	points( (1:NG)[ord], log10pvDown[ord], pch=19, col=myColorDown[ord], cex=1)
	lines( c(-100,NG+100), c(0,0), lty=1, col=1, lwd=1)

	whoShow <- order( log10pvUp, decreasing=T)[1:nGenesToLabel]
	if ( length(whoShow)) text( whoShow, log10pvUp[whoShow], gName[whoShow], cex=gene.cex, col=1, pos=gene.pos, offset=0.35)
	whoShow <- order( log10pvDown, decreasing=F)[1:nGenesToLabel]
	if ( length(whoShow)) text( whoShow, log10pvDown[whoShow], gName[whoShow], cex=gene.cex, col=1, pos=4-gene.pos, offset=0.35)

	legend( 'topleft', paste( "Kmers UP in", grp2Name), pch=19, pt.cex=1.25, col='red', bty="n", cex=0.95)
	legend( 'bottomleft', paste( "Kmers UP in", grp1Name), pch=19, pt.cex=1.25, col='blue', bty="n", cex=0.95)

	# send back the summary by gene
	out <- data.frame( "GENE_ID"=gName, "PRODUCT"=gene2Product(gName), "LifeStage"=gene2CellType(gName), 
			"Avg.Up.Log10.Pvalue"=log10pvUp, "Avg.Down.Log10.Pvalue"=log10pvDown, stringsAsFactors=F)
	return( invisible( out))
}


`plotSignificantKmersOneGene` <- function( kmerTbl, geneID, speciesID=getCurrentSpecies(), domainMap=NULL, ...) {

	neededColumns <- c( "SEQ_ID", "POSITION", "GENE_ID", "Log2.Fold", "P.Value", "DIF_FROM_REF")
	if ( ! all( neededColumns %in% colnames(kmerTbl))) {
		cat( "\nError: Kmer table is missing some of these needed columns: ", neededColumns)
		return( NULL)
	}

	# make sure the table only has this one gene
	keep <- which( kmerTbl$GENE_ID == geneID)
	kmerTbl <- kmerTbl[ keep, ]
	# force the data into protein and chromosomal order
	ord <- order( as.numeric(kmerTbl$AA_POSITION), as.numeric(kmerTbl$POSITION))
	kmerTbl <- kmerTbl[ ord, ]
	NK <- nrow(kmerTbl)
	if ( ! NK) {
		cat( "\nError: no Kmers for gene: ", geneID)
		return( NULL)
	}
	NAA <- max( as.numeric( kmerTbl$AA_POSITION), na.rm=T)

	if ( speciesID != getCurrentSpecies()) setCurrentSpecies( speciesID)

	# force the data to be numeric
	kmerTbl$P.Value <- as.numeric( kmerTbl$P.Value)
	kmerTbl$Log2.Fold <- as.numeric( kmerTbl$Log2.Fold)

	# get the names of the two comparison groups. Take it from the 'N.Samples' columns
	groupColumns <- grep( "N.Samples", colnames(kmerTbl), value=T)
	groupColumns <- sub( "^N.Samples.", "", groupColumns)
	grp1Name <- groupColumns[1]
	grp2Name <- groupColumns[2]

	# for each AA, we will gather a few metrics that summarize how big the fold and strong the P-value
	# trap any zero P-values at something small but not zero
	smallP <- min( kmerTbl$P.Value[ kmerTbl$P.Value > 0])
	if (smallP < 1e-100) smallP <- 1e-100
	kmerTbl$P.Value[ kmerTbl$P.Value < smallP] <- smallP
	aaUpFold <- aaDownFold <- aaUpPval <- aaDownPval <- rep.int( NA, NK)
	kmerAAloc <- rep.int( NA, NK)
	tapply( 1:nrow(kmerTbl), factor(kmerTbl$AA_POSITION), function(x) {
		# given all the kmers that map to one AA
		myAA <- as.numeric( kmerTbl$AA_POSITION[x[1]])
		isUP <- x[ which( kmerTbl$Log2.Fold[x] > 0)]
		isDOWN <- x[ which( kmerTbl$Log2.Fold[x] < 0)]
		aaUpFold[myAA] <<- if (length(isUP)) max( kmerTbl$Log2.Fold[isUP], na.rm=T) else NA
		aaDownFold[myAA] <<- if (length(isDOWN)) min( kmerTbl$Log2.Fold[isDOWN], na.rm=T) else NA
		aaUpPval[myAA] <<- if (length(isUP)) -log10( min( kmerTbl$P.Value[isUP], na.rm=T)) else NA
		aaDownPval[myAA] <<- if (length(isDOWN)) -log10( min( kmerTbl$P.Value[isDOWN], na.rm=T)) else NA
		# stagger the X limits based on how many there are
		kmerAAloc[x] <<- seq( myAA, myAA+0.9, length.out=length(x))
	})

	# some Kmer metrics are direction specific, others are not
	isUP <- which( kmerTbl$Log2.Fold > 0)
	isDOWN <- which( kmerTbl$Log2.Fold < 0)
	bigFC <- max( aaUpFold, abs(aaDownFold), na.rm=T)
	kmerYloc <- log10( -log10( kmerTbl$P.Value))
	bigPV <- max( kmerYloc, na.rm=T)
	kmerYloc[ isDOWN] <- -(kmerYloc[isDOWN])

	# have the color show how far from NF54 reference the Kmer is, so
	# first color is NF54 match (white) and gets more red with more mutations
	colorRamp <- heatMapColors( 13, palette="red-white-blue", inflex=0.5, rampExponent=1.0)[8:13]

	# add some space for legends
	bigPVshow <- bigPV * 1.1
	mainText <- paste( "Gene Manhattan plot:  Significant Kmer Differences in: ", geneID)
	plot( 1,1, type="n", main=mainText, xlab="Kmer Location along Protein AA Sequence", xaxt="n", yaxt="n",
		ylab=paste( "Log10( Log10( EdgeR P-Value))"), xlim=c(-(NAA*0.1),NAA*1.01), ylim=c(-bigPVshow,bigPVshow))
	xTicks <- pretty( c(1,NAA))
	if (xTicks[1] < 1) xTicks[1] <- 1
	axis( side=1, at=xTicks)
	yTicks <- pretty( c(-bigPVshow,bigPVshow))
	axis( side=2, at=yTicks, abs(yTicks))
	
	# point scaling is by Fold Change
	# point color is by how different from Reference AA
	pt.cex <- (abs(kmerTbl$Log2.Fold) / bigFC) * 4
	pt.dist <- sapply( kmerTbl$DIF_FROM_REF, function(x) { chs <- strsplit(x,split="")[[1]]; return(sum(chs != "."))})
	pt.col <- colorRamp[ pmin( (pt.dist+1), length(colorRamp))]

	# draw all the points that are big enough to see
	toShow <- which( pt.cex > 0.2 & kmerTbl$P.Value <= 0.05) 
	ord <- order( -pt.dist[toShow], kmerTbl$P.Value[toShow], decreasing=T)
	#ord <- order( kmerTbl$P.Value[toShow], decreasing=T)
	toShow <- toShow[ord]
	points( kmerAAloc[toShow], kmerYloc[toShow], pch=21, bg=pt.col[toShow], cex=pt.cex[toShow], lwd=0.5)

	# if given a domain map show it
	ylo <- -bigPVshow * 0.075;  yhi <- bigPVshow * 0.075
	if ( is.null( domainMap)) {
		rect( 1, ylo, NAA+1, yhi, border=1, col='grey80', lwd=2)
		text( NAA/2, (ylo+yhi)/2, geneID, cex=1)
	} else {
		rect( 1, ylo, NAA+1, yhi, border=1, col='grey80', lwd=2)
		if ( ! ("COLOR" %in% colnames(domainMap))) domainMap$COLOR <- rainbow( nrow(domainMap), end=0.85)
		for (k in 1:nrow(domainMap)) {
			myXlo <- max( 1, domainMap$REF_START[k])
			myXhi <- min( NAA, domainMap$REF_STOP[k]) + 1
			xWide <- myXhi - myXlo + 1
			myDID <- domainMap$DOMAIN_ID[k]
			myColor <- domainMap$COLOR[k]
			rect( myXlo, ylo, myXhi, yhi, border=1, col=myColor, lwd=1.5)
			if (xWide > 10) {
				if (nchar(myDID)) text( (myXlo+myXhi)/2, (ylo+yhi)/2, myDID, cex=0.75)
			} else {
				if (nchar(myDID)) text( (myXlo+myXhi)/2, yhi, myDID, cex=0.75, pos=3, offset=0.2)
			}
		}
	}
	legend( 'top', paste( "Kmers more abundant in", grp2Name), pch=1, pt.cex=1.25, col=1, bty="n", cex=0.95)
	legend( 'bottom', paste( "Kmers more abundant in", grp1Name), pch=1, pt.cex=1.25, col=1, bty="n", cex=0.95)
	legend( 'topleft', c( "Match Ref NF54", paste( c('1','2','3','4','5+'), "AA mutations")), fill=colorRamp, bty="o", cex=0.95)
	legend( 'bottomleft', paste( "Log2 FC =", 1:4), pch=21, pt.cex=1:4, bty="o", cex=0.95, title="Abundance Difference")
}


`debugBigKmerTable` <- function( kmer.debug) {

	# mine the big table for the presence of kmers
	kmers <- DNAStringSet( kmer.debug)
	rcKmers <- reverseComplement( kmers)
	targets <- c( kmers, rcKmers)

	Nbins <- length( bigKmerStrings)
	for ( i in 1:Nbins) {
		thisStrs <- bigKmerStrings[[i]]
		thisCnts <- bigKmerCounts[[i]]
		wh <- match( targets, thisStrs, nomatch=0)
		myCounts <- thisCnts[wh]
		if ( length( myCounts)) {
			sml <- data.frame( "Kmer"=as.character(targets[wh > 0]), "Count"=myCounts, stringsAsFactors=F)
			cat( "\nKmer Debug hits in Chunk: ", i, "\n")
			print( sml)
		} else {
			cat( "\nKmer Debug: No hits in Chunk: ", i)
		}
	}
}


`kmerGrep` <- function( kmerTbl, dnaPattern, fixed=TRUE, value=FALSE) {

	# try to find Kmer substrings quickly
	neededColumns <- c( "Kmer")
	if ( ! all( neededColumns %in% colnames(kmerTbl))) {
		cat( "\nError: Kmer table is missing some needed columns: ", neededColumns)
		return( NULL)
	}
	
	rcPattern <- myReverseComplement( dnaPattern)
	
	ans1 <- grep( dnaPattern, kmerTbl$Kmer, fixed=fixed)
	ans2 <- grep( rcPattern, kmerTbl$Kmer, fixed=fixed)
	
	out <- unique( sort( c( ans1, ans2)))
	if (value) out <- kmerTbl$Kmer[out]
	out
}
