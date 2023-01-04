# pipe.SNPmatrix.R


# wrapper function to do all the steps we need
pipe.BuildSNP.FreqMatrix <- function( sampleIDset, outfileKeyword="AllSamples", optionsFile="Options.txt", 
					results.path=NULL, speciesID=getCurrentSpecies(), missingOnly=TRUE, 
					# args for the Variant Caller
					prob.variant=0.95, snpCallMode=c("consensus","multiallelic"), min.depth=1, 
					max.depth=10000, exonOnly=FALSE, snpOnly=FALSE,
					# args for the SNP Merge step
					min.qual=5, min.score=40, dropIndels=snpOnly,
					# args for Base Depth step
					SNPtablePath="~/SNPs/", keepIntergenics=TRUE,
					# args for the base matrix creation
					na.rm="half", min.freq=1.0, min.diff=5, min.reads=1) {

	if ( speciesID != getCurrentSpecies()) setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	}
	variants.path <- file.path( results.path, "VariantCalls")

	# part 1:  each sample should have its SNPs called.
	snpFiles <- paste( sampleIDset, prefix, "Summary.VCF.txt", sep=".")
	snpFiles <- file.path( variants.path, sampleIDset, snpFiles)
	need <- which( ! file.exists( snpFiles))
	if ( ! missingOnly) need <- 1:length(snpFiles)
	if ( length(need)) {
		cat( "\nStep 1:  Finding Variant SNP sites for", length(need), "samples..")
		snpCallMode <- match.arg( snpCallMode)
		for (s in sampleIDset[need]) pipe.VariantCalls( s, optionsFile=optionsFile, results.path=results.path, 
							speciesID=speciesID, prob.variant=prob.variant, 
							snpCallMode=snpCallMode, min.depth=min.depth, max.depth=max.depth,
							exonOnly=exonOnly, snpOnly=snpOnly)
	} else {
		cat( "\nStep 1:  Using existing 'VariantCall' summary files..")
	}

	# part 2:  combine them to find the set of SNPs seen in any sample
	# current species 'could' get altered by each step, in mixed target mode
	if ( speciesID != getCurrentSpecies()) setCurrentSpecies( speciesID)
	vcfFile <- paste( outfileKeyword, prefix, "AllVCFs.txt", sep=".")
	vcfFile <- file.path( variants.path, vcfFile)
	if ( ! file.exists( vcfFile)) {
		cat( "\n\nStep 2:  Merging all SNPs sites..")
		vcfSet <- pipe.VariantMerge( sampleIDset, outfile=vcfFile, optionsFile=optionsFile, results.path=results.path,
						speciesID=speciesID, min.qual=min.qual, min.depth=min.depth, min.score=min.score,
						dropIndels=dropIndels)
	} else {
		cat( "\nStep 2:  Using existing 'xxxx.AllVCFs.txt' SNP file..")
		vcfSet <- read.delim( vcfFile, as.is=T)
	}

	# part 3:  measure actual pileup depths and base calls at all those sites
	# current species 'could' get altered by each step, in mixed target mode
	if ( speciesID != getCurrentSpecies()) setCurrentSpecies( speciesID)
	baseDepthFiles <- paste( sampleIDset, prefix, "SNP.BaseDepth.txt", sep=".")
	baseDepthFiles <- file.path( variants.path, sampleIDset, baseDepthFiles)
	need <- which( ! file.exists( baseDepthFiles))
	if ( ! missingOnly) need <- 1:length(baseDepthFiles)
	if ( length(need)) {
		cat( "\n\nStep 3:  Measure Base Depth at all SNP sites..   Slow...")
		multicore.lapply( sampleIDset[need], pipe.SNP.BaseDepth, otherSNPs=vcfSet, indelsToo=(!dropIndels),
				optionsFile=optionsFile, results.path=results.path, keepIntergenics=keepIntergenics)
	} else {
		cat( "\nStep 3:  Using existing 'BaseDepth' files..")
	}

	# part 4:  combine into the giant table
	# current species 'could' get altered by each step, in mixed target mode
	if ( speciesID != getCurrentSpecies()) setCurrentSpecies( speciesID)
	cat( "\n\nStep 4:  Turn all Base Depths into one Frequency Matrix of all samples..")
	freqM <- pipe.SNP.FreqMatrix( sampleIDset, optionsFile=optionsFile, results.path=results.path, 
					na.rm=na.rm, min.freq=min.freq, min.diff=min.diff, min.reads=min.reads,
					exonOnly=exonOnly, indelsToo=(!dropIndels))

	# write out the results
	outfile <- paste( outfileKeyword, prefix, "BaseFreqMatrix.txt", sep=".")
	outfile <- file.path( variants.path, outfile)
	cat( "\nWriting final results: ", outfile)
	write.table( freqM, outfile, sep="\t", quote=F, row.names=T)

	#done
	return( nrow( freqM))
}


`pipe.SNP.BaseDepth` <- function( sampleID, optionsFile="Options.txt", 
				results.path=NULL, seqIDset=NULL, SNPtablePath="~/SNPs/",
				otherSNPs=NULL, indelsToo=TRUE, keepIntergenics=TRUE,
				baseDepthTableToMatch=NULL) {

	BASES <- c("A","C","G","T")
	N_BASES <- 4
	if (indelsToo) {
		BASES <- c("A","C","G","T","Indel")
		N_BASES <- 5
	}

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	curSpecies <- getCurrentSpecies()
	geneMap <- getCurrentGeneMap()
	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	fastaFile <- getOptionValue( optT, "genomicFastaFile", verbose=F)

	if ( ! is.null(SNPtablePath)) SNP_curSNPpath <<- SNPtablePath

	bamfile <- file.path( results.path, "align", paste( sampleID, "genomic.bam", sep="."))
	bamfile <- BAM.verifySorted( bamfile)

	snpfolder <- file.path( results.path, "VariantCalls", sampleID)
	if ( ! file.exists( snpfolder)) dir.create( snpfolder, recursive=T)
	snpfile <- file.path( snpfolder, paste( sampleID, prefix, "SNP.BaseDepth.txt", sep="."))

	if ( is.null( seqIDset)) seqIDset <- getCurrentSeqMap()$SEQ_ID
	allSeqID <- seqIDset

	# we can be given a data.frame of other locations to query
	otherSNPsToo <- FALSE
	if ( ! is.null( otherSNPs)) {
		if ( ! all( c("SEQ_ID","POSITION","GENE_ID") %in% colnames(otherSNPs))) {
			stop( "'otherSNPs' data.frame must have columns 'SEQ_ID','POSITION','GENE_ID'")
		}
		otherSNPsToo <- TRUE
	}

	# default behavior is to query all the SNP sites explicitly seen in this sample, with
	# all 'known' SNP sites from database object for this species.

	# new optional behavior:  if given a base depth table object, then use that as the set
	# of sites.  Assuring that the new created base depth for this sample is the same SNP
	# sites as the given one
	ForceSnpSites <- FALSE
	if ( ! is.null( baseDepthTableToMatch)) {
		rnames <- rownames( baseDepthTableToMatch)
		snpTerms <- strsplit( rnames, split="::", fixed=T)
		forceSnpSeq <- sapply( snpTerms, `[`, 1)
		forceSnpGene <- sapply( snpTerms, `[`, 2)
		forceSnpPos <- as.integer( sapply( snpTerms, `[`, 3))
		ForceSnpSites <- TRUE
	}


	# local functoin to process one SeqID
	getSNPmatrixOneSeq <- function( seqID) {

		# gmap <- subset.data.frame( geneMap, REAL_G == TRUE & SEQ_ID %in% allSeqID)
		# do not throw away the intergenics at this point, unless asked to
		gmap <- subset.data.frame( geneMap, SEQ_ID == seqID)
		if ( ! keepIntergenics) gmap <- subset( gmap, REAL_G == TRUE)
		allGenes <- gmap$GENE_ID

		# we can pre-build the storage for the answer
		snpSeq <- snpPos <- snpGene <- vector()
		nSNP <- 0

		if ( ! ForceSnpSites) {
			loadKnownSNPtable( seqID)
			if ( nrow(SNP_curSNPtable)) {
				keepers <- which( SNP_curSNPtable$GENE_ID %in% allGenes)
				snpSeq <- SNP_curSNPtable$SEQ_ID[ keepers]
				snpPos <- SNP_curSNPtable$POSITION[ keepers]
				snpGene <- SNP_curSNPtable$GENE_ID[ keepers]
				nSNP <- length( snpPos)	
				cat( "\nKnown SNPs table:    ", length(keepers))
			}
			if (otherSNPsToo) {
				smlOtherSNPs <- subset.data.frame( otherSNPs, SEQ_ID == seqID)
				if ( nrow(smlOtherSNPs)) {
					keepers <- which( smlOtherSNPs$GENE_ID %in% allGenes)
					snpSeq <- c( snpSeq, smlOtherSNPs$SEQ_ID[ keepers])
					snpPos <- c( snpPos, smlOtherSNPs$POSITION[ keepers])
					snpGene <- c( snpGene, smlOtherSNPs$GENE_ID[ keepers])
					nSNP <- length( snpPos)	
					cat( "\nOther SNPs table:    ", length(keepers))
				}
			}
			# catch any duplicates!
			drops <- which( duplicated( snpPos))
			if ( length(drops)) {
				snpSeq <- snpSeq[ -drops]
				snpPos <- snpPos[ -drops]
				snpGene <- snpGene[ -drops]
				nSNP <- length( snpPos)	
				cat( "\nDrop SNP duplicates: ", length(drops))
			}
		} else {
			use <- which( forceSnpSeq == seqID)
			snpSeq <- forceSnpSeq[use]
			snpGene <- forceSnpGene[use]
			snpPos <- forceSnpPos[use]
			nSNP <- length( snpPos)	
		}
		if ( ! nSNP) return( NULL)

		# turn into chromosomal order
		ord <- order( snpSeq, snpPos)
		snpSeq <- snpSeq[ ord]
		snpPos <- snpPos[ ord]
		snpGene <- snpGene[ ord]

		if ( ! ForceSnpSites) {
			cat( "\nMeasuring Base Depth @ Known SNP sites:  ", sampleID, "\t", seqID, "\tN_SNP: ", nSNP, "\n")
		} else {
			cat( "\nMeasuring Base Depth @ Match SNP sites:  ", sampleID, "\t", seqID, "\tN_SNP: ", nSNP, "\n")
		}
		out <- matrix( NA, nrow=N_BASES, ncol=nSNP)
		rownames(out) <- BASES
		colnames(out) <- paste( seqID, snpGene, snpPos, sep="::")

		# visit every gene, and build it up
		for ( i in 1:nrow(gmap)) {
	
			geneID <- gmap$GENE_ID[i]

			# get all the SNP sites for this gene
			who <- which( snpGene == geneID)
			# if no SNPs to a gene, skip it now
			if ( ! length(who)) next

			# grab the reads in this gene, that land in the range of these SNPs
			genePosRange <- range( snpPos[who])
			xLo <- genePosRange[1]
			xHi <- genePosRange[2]
			curMPU <- BAM.mpileup( bamfile, seqID, fastaFile, start=xLo, stop=xHi, summarize.calls=FALSE,
					verbose=FALSE)
			if ( is.null( curMPU)) return(NULL)
			if (nrow( curMPU) < 1) next

			# only keep the part at the defined SNPs
			curMPU <- subset.data.frame( curMPU, POSITION %in% snpPos[ who])
			if (nrow( curMPU) < 1) next

			# turn these to the tabular form
			calls <- MPU.callBases( curMPU$CALL_BASE, curMPU$REF_BASE)
			callStrs <- MPU.callTableToString( calls$depth.table)
			matrx <- MPU.callStringsToMatrix( callStrs)
			baseCounts <- MPU.callMatrixToBaseCounts( matrx, curMPU$REF_BASE, indelsToo=indelsToo, normalize=FALSE)

			# now we can stash this data where it goes
			theseSNPs <- match( curMPU$POSITION, snpPos)
			out[ , theseSNPs] <- t( baseCounts )
			if ( i %% 100 == 0) cat( "\r", i, geneID, substr( gmap$PRODUCT[i],1,50), 
				"\tN_SNP: ", nrow(curMPU), "   ")
		}

		# send back the bases as the column
		t( out)
	}


	# do each seqID
	cat( "\n\nSampleID:  ", sampleID)
	doSeqID <- multicore.seqID.order( allSeqID)
	if ( length( doSeqID) > 1) {
		ans <- multicore.lapply( doSeqID, getSNPmatrixOneSeq)
		names(ans) <- doSeqID

		# repackage the results
		out <- NULL
		cat( "\nCombining by chromosome..\n")
		combineOrder <- if ( all( allSeqID == doSeqID)) 1:length(allSeqID) else length(allSeqID):1
		for ( i in combineOrder) {
			thisM <- ans[[i]]
			if ( is.null( thisM)) next
			if ( ! is.matrix( thisM)) next
	
			if ( is.null( out)) {
				out <- thisM
			} else {
				out <- rbind( out, thisM)
			}
			cat( "\n", names(ans)[i], nrow(out))
		}
	} else {
		out <- getSNPmatrixOneSeq( doSeqID[1])
	}

	# write it out
	write.table( out, snpfile, sep="\t", quote=F, row.names=TRUE)
	cat( "\nWrote SNP file:  ", snpfile, "\n")
}


`exonSNPsOnly` <- function( tbl, seqID=tbl$SEQ_ID, pos=tbl$POSITION, exonMap=getCurrentExonMap()) {

	if ( nrow(tbl) < 1) return(tbl)

	mySID <- unique.default( seqID) 
	if ( length( mySID) > 1) stop( "Error in 'exonSNPsOnly':  Expected exactly one chromosome of SNP data")
	eSID <- unique.default( exonMap$SEQ_ID) 
	if ( length( eSID) > 1) {
		exonMap <- subset.data.frame( exonMap, SEQ_ID == mySID)
		if ( nrow( exonMap) < 2) return(tbl)
	}
	
	# build a lookup table of exon starts and stops
	NE2 <- nrow( exonMap) * 2
	posVec <- rep.int( 0, NE2)
	geneVec <- rep.int( "", NE2)
	typeVec <- rep.int( "", NE2)
	i <- seq.int( 1, NE2, by=2)
	posVec[ i] <- exonMap$POSITION
	posVec[ i+1] <- exonMap$END      #+ 1
	typeVec[ i] <- "E"
	typeVec[ i+1] <- "I"
	geneVec[ i] <- exonMap$GENE_ID
	geneVec[ i+1] <- exonMap$GENE_ID

	# force to increasing order for 'findInterval'
	ord <- order( posVec)
	posVec <- posVec[ ord]
	typeVec <- typeVec[ ord]
	geneVec <- geneVec[ ord]

	# find the location of all the SNPs
	snpPos <- findInterval( pos, posVec, all.inside=FALSE)

	# now call each SNP position as being in an exon or not
	isExon <- sapply( snpPos, function(x) {
		
		# outside the boundaries is instant NO
		if ( x < 1 || x >= NE2) return( FALSE)
		# in an exon is YES
		if ( typeVec[ x] == "E") return( TRUE)
		FALSE
	})

	isIntron <- ! isExon
	if ( any( isIntron)) {
		return( tbl[ isExon, ])
	} else {
		return( tbl)
	}
}


`pipe.SNP.FreqMatrix` <- function( sampleIDset, optionsFile="Options.txt", 
				speciesID=getCurrentSpecies(), results.path=NULL, 
				na.rm=c('all','any','half','none'), min.freq=1.0, min.diff=NULL, min.reads=NULL,
				indelsToo=TRUE, depthToo=TRUE, exonOnly=FALSE) {

	BASES <- c("A","C","G","T")
	N_BASES <- 4
	if (indelsToo) {
		BASES <- c("A","C","G","T","Indel")
		N_BASES <- 5
	}
	na.rm <- match.arg( na.rm)

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	if ( speciesID != getCurrentSpecies()) setCurrentSpecies( speciesID)
	curSpecies <- getCurrentSpecies()
	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	snpPath <- file.path( results.path, "VariantCalls")

	snpfiles <- file.path( snpPath, sampleIDset, paste( sampleIDset, prefix, "SNP.BaseDepth.txt", sep="."))
	finfo <- file.info( snpfiles)
	missing <- which( is.na( finfo$size))
	if ( length(missing)) {
		cat( "\nSome 'BaseDepth' files not found: ", sampleIDset[missing],  "  Dropping..")
		snpfiles <- snpfiles[ -missing]
		sampleIDset <- sampleIDset[ -missing]
	}
	nFiles <- length( snpfiles)

	# get each file, and accumulate the data
	out <- NULL

	cat( "\nLoading..")
	for ( i in 1:nFiles) {

		cat( " ", sampleIDset[i])
		# read in the file of base counts
		tbl <- read.delim( snpfiles[i], as.is=T)
		if ( ! all( BASES %in% colnames(tbl))) {
			cat( "\nUnexpected column names in file..  Skipping: ", basename(snpfiles[i]))
			next
		}

		# if we have but don't want InDels, remove them now
		nBases <- N_BASES
		if ( !indelsToo && ncol(tbl) > 4) {
			tbl <- tbl[ , 1:4, drop=FALSE]
			nBases <- 4
		}
		# if we asked for Indels, and they're not present, ignore
		if (nBases > ncol(tbl)) nBases <- ncol(tbl)

		tblFreq <- tblCnt <- as.matrix( tbl)

		# convert to frequency
		rsum <- apply( tblCnt, 1, sum, na.rm=T)
		for ( j in 1:nBases) {
			tblFreq[ ,j] <- tblCnt[ ,j] * 100 / rsum
		}
		tblFreq[ is.nan(tblFreq)] <- NA

		# if too few reads, flag it for later NA removal
		if ( ! is.null( min.reads)) {
			tooFew <- which( rsum < min.reads)
			if ( length( tooFew)) {
				tblCnt[ tooFew, ] <- NA
				tblFreq[ tooFew, ] <- NA
			}
		}

		# transform to a row for each base
		v <- as.vector( t( tblFreq))
		names(v) <- paste( rep( rownames(tblFreq), each=nBases), rep( BASES[1:nBases], times=nrow(tblFreq)), sep="::")
		nSNP <- length(v)
		depth <- rep( rsum, each=nBases)

		# stash the data, either making new storage or adding to 
		if ( is.null( out)) {
			if (depthToo) {
				out <- matrix( NA, nrow=nSNP, ncol=nFiles*2)
				colnames(out) <- paste( rep(sampleIDset,times=2), rep(c("Freq","Depth"),each=nFiles), sep="_")
			} else {
				out <- matrix( NA, nrow=nSNP, ncol=nFiles)
				colnames(out) <- paste( sampleIDset, "Freq", sep="_")
			}
			rownames(out) <- names(v)
		} else {
			if ( nSNP != nrow(out)) {
				cat( "\nSNP Data is wrong size!   File: ", basename(snpfiles[i]))
				next
			}
			if ( any( names(v) != rownames(out))) {
				bad <- which( names(v) != rownames(out))
				cat( "\nSNP names do not match!  ", length(bad), " of ", nrow(out), "\n")
				toShow <- data.frame( "Current"=rownames(out)[bad], "Now"=names(v)[bad])
				print( head( toShow), min( 12, length(bad)))
				next
			}
		}

		# data is OK, add it
		out[ , i] <- round(v,digits=2)
		if (depthToo) out[ ,i+nFiles] <- depth
	}

	# we have it all

	# test for reasons to drop data, but base the tests on just the frequency data
	cat( "\nInitial Table of SNP site frequencies:   N_Samples: ", nFiles, "  N_Rows: ", nrow(out))
	# drop any empty columns
	drops <- vector()
	for ( i in 1:nFiles) if ( all( is.na( out[ ,i]))) drops <- c( drops, i)
	if ( length(drops)) {
		dropFreqs <- drops
		if (depthToo) {
			dropDepths <- (drops+nFiles)
			drops <- c( dropFreqs, dropDepths)
		}
		out <- out[ , -drops, drop=FALSE]
		nFiles <- nFiles - length(dropFreqs)
		cat( "\nDropped entire samples for missing/invalid data: ", length(dropFreqs))
	}

	# drop any empty rows
	if ( nrow(out) && na.rm != "none") {
		nNA <- apply( out[ ,1:nFiles, drop=FALSE], 1, function(x) sum( is.na(x)))
		if ( na.rm == "all") drops <- which( nNA == nFiles)
		if ( na.rm == "any") drops <- which( nNA > 0)
		if ( na.rm == "half") drops <- which( nNA >= nFiles/2)
		if ( length(drops)) {
			out <- out[ -drops, , drop=FALSE]
			cat( "\nDropped rows with at least '",na.rm, "' NA:  ", length(drops), sep="")
		}
	}

	# drop any zero rows, or any rows with the 'too small' a frequency value
	if ( nrow(out) && ! is.null( min.freq)) {
		rmaxs <- apply( out[ ,1:nFiles, drop=FALSE], 1, max, na.rm=T)
		drops <- which( rmaxs < min.freq)
		if ( length(drops)) {
			out <- out[ -drops, , drop=FALSE]
			cat( "\nDropped rows with no value above 'min.freq' =", min.freq,"   N:  ", length(drops))
		}
	}

	# drop any minimal difference rows
	if ( nrow(out) && ! is.null( min.diff)) {
		rdiffs <- apply( out[ ,1:nFiles, drop=FALSE], 1, function(x) diff( range(x)))
		drops <- which( rdiffs < min.diff)
		if ( length(drops)) {
			out <- out[ -drops, , drop=FALSE]
			cat( "\nDropped rows with no values differing by at least 'min.diff' =", 
					min.diff, "   N: ", length(drops))
		}
	}
	
	# drop any SNPs that land outside Exons?
	if ( nrow(out) && exonOnly) {
		# exon search is one chromosome at a time, so do it in chuncks
		cat( "\nAbout to remove non-exon SNPs..")
		nameTerms <- strsplit( rownames(out), split="::")
		seqVec <- sapply( nameTerms, FUN=`[`, 1)
		posVec <- as.numeric( sapply( nameTerms, FUN=`[`, 3))
		ans <- tapply( 1:nrow(out), factor(seqVec), function(x) {
						smlDF <- out[ x, , drop=F]
						smlSeq <- seqVec[x]
						smlPos <- posVec[x]
						return( exonSNPsOnly( smlDF, seqID=smlSeq, pos=smlPos))
					})
		newOut <- data.frame()
		for ( i in 1:length(ans)) {
			smlAns <- ans[[i]]
			newOut <- rbind( newOut, smlAns)
		}
		cat( "  Done.  N_Dropped =", nrow(out) - nrow(newOut))
		out <- newOut
	}
	
	cat( "\nFinal Table of SNP site frequencies:   N_Samples: ", nFiles, "  N_Rows: ", nrow(out), "\n")

	# we are really done
	return( out)
}


`freqMatrix.combineReplicates` <- function( m, ids ) {

	# given a Freq Matrix with multiple columns per biological sample, 
	# reduce technical replicates by averaging
	m <- as.matrix(m)

	freqCols <- grep( "_Freq$", colnames(m))
	depthCols <- grep( "_Depth$", colnames(m))
	freqM <- m[ , freqCols]
	depthToo <- FALSE
	if ( length(depthCols)) {
		depthToo <- TRUE
		depthM <- m[ , depthCols]
	}

	# the set of IDs must be the right length
	if ( length(ids) != ncol(freqM)) {
		stop( paste( "Length of 'ids' must match size of Freq Matrix: ", length(ids), ncol(freqM)))
	}

	idFac <- factor( ids)
	idLevels <- levels(idFac)
	NIDS <- nlevels(idFac)
	NR <- nrow(freqM)
	NC <- ncol(freqM)

	fm <- matrix( NA, nrow=NR, ncol=NIDS)
	rownames(fm) <- rownames(freqM)
	colnames(fm) <- paste( idLevels, "Freq", sep="_")
	k <- 0
	tapply( 1:NC, idFac, function(x) {
			smlM <- freqM[ ,x]
			k <<- k + 1
			cat( "\rFreq:  ", k, "of", NIDS)
			if ( length(x) == 1) {
				fm[ ,k] <<- smlM
			} else {
				fm[ ,k] <<- apply( smlM, 1, mean, na.rm=T)
			}
		})
	#for ( k in 1:NR) fm[k, ] <- tapply( freqM[k, ], idFac, mean, na.rm=T)
	fm <- round( fm, digits=2)

	if (depthToo) {
		dm <- matrix( NA, nrow=NR, ncol=NIDS)
		rownames(dm) <- rownames(freqM)
		colnames(dm) <- paste( idLevels, "Depth", sep="_")
		k <- 0
		tapply( 1:NC, idFac, function(x) {
				smlM <- depthM[ ,x]
				k <<- k + 1
				cat( "\rDepth: ", k, "of", NIDS)
				if ( length(x) == 1) {
					dm[ ,k] <<- smlM
				} else {
					dm[ ,k] <<- apply( smlM, 1, mean, na.rm=T)
				}
			})
		dm <- round( dm)
	}

	if (depthToo) {
		out <- cbind( fm, dm)
	} else {
		out <- fm
	}
	return(out)
}


`freqMatrix.PCAplot` <- function( m, min.depth=1, min.diff=1, dropNA=TRUE, ... ) {

	# given a SNP Freq Matrix
	m <- as.matrix(m)

	freqCols <- grep( "_Freq$", colnames(m))
	depthCols <- grep( "_Depth$", colnames(m))
	freqM <- m[ , freqCols]
	colnames(freqM) <- sub( "_Freq$", "", colnames(freqM))
	depthToo <- FALSE
	if ( length(depthCols)) {
		depthToo <- TRUE
		depthM <- m[ , depthCols]
		colnames(depthM) <- sub( "_Depth$", "", colnames(depthM))
	}

	# use the depth and diff and NA terms to reduce rows before PCA.
	# to narrow down to the most relevant features
	if ( dropNA) {
		# find all rows that have missing data
		nNAf <- apply( freqM, 1, function(x) sum( is.na(x)))
		dropF <- which( nNAf > 0)
		if (depthToo) {
			nNAd <- apply( depthM, 1, function(x) sum( is.na(x)))
			dropD <- which( nNAd > 0)
			dropF <- sort( union( dropF, dropD))
		}
		if ( length( dropF)) {
			freqM <- freqM[ -dropF, ]
			depthM <- depthM[ -dropF, ]
		}
	}
	if ( depthToo) {
		minDeep <- apply( depthM, 1, min, na.rm=T)
		dropD <- which( minDeep < min.depth)
		if ( length( dropD)) {
			freqM <- freqM[ -dropD, ]
			depthM <- depthM[ -dropD, ]
		}
	}
	minDiff <- apply( freqM, 1, function(x) diff( range( x, na.rm=T)))
	dropD <- which( minDiff < min.diff)
	if ( length( dropD)) {
		freqM <- freqM[ -dropD, ]
		depthM <- depthM[ -dropD, ]
	}

	# OK, ready to PCA and draw
	ans <- matrix.PCAplot( freqM, ...)
	return( ans)
}


`freqMatrix.PhyloTree` <- function( m, min.depth=1, min.diff=1, dropNA=TRUE, 
					keepIntergenics=TRUE, keepIndels=TRUE, ... ) {

	# given a SNP Freq Matrix
	m <- as.matrix(m)

	freqCols <- grep( "_Freq$", colnames(m))
	depthCols <- grep( "_Depth$", colnames(m))
	freqM <- m[ , freqCols]
	colnames(freqM) <- sub( "_Freq$", "", colnames(freqM))
	depthToo <- FALSE
	if ( length(depthCols)) {
		depthToo <- TRUE
		depthM <- m[ , depthCols]
		colnames(depthM) <- sub( "_Depth$", "", colnames(depthM))
	}

	# use the depth and diff and NA terms to reduce rows before PCA.
	# to narrow down to the most relevant features
	if ( dropNA) {
		# find all rows that have missing data
		nNAf <- apply( freqM, 1, function(x) sum( is.na(x)))
		dropF <- which( nNAf > 0)
		if (depthToo) {
			nNAd <- apply( depthM, 1, function(x) sum( is.na(x)))
			dropD <- which( nNAd > 0)
			dropF <- sort( union( dropF, dropD))
		}
		if ( length( dropF)) {
			freqM <- freqM[ -dropF, ]
			depthM <- depthM[ -dropF, ]
		}
	}
	if ( depthToo) {
		minDeep <- apply( depthM, 1, min, na.rm=T)
		dropD <- which( minDeep < min.depth)
		if ( length( dropD)) {
			freqM <- freqM[ -dropD, ]
			depthM <- depthM[ -dropD, ]
		}
	}
	minDiff <- apply( freqM, 1, function(x) diff( range( x, na.rm=T)))
	dropD <- which( minDiff < min.diff)
	if ( length( dropD)) {
		freqM <- freqM[ -dropD, ]
		depthM <- depthM[ -dropD, ]
	}

	# OK, ready to see the actual base calls at those lacations
	baseCalls <- freqMatrix.BaseCalls( freqM, keepIntergenics=keepIntergenics, keepIndels=keepIndels)

	# now decide which rows to use for the final base 'string'
	nN <- apply( baseCalls, MARGIN=1, function(x) sum( x == "N"))
	drops <- which( nN > 0)
	if ( length( drops)) {
		cat( "\nDropping sites with no call in at least one sample: ", length(drops))
		baseCalls <- baseCalls[ -drops, ]
	}

	# Only keep those base calls where 'something' is different
	nCh <- apply( baseCalls, MARGIN=1, function(x) length( unique.default( x)))
	drops <- which( nCh == 1)
	if ( length( drops)) {
		cat( "\nDropping sites with all same one base called: ", length(drops))
		baseCalls <- baseCalls[ -drops, ]
	}

	# since we want each string to be exactly the same length, and any indels to be 1 char, 
	# force that now
	baseCalls <- substr( baseCalls, 1, 1)

	# good, now turn these into strings
	baseString <- apply( baseCalls, MARGIN=2, paste, collapse="")
	cat( "\nFinal Base Call String length: ", nrow(baseCalls))
	
	# finally make this into a distance matrix
	# because we have exact same length strings, just count up mismatches.
	# No need for distance tools
	cat( "\nBuilding String Distance Matrix..")
	NS <- ncol(baseCalls)
	dm <- matrix( 0, NS, NS)
	colnames(dm) <- row.names(dm) <- colnames(baseCalls)
	for ( i in 1:(NS-1)) for ( j in (i+1):NS) dm[i,j] <- dm[j,i] <- sum( baseCalls[,i] != baseCalls[,j])
	
	# and show it
	plotPhyloTree( baseString, dm=dm, ...)

	out <- list( "dm"=dm, "base.calls"=baseCalls, "base.strings"=baseString)
	return( invisible(out))
}


`freqMatrix.BaseCalls` <- function( tbl, keepIntergenics= TRUE, keepIndels=TRUE) {

	# the freq matrix may have the "Freq" and "Depth" columns both
	if ( length( grep( "(Freq|Depth)$", colnames(tbl)))) {
		keep <- grep( "_Freq$", colnames(tbl))
		tbl <- tbl[ , keep]
		colnames(tbl) <- sub( "_Freq$", "", colnames(tbl))
	}

	# do we want the intergenic stuff removed?
	if ( ! keepIntergenics) {
		drops <- grep( "(ng)", rownames(tbl), fixed=T)
		if ( length(drops)) tbl <- tbl[ -drops, ]
		cat( "\nDropped intergenic SNP sites: ", length(drops))
	}

	# do we want the Indels?
	if ( ! keepIndels) {
		drops <- grep( "Indel$", rownames(tbl))
		if ( length(drops)) tbl <- tbl[ -drops, ]
		cat( "\nDropped Indel SNP calls: ", length(drops))
	}

	# extract the gene/position and the base call for each row
	snpID <- rownames(tbl)
	snpLoc <- sub( "::[ACGT]$", "", snpID)
	if (keepIndels) snpLoc <- sub( "::Indel$", "", snpLoc)
	snpBase <- sub( ".+::[0-9]+::", "", snpID)
	snpFac <- factor( snpLoc)

	# build storage for the result
	out <- matrix( "N", nrow=nlevels(snpFac), ncol=ncol(tbl))
	colnames(out) <- colnames(tbl)
	rownames(out) <- levels(snpFac)
	rowOut <- 0

	# visit each base location and see what's the top call
	cat( "\n")
	tapply( 1:nrow(tbl), snpFac, function(x) {

			rowsIn <- x
			rowOut <<- rowOut + 1
			x1 <- x[1]
			for ( j in 1:ncol(tbl)) {
				cnts <- tbl[ rowsIn, j]
				if ( all( is.na( cnts))) next
				best <- which.max( cnts)
				mybase <- snpBase[ x1 + best - 1]
				out[ rowOut, j] <<- mybase
			}
			if ( rowOut %% 1000 == 0) cat( "\r", rowOut, x1, snpLoc[x1])
		})

	return( out)
}


`SNP.Calls2YesNo` <- function( charM, optionsFile="Options.txt",
				fastaPath=getOptionValue( optionsFile, "genomicFastaFile")) {

	snpID <- rownames(charM)

	snpLoc <- as.integer( sub( "(.+::)([0-9]+$)", "\\2", snpID))
	snpSeq <- sub( "::.+", "", snpID)

	seqFac <- factor( snpSeq)

	out <- matrix( 0, nrow=nrow(charM), ncol=ncol(charM))
	colnames(out) <- colnames(charM)
	rownames(out) <- rownames(charM)

	tapply( 1:nrow(charM), seqFac, function(x) {

			seqID <- snpSeq[x[1]]
			fa <- getFastaSeqFromFilePath( fastaPath, seqID)

			myLocs <- snpLoc[ x]
			refBase <- strsplit( fa, split="")[[1]][ myLocs]

			for ( j in 1:ncol(charM)) {
				myBase <- charM[ x, j]
				isRef <- (refBase == myBase)
				out[ x[ isRef], j] <<- 1
			}
		})

	return( out)
}


`pipe.SNP.MOIcall` <- function( sampleID, optionsFile="Options.txt", 
				results.path=NULL, min.total.depth=50, min.allele.freq=0.10, 
				drop.genes=NULL, ignore.introns=TRUE, verbose=TRUE) {

	BASES <- c("A","C","G","T")
	N_BASES <- 4

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	curSpecies <- getCurrentSpecies()
	geneMap <- getCurrentGeneMap()
	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}

	snpfile <- paste( sampleID, prefix, "SNP.BaseDepth.txt", sep=".")
	snpfile <- file.path( results.path, "VariantCalls", sampleID, snpfile)
	if (verbose) cat( "\nReading SNP BaseDepth table..")
	snpTbl <- as.matrix( read.delim( snpfile, as.is=T))
	if (verbose) cat( "\nTotal SNP sites: ", nrow(snpTbl))

	moifile <- paste( sampleID, prefix, "SNP.MOIcall.txt", sep=".")
	moifile <- file.path( results.path, "VariantCalls", sampleID, moifile)

	# drop any with NA no data
	drops <- which( is.na( snpTbl[ ,1]))
	if ( length(drops)) {
		snpTbl <- snpTbl[ -drops, ]
		if (verbose) cat( "\nDropped for no read data: ", length(drops))
	}

	# drop any with not enough reads
	totalReads <- apply( snpTbl, MARGIN=1, sum)
	drops <- which( totalReads < min.total.depth)
	if ( length(drops)) {
		snpTbl <- snpTbl[ -drops, ]
		totalReads <- totalReads[ -drops]
		if (verbose) cat( "\nDropped for insufficient depth ( <",min.total.depth,"): ", length(drops))
	}

	snpSeq <- sub( ":.+:[0-9]+$", "", rownames(snpTbl))
	snpBase <- suppressWarnings( as.integer( sub( "(^.+:)([0-9]+$)", "\\2", rownames(snpTbl))))
	snpGene <- sub( "(^.+:)(.+)(:[0-9]+$)", "\\2", rownames(snpTbl))

	# do we ignore any genes?
	if ( ! is.null( drop.genes)) {
		drops <- which( snpGene %in% drop.genes)
		if ( length(drops)) {
			snpTbl <- snpTbl[ -drops, ]
			totalReads <- totalReads[ -drops]
			snpSeq <- snpSeq[ -drops]
			snpBase <- snpBase[ -drops]
			snpGene <- snpGene[ -drops]
			if (verbose) cat( "\nDropped as 'drop.genes' loci:   ", length(drops))
		}
	}

	# do we ignore the introns?
	if ( ignore.introns) {
		oldN <- nrow( snpTbl)
		emap <- getCurrentExonMap()
		newTbl <- NULL
		for (sid in sort( unique( snpSeq))) {
			who <- which( snpSeq == sid)
			if ( ! length(who)) next
			smlTbl <- snpTbl[ who, , drop=F]
			pos <- snpBase[ who]
			emp <- subset.data.frame( emap, SEQ_ID == sid)
			smlNew <- exonSNPsOnly( smlTbl, seqID=sid, pos=pos, exonMap=emp)
			if ( is.null( newTbl)) {
				newTbl <- smlNew
			} else {
				newTbl <- rbind( newTbl, smlNew)
			}
		}
		snpTbl <- newTbl
		snpSeq <- sub( ":[0-9]+:.+", "", rownames(snpTbl))
		snpBase <- suppressWarnings( as.integer( sub( "(^.+:)([0-9]+)(:.+)", "\\2", rownames(snpTbl))))
		totalReads <- apply( snpTbl, MARGIN=1, sum)
		snpGene <- sub( "^.+:[0-9]+:", "", rownames(snpTbl))
		if (verbose) cat( "\nDropped as Intron loci:   ", oldN - nrow(snpTbl))
	}

	# OK, count how many alleles are deep enough
	nDeep <- alleles <- vector( length=nrow(snpTbl))
	if ( is.null(snpTbl) || nrow(snpTbl) < 1) {
		cat( "\nNo SNPs found.  Perhaps there was errors...")
		file.delete( moifile)
		return(NULL)
	}

	sapply( 1:nrow(snpTbl), FUN=function(x) {
			mycut <- totalReads[x] * min.allele.freq
			hits <- which( snpTbl[ x, ] >= mycut)
			hits <- intersect( hits, 1:4)
			nDeep[x] <<- length(hits)
			alleles[x] <<- paste( BASES[hits], collapse="/")
			return( NULL)
		})

	freqTbl <- snpTbl
	for ( i in 1:nrow(snpTbl)) {
		freqTbl[ i, ] <- snpTbl[ i, ] / totalReads[i]
	}
	freqTbl <- round( freqTbl, digits=2)
	colnames(freqTbl) <- paste( "PCT", colnames(freqTbl), sep="_")

	moiSNPs <- data.frame( "SNP_ID"=rownames(snpTbl), "N_ALLELES"=nDeep, "ALLELES"=alleles, 
				snpTbl, freqTbl, stringsAsFactors=F)
	rownames(moiSNPs) <- 1:nrow(moiSNPs)

	if (verbose) cat( "\nWriting MOI summary file:  ", moifile)
	write.table( moiSNPs, moifile, sep="\t", quote=F, row.names=F)
	if (verbose) cat( "\nN_SNPS: ", nrow(moiSNPs))

	if (verbose) cat( "\n\nTable of Allele Frequencies:\n")
	if (verbose) print( tbl <- table( moiSNPs$ALLELES))
	if (verbose) print( as.percent( tbl, big.value=nrow(moiSNPs), digits=2))

	return( invisible(moiSNPs))
}


`pipe.SNP.MOIscore` <- function( sampleID, optionsFile="Options.txt", 
				results.path=NULL, falseMOIsites=NULL, verbose=TRUE, ...) {

	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	}
	prefix <- getCurrentSpeciesFilePrefix()

	f <- paste( sampleID, prefix, "SNP.MOIcall.txt", sep=".")
	f <- file.path( results.path, "VariantCalls", sampleID, f)
	# always rebuild the MOI calls..  not that slow..
	if (verbose) cat( "\nCalling function 'pipe.SNP.MOIcall()'...")
	pipe.SNP.MOIcall( sampleID, optionsFile=optionsFile, results.path=results.path, verbose=verbose, ...)
	if ( ! file.exists( f)) {
		cat( "\nCall failed..  Perhaps run 'pipe.SNP.BaseDepth()' first..")
		return( list( "Score"=NA, "Counts"=NA, "Weights"=NA))
	}

	moiTbl <- read.delim( f, as.is=T)
	if (verbose) {
		cat( "\nLoaded SNP MOI call table.  N_SNPs:      ", nrow( moiTbl))
		cat( "\nNumber of multi-allelic candidate sites: ", sum( moiTbl$N_ALLELES > 1))
	}

	if ( ! is.null( falseMOIsites)) {
		if ( is.character(falseMOIsites)) {
			falseMOIsites <- read.delim( falseMOIsites, as.is=T)
		}
		if (verbose) cat( "\n\nGiven a set of false MOI sites to drop from consideration. N_Row: ", nrow(falseMOIsites))
	
		drops <- which( moiTbl$SNP_ID %in% falseMOIsites$SNP_ID)
		if ( N <- length(drops)) {
			if (verbose) cat( "\nRemoving some MOI as false sites.  N: ", N)
			moiTbl <- moiTbl[ -drops, ]
		}
	}
	cat( "\n")

	# score it...
	nSingle <- sum( moiTbl$N_ALLELES == 1)
	nDouble <- sum( moiTbl$N_ALLELES == 2)
	nTriple <- sum( moiTbl$N_ALLELES == 3)
	nMore <- sum( moiTbl$N_ALLELES > 3)
	counts <- c( "1"=nSingle, "2"=nDouble, "3"=nTriple, "4"=nMore)
	weights <- c( "1"=(nSingle*1), "2"=(nDouble*10^1), "3"=(nTriple*10^2), "4"=(nMore*10^3))
	score <- round( sum(weights) / sum( counts), digits=2)

	# try making a histogram of the multi-allelic site percentages
	sml <- as.matrix( subset( moiTbl, N_ALLELES %in% 2:3)[ , c("PCT_A","PCT_C","PCT_G","PCT_T","PCT_Indel")])
	# only use the rows that do not have any Indel component, as that distorts the relativer percentages 
	sml <- sml[ sml[ , "PCT_Indel"] == 0, ]
	# remove the 0's, then turn them into 5% bins
	pctV <- as.numeric( sml) * 100
	pctV <- round( pctV/5) * 5
	pctV <- pctV[ pctV > 0]
	pctHist <- sort( table( factor( as.integer(pctV), levels=seq( 5, 95, by=5))), decreasing=T)
	names(pctHist) <- paste( names(pctHist), "%", sep="")

	# make the base call percentags table
	pctTable <- table( moiTbl$ALLELES)
	pctTable <- round( pctTable * 100 / sum(pctTable),, digits=2)
	return( list( "Score"=score, "Counts"=counts, "Weights"=weights, "Histogram"=pctHist, 
			"Percentages"=pctTable))
}


`gatherFalseMOIsites` <- function( sampleIDset, optionsFile="Options.txt", results.path=NULL,
				minTimesSeen=2) {

	# grab MOI sites from samples that we know to be monoclonal, for use a a False set to ignore in field samples

	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	}
	prefix <- getCurrentSpeciesFilePrefix()

	bigDF <- data.frame()

	cat( "\nGathering Multiple Allele Sites..")
	for (s in sampleIDset) {
		f <- paste( s, prefix, "SNP.MOIcall.txt", sep=".")
		f <- file.path( results.path, "VariantCalls", s, f)
		if ( ! file.exists( f)) {
			cat( "\nMOI calls file not found:  ", f)
			next
		}
		tbl <- read.delim( f, as.is=T)
		smlDF <- subset( tbl, N_ALLELES > 1)
		if ( (N <- nrow(smlDF)) > 0) {
			smlDF$SampleID <- s
			bigDF <- rbind( bigDF, smlDF)
		}
		cat( "\n", s, N, nrow(bigDF))
	}

	# let's require multiple sightings to be a true 'false'
	nSeen <- tapply( 1:nrow(bigDF), factor( bigDF$SNP_ID), length)

	cat( "\n\nTable of frequencies:  How often was each site called an MOI site:\n")
	print( table( nSeen))
	keepers <- which( nSeen >= minTimesSeen)
	if (length( keepers) < 1) {
		cat( "\nNo SNPs were seen at least ", minTimesSeen, " times.")
		return( NULL)
	}

	# tally up the average percentages for these 'true' false MOI sites
	mySNPs <- names(nSeen)[ keepers]
	nSNPs <- length(mySNPs)
	myNseen <- nSeen[ keepers]
	cat( "\nAveraging observed allele percentages for ", nSNPs, " keepers..")
	avgPcts <- matrix( NA, nrow=nSNPs, ncol=4)
	colnames(avgPcts) <- myColumns <- c( "PCT_A", "PCT_C", "PCT_G", "PCT_T")
	for ( i in 1:nSNPs) {
		wh <- which( bigDF$SNP_ID == mySNPs[i])
		avgs <- apply( as.matrix( bigDF[ wh, myColumns]), MARGIN=2, mean)
		avgPcts[ i, ] <- avgs
	}

	cat( "  Done.\n")
	out <- data.frame( "SNP_ID"=mySNPs, "N_Samples"= myNseen, round( avgPcts, 2), stringsAsFactors=F)
	return( out)
}

