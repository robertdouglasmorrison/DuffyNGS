# pipe.SNPmatrix.R


`pipe.SNP.BaseDepth` <- function( sampleID, optionsFile="Options.txt", 
				results.path=NULL, seqIDset=NULL, SNPtablePath="~/SNPs/",
				otherSNPs=NULL, indelsToo=TRUE, keepIntergenics=TRUE) {

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
	snpfile <- file.path( snpfolder, paste( sampleID, "SNP.BaseDepth.txt", sep="."))

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


	getSNPmatrixOneSeq <- function( seqID) {

		# gmap <- subset.data.frame( geneMap, REAL_G == TRUE & SEQ_ID %in% allSeqID)
		# do not throw away the intergenics at this point, unless asked to
		gmap <- subset.data.frame( geneMap, SEQ_ID == seqID)
		if ( ! keepIntergenics) gmap <- subset( gmap, REAL_G == TRUE)
		allGenes <- gmap$GENE_ID

		# we can pre-build the storage for the answer
		snpSeq <- snpPos <- snpGene <- vector()
		nSNP <- 0

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
		if ( ! nSNP) return( NULL)

		# turn into chromosomal order
		ord <- order( snpSeq, snpPos)
		snpSeq <- snpSeq[ ord]
		snpPos <- snpPos[ ord]
		snpGene <- snpGene[ ord]

		cat( "\nMeasuring Base Depth @ Known SNPs:  ", seqID, "\tN_SNP: ", nSNP, "\n")
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
	cat( "\nWrote SNP file:  ", snpfile)
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


`pipe.SNP.FreqMatrix` <- function( sampleIDset, optionsFile="Options.txt", results.path=NULL, 
				na.rm=c('all','any','half','none'), min.freq=1.0, min.diff=NULL, min.reads=NULL,
				indelsToo=TRUE, depthToo=TRUE) {

	BASES <- c("A","C","G","T")
	N_BASES <- 4
	if (indelsToo) {
		BASES <- c("A","C","G","T","Indel")
		N_BASES <- 5
	}
	na.rm <- match.arg( na.rm)

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	curSpecies <- getCurrentSpecies()
	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	snpPath <- file.path( results.path, "VariantCalls")

	snpfiles <- file.path( snpPath, sampleIDset, paste( sampleIDset, "SNP.BaseDepth.txt", sep="."))
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
		if ( ! all( colnames(tbl) %in% BASES)) {
			cat( "\nUnexpected column names in file..  Skipping: ", basename(snpfiles[i]))
			next
		}

		# if we have but don't want InDels, remove them now
		nBases <- N_BASES
		if ( !indelsToo && ncol(tbl) > 4) {
			tbl <- tbl[ , 1:4]
			nBases <- 4
		}

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
	}

	# drop any empty rows
	if ( na.rm != "none") {
		nNA <- apply( out[ ,1:nFiles], 1, function(x) sum( is.na(x)))
		if ( na.rm == "all") drops <- which( nNA == nFiles)
		if ( na.rm == "any") drops <- which( nNA > 0)
		if ( na.rm == "half") drops <- which( nNA >= nFiles/2)
		if ( length(drops)) {
			out <- out[ -drops, , drop=FALSE]
			cat( "\nDropped rows with '",na.rm, "' NA:  ", length(drops), sep="")
		}
	}

	# drop any zero rows, or any rows with the 'too small' a frequency value
	if ( ! is.null( min.freq)) {
		rmaxs <- apply( out[ ,1:nFiles], 1, max, na.rm=T)
		drops <- which( rmaxs < min.freq)
		if ( length(drops)) {
			out <- out[ -drops, , drop=FALSE]
			cat( "\nDropped rows below 'min.freq'=", min.freq,"  N_rows:  ", length(drops))
		}
	}

	# drop any minimal difference rows
	if ( ! is.null( min.diff)) {
		rdiffs <- apply( out[ ,1:nFiles], 1, function(x) diff( range(x)))
		drops <- which( rdiffs < min.diff)
		if ( length(drops)) {
			out <- out[ -drops, , drop=FALSE]
			cat( "\nDropped rows < min.diff:   ", length(drops))
		}
	}

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


`SNP.FreqMatrix2Calls` <- function( tbl, indelsToo=TRUE) {

	BASES <- c("A","C","G","T")
	N_BASES <- 4
	if (indelsToo) {
		BASES <- c("A","C","G","T","Indel")
		N_BASES <- 5
	}

	snpID <- rownames(tbl)
	snpLoc <- sub( ":[ACGT]$", "", snpID)
	snpBase <- sub( ".+:[0-9]+:", "", snpID)

	snpFac <- factor( snpLoc)

	out <- matrix( NA, nrow=nlevels(snpFac), ncol=ncol(tbl))
	colnames(out) <- colnames(tbl)
	rownames(out) <- levels(snpFac)
	rowOut <- 0

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
			if ( rowOut %% 100 == 0) cat( "\r", rowOut, x1, snpLoc[x1])
		})

	return( out)
}


`SNP.Calls2YesNo` <- function( charM, optionsFile="Options.txt",
				fastaPath=getOptionValue( optionsFile, "genomicFastaFile")) {

	snpID <- rownames(charM)

	snpLoc <- as.integer( sub( "(.+:)([0-9]+$)", "\\2", snpID))
	snpSeq <- sub( ":.+", "", snpID)

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
				results.path=NULL, min.total.depth=15, min.allele.freq=0.15, 
				ignore.vargenes=TRUE, ignore.rifins=TRUE, ignore.introns=TRUE) {

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

	snpfile <- file.path( results.path, "VariantCalls", sampleID, paste( sampleID, "SNP.BaseDepth.txt", sep="."))
	cat( "\nReading SNP BaseDepth table..")
	snpTbl <- as.matrix( read.delim( snpfile, as.is=T))
	cat( "\nTotal SNP sites: ", nrow(snpTbl))

	moifile <- file.path( results.path, "VariantCalls", sampleID, paste( sampleID, "SNP.MOIcall.txt", sep="."))

	# drop any with NA no data
	drops <- which( is.na( snpTbl[ ,1]))
	if ( length(drops)) {
		snpTbl <- snpTbl[ -drops, ]
		cat( "\nDropped for no read data: ", length(drops))
	}

	# drop any with not enough reads
	totalReads <- apply( snpTbl, MARGIN=1, sum)
	drops <- which( totalReads < min.total.depth)
	if ( length(drops)) {
		snpTbl <- snpTbl[ -drops, ]
		totalReads <- totalReads[ -drops]
		cat( "\nDropped for insufficient depth ( <",min.total.depth,"): ", length(drops))
	}

	snpSeq <- sub( ":.+:[0-9]+$", "", rownames(snpTbl))
	snpBase <- as.integer( sub( "(^.+:)([0-9]+$)", "\\2", rownames(snpTbl)))
	snpGene <- sub( "(^.+:)(.+)(:[0-9]+$)", "\\2", rownames(snpTbl))

	# do we ignore the var genes?
	if ( ignore.vargenes) {
		vmap <- getVargeneDomainMap()
		vgenes <- sort( unique( vmap$GENE_NAME))
		drops <- which( snpGene %in% vgenes)
		if ( length(drops)) {
			snpTbl <- snpTbl[ -drops, ]
			totalReads <- totalReads[ -drops]
			snpSeq <- snpSeq[ -drops]
			snpBase <- snpBase[ -drops]
			snpGene <- snpGene[ -drops]
			cat( "\nDropped as VarGene loci:   ", length(drops))
		}
	}
	# do we ignore the rifin genes?
	if ( ignore.rifins) {
		gmap <- getCurrentGeneMap()
		rifgenes <- gmap$GENE_ID[ grep( "rifin|RIF", gmap$PRODUCT)]
		drops <- which( snpGene %in% rifgenes)
		if ( length(drops)) {
			snpTbl <- snpTbl[ -drops, ]
			totalReads <- totalReads[ -drops]
			snpSeq <- snpSeq[ -drops]
			snpBase <- snpBase[ -drops]
			snpGene <- snpGene[ -drops]
			cat( "\nDropped as Rifin loci:   ", length(drops))
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
		snpBase <- as.integer( sub( "(^.+:)([0-9]+)(:.+)", "\\2", rownames(snpTbl)))
		totalReads <- apply( snpTbl, MARGIN=1, sum)
		snpGene <- sub( "^.+:[0-9]+:", "", rownames(snpTbl))
		cat( "\nDropped as Intron loci:   ", oldN - nrow(snpTbl))
	}

	# OK, count how many alleles are deep enough
	nDeep <- alleles <- vector( length=nrow(snpTbl))
	sapply( 1:nrow(snpTbl), FUN=function(x) {
			mycut <- totalReads[x] * min.allele.freq
			hits <- which( snpTbl[ x, ] >= mycut)
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

	cat( "\nWriting MOI summary file:  ", moifile)
	write.table( moiSNPs, moifile, sep="\t", quote=F, row.names=F)
	cat( "\nN_SNPS: ", nrow(moiSNPs))

	cat( "\n\nTable of Allele Frequencies:\n")
	print( tbl <- table( moiSNPs$ALLELES))
	print( as.percent( tbl, big.value=nrow(moiSNPs), digits=2))

	return( invisible(moiSNPs))
}


`pipe.SNP.MOIscore` <- function( sampleID, optionsFile="Options.txt", 
				results.path=NULL, falseMOIsites=NULL, verbose=TRUE) {

	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	}

	f <- paste( sampleID, "SNP.MOIcall.txt", sep=".")
	f <- file.path( results.path, "VariantCalls", sampleID, f)
	if ( ! file.exists( f)) {
		cat( "\nMOI calls file not found:  ", f)
		cat( "\nCalling function 'pipe.SNP.MOIcall()'...")
		pipe.SNP.MOIcall( sampleID, optionsFile=optionsFile, results.path=results.path)
		if ( ! file.exists( f)) stop( "Call failed..  Build MOIcall results first...")
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

	score <- sum(weights) / sum( counts)
	score <- round(score,digits=2)

	return( list( "Score"=score, "Counts"=counts, "Weights"=weights))
}


`gatherFalseMOIsites` <- function( sampleIDset, optionsFile="Options.txt", results.path=NULL,
				minTimesSeen=2) {

	# grab MOI sites from samples that we know to be monoclonal, for use a a False set to ignore in field samples

	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	}

	bigDF <- data.frame()

	cat( "\nGathering Multiple Allele Sites..")
	for (s in sampleIDset) {
		f <- paste( s, "SNP.MOIcall.txt", sep=".")
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
