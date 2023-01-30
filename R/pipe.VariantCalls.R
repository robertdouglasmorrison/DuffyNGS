# pipe.VariantCalls.R

`pipe.VariantCalls` <- function( sampleIDset, annotationFile="Annotation.txt",
				optionsFile="Options.txt", speciesID=getCurrentSpecies(), results.path=NULL,
				seqIDset=NULL, start=NULL, stop=NULL, geneID=NULL, prob.variant=0.95, min.qual=10,
				snpCallMode=c("consensus", "multiallelic"), min.depth=1, max.depth=10000, 
				ploidy=if (speciesID %in% MAMMAL_SPECIES) "" else "1",
				mpileupArgs="", vcfArgs="", comboSamplesName="Combined", 
				# args passed to 'Summary' step
				exonOnly=FALSE, snpOnly=FALSE, 
				verbose=TRUE) {

	# the SNP calling probability threshold is tuned for genomic DNA of uniform depth, and diploid organism.
	# Use much higher cutoff for RNA-seq with hihgly variant read depth

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	target <- getAndSetTarget( optionsFile, sampleID=sampleIDset[1], annotationFile=annotationFile)

	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	} else {
		results.path <- results.path
	}
	fastaFile <- getOptionValue( optT, "genomicFastaFile")

	finalName <- sampleIDset[1]
	if (length(sampleIDset) > 1) finalName <- comboSamplesName
	vcfPath <- file.path( results.path, "VariantCalls", finalName)
	if ( ! file.exists( vcfPath)) dir.create( vcfPath, recursive=TRUE)

	bamfilelist <- vector()
	for ( sampleID in sampleIDset) {
		# make sure we have the BAM file already sorted
		bamfile <- paste( sampleID, "genomic.bam", sep=".")
		bamfile <- file.path( results.path, "align", bamfile)
		sortedbamfile <- BAM.verifySorted( bamfile, index=TRUE)
		bamfilelist <- c( bamfilelist, sortedbamfile)
	}

	if ( is.null( speciesID)) {
		speciesSet <- getCurrentTargetSpecies()
	} else { 
		speciesSet <- speciesID
	}
	allIDs <- allCounts <- vector()
	snpCallMode <- match.arg( snpCallMode)

	
	`variantCallOneSeq` <- function( sid, ploidy="1") {
		ans <- BAM.variantCalls( bamfilelist, seqID=sid, fastaFile=fastaFile, 
					start=start, stop=stop, min.depth=min.depth, max.depth=max.depth, 
					prob.variant=prob.variant, min.qual=min.qual, mpileupArgs=mpileupArgs, 
					vcfArgs=vcfArgs, ploidy=ploidy, snpCallMode=snpCallMode,
					verbose=verbose)
		cat( "\n", sid, "\tN_Variants: ", nrow(ans),"\n")
		N <- nrow(ans)
		outfile <- paste( finalName, sid, "VCF.txt", sep=".")
		if ( ! is.null(geneID)) {
			outfile <- paste( finalName, geneID, "VCF.txt", sep=".")
		}
		outfile <- file.path( vcfPath, outfile)
		file.delete( outfile)
		if ( nrow(ans) > 0) write.table( ans, outfile, sep="\t", quote=FALSE, row.names=FALSE)
		return( N)
	}
	

	for ( speciesID in speciesSet) {
		setCurrentSpecies( speciesID)
		seqMap <- getCurrentSeqMap()
		if ( length(speciesSet) > 1) {
			ploidy <- if ( speciesID %in% MAMMAL_SPECIES) "" else "1"
		}

		# order to speed up the parallel computation
		seqMap <- seqMap[ order( seqMap$LENGTH, decreasing=TRUE), ]
		seqIDs <- seqMap$SEQ_ID
		if ( ! is.null( seqIDset)) seqIDs <- intersect( seqIDs, seqIDset)

		# allow to be given a single gene name
		if ( ! is.null( geneID)) {
			geneMap <- subset( getCurrentGeneMap(), GENE_ID == geneID)
			if ( nrow(geneMap)) {
				seqIDs <- geneMap$SEQ_ID[1]
				start <- geneMap$POSITION[1]
				stop <- geneMap$END[1]
			} else {
				seqIDs <- character(0)
				start <- stop <- NULL
			}
		}

		N <- length(seqIDs)
		if ( N < 1) next
		vCounts <- rep.int( 0, N)

		cat( "\n")
		ans <- multicore.lapply( seqIDs, FUN=variantCallOneSeq, ploidy=ploidy)
		vCounts <- as.integer( unlist( ans))
		
		allIDs <- c( allIDs, seqIDs)
		allCounts <- c( allCounts, vCounts)

		# if all the chromosomes were done, go ahead and summarize too
		if ( length( sampleIDset) == 1 && is.null( seqIDset) && is.null(geneID)) {
			pipe.VariantSummary( sampleIDset[1], speciesID, annotationFile=annotationFile,
					optionsFile=optionsFile, results.path=results.path,
					exonOnly=exonOnly, snpOnly=snpOnly)
		} else if ( length( sampleIDset) == 1 && ! is.null( geneID)) {
			out <- pipe.VariantSummary( sampleIDset[1], speciesID, annotationFile=annotationFile,
					optionsFile=optionsFile, results.path=results.path, seqIDset=seqIDs,
					exonOnly=exonOnly, snpOnly=snpOnly, geneID=geneID)
			return( invisible(out))
		}
	}
	if ( length( allIDs) < 1) return( NULL)

	out <- data.frame( "SEQ_ID"=allIDs, "N_VARIANTS"=allCounts, stringsAsFactors=F)
	return( out)
}


pipe.VariantSummary <- function( sampleID, speciesID=getCurrentSpecies(), annotationFile="Annotation.txt", 
				optionsFile="Options.txt", results.path=NULL, seqIDset=NULL, geneID=NULL,
				min.depth=1, min.score=5, exonOnly=FALSE, snpOnly=FALSE) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	vcfPath <- file.path( results.path, "VariantCalls", sampleID)
	genomeFastaFile <- getOptionValue( optT, "genomicFastaFile")

	setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.null( seqIDset)) seqIDset <- getCurrentSeqMap()$SEQ_ID
	exonMap <- getCurrentExonMap()

	out <- data.frame()

	for ( sid in seqIDset) {
		infile <- paste( sampleID, sid, "VCF.txt", sep=".")
		# allow that we may have called for a single gene
		if ( ! is.null( geneID)) {
			infile <- paste( sampleID, geneID, "VCF.txt", sep=".")
		}
		infile <- file.path( vcfPath, infile)
		if ( ! file.exists(infile)) {
			cat( "\rVariant Calls file not found:  ", infile, "  Skip..\n")
			next
		}
		tbl <- read.delim( infile, as.is=T)
		if ( ! all( tbl$FORMAT %in% c( "GT:PL:DP", "GT:PL"))) stop( "Required VCF format of GT:PL:DP not present...")
		tbl <- tbl[ , -match( c("FILTER","FORMAT"), colnames(tbl))]
		colnames(tbl)[ ncol(tbl)] <- "GENOTYPE_CALL"

		if (exonOnly) {
			emap <- subset.data.frame( exonMap, SEQ_ID == sid)
			tbl <- exonVariantsOnly( tbl, exonMap=emap)
		}
		if (snpOnly) {
			tbl <- snpVariantsOnly( tbl, mode="single")
		}
		if ( ! nrow(tbl)) next

		depth <- VCF.TotalDepth( tbl$INFO)
		score <- VCF.Score( tbl$GENOTYPE_CALL)
		pctRef <- VCF.PctReference( tbl$INFO)
		tbl$INFO <- depth
		colnames(tbl)[ match( "INFO",colnames(tbl))] <- "DEPTH"
		tbl$SCORE <- score
		tbl$PCT_REF <- pctRef
		drops <- which( score < min.score | depth < min.depth)
		if ( length( drops) > 0) tbl <- tbl [ -drops, ]
		if ( ! nrow(tbl)) next

		# let's augment with the AA SNP if at all possible
		tbl <- DNAvariantsToAAvariants( tbl, seqID=sid, genomeFastaFile=genomeFastaFile)

		out <- rbind( out, tbl)
		cat( "\r", sid, nrow(tbl))
	}
	

	# we now allow the tool to have been called on a single gene
	# if so, combine all the actual base call counts with the SNP calling
	snpOut <- out
	if ( ! is.null( geneID) && nrow(out)) {
		geneAns1 <- pipe.BAMpileup( sampleID, geneID=geneID, max.depth=0)
		geneAns2 <- MPU.callStringsToMatrix( geneAns1$BASE_TABLE)
		# merge and recombine all these perspectives of the data
		out1 <- cbind( geneAns1[ , 1:2], "GENE_ID"=geneID, geneAns1[ , c(3,5)], "RAW_DEPTH"=geneAns1$DEPTH, stringsAsFactors=F)
		colnames(geneAns2) <- paste( sub( "^,$","REF",colnames(geneAns2)), "DEPTH", sep="_")
		# see which positions had a SNP result, and use those too
		where <- match( out1$POSITION, snpOut$POSITION, nomatch=NA)
		noData <- which( is.na(where))
		out2 <- snpOut[ where, 5:11]
		noData <- which( is.na(where))
		if ( length(noData)) {
			out2$ALT_BASE[noData] <- ""
			out2$QUAL[noData] <- 0
			out2$SCORE[noData] <- 0
			out2$PCT_REF[noData] <- 100
			out2$ALT_AA[noData] <- ""
			out2$GENOTYPE_CALL[noData] <- ""
		}
		colnames(out2) <- sub( "^DEPTH$", "VCF_DEPTH", colnames(out2))
		out <- cbind( out1, geneAns2, out2, stringsAsFactors=F)
		rownames(out) <- 1:nrow(out)

		# the calculation about "PCT_REF" that the SNP caller returns only uses the depths as seen by VCFTOOLS, that
		# is often much lower that the true depths.  So re-calc from first principles
		hasData <- which( ! is.na(where))
		if ( length( hasData)) {
			myRefCnt <- geneAns2[ hasData, "REF_DEPTH"]
			myTotalCnt <- apply( geneAns2[ hasData, , drop=F], 1, sum, na.rm=T)
			out$PCT_REF[ hasData] <- round( myRefCnt * 100 / myTotalCnt, digits=2)
		}
	}

	# tiny chance of no SNPs at all....
	if ( ! nrow(out)) {
		null <- vector( mode="character", length=0)
		out <- data.frame( "SEQ_ID"=null, "POSITION"=null, "GENE_ID"=null, "REF_BASE"=null,
				"ALT_BASE"=null, "QUAL"=null, "DEPTH"=null, "GENOTYPE_CALL"=null,
				"SCORE"=null, "PCT_REF"=null, "ALT_AA"=null)
	}
	outfile <- paste( sampleID, prefix, "Summary.VCF.txt", sep=".")
	if ( ! is.null( geneID)) {
		outfile <- paste( sampleID, prefix, geneID, "Summary.VCF.txt", sep=".")
	}
	outfile <- file.path( vcfPath, outfile)
	write.table( out, outfile, sep="\t", quote=F, na="", row.names=F)
	cat( "\nWrote variant call summary: ", outfile, "\nN_Calls: ", nrow(snpOut), "\n")

	if ( ! is.null( geneID)) {
		return( invisible( out))
	}
	return()
}


`pipe.VariantCalls2html` <- function( sampleID, annotationFile="Annotation.txt", optionsFile="Options.txt", 
					speciesID=getCurrentSpecies(), results.path=NULL, max.plots=200, min.score=100,
					tailWidth=46, altAA.only=TRUE, label="", ...) {

	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	}
	if ( speciesID != getCurrentSpecies()) setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()

	snp.path <- file.path( results.path, "VariantCalls", sampleID)
	snpFile <- file.path( snp.path, paste( sampleID, prefix, "Summary.VCF.txt", sep="."))
	if ( ! file.exists( snpFile)) {
		cat( "\nUnable to find SNP summary file: ", snpFile)
		cat( "\nPerhaps run 'pipe.VariantCalls()' first...")
		return( NULL)
	}
	tbl <- read.delim( snpFile, as.is=T)
	cat( "\nTotal SNP Count: ", nrow(tbl))
	drops <- which( tbl$SCORE < min.score)
	if ( length( drops)) {
		tbl <- tbl[ -drops, ]
		cat( "\nDropped for SCORE < ", min.score, ":  ", length(drops))
	}
	if ( altAA.only) {
		drops <- which( tbl$ALT_AA == "")
		if ( length( drops)) {
			tbl <- tbl[ -drops, ]
			cat( "\nDropped as synonymous: ", length(drops))
		}
	}
	# order by SCORE, then by coordinates
	ord <- order( -tbl$SCORE, tbl$SEQ_ID, tbl$POSITION)
	tbl <- tbl[ ord, ]
	rownames(tbl) <- 1:nrow(tbl)

	N <- min( max.plots, nrow( tbl))
	if ( N < 1) return(NULL)

	localPlotPath <- "SNP_Plots"
	plotPath <- file.path( snp.path, localPlotPath)
	if ( ! file.exists(plotPath)) dir.create( plotPath, recursive=T)

	# format up the table
	out <- tbl[ 1:N, ]
	# convert the GeneIDs to have the genomic coord...
	origID <- shortGeneName( out$GENE_ID, keep=1)
	newID <- paste( origID, out$POSITION, sep=".")
	out$GENE_ID <- newID

	# the "GENOTYPE" field is not wanted
	isGENO <- match( "GENOTYPE_CALL", colnames(out), nomatch=0)
	if (isGENO) out <- out[ , -isGENO]
	# modify the SNP column names
	colnames(out)[4:ncol(out)] <- gsub( "_", " ", colnames(out)[4:ncol(out)])

	globalOutfile <- file.path( snp.path, sub( "txt$", "html", basename(snpFile)))
	table2html( out, globalOutfile, title=paste( "Top SNPs:  &nbsp; ", sampleID, " &nbsp; ", label), 
			maxRows=N, linkPaths=localPlotPath)

	for ( i in 1:N) {
		myPos <- out$POSITION[i]
		mySeq <- out$SEQ_ID[i]
		myGene <- origID[i]
		myFileName <- paste( myGene, myPos, "png", sep=".")
		pipe.PlotSNP( sampleID, seqID=mySeq, pos=myPos, tailWidth=tailWidth, 
				annotationFile=annotationFile, optionsFile=optionsFile, results.path=results.path,
				plotFormat="png", plotFileName=myFileName, plot.path=plotPath, label=label, ...)
		cat( "\r", i, mySeq, myGene, myPos)
	}
}


`DNAvariantsToAAvariants` <- function( tbl, seqID, genomeFastaFile, altBaseColumn="ALT_BASE") {

	# given the SNPs to one chromosome, use the genomic sequence and current maps to
	# find which are no-synonymous
	chromoDNA <- getFastaSeqFromFilePath( genomeFastaFile, seqID, verbose=F)
	chromoDNA <- base::strsplit( chromoDNA, split="")[[1]]

	# get the set of coding genes that we need to assay
	genes <- sort( unique( tbl$GENE_ID))
	gmap <- subset( getCurrentGeneMap(), REAL_G == TRUE & SEQ_ID == seqID)
	cmap <- subset( getCurrentCdsMap(), SEQ_ID == seqID)
	genes <- intersect( genes, gmap$GENE_ID)
	NG <- length(genes)

	# get the set of SNP calls.  If multiple variants, just use the first
	pos <- tbl$POSITION
	altBase <- tbl[[ altBaseColumn]]
	altBase <- sub( ",.+", "", altBase)
	# slight chance of more than one call at same loci
	isDupPos <- which( duplicated( pos))
	if ( length( isDupPos)) {
		pos <- pos[ -isDupPos]
		altBase <- altBase[ -isDupPos]
		tbl <- tbl[ -isDupPos, ]
	}

	# OK, we are ready to visit each gene, modify its DNA and see what the AA calls are
	aaOut <- rep.int( "", nrow(tbl))
	whereInGeneMap <- match( genes, gmap$GENE_ID)
	if ( NG) {
	    for ( i in 1:NG) {
		thisGene <- genes[i]
		thisStart <- gmap$POSITION[ whereInGeneMap[i]]
		thisStop <- gmap$END[ whereInGeneMap[i]]
		thisStrand <- gmap$STRAND[ whereInGeneMap[i]]
		myRows <- which( tbl$GENE_ID == thisGene)
		myPositions <- pos[ myRows]
		# for simplicity, just treat all as 1-base SNPs, indels will kill reading frame
		myAltBases <- substr( altBase[ myRows], 1, 1)
		# also convert the chromosome locations to intra-gene, and make sure all are inside
		myGenePositions <- myPositions - thisStart + 1
		drops <- which( myPositions < thisStart | myPositions > thisStop)
		if ( length( drops)) {
			myRows <- myRows[ -drops]
			myPositions <- myPositions[ -drops]
			myGenePositions <- myGenePositions[ -drops]
			myAltBases <- myAltBases[ -drops]
		}
		# fill our query sequence
		dnaQuery <- chromoDNA[ thisStart:thisStop]
		for (j in 1:length(myPositions)) {
			where <- myGenePositions[j]
			dnaQuery[ where] <- myAltBases[j]
		}
		# now see how that sequence translates!
		#cat( "\n", i, thisGene, thisStart, thisStop, "|", myPositions)
		ans <- convertGenomicBasesToCodingAA( seqID, thisStart, thisStop, strand=thisStrand,
					dnaQuery=dnaQuery, genomeDNA=chromoDNA, geneMap=gmap, 
					cdsMap=cmap, geneID=thisGene)

		# label the locations so we can keep track of who's who
		names(ans$genomic) <- names(ans$query) <- thisStart : thisStop
		#SAVE_ANS <<- ans

		# now back-deduce the AA call and location
		genomeAA <- ans$genomic[ seq( 2, length(ans$genomic), by=3)]
		queryAA <- ans$query[ seq( 2, length(ans$query), by=3)]
		if ( thisStrand == "-") {
			aaCalls <- paste( genomeAA, length(genomeAA):1, queryAA, sep="")
		} else {
			aaCalls <- paste( genomeAA, 1:length(genomeAA), queryAA, sep="")
		}
		genomeAA.long <- rep( genomeAA, each=3)
		queryAA.long <- rep( queryAA, each=3)
		aaCalls.long <- rep( aaCalls, each=3)
		posOut <- thisStart:thisStop
		#n <- length(posOut)
		#if ( length( posOut) != length(genomeAA.long)) {
			#n <- min( length(posOut), length( genomeAA.long))
		#}
		#smallSAVE <<- data.frame( "POSITION"=posOut[1:n], "Genome"=genomeAA.long[1:n], 
					#"Query"=queryAA.long[1:n], "AA_Call"=aaCalls.long[1:n])

		whoSNP <- which( genomeAA.long != queryAA.long)
		realSNPs <- which( myGenePositions %in% whoSNP)

		if ( length( whoSNP)) {
			#cat( "\nDebug: ", thisGene, "|", whoSNP, "|", realSNPs, "\n")
			#print( smallSAVE[ whoSNP, ])

			aaOut[ myRows[ realSNPs]] <- aaCalls.long[ myGenePositions[realSNPs]]
		}
	    }
	}

	tbl$ALT_AA <- aaOut
	return( tbl)
}


`exonVariantsOnly` <- function( tbl, exonMap=getCurrentExonMap()) {

	if ( nrow(tbl) < 1) return(tbl)

	mySID <- unique.default( tbl$SEQ_ID) 
	if ( length( mySID) > 1) stop( "Error in 'exonVariantsOnly':  Expected exactly one chromosome of SNP data")
	eSID <- unique.default( exonMap$SEQ_ID) 
	if ( length( eSID) > 1) {
		exonMap <- subset.data.frame( exonMap, SEQ_ID == mySID)
	}
	if ( nrow( exonMap) < 2) return(tbl)
	
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
	snpPos <- findInterval( tbl$POSITION, posVec, all.inside=FALSE)

	# now call each SNP position as being in an exon or not
	isExon <- sapply( snpPos, function(x) {
		
		# outside the boundaries is instant NO
		if ( x < 1 || x >= NE2) return( FALSE)
		# in an exon is NO
		if ( typeVec[ x] == "E") return( TRUE)
		FALSE
	})

	isIntron <- ! isExon
	if ( any( isIntron)) {
		cat( "  N_NonExon: ", sum( isIntron))
		return( tbl[ isExon, ])
	} else {
		return( tbl)
	}
}


`pipe.VariantComparison` <- function( sampleIDset, groupSet=sampleIDset, annotationFile="Annotation.txt",
				optionsFile="Options.txt", speciesID=getCurrentSpecies(), results.path=NULL, 
				min.deltaScore=40, exonOnly=TRUE, snpOnly=TRUE, AAsnpOnly=exonOnly, 
				capScore=200, min.depth=1) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	vcfPath <- file.path( results.path, "VariantCalls")

	setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()

	groupSet <- checkGroupNames( groupSet)

	# gather all the VCF calls from all these samples
	seqV <- posV <- geneV <- refV <- altV <- altAA <- scoreV <- qualV <- sidV <- gidV <- depthV <- vector()
	nall <- 0

	NS <- length(sampleIDset)
	cat( "\nLoading all VCF summary data..\n")
	for ( i in 1:NS) {
		sampleID <- sampleIDset[i]
		groupID <- groupSet[i]

		vcffile <- file.path( vcfPath, sampleID, paste( sampleID, prefix, "Summary.VCF.txt", sep="."))
		if ( ! file.exists( vcffile)) {
			cat( "VCF Summary file not found: ", vcffile)
			return(NULL)
		}

		tbl <- read.delim( vcffile, as.is=T)
		NC <- nrow(tbl)
		cat( "\n", i, sampleID, "\tN_SNP: ", NC)
		if ( NC < 1) next

		if (snpOnly) {
			tbl <- snpVariantsOnly( tbl, mode="single")
			NC <- nrow(tbl)
			if ( NC < 1) next
		}

		if (exonOnly) {
			# clean it by SeqID
			sids <- sort( unique( tbl$SEQ_ID))
			newtbl <- data.frame()
			for ( s in sids) {
				isS <- which( tbl$SEQ_ID == s)
				newsml <- exonVariantsOnly( tbl[ isS, ])
				newtbl <- rbind( newtbl, newsml)
			}
			tbl <- newtbl
			NC <- nrow(tbl)
			if ( NC < 1) next
		}

		# don't delete them here...  wait until we combine...
		#if (AAsnpOnly) {
			#if ( "ALT_AA" %in% colnames(tbl)) {
				#tbl <- subset( tbl, ALT_AA != "")
				#NC <- nrow(tbl)
				#if ( NC < 1) next
			#}
		#}

		now <- (nall+1) : (nall+NC)
		seqV[ now] <- tbl$SEQ_ID
		posV[ now] <- tbl$POSITION
		geneV[ now] <- tbl$GENE_ID
		refV[ now] <- tbl$REF_BASE
		altV[ now] <- tbl$ALT_BASE
		altAA[ now] <- tbl$ALT_AA
		qualV[ now] <- tbl$QUAL
		depthV[ now] <- tbl$DEPTH
		scoreV[ now] <- tbl$SCORE
		sidV[ now] <- sampleID
		gidV[ now] <- groupID
		nall <- max( now)
		cat( "\tTotal SNPs: ", nall)
	}
	cat( "\nDone.\n")

	# any tweaking of fields...
	# 1. the ALT base may be comma separated...  keep only the first
	altV <- sub( ",.+","", altV)

	# now that we have all the VCF data, we need to visit each loci and summarize any 
	# differences found between groups
	cat( "\nOranizing all loci..")
	posKey <- paste( seqV, posV, sep="::")
	keyFac <- factor( posKey)
	nout <- nlevels( keyFac)
	groupFac <- factor( groupSet)
	NG <- nlevels(groupFac)
	groupNames <- levels(groupFac)
	gBase <- gAA <- gScore <- gDepth <- vector( length=NG)

	# storage for the final answers
	seqOut <- posOut <- geneOut <- refOut <- depthOut <- deltaOut <- vector( length=nout)
	grpBaseOut <- matrix( "", nrow=nout, ncol=NG)
	grpAAout <- matrix( "", nrow=nout, ncol=NG)
	grpScoreOut <- matrix( 1, nrow=nout, ncol=NG)
	grpDepthOut <- matrix( 1, nrow=nout, ncol=NG)
	colnames(grpBaseOut) <- colnames(grpAAout) <- colnames(grpScoreOut) <- colnames(grpDepthOut) <- groupNames

	# apply a cap on scores to temper the effect of read depth variation
	cat( "\nCapping high SNP scores at:  ", capScore)
	scoreV[ scoreV > capScore] <- capScore

	cat( "\nEvaluating ", nout, " VCF loci..\n")
	iout <- 0
	tapply( 1:nall, keyFac, function(x) {

		# start from the null state that all are the reference
		# we will give the reference a high score of being the reference, 
		# and then revert it back to '1 = 'not a SNP' later
		trueRefBase <- refV[x[1]]

		# and estimate the depth and scoring from what SNPs were seen here
		thisAvgDepth <- mean( depthV[x])
		thisAvgScore <- mean( scoreV[x])
		myBase <- rep.int( trueRefBase, NS)
		myAA <- rep.int( "", NS)
		myScore <- rep.int( thisAvgScore, NS)   # used to use '1'
		myDepth <- rep.int( thisAvgDepth, NS)   # used to use '0'
		myScore <- rep.int( NA, NS)   # used to use '1'
		myDepth <- rep.int( NA, NS)   # used to use '0'

		# now fill in the ALTs that we saw
		who <- match( sidV[x], sampleIDset)
		myBase[who] <- altV[x]
		myAA[who] <- altAA[x]
		myScore[who] <- scoreV[x]
		myDepth[who] <- depthV[x]

		# now we combine by group to get the average of each
		igrp <- 0
		tapply( 1:NS, groupFac, function(j) {
			allBases <- myBase[j]
			allAA <- myAA[j]
			allScores <- myScore[j]
			allDepths <- myDepth[j]
			# catch those samples that had no SNP info
			#allScores[ is.na(allScores)] <- 1
			allDepths[ is.na(allDepths)] <- 1
			igrp <<- igrp + 1
			thisGroupsBases <- sort.int( table( allBases), decreasing=T)
			thisGroupsAAs <- sort.int( table( allAA), decreasing=T)
			# when they all agree, easy case!
			if ( length( thisGroupsBases) == 1) {
				gBase[igrp] <<- names( thisGroupsBases)[1]
				gAA[igrp] <<- names( thisGroupsAAs)[1]
				gScore[igrp] <<- mean( allScores, na.rm=T)
				gDepth[igrp] <<- mean( allDepths, na.rm=T)
			} else {
				# more than one called, so try to do a weigthed average 
				grpWeightedScores <- tapply( allScores, factor(allBases), function(k) {
							tmpScore <- mean( k, na.rm=T)
							tmpPct <- length(k) / length(allBases)
							return( tmpScore * tmpPct)
						})
				bestOne <- which.max(grpWeightedScores)
				gBase[igrp] <<- bestBase <- names( grpWeightedScores)[ bestOne]
				# for the AA call, take the most common non-blank
				hasAA <- which( allAA != "")
				if( length(hasAA)) {
					aaWtScores <- tapply( allScores[hasAA], factor(allAA[hasAA]), function(k) {
							tmpScore <- mean( k, na.rm=T)
							tmpPct <- length(k) / length(hasAA)
							return( tmpScore * tmpPct)
						})
					bestOneAA <- which.max(aaWtScores)
					gAA[igrp] <<- names( aaWtScores)[ bestOneAA]
				} else {
					gAA[igrp] <<- ""
				}
				gAA[igrp] <<- allAA[ which( allBases == bestBase)[1]]
				gScore[igrp] <<- grpWeightedScores[ bestOne]
				gDepth[igrp] <<- mean( allDepths[ allBases == bestBase], na.rm=T)
			}
		})

		# OK, the reference's score weight was used in group averages, now revert it to its original
		# meaning (i.e.  100% sure reference == 1)
		isREF <- which( gBase == trueRefBase)
		gScore[ isREF] <- capScore - gScore[isREF] + 1
		gScore[ is.na(gScore)] <- 1
		#gDepth[ is.na(gDepth)] <- NA

		# find the biggest differences between the groups
		delta <- diff( range( gScore))
		if (NG == 1) delta <- gScore[1]

		# penalize when some groups are too low depth
		pctTooLow <- (min.depth - gDepth) / min.depth
		pctTooLow[ pctTooLow <= 0] <- 0
		if ( any( pctTooLow > 0)) {
			penalty <- mean( pctTooLow, na.rm=T)
			delta <- delta * ( 1 - penalty)
		}
		# penalize when all the base calls are the same, and all the AA calls are the same
		if ( length( unique( gBase)) == 1) delta <- delta / 2
		if ( length( unique( gAA)) == 1) delta <- delta / 2

		# fill the results data for this loci
		iout <<- iout + 1
		seqOut[iout] <<- seqV[ x[1]]
		posOut[iout] <<- posV[ x[1]]
		geneOut[iout] <<- geneV[ x[1]]
		refOut[iout] <<- refV[ x[1]]
		deltaOut[iout] <<- delta
		grpBaseOut[ iout, ] <<- gBase
		grpAAout[ iout, ] <<- gAA
		grpScoreOut[ iout, ] <<- gScore
		grpDepthOut[ iout, ] <<- gDepth
		if ( iout %% 1000 == 0) cat( "\r", iout, seqOut[iout], shortGeneName( geneOut[iout], keep=1), 
						posOut[iout], gBase, gAA, gScore, gDepth, "      ")
	})

	# now build the final result
	out <- data.frame( "SEQ_ID"=seqOut, "POSITION"=posOut, "GENE_ID"=geneOut, "PRODUCT"=gene2Product( geneOut),
			"REF_BASE"=refOut, "DELTA_SCORE"=round(deltaOut, digits=2), stringsAsFactors=FALSE)
	for ( j in 1:NG) {
		small <- data.frame( grpBaseOut[ ,j], grpAAout[ ,j], round( grpScoreOut[ ,j], digits=2), 
				round( grpDepthOut[ ,j]), stringsAsFactors=FALSE)
		colnames(small) <- paste( groupNames[j], c( "BASE","AA","SCORE","DEPTH"), sep="_")
		out <- cbind( out, small, stringsAsFactors=FALSE)
	}

	# give more detail to VAR gene SNPs
	if ( speciesID == "Pf3D7") {
		vgmap <- getVargeneDomainMap()
		isVAR <- which( out$GENE_ID %in% vgmap$GENE_NAME)
		for ( j in isVAR) {
			smlvgmap <- subset( vgmap, GENE_NAME == out$GENE_ID[j] & POSITION <= out$POSITION[j] &
						END >= out$POSITION[j])
			if ( nrow(smlvgmap)) {
				domainID <- smlvgmap$DOMAIN_ID[1]
				cassette <- smlvgmap$CASSETTE[1]
				grp <- smlvgmap$GENE_GROUP[1]
				extraText <- paste( "{", cassette, ", ", domainID, ", ", "Grp:", grp, "}", sep="")
				out$PRODUCT[j] <- paste( out$PRODUCT[j], extraText, sep="  ")
			}
		}
	}

	# when there are more than one different SNP base at the same loci, that should get a higher score 
	# than returned by 'diff(range())'.
	cat( "\nRescoring multiple alternate alleles..")
	firstAltColumn <- 7
	altColumnStepSize <- 4
	if ( NG > 1) {
		altcolumns <- seq( firstAltColumn, ncol(out), by=altColumnStepSize)
		for ( i in 1:nrow(out)) {
			ref <- out$REF_BASE[i]
			alts <- substr( as.character(out[ i, altcolumns]), 1,1)
			uniqueAlts <- setdiff( unique( alts), ref)
			if ( length( uniqueAlts) > 1) {
				isAlt <- which( alts %in% uniqueAlts)
				altScores <- as.numeric( out[ i, altcolumns+1])
				altFac <- factor( alts[ isAlt])
				altScores <- tapply( altScores[isAlt], altFac, mean)
				out$DELTA_SCORE[i] <- sum( altScores)
			}
		}
	}

	# now at this 'almost final' point, before re-ordering, see if we will drop rows that 
	# did not cause an AA mutation
	if (AAsnpOnly) {
		# rows that are all blank mean never saw a AA SNP
		nAAsnp <- apply( grpAAout, MARGIN=1, FUN=function(x) sum( x != ""))
		drops <- which( nAAsnp < 1)
		if ( length( drops)) {
			cat( "\nDropping fully synonymous SNP sites.  N_Given: ", nrow(out), "  N_Drop: ", length(drops))
			out <- out[ -drops, ]
		}
	}

	# final order is those most different bwtween the groups
	keep <- which( out$DELTA_SCORE >= min.deltaScore)
	if ( length(keep)) {
		out <- out[ keep, ]
		ord <- order( -(out$DELTA_SCORE), out$SEQ_ID, out$POSITION)
		out <- out[ ord, ]
		rownames(out) <- 1:nrow(out)
		cat( "\nDone.  N_Differential_SNPs: ", nrow(out), "\n")
		return(out)
	} else {
		cat( "\nNo Differential SNPs.  Perhap loosen 'min.deltaScore' or other arguments..? \n")
		return( data.frame())
	}
}


`snpVariantsOnly` <- function( tbl, mode=c("compare", "single"), min.depth=1) {

	mode <- match.arg( mode)
	if ( nrow(tbl) < 1) return( tbl)

	ref <- tbl$REF_BASE
	if ( mode == "compare") {
		alt <- alt1 <- tbl[[ 7]]
		depth1 <- as.numeric( tbl[[9]])
		alt2 <- tbl[[ 10]]
		depth2 <- as.numeric( tbl[[12]])
	} else {
		alt <- sub( ",.+", "", tbl$ALT_BASE)
	}

	# we want to leave out any indels, and keep only simple SNPs
	isSNPref <- which( nchar( ref) == 1)
	isSNPalt <- which( nchar( alt) == 1)
	if ( mode == "compare") {
		isSNPalt2 <- which( nchar( alt2) == 1)
		isSNPalt <- intersect( isSNPalt, isSNPalt2)
		# secondly, are these 2 alts really different from each other
		base1 <- alt1[ isSNPalt]
		base2 <- alt2[ isSNPalt]
		good <- which( base1 != base2)

		# secondly, if a sample has too few reads, it not proof of a SNP
		#deep1 <- which( depth1[isSNPalt] >= min.depth)
		#deep2 <- which( depth2[isSNPalt] >= min.depth)
		#deepEnough <- intersect( deep1, deep2)
		#good <- intersect( good, deepEnough)

		# recall that a location matching the reference was not called a SNP, so we have no depth info...!!
		mostDeep <- pmax( depth1[isSNPalt], depth2[isSNPalt])
		deepEnough <- which( mostDeep >= min.depth)
		good <- intersect( good, deepEnough)

		# OK, its a real SNP
		isSNPalt <- isSNPalt[ good]
	} else {
		# secondly, is this alt really different from reference
		base1 <- tbl$REF_BASE[ isSNPalt]
		base2 <- tbl$ALT_BASE[ isSNPalt]
		good <- which( base1 != base2)
		isSNPalt <- isSNPalt[ good]
	}

	# those that satisfy all are the ones we want
	keep <- intersect( isSNPref, isSNPalt)
	out <- tbl[ keep, ]
	if ( nrow(out)) rownames(out) <- 1:nrow(out)
	return( out)
}


`pipe.VariantCompare2html` <- function( tbl, outfile="VariantCompare.html", 
					sampleIDset, groupSet=sampleIDset, 
					annotationFile="Annotation.txt", optionsFile="Options.txt", 
					speciesID=getCurrentSpecies(), results.path=NULL, Ngenes=200,
					tailWidth=26, indelCharLen=8, out.path=sub( ".html", "", outfile), 
					showDepth=TRUE, dropNonSNPs=TRUE, label="", ...) {

	REF_column <- 5
	SNP_columns <- seq( 7, ncol(tbl), by=3)

	# some SNPs will be the same call, just different depths... ignore those
	if (dropNonSNPs) {
		ndiff <- apply( as.matrix( tbl[ ,SNP_columns]), MARGIN=1, function(x) length( unique.default( x)))
		same <- which( ndiff == 1)
		if ( length( same) > 0) tbl <- tbl[ -same, ]
	}

	N <- min( Ngenes, nrow( tbl))
	if ( N < 1) return()

	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	outPath <- file.path( results.path, "VariantCalls", out.path)
	if ( ! file.exists(outPath)) dir.create( outPath, recursive=T)
	localPlotPath <- paste( out.path, "pngPlots", sep=".")
	plotPath <- file.path( outPath, localPlotPath)
	if ( ! file.exists(plotPath)) dir.create( plotPath, recursive=T)

	fastaFile <- getOptionValue( optT, "genomicFastaFile")

	out <- tbl[ 1:N, ]
	rownames(out) <- 1:N
	toplot <- out[ order( out$SEQ_ID, out$POSITION), ]

	# format up the table
	# convert the GeneIDs to have the genomic coord...
	newID <- paste( shortGeneName( out$GENE_ID, keep=1), out$POSITION, sep=".")
	out$GENE_ID <- newID
	# trim long indels
	for ( j in c( REF_column, SNP_columns)) {
		txt <- out[[j]]
		islong <- which( nchar(txt) > indelCharLen)
		if ( length(islong) > 0) {
			newtxt <- paste( substr( txt[islong], 1, indelCharLen), "...", sep="")
			out[[j]][ islong] <- newtxt
		}
	}
	# modify the SNP column names
	colnames(out)[5:ncol(out)] <- gsub( "_", " ", colnames(out)[5:ncol(out)])
	colnames(out)[5:ncol(out)] <- gsub( "\\.", " ", colnames(out)[5:ncol(out)])

	if ( ! showDepth) {
		drops <- grep( "DEPTH", colnames(out))
		if ( length(drops)) out <- out[ , -drops]
	}

	globalOutfile <- file.path( outPath, outfile)
	table2html( out, globalOutfile, title=paste( "Variant Comparison: &nbsp; ", label), maxRows=N,
			linkPaths=localPlotPath)
	# write a text version too
	globalOutfile <- file.path( outPath, sub( "html$", "txt", outfile))
	write.table( tbl, globalOutfile, sep="\t", quote=F, row.names=F)

	for ( i in 1:N) {
		myPos <- toplot$POSITION[i]
		mySeq <- toplot$SEQ_ID[i]
		myGene <- shortGeneName( toplot$GENE_ID[i], keep=1)
		myFileName <- paste( myGene, myPos, "png", sep=".")
		pipe.PlotSNP( sampleIDset, seqID=mySeq, position= myPos, geneID=myGene, 
				tailWidth=tailWidth, groups=groupSet,
				annotationFile=annotationFile, optionsFile=optionsFile, results.path=results.path,
				plotFormat="png", plotFileName=myFileName, plot.path=plotPath, label=label, ...)
		cat( "\r", i, mySeq, myGene, myPos)
	}
}


VCF.Score <- function( PLterms) {

	# we 'are back to' doing just the homozygous reference value

	# now the genotype field is GT:PL:DP, so we need the middle one
	#scoreStr <- sub( ":.+", "", PLterms)
	myPLterms <- base::strsplit( PLterms, split=":", fixed=T)
	scoreStr <- sapply( myPLterms, function(x) x[2])
	scoreTerms <- base::strsplit( scoreStr, split=",", fixed=T)
	scores <- sapply( scoreTerms, function(x) {
		max( as.integer(x[1]), na.rm=T)
	})
	scores[ is.na(scores)] <- 0
	return( scores)
}


VCF.Pvalue <- function( PLterms) {

	scores <- VCF.Score( PLterms)
	pval <- 10 ^ (scores/-10)
	return( pval)
}


VCF.HighScoreDepth <- function( PLterms) {

	highScoringDepth <- sub( ".+:", "", PLterms)
	highScoringDepth <- as.integer( highScoringDepth)
	highScoringDepth[ is.na(highScoringDepth)] <- 0
	return( highScoringDepth)
}


VCF.TotalDepth <- function( info) {

	totalDepth <- sub( "(^.*)(DP=)([0-9]+)(;.*$)", "\\3", info)
	totalDepth <- as.integer( totalDepth)
	totalDepth[ is.na(totalDepth)] <- 0
	return( totalDepth)
}


VCF.PctReference <- function( info) {

	pctRef <- rep.int( NA, length(info))

	dp4Str <- sub( "(^.*)(DP4=)([0-9]+,[0-9]+,[0-9]+,[0-9]+)(;.*$)", "\\3", info)
	dp4 <- base::strsplit( dp4Str, split=",")
	# these should be exactly 4 integers
	dp4Len <- sapply( dp4, length)
	good <- which( dp4Len == 4)
	if ( length( good)) {
		dp4M <- matrix( as.numeric( unlist( dp4[good])), nrow=4, ncol=length(good))
		nRef <- apply( dp4M[ 1:2, , drop=F], MARGIN=2, sum, na.rm=T)
		nAlt <- apply( dp4M[ 3:4, , drop=F], MARGIN=2, sum, na.rm=T)
		myPct <- 100 * nRef / (nRef + nAlt)
		pctRef[ good] <- round( myPct, digits=2)
	}
	return( pctRef)
}


`pipe.VariantTable` <- function( sampleIDset, seqID, annotationFile="Annotation.txt",
				optionsFile="Options.txt", speciesID=getCurrentSpecies(), results.path=NULL,
				min.depth=1, max.depth=max.depth, mpileupArgs="", 
				explicit.reference=TRUE, call.mode=c("percents","counts"), findVariants=TRUE) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	target <- getAndSetTarget( optionsFile, sampleID=sampleIDset[1], annotationFile=annotationFile)
	call.mode <- match.arg( call.mode)

	fastaFile <- getOptionValue( optT, "genomicFastaFile")
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	} else {
		results.path <- results.path
	}
	vcfPath <- file.path( results.path, "VariantCalls")
	if ( ! file.exists( vcfPath)) dir.create( vcfPath, recursive=TRUE)

	# make sure we have the BAM file already sorted
	bamfilelist <- vector()
	for ( sampleID in sampleIDset) {
		bamfile <- paste( sampleID, "genomic.bam", sep=".")
		bamfile <- file.path( results.path, "align", bamfile)
		sortedbamfile <- BAM.verifySorted( bamfile, index=TRUE)
		bamfilelist <- c( bamfilelist, sortedbamfile)
	}

	
	# the function that does one chromosome of one sample
	`variantTableOneSeq` <- function( i) {

		sampleID <- sampleIDset[i]
		bamfile <- bamfilelist[i]

		outPath <- file.path( vcfPath, sampleID)
		if ( ! file.exists( outPath)) dir.create( outPath, recursive=T)
		outfile <- file.path( outPath, paste( sampleID, seqID, "VariantTable.txt", sep="."))

		# if we don't need to calc, and its already there, just grab it
		if ( file.exists( outfile) && ! findVariants) {
			out <- read.delim( outfile, as.is=T)
			cat( "\nTable exists: ", sampleID, seqID, "  N_Bases:", nrow(out))
			return( dim(out))
		}

		ans <- BAM.mpileup( bamfile, seqID=seqID, fastaFile=fastaFile, 
					min.depth=min.depth, max.depth=max.depth, mpileupArgs=mpileupArgs, 
					summarize.calls=FALSE)
		if ( is.null(ans)) return(NULL)
		N <- nrow(ans)
		cat( "\nMpileups: ", sampleID, seqID, "  N_Bases:", N)

		# now convert the calls to table format
		cat( "\tTablify base calls..")
		rawBases <- ans$CALL_BASE
		callAns <- MPU.callBases( rawBases, ans$REF_BASE, mode=call.mode)
		ans$CALL_BASE <- callAns$call
		baseTbl <- MPU.callTableToString( callAns$depth.table)
		if (explicit.reference) {
			rbase <- ans$REF_BASE
			baseTbl <- mapply( rbase, baseTbl, SIMPLIFY=TRUE, 
					FUN=function(x,y) gsub( ",", x, y, fixed=T))
		}
		# over write the call_score field
		colnames(ans) <- sub( "CALL_SCORE", "BASE_TABLE", colnames(ans))
		ans$BASE_TABLE <- baseTbl
		
		# lastly, stick some GeneIDs in the result
		gmap <- subset( getCurrentGeneMap(), SEQ_ID == seqID)
		where <- fastFindInterval( ans$POSITION, gmap$POSITION)
		gid <- rep.int( "", N)
		gid[ where <= N] <- gmap$GENE_ID[where]
		ans <- data.frame( "GENE_ID"=gid, ans, stringsAsFactors=F)

		# and now we can append any AA substitutions
		ans <- DNAvariantsToAAvariants( ans, seqID=seqID, genomeFastaFile=fastaFile, altBaseColumn="CALL_BASE")


		# save it and return it
		write.table( ans, outfile, sep="\t", quote=F, row.names=F)
		return(  dim(ans))
	}
	

	setCurrentSpecies( speciesID)
	seqMap <- getCurrentSeqMap()
	if ( ! seqID %in% seqMap$SEQ_ID) stop( paste( "'seqID' not in current seq map: ", seqIDset))

	N <- length( sampleIDset)
	if ( N > 1) {
		ans <- multicore.lapply( 1:N, FUN=variantTableOneSeq)
	} else { 
		ans <- variantTableOneSeq( 1)
	}
	
	# now we read up all these samples, and find those rows that are variants in 1+ samples
	varfileset <- file.path( vcfPath, sampleIDset, paste( sampleIDset, seqID, "VariantTable.txt", sep="."))

	# premake the final storage
	Ndetails <- 4
	NS <- length( sampleIDset)
	chromoLen <- subset( seqMap, SEQ_ID == seqID)$LENGTH
	outPos <- 1:chromoLen
	outGID <- outRef <- rep.int( "", chromoLen)
	outDetails <- matrix( "", nrow=chromoLen, ncol=NS*Ndetails)
	outDetailNames <- 1:ncol(outDetails)

	cat( "\nCombining..")
	for ( i in 1:NS) {
		thisFile <- varfileset[i]
		thisSample <- sampleIDset[i]
		tbl <- read.delim( thisFile, as.is=T)

		# some results may not yet have the "ALT_AA" column
		if ( ! ("ALT_AA" %in% colnames(tbl))) {
			tbl <- DNAvariantsToAAvariants( tbl, seqID=seqID, genomeFastaFile=fastaFile, altBaseColumn="CALL_BASE")
			write.table( tbl, thisFile, sep="\t", quote=F, row.names=F)
		}

		# map what we found to the full chromosome
		wh <- match( tbl$POSITION, outPos)
		outGID[ wh] <- tbl$GENE_ID
		outRef[ wh] <- tbl$REF_BASE

		myCols <- (i-1) * Ndetails + (1:Ndetails)
		outDetails[ wh, myCols[1]] <- tbl$DEPTH
		outDetails[ wh, myCols[2]] <- tbl$CALL_BASE
		outDetails[ wh, myCols[3]] <- tbl$BASE_TABLE
		outDetails[ wh, myCols[4]] <- tbl$ALT_AA
		outDetailNames[myCols] <- make.names( paste( c("DEPTH", "CALL_BASE", "BASE_TABLE", "ALT_AA"), thisSample, sep="_"))
		cat( "  ", thisSample)
	}
	colnames(outDetails) <- outDetailNames

	# now we can trim down to just the variant rows
	cat( "\nDropping non-variant rows..")
	# case 1:  base was never seen at all
	drops1 <- which( outRef == "")
	cat( "  N_Blank:", length(drops1))
	# case 2:  base was never seen deeply enough
	depths <- matrix( 0, nrow=chromoLen, ncol=NS)
	for ( i in 1:NS) depths[,i] <- as.numeric( outDetails[ , (i-1)*Ndetails+1])
	depths[ is.na(depths)] <- 0
	maxDepth <- apply( depths, 1, max)
	drops2 <- which( maxDepth < min.depth)
	cat( "  N_NotDeep:", length(drops2))
	# case 3:  called base was always the reference
	calls <- outDetails[ , seq( 2, ncol(outDetails), by=Ndetails), drop=FALSE]
	allCalls <- apply( calls, 1, function(x) { xx <- unique(x);  if ( length(xx) < 2) xx[1] else paste( xx, collapse="")})
	drops3 <- which( outRef == allCalls)
	cat( "  N_NotVariant:", length(drops3))

	# OK, we have what we need
	drops <- sort( unique( c( drops1, drops2, drops3)))
	cat( "  N_To_Drop: ", length( drops))
	outGID <- outGID[ -drops]
	outRef <- outRef[ -drops]
	outPos <- outPos[ -drops]
	outDetails <- outDetails[ -drops, ]

	out <- data.frame( "GENE_ID"=outGID, "SEQ_ID"=seqID, "POSITION"=outPos, "REF_BASE"=outRef,
				outDetails, stringsAsFactors=FALSE)
	rownames(out) <- 1:nrow(out)
	cat( "\nFinal Table of Variants: ", nrow( out), "\n")

	return( out)
}


`pipe.VariantMerge` <- function( sampleIDset, annotationFile="Annotation.txt",
				optionsFile="Options.txt", speciesID=getCurrentSpecies(), results.path=NULL,
				min.qual=5,  min.depth=1, min.score=40, dropDuplicates=FALSE, 
				dropIndels=FALSE, verbose=TRUE, outfile=NULL) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	target <- getAndSetTarget( optionsFile, sampleID=sampleIDset[1], annotationFile=annotationFile)

	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	} else {
		results.path <- results.path
	}
	if ( getCurrentSpecies() != speciesID) setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()

	# go gather all the VCF calls over a set of samples
	vcfSet <- data.frame()
	for ( sid in sampleIDset) { 
		vcf.path <- file.path( results.path, "VariantCalls", sid)
		f <- file.path( vcf.path, paste( sid, prefix, "Summary.VCF.txt", sep=".")) 
		if ( ! file.exists( f)) {
			cat( "\nVariant Call file not found.  Skip:  ", f)
			next
		}
		sml <- read.delim( f, as.is=T) 
		sml$SampleID <- sid

		# small chance of having more than one call at a location
		if ( any( duplicated( sml$POSITION))) {
			ord <- order( sml$SEQ_ID, sml$POSITION, -sml$SCORE)
			sml <- sml[ ord, ]
			key <- paste( sml$GENE_ID, sml$POSITION, sep="::")
			keep <- which( ! duplicated(key))
			sml <- sml[ keep, ]
		}

		# perhaps drop low quality/score
		if ( ! is.null( min.qual)) {
			drops <- which( sml$QUAL < min.qual)
			if ( length( drops)) sml <- sml[ -drops, ]
		}
		if ( ! is.null( min.depth)) {
			drops <- which( sml$DEPTH < min.depth)
			if ( length( drops)) sml <- sml[ -drops, ]
		}
		if ( ! is.null( min.score)) {
			drops <- which( sml$SCORE < min.score)
			if ( length( drops)) sml <- sml[ -drops, ]
		}
		if ( dropIndels) {
			drops1 <- which( nchar( sml$REF_BASE) > 1)
			drops2 <- which( nchar( sml$ALT_BASE) > 1)
			drops <- sort( union( drops1, drops2))
			if ( length( drops)) sml <- sml[ -drops, ]
		}

		vcfSet <- rbind( vcfSet, sml); 
		if (verbose) cat( "\n", sid, nrow( sml), nrow(vcfSet))
	}

	# keep any site seen in at least 1 location, and keep the highest scoring of each
	if (dropDuplicates) {
		cat( "\nRemoving duplicates..\n")
		key <- paste( vcfSet$GENE_ID, vcfSet$POSITION, sep="::")
		keyFac <- factor( key)
		keyPtrs <- tapply( 1:length(key), keyFac, FUN=NULL)
		nKeys <- nlevels(keyFac)
		drops <- vector()
		for ( i in 1:nKeys) {
			thisSet <- which( keyPtrs == i)
			if (length(thisSet) < 2) next
			keeper <- thisSet[ which.max( vcfSet$SCORE[thisSet])]
			drops <- c( drops, setdiff( thisSet, keeper))
			cat( "\r", i, levels(keyFac)[i], length(drops))
		}
		if ( length(drops)) {
			drops <- sort( drops)
			vcfSet <- vcfSet[ -drops, ]
		}
	}

	ord <- order( vcfSet$SEQ_ID, vcfSet$POSITION, -vcfSet$SCORE, vcfSet$SampleID)
	vcfSet <- vcfSet[ ord, ]
	if ( nrow(vcfSet)) rownames(vcfSet) <- 1:nrow(vcfSet)

	if ( ! is.null( outfile)) write.table( vcfSet, outfile, sep="\t", quote=F, row.names=F)

	return( vcfSet)
}


`pipe.VariantSimilarity` <- function( sampleIDset, speciesID=getCurrentSpecies(), 
					annotationFile="Annotation.txt", optionsFile="Options.txt", 
					results.path=NULL, min.qual=1,  min.depth=1, min.score=1, 
					max.pctRef=NULL, max.SNPperGene=NULL,
					dropIntergenics=FALSE, dropIndels=TRUE, dropIntrons=TRUE, geneSubset=NULL, 
					debugPattern=NULL, freqM=NULL, min.freqDifference=60, verbose=TRUE) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	#target <- getAndSetTarget( optionsFile, sampleID=sampleIDset[1], annotationFile=annotationFile)

	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	if ( getCurrentSpecies() != speciesID) setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()

	# if we trim to just SNPs in Exons, we need to do some setup
	if (dropIntrons) {
		if (verbose) cat( "\nSetting up for intron SNP removal..")
		emap <- getCurrentExonMap()
		seqIDs <- sort( unique( emap$SEQ_ID))
		NSEQ <- length(seqIDs)
		exonPosSets <- vector( mode="list", length=NSEQ)
		names(exonPosSets) <- seqIDs
		for( i in 1:NSEQ) {
			sml <- subset.data.frame( emap, SEQ_ID == seqIDs[i])
			ans <- mapply( sml$POSITION, sml$END, FUN=`:`, SIMPLIFY=F, USE.NAMES=F)
			exonPosSets[[i]] <- unlist( ans)
			if (verbose) cat( " ", seqIDs[i])
		}
		if (verbose) cat( "  Done.")
	}

	useFM <- FALSE
	fm <- NULL
	if ( ! is.null( freqM)) {
		# may be a filename
		if ( is.character(freqM) && file.exists(freqM)) {
			freqM <- read.delim( freqM, as.is=T)
		}
		# make sure we don't have the depth info
		fm <- freqM[ , grep( "_Freq$", colnames(freqM))]
		dm <- freqM[ , grep( "_Depth$", colnames(freqM))]
		colnames(fm) <- sub( "_Freq$", "", colnames(fm))
		colnames(dm) <- sub( "_Depth$", "", colnames(dm))
		# set up the location info
		fmKey <- rownames(fm)
		fmTerms <- base::strsplit( fmKey, split="::")
		fmSID <- sapply( fmTerms, `[`, 1)
		fmGID <- sapply( fmTerms, `[`, 2)
		fmPOS <- sapply( fmTerms, `[`, 3)
		#fmBASE <- sapply( fmTerms, `[`, 4)
		sidToFMptr <- match( sampleIDset, colnames(fm), nomatch=0)
		if ( any( sidToFMptr == 0)) {
			cat( "\nError:  given BaseFreqMatrix does not contain all SampleIDs.")
			cat( "\nMissing:  ", sampleIDset[ sidToFMptr == 0])
			cat( "\nIgnoring 'freqM' argument..")
			useFM <- FALSE
		} else {
			useFM <- TRUE
		}
	}

	# go gather all the VCF summary calls over a set of samples
	vcfSet <- list()
	N <- 0
	badGenes <- vector()

	for ( sid in sampleIDset) { 
		vcf.path <- file.path( results.path, "VariantCalls", sid)
		f <- file.path( vcf.path, paste( sid, prefix, "Summary.VCF.txt", sep=".")) 
		if ( ! file.exists( f)) {
			cat( "\nVariant Call file not found.  Skip:  ", f)
			next
		}
		sml <- read.delim( f, as.is=T) 
		if (verbose) cat( "\nRead in: ", sid, "  N_SNP:", nrow(sml))

		# small chance of having more than one call at a location
		if ( any( duplicated( sml$POSITION))) {
			ord <- order( sml$SEQ_ID, sml$POSITION, -sml$SCORE)
			sml <- sml[ ord, ]
			key <- paste( sml$GENE_ID, sml$POSITION, sep="::")
			drops <- which( duplicated(key))
			if ( length( drops)) {
				sml <- sml[ -drops, ]
				if (verbose) cat( "  Drop Duplicates:", length(drops))
			}
		}

		# new idea:  don't frop them now, wait until we pair them up?
		if (DROP_AT_READ <- FALSE) {
		# perhaps drop low quality/score
		if ( ! is.null( min.qual)) {
			drops <- which( sml$QUAL < min.qual)
			if ( length( drops)) {
				sml <- sml[ -drops, ]
				if (verbose) cat( "  Drop Low Qual:", length(drops))
			}
		}
		if ( ! is.null( min.depth)) {
			drops <- which( sml$DEPTH < min.depth)
			if ( length( drops)) {
				sml <- sml[ -drops, ]
				if (verbose) cat( "  Drop Low Depth:", length(drops))
			}
		}
		if ( ! is.null( min.score)) {
			drops <- which( sml$SCORE < min.score)
			if ( length( drops)) {
				sml <- sml[ -drops, ]
				if (verbose) cat( "  Drop Low Score:", length(drops))
			}
		}
		if ( ! is.null( max.pctRef) && "PCT_REF" %in% colnames(sml)) {
			drops <- which( sml$PCT_REF > max.pctRef)
			if ( length( drops)) {
				sml <- sml[ -drops, ]
				if (verbose) cat( "  Drop High PctRef:", length(drops))
			}
		}
		}
		# there is a format of Alternate base like "A,G" if more than one seen
		# It's not an indel.  So catch here
		sml$ALT_BASE <- sub( "\\,[ACGT,]+", "", sml$ALT_BASE)
		if ( dropIndels) {
			drops1 <- which( nchar( sml$REF_BASE) > 1)
			drops2 <- which( nchar( sml$ALT_BASE) > 1)
			drops <- sort( union( drops1, drops2))
			if ( length( drops)) {
				sml <- sml[ -drops, ]
				if (verbose) cat( "  Drop Indels:", length(drops))
			}
		}
		if ( dropIntergenics) {
			drops <- grep( "(ng)", sml$GENE_ID, fixed=T)
			if ( length( drops)) {
				sml <- sml[ -drops, ]
				if (verbose) cat( "  Drop Intergenics:", length(drops))
			}
		}
		if ( dropIntrons) {
			smlSIDs <- sort( unique( sml$SEQ_ID))
			drops <- vector()
			for ( s in smlSIDs) {
				myrows <- which( sml$SEQ_ID == s)
				myhits <- which( sml$POSITION[ myrows] %in% exonPosSets[[ s]])
				mydrops <- setdiff( myrows, myrows[myhits])
				drops <- c( drops, mydrops)
			}
			if ( length( drops)) {
				sml <- sml[ -drops, ]
				if (verbose) cat( "  Drop Introns:", length(drops))
			}
		}
		if ( ! is.null( geneSubset)) {
			drops <- which( ! ( sml$GENE_ID %in% geneSubset))
			if ( length( drops)) {
				sml <- sml[ -drops, ]
				if (verbose) cat( "  Drop Unwanted Genes:", length(drops))
			}
		}

		# perhaps exclude whole genes for being too variant?
		if ( ! is.null( max.SNPperGene)) {
			snpTbl <- table( sml$GENE_ID)
			flag <- which( snpTbl > max.SNPperGene)
			if ( length(flag)) {
				badGenes <- sort( unique( c( badGenes, names(snpTbl)[flag])))
				if (verbose) cat( "  Flag HighSNP Genes:", length(flag))
			}
		}
		if (verbose) cat( "  Final N_Variants:", nrow(sml))

		N <- N + 1
		vcfSet[[N]] <- sml
		names(vcfSet)[N] <- sid
	}
	if (verbose) cat( "\n")

	# report about bad genes
	if ( length( badGenes)) {
		badProds <- gene2Product( badGenes)
		cat( "\n\nIgnoring genes with excessive SNP sites: ", length( badGenes), "\nTop Culprits:\n")
		print( head( sort( table( badProds), decreasing=T), 10))
	}

	# OK, we have all the wanted SNP sites, get ready to compare them pairwise
	outDist <- outJaccard <- matrix( NA, nrow=N, ncol=N)
	colnames(outDist) <- colnames(outJaccard) <- names(vcfSet)
	rownames(outDist) <- rownames(outJaccard) <- names(vcfSet)

	# make a unique key string that fully captures the base call, so we can easily count matches & mismatches
	smlDiffDF <- NULL
	cat( "\nComparing ")
	for ( i in 1:N) {
		sml1 <- vcfSet[[i]]
		if ( length( badGenes)) {
			drops <- which( sml1$GENE_ID %in% badGenes)
			if ( length( drops)) sml1 <- sml1[ -drops, ]
		}
		key1 <- paste( sml1$SEQ_ID, sml1$GENE_ID, sml1$POSITION, sml1$ALT_BASE, sep="::")
		cat( " ", sampleIDset[i])
		for ( j in i:N) {
			sml2 <- vcfSet[[j]]
			if ( length( badGenes)) {
				drops <- which( sml2$GENE_ID %in% badGenes)
				if ( length( drops)) sml2 <- sml2[ -drops, ]
			}
			key2 <- paste( sml2$SEQ_ID, sml2$GENE_ID, sml2$POSITION, sml2$ALT_BASE, sep="::")
			myKey1 <- key1

			# lets use the FreqM if we have it, instead of the simplistic set intersection method
			either <- union( key1, key2)
			nEither <- length( either)
			# find where they live in the freqM
			if (useFM) {
				whereFM <- match( either, fmKey, nomatch=0)
				keep <- which( whereFM > 0)
				either <- either[ keep]
				whereFM <- whereFM[ keep]
				nEither <- length(keep)
				myKey1 <- key1[ key1 %in% either]
				key2 <- key2[ key2 %in% either]
				column1 <- sidToFMptr[i]
				column2 <- sidToFMptr[j]
				freq1 <- fm[ whereFM, column1]
				freq2 <- fm[ whereFM, column2]
				deep1 <- dm[ whereFM, column1]
				deep2 <- dm[ whereFM, column2]
				minDeep <- pmin( deep1, deep2)
				# given all freq values are in the range 0..100, let's require just over 50% difference
				# to call the site truly a different allele
				deltaF <- abs( freq1 - freq2)
				# only use sites where both have data
				keep1 <- which( ! is.na(deltaF))
				deltaF <- deltaF[ keep1]
				minDeep <- minDeep[ keep1]
				either <- either[ keep1]
				whereFM <- whereFM[ keep1]
				nEither <- length(keep1)
				# only keep those with enough read depth
				keep2 <- which( minDeep >= min.depth)
				deltaF <- deltaF[ keep2]
				minDeep <- minDeep[ keep2]
				either <- either[ keep2]
				whereFM <- whereFM[ keep2]
				nEither <- length(keep2)
				# now we are ready to do the match / no match question
				whoDiff <- which( deltaF >= min.freqDifference)
				nDiff <- length(whoDiff)
				nSame <- nEither - nDiff
				# lastly, use these values for the Jaccard, Distance
				nBoth <- nSame
				nOnlyOne <- nDiff
				jaccard <- nBoth / nEither
				# make a data.frame of who are the diffs
				if ( N == 2 && nDiff > 0) {
					diffPtr <- whereFM[ whoDiff]
					smlDiffDF <- data.frame( "SEQ_ID"=fmSID[diffPtr], "GENE_ID"=fmGID[diffPtr],
								"POSITION"=as.numeric(fmPOS[diffPtr]), "FREQ_Diff"=round(deltaF[whoDiff]),
								stringsAsFactors=F)
					ord <- order( -smlDIffDF$FREQ_Diff, smlDiffDF$SEQ_ID, smlDiffDF$POSITION)
					smlDiffDF <- smlDiffDF[ ord, ]
					rownames(smlDiffDF) <- 1:nrow(smlDiffDF)
				}
			} else {
				# Jaccard is intersection over union
				both <- intersect( key1, key2)
				nBoth <- length( both)
				onlyOne <- setdiff( either, both)
				nOnlyOne <- length(onlyOne)
				jaccard <- nBoth / nEither
			}
			# special case if no SNPs at all
			if (nEither == 0) jaccard <- 1

			outDist[i,j] <- outDist[j,i] <- nOnlyOne
			outJaccard[i,j] <- outJaccard[j,i] <- jaccard

			if ( ! is.null( debugPattern)) {
				if ( i != j && grepl(debugPattern, names(vcfSet)[i]) && grepl(debugPattern, names(vcfSet)[j])) {
					cat( "\nDebug:   ", i, j, "|", length(key1), length(key2), "|", nBoth, nEither, nOnlyOne)
					cat( "\nDebug I: ", head( key1, 100))
					cat( "\nDebug J: ", head( key2, 100))
					cat( "\nDebug ?: ", head( onlyOne, 50))
				}
			}
		}
	}

	if ( is.null(smlDiffDF)) {
		return( list( "Distance"=outDist, "Jaccard"=round(outJaccard, digits=3)))
	} else {
		return( list( "Distance"=outDist, "Jaccard"=round(outJaccard, digits=3), "Differences"=smlDiffDF))
	}
}

