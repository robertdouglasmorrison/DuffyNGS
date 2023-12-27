# pipe.BAMproteins.R -- get the consensus proteins directly from the BAM file results


# convert a BAM file into a FASTA of all proteins

`pipe.BAMproteins.CreateFasta` <- function( sampleID, geneIDset=NULL, annotationFile="Annotation.txt",
				optionsFile="Options.txt", speciesID=getCurrentSpecies(), results.path=NULL,
				noReadCalls=NULL, SNP.only=FALSE, minReadCalls=NULL, minPercentSNP=NULL, 
				makeFastaFile=TRUE, verbose=TRUE) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	if ( is.null( noReadCalls)) {
		dataType <- getAnnotationValue( annotationFile, key=sampleID, columnArg="DataType", notfound="DNA-seq", verbose=verbose)
		if (dataType == "DNA-seq") noReadCalls <- "blank"
		if (dataType == "RNA-seq") noReadCalls <- "genomic"
	}
	if ( is.null(noReadCalls) || !(noReadCalls %in% c("blank","genomic"))) {
		cat( "\nArgument 'noReadCalls' must be one of 'blank' or 'genomic'")
		stop( "See command 'pipe.BAMproteins.CreateFasta()'")
	}

	# make sure we have the BAM file already sorted
	bamfile <- paste( sampleID, "genomic.bam", sep=".")
	bamfile <- file.path( results.path, "align", bamfile)
	sortedbamfile <- BAM.verifySorted( bamfile, index=TRUE)
	if ( is.null( sortedbamfile)) return(NULL)

	if ( speciesID != getCurrentSpecies()) setCurrentSpecies( speciesID)
	geneMap <- getCurrentGeneMap()
	cdsMap <- getCurrentCdsMap()
	allGenes <- sort( unique( cdsMap$GENE_ID))
	if ( is.null( geneIDset)) {
		geneIDset <- allGenes
	} else {
		geneIDset <- intersect( geneIDset, allGenes)
	}
	if ( length( geneIDset) < 1) {
		cat( "\nNo Genes selected")
		return(NULL)
	}

	require( Biostrings)

	# get the genome into vectors of bases
	genomicFastaFile <- getOptionValue( optionsFile, "genomicFastaFile", verbose=verbose)
	fa <- loadFasta( genomicFastaFile, verbose=verbose)
	baseVectors <- strsplit( fa$seq, split="")
	names(baseVectors) <- fa$desc

	MATCH <- base::match
	NCHAR <- base::nchar
	PASTE <- base::paste
	WHICH <- base::which

	`getProteinOneGene` <- function( gid) {

		where <- MATCH( gid, geneMap$GENE_ID)
		myStart <- geneMap$POSITION[where]
		myStop <- geneMap$END[where]
		myStrand <- geneMap$STRAND[where]
		mySID <- geneMap$SEQ_ID[where]
		
		myBaseVecPt <- MATCH( mySID, names(baseVectors))

		ans <- pipe.ConsensusBaseCalls( sampleID, geneID=gid, seqID=mySID, start=myStart, stop=myStop, 
					annotationFile=annotationFile, genomicFastaFile=genomicFastaFile,
					genomicVector=baseVectors[[myBaseVecPt]],
					optionsFile=optionsFile, results.path=results.path, noReadCalls=noReadCalls,
					aaToo=TRUE, as.cDNA=TRUE, best.frame=FALSE, SNP.only=SNP.only, 
					minReadCalls=minReadCalls, minPercentSNP=minPercentSNP, verbose=FALSE)

		# by default, indels can hose the translation into proteins
		# do we want to prevent that?
		if ( SNP.only) {
			# indels will be detectable as elements that are not a single character long.
			# to prevent them from causing trouble, let's find them and replace with the reference
			dna <- ans$dna.consensus
			isIndel <- WHICH( nchar(dna) != 1)
			if ( length( isIndel)) {
				refDNA <- ans$ref[ isIndel]
				refDNA[ nchar(refDNA) != 1] <- "N"
				dna[ isIndel] <- refDNA
			}
			dna <- PASTE( dna, collapse="")
			myProt <- DNAtoAA( dna, clipAtStop=F, readingFrame=1)
		} else {
			myAA <- ans$aa.consensus
			myProt <- PASTE( myAA, collapse="")
		}

		# the consensus now returns a confidence score
		conf <- if ( is.null( ans$aa.confidence)) 0 else round( as.numeric( ans$aa.confidence) * 100, digits=2)

		if (verbose) cat( "\r", gid, NCHAR(myProt), conf, parallel:::isChild())
		return( list( "seq"=myProt, "confidence"=conf))
	}
	

	if (verbose) cat( "\nExtracting proteins from BAM consensus:  N_Genes =", length(geneIDset), "\n")
	ans <- multicore.lapply( geneIDset, FUN=getProteinOneGene)
	#if ( exists( "MCLAPPLY_DEBUG")) rm( MCLAPPLY_DEBUG, inherits=T)
	
	if ( length(geneIDset) > 1) {
		allProts <- sapply( ans, function(x) if (is.null(x)) NA else x[[1]])
		allConf <- sapply( ans, function(x) if (is.null(x) || length(x) < 2) NA else x[[2]])
	} else {
		allProts <- ans[[1]]
		allConf <- ans[[2]]
	}

	out <- data.frame( "GENE_ID"=geneIDset, "Confidence"=allConf, "Protein"=allProts, stringsAsFactors=F)
	rownames(out) <- 1:nrow(out)

	if (makeFastaFile) {
		genes <- out$GENE_ID
		prods <- gene2Product( genes)
		confs <- out$Confidence
		seqs <- out$Protein
		desc <- paste( genes, " | ", "Confidence=", confs, " | ", prods, sep="")
		fa <- as.Fasta( desc, seqs)
		fa.path <- file.path( results.path, "ConsensusProteins", sampleID)
		if ( ! file.exists( fa.path)) dir.create( fa.path, recursive=T)
		fa.file <- file.path( fa.path, paste( "All.BAM.Proteins", sampleID, "fasta", sep="."))
		writeFasta( fa, fa.file, line.width=100)
		if (verbose) cat( "\nWrote Fasta Protein file: ", fa.file, "\n")
		return( invisible( out))
	} else {
		return( out)
	}
}


# summarize BAM protein differences between 2 samples

`pipe.BAMproteins.SampleCompare` <- function( sampleID1, sampleID2=NULL, geneIDset=NULL, annotationFile="Annotation.txt",
				optionsFile="Options.txt", speciesID=getCurrentSpecies(), results.path=NULL, dropGenes=NULL,
				min.confidence=40, min.editDist=1, nWorst=20, show.details=T, nMutations=10, 
				folder=NULL, otherIDs=NULL, verbose=T, ...) {

	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	if (speciesID != getCurrentSpecies()) setCurrentSpecies( speciesID)

	# grab 2 FASTA of BAM file consensus proteins
	f1 <- paste( "All.BAM.Proteins", sampleID1, "fasta", sep=".")
	f1 <- file.path( results.path, "ConsensusProteins", sampleID1, f1)
	if ( ! file.exists( f1)) {
		cat( "\nExpected FASTA file not found: ", f1)
		cat( "\nPerhaps run 'pipe.BAMproteins.CreateFasta()' first..")
		return( NULL)
	}

	# the second file can be a sample, or NULL means the reference genome
	if ( is.null( sampleID2)) {
		sampleID2 <- "Reference"
		f2 <- paste( speciesID, "ReferenceProteins.fasta", sep=".")
		hasConfidence2 <- FALSE
	} else {
		f2 <- paste( "All.BAM.Proteins", sampleID2, "fasta", sep=".")
		f2 <- file.path( results.path, "ConsensusProteins", sampleID2, f2)
		hasConfidence2 <- TRUE
	}
	if ( ! file.exists( f2)) {
		if ( sampleID2 == "Reference") {
			cat( "\nMaking FASTA of Reference Proteins..\n")
			genomicFastaFile <- getOptionValue( optT, "genomicFastaFile", notfound="Pf_genomicDNA.fasta", verbose=F)
			allGenes <- unique( getCurrentCdsMap()$GENE_ID)
			fa <- gene2Fasta( allGenes, genomicFastaFile, mode="aa", verbose=T)
			writeFasta( fa, f2, line.width=100)
			cat( "\nDone.\n")
		} else {
			cat( "\nExpected FASTA file not found: ", f2)
			cat( "\nPerhaps run 'pipe.BAMproteins.CreateFasta()' first..")
			return( NULL)
		}
	}
	fa1 <- loadFasta( f1, short=F, verbose=F)
	fa2 <- loadFasta( f2, short=F, verbose=F)

	NCHAR <- base::nchar
	PASTE <- base::paste
	STRSPLIT <- base::strsplit
	SAPPLY <- base::sapply

	# the descriptors may have confidence scores and products
	desc1Terms <- STRSPLIT( fa1$desc, split=" | ", fixed=T)
	desc1 <- SAPPLY( desc1Terms, `[`, 1)
	conf1 <- as.numeric( sub( "Confidence=", "", SAPPLY( desc1Terms, `[`, 2)))
	if (hasConfidence2) {
		desc2Terms <- STRSPLIT( fa2$desc, split=" | ", fixed=T)
		desc2 <- SAPPLY( desc2Terms, `[`, 1)
		conf2 <- as.numeric( sub( "Confidence=", "", SAPPLY( desc2Terms, `[`, 2)))
	} else {
		desc2Terms <- STRSPLIT( fa2$desc, split=" | ", fixed=T)
		desc2 <- SAPPLY( desc2Terms, `[`, 1)
		conf2 <- rep.int( 100, length(desc2))
	}

	# decide which genes can be compared
	gids <- sort( intersect( desc1, desc2))
	if ( ! is.null( geneIDset)) gids <- intersect( gids, geneIDset)
	NG <- length( gids)
	if (verbose) cat( "\nStarting difference search   N_Genes: ", NG)

	# allow use of confidence scores too
	good1 <- desc1[ conf1 >= min.confidence]
	good2 <- desc2[ conf2 >= min.confidence]
	goodConf <- intersect( good1, good2)
	gids <- intersect( gids, goodConf)
	NG <- length( gids)
	if (verbose) cat( "\nAfter Confidence cutoff=", min.confidence, "   N_Genes:  ", NG, sep="")

	if ( ! is.null( dropGenes)) {
		gids <- setdiff( gids, dropGenes)
		NG <- length( gids)
		if (verbose) cat( "\nAfter explicit exclusions    N_Genes: ", NG)
	}
	wh1 <- match( gids, desc1)
	wh2 <- match( gids, desc2)
	prot1 <- fa1$seq[ wh1]
	prot2 <- fa2$seq[ wh2]

	# to prevent any runtime erorrs, catch/fix any invalid AA calls here
	prot1 <- gsub( "?", "X", prot1, fixed=T)
	prot2 <- gsub( "?", "X", prot2, fixed=T)
	prot1 <- gsub( "U", "X", prot1, fixed=T)
	prot2 <- gsub( "U", "X", prot2, fixed=T)

	# find out how many differences
	ed <- td <- conf <- nistop <- dSize <- pctDiff <- ngap <- vector( length=NG)
	details <- rep.int( "", NG)

	require( Biostrings)
	data( BLOSUM62)

	if (verbose) cat( "\n")
	for( i in 1:NG) {
		s1 <- prot1[i]
		s2 <- prot2[i]
		# gracefully catch completed deleted proteins
		if ( s1 == "") s1 <- "*"
		if ( s2 == "") s2 <- "*"
		nc <- NCHAR( c( s1, s2))
		td[i] <- l <- max( nc)
		if ( s1 == s2) {
			ed[i] <- d <- 0
		} else {
			ed[i] <- d <- adist( s1, s2)[1]
		}
		conf[i] <- min( conf1[wh1[i]], conf2[wh2[i]])
		dSize[i] <- nc[1] - nc[2]
		pctDiff[i] <- d * 100 / l
		if (show.details && d > 0) {
			pa <- pairwiseAlignment( s1, s2, type="global", substitutionMatrix=BLOSUM62)
			ch1 <- STRSPLIT( as.character( alignedPattern(pa)), split="")[[1]]
			ch2 <- STRSPLIT( as.character( alignedSubject(pa)), split="")[[1]]
			diffs <- which( ch1 != ch2)
			if ( length(diffs)) {
				deltas <- PASTE( ch2[diffs], as.character(diffs), ch1[diffs], sep="")
				if ( length( deltas) > nMutations) deltas <- c( deltas[1:nMutations], "...")
				details[i] <- PASTE( deltas, collapse="; ")
				nistop[i] <- sum( ch1[diffs] == "*")
				ngap[i] <- sum( ch1[diffs] == "-")
			}
		} else {
			details[i] <- ""
			# second method for counting internal stops when not doing PA
			nistop[i] <- 0
			ngap[i] <- 0
			if ( d > 0) {
				chV <- STRSPLIT( s1, split="")[[1]]
				lm1 <- length(chV) - 1
				if ( lm1) {
					nistop[i] <- sum( chV[1:lm1] == "*")
					ngap[i] <- sum( chV[1:lm1] == "-")
				}
			}
		}
		if ( verbose && i %% 100 == 0) cat( "\r", i, gids[i], l, d, round(pctDiff[i],digits=2))
	}
	if (verbose) cat( "\nDone.\n") else cat( " ", sampleID1)

	# total it up
	ted <- sum( ed, na.rm=T)
	ttd <- sum( td, na.rm=T)

	nWorst <- min( nWorst, NG)
	worst <- order( pctDiff, ed, decreasing=T)[1:nWorst]
	worstDF <- data.frame( "GENE_ID"=gids[worst], "Length"=td[worst], "EditDistance"=ed[worst], 
			"DeltaLength"=formatC( dSize[worst], format="d", flag="+"), 
			"Pct_Different"=round( pctDiff[worst], digits=2), "InternalStopCodons"=nistop[worst], 
			"Indel_Gaps"=ngap[worst], "Confidence"=conf[worst], "PRODUCT"=gene2Product(gids[worst]), 
			stringsAsFactors=F)
	if (show.details) worstDF <- cbind( worstDF, "MutationDetails"=details[worst], stringsAsFactors=F)
	worstDF <- subset( worstDF, EditDistance >= min.editDist)
	if (nrow(worstDF)) rownames(worstDF) <- 1:nrow(worstDF)

	out <- list( "SampleID"=sampleID1, "Comparitor"=sampleID2, "N_Genes"=NG, 
			"N_ExactMatch"=sum( ed == 0, na.rm=T), "N_Different"=sum( ed > 0, na.rm=T), 
			"TotalEditDist"=ted, "AvgEditDist"=round( mean(ed,na.rm=T), digits=2), 
			"Overall_PctDifferent"=round( ted*100/ttd, digits=3), "Worst.Genes"=worstDF)

	# call the plotter function too?
	if ( ! is.null( folder)) {
		plotBAMproteinDifferences( out, folder=folder, otherIDs=otherIDs,
				results.path=results.path, optionsFile=optionsFile, ...)
	}

	return( out)
}


`plotBAMproteinDifferences` <- function( differences, folder=paste("Top.BAM.Protein.Differences", differences$SampleID, sep="_"), 
					otherIDs=NULL, results.path=NULL, optionsFile="Options.txt", ...) {

	# create a HTML file with plot links to images of the biggest protein differences

	# given the list of results from a call to 'pipe.BAMproteinDifference()'
	if ( ! ("Worst.Genes" %in% names(differences))) {
		cat( "\nError: Difference object does not contain expected 'Worst.Genes' data frame.")
		return( NULL)
	}
	sampleID1 <- differences$SampleID
	sampleID2 <- differences$Comparitor

	# create the folders to hold these results
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	}
	out.path <- file.path( results.path, "ConsensusProteins", folder)
	if ( ! file.exists( out.path)) dir.create( out.path, recursive=T)
	png.path <- file.path( out.path, "SNP_Plots")
	if ( ! file.exists( png.path)) dir.create( png.path, recursive=T)
	htmlFile <- paste( "BAM.Protein.Differences_", sampleID1, ".v.", sampleID2, ".html", sep="")
	htmlFile <- file.path( out.path, htmlFile)
	textFile <- paste( "BAM.Protein.Differences_", sampleID1, ".v.", sampleID2, ".txt", sep="")
	textFile <- file.path( out.path, textFile)

	# extract what we want to make the HTML table
	tbl <- differences$Worst.Genes
	if ( ! nrow(tbl)) {
		cat( "\nNo Proteins flagged as different...")
		return(NULL)
	}

	# re-arrange the columns and rename a bit
	tbl <- tbl[ , c(1,9,2,5,3,6,7,10)]
	colnames(tbl) <- c( "GENE_ID", "PRODUCT", "N_AA", "Pct Different", "Edit Dist", "Internal Stop Codons",
				"Indel Gaps", "MutationDetails")

	# put some of the other metrics into the title
	titleString <- paste( "Top BAM Protein Differences <br> Sample: &nbsp; ", sampleID1, "<br>",
				"Comparitor Genome: &nbsp; ", sampleID2, "<br>",
				"N_Proteins that Exactly Match: &nbsp; ", differences$N_ExactMatch, "<br>",
				"N_Proteins with Mis-Matches: &nbsp; ", differences$N_Different, "<br>",
				"Average A.A. Edits per Protein: &nbsp; ", differences$AvgEditDist)


	cat( "\nWriting BAM Protein Difference result files to: ", out.path)
	write.table( tbl, textFile, sep="\t", quote=F, row.names=F)

	# tweak the gene names to see the mutation text too
	htmlTbl <- tbl
	firstMutation <- sub( "; .+", "", tbl$MutationDetails)
	firstMutation <- sub( "*", "stop", firstMutation, fixed=T)
	firstMutation <- sub( "-", "minus", firstMutation, fixed=T)
	htmlTbl$GENE_ID <- paste( tbl$GENE_ID, firstMutation, sep=".")
	table2html( htmlTbl, htmlFile, title=titleString, linkPaths="SNP_Plots")

	# now make SNP plots for the first mutation in each gene
	# get the set of IDs to plot, we may have been given 'others'...
	sidSet <- sampleID1
	if ( ! is.null( otherIDs)) sidSet <- unique( c( sidSet, otherIDs))

	# visit each entry and make one plot
	checkX11()
	gmap <- getCurrentGeneMap()
	where <- match( tbl$GENE_ID, gmap$GENE_ID)
	visitOrder <- order( where)
	cat( "\n")
	for ( i in 1:nrow(tbl)){
		ii <- visitOrder[i]
		thisGene <- tbl$GENE_ID[ii]
		# grab the AA position of the first mutation
		firstMutation <- sub( "; .+", "", tbl$MutationDetails[ii])
		aaPosition <- as.integer( gsub( "[\\*\\-]", "", gsub( "[A-Z]", "", firstMutation)))
		#cat( "\nDebug: ", i, thisGene, firstMutation, "|", aaPosition)
		if ( is.na( aaPosition)) {
			cat( "\nError:  unable to parse integer location from mutation string: ", firstMutation, "  Skip..")
			next
		}
		# turn this into into genomic location
		genomeAns <- convertAApositionToGenomicDNAposition( thisGene, aaPosition)
		thisSeqID <- genomeAns$SEQ_ID
		thisPos <- genomeAns$SEQ_POSITION + 1
		# OK, plot that SNP site
		plotFile <- paste( thisGene, firstMutation, "png", sep=".")
		# trap stop codons
		plotFile <- sub( "*", "stop", plotFile, fixed=T)
		plotFile <- sub( "-", "minus", plotFile, fixed=T)
		pipe.PlotSNP( sidSet, seqID=thisSeqID, pos=thisPos, geneID=thisGene,
				results.path=results.path, optionsFile=optionsFile,
				plotFormat="png", plotFileName=plotFile, plot.path=png.path, ...)
		cat( "\r", i, thisGene, firstMutation)
	}
}


`pipe.BAMproteins.GroupCompare` <- function( sampleIDset, groupSet, geneIDset=NULL, annotationFile="Annotation.txt",
				optionsFile="Options.txt", speciesID=getCurrentSpecies(), results.path=NULL,
				comparison=c("reference","pairwise"), tbl=NULL, dropGenes=NULL, min.confidence=40,
				wt.estimate=1, wt.pvalue=1, folder=NULL, Ngenes=50, verbose=T) {

	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	if (speciesID != getCurrentSpecies()) setCurrentSpecies( speciesID)
	gmap <- getCurrentGeneMap()
	NG <- nrow(gmap)
	NS <- length( sampleIDset)

	comparison <- match.arg( comparison)

	# we (for now) only handle 2 groups
	if ( length( groupSet) != NS) stop( paste( "Must provide group for every sample"))
	grpLevels <- sort( unique( groupSet))
	if ( length( grpLevels) != 2) stop( paste( "Must provide samples from exactly 2 groups"))

	# in most cases, we need to build the giant table of all comparisons
	if ( is.null( tbl)) {
		bigDF <- data.frame()
		if ( comparison == "reference") {
			cat( "\nBuilding protein Difference data for", length(sampleIDset), "samples against reference proteome..")
			cat( "\nUsing multicore..")
			mcAns <- multicore.lapply( sampleIDset, FUN=pipe.BAMprotein.Difference, sampleID2=NULL, optionsFile=optionsFile,
							results.path=results.path, geneIDset=geneIDset, dropGenes=dropGenes, 
							min.confidence=0, min.editDist=0, nWorst=NG, show.details=T, verbose=F)
			cat( "\nMerging..")
			for (i in 1:NS) {
				s <- sampleIDset[i]
				g <- groupSet[i]
				ans <- mcAns[[i]]
				sml <- data.frame( "SampleID"=s, "Group"=g, ans$Worst.Genes, stringsAsFactors=F)
				bigDF <- rbind( bigDF, sml)
				cat( " ", s, sep="")
			}
		} else {
			cat( "\nBuilding protein Difference data between", length(sampleIDset), "samples against each other..")
			visitOrd <- order( groupSet)
			for (i in 1:(NS-1)) {
				ii <- visitOrd[i]
				s1 <- sampleIDset[ii]
				g1 <- groupSet[ii]
				for (j in (i+1):NS) {
					jj <- visitOrd[j]
					s2 <- sampleIDset[jj]
					g2 <- groupSet[jj]
					ans <- pipe.BAMprotein.Difference( s1, sampleID2=s2, optionsFile=optionsFile,
							results.path=results.path, geneIDset=geneIDset, dropGenes=dropGenes, 
							min.confidence=0, min.editDist=0, nWorst=NG, show.details=T, verbose=F)
					sml <- data.frame( "SampleID"=s1, "Group"=g1, "Sample2"=s2, "Group2"=g2, 
							ans$Worst.Genes, stringsAsFactors=F)
					bigDF <- rbind( bigDF, sml)
					cat( ".")
				}
			}
		} 
		cat( "  Done.\n")
		if ( "Group2" %in% colnames(bigDF)) {
			bigDF$Group <- paste( bigDF$Group, bigDF$Group2, sep=".")
		}
		tbl <- bigDF
	}

	# with all the details in place, before doing the comparison steps, see what genes should get dropped
	# due to the confidence cutoff
	geneMinConf <- tapply( tbl$Confidence, factor( tbl$GENE_ID), min, na.rm=T)
	drops <- names(geneMinConf)[ as.numeric(geneMinConf) < min.confidence]
	if ( length(drops)) {
		dropRows <- which( tbl$GENE_ID %in% drops)
		tbl <- tbl[ -dropRows, ]
		if (verbose) cat( "\nDropping", length(drops), "genes due to minimum confidence cutoff of", min.confidence)
	}

	# with giant table in hand, now do all the comparisons
	gFac <- factor( tbl$GENE_ID)
	NG <- nlevels(gFac)
	N <- nrow(tbl)
	grpLevels <- sort( unique( tbl$Group))

	# storage for what we find
	gOut <- tOut <- eOut <- pOut <- v1Out <- v2Out <- vector( length=NG*6)
	nout <- 0

	# visit every gene
	cat( "\nModeling", NG, "genes against 6 protein difference metrics..\n")
	tapply( 1:N, gFac, function(x) {
		myG <- tbl$GENE_ID[ x[1]]
		myGrps <- factor(tbl$Group[x])
		if ( length(x) < 2) return()
		if ( nlevels(myGrps) < 2) return()
	
		# do the tests we can
		ed <- as.numeric( tbl$EditDistance[x])
		if ( ! all( ed == 0)) {
			lmAns <- summary( lm( ed ~ myGrps))
			cf <- coefficients(lmAns)
			est1 <- cf[ 1, 1]
			est2 <- cf[ 2, 1]
			pv <- cf[ 2, 4]
			nout <<- nout + 1
			gOut[nout] <<- myG
			tOut[nout] <<- "EditDistance"
			eOut[nout] <<- est2
			pOut[nout] <<- pv
			v1Out[nout] <<- est1
			v2Out[nout] <<- est1 + est2
		}

		dl <- as.numeric( tbl$DeltaLength[x])
		if ( ! all( dl == 0)) {
			lmAns <- summary( lm( dl ~ myGrps))
			cf <- coefficients(lmAns)
			est1 <- cf[ 1, 1]
			est2 <- cf[ 2, 1]
			pv <- cf[ 2, 4]
			nout <<- nout + 1
			gOut[nout] <<- myG
			tOut[nout] <<- "DeltaLength"
			eOut[nout] <<- est2
			pOut[nout] <<- pv
			v1Out[nout] <<- est1
			v2Out[nout] <<- est1 + est2
		}

		pd <- as.numeric( tbl$Pct_Different[x])
		if ( ! all( pd == 0)) {
			lmAns <- summary( lm( pd ~ myGrps))
			cf <- coefficients(lmAns)
			est1 <- cf[ 1, 1]
			est2 <- cf[ 2, 1]
			pv <- cf[ 2, 4]
			nout <<- nout + 1
			gOut[nout] <<- myG
			tOut[nout] <<- "PctDifferent"
			eOut[nout] <<- est2
			pOut[nout] <<- pv
			v1Out[nout] <<- est1
			v2Out[nout] <<- est1 + est2
		}

		sc <- as.numeric( tbl$InternalStopCodons[x])
		if ( ! all( sc == 0)) {
			lmAns <- summary( lm( sc ~ myGrps))
			cf <- coefficients(lmAns)
			est1 <- cf[ 1, 1]
			est2 <- cf[ 2, 1]
			pv <- cf[ 2, 4]
			nout <<- nout + 1
			gOut[nout] <<- myG
			tOut[nout] <<- "InternalStopCodons"
			eOut[nout] <<- est2
			pOut[nout] <<- pv
			v1Out[nout] <<- est1
			v2Out[nout] <<- est1 + est2
		}
		
		ig <- as.numeric( tbl$Indel_Gaps[x])
		if ( ! all( ig == 100)) {
			lmAns <- summary( lm( ig ~ myGrps))
			cf <- coefficients(lmAns)
			est1 <- cf[ 1, 1]
			est2 <- cf[ 2, 1]
			pv <- cf[ 2, 4]
			nout <<- nout + 1
			gOut[nout] <<- myG
			tOut[nout] <<- "Indel_Gaps"
			eOut[nout] <<- est2
			pOut[nout] <<- pv
			v1Out[nout] <<- est1
			v2Out[nout] <<- est1 + est2
		}
		
		cd <- as.numeric( tbl$Confidence[x])
		if ( ! all( cd == 100)) {
			lmAns <- summary( lm( cd ~ myGrps))
			cf <- coefficients(lmAns)
			est1 <- cf[ 1, 1]
			est2 <- cf[ 2, 1]
			pv <- cf[ 2, 4]
			nout <<- nout + 1
			gOut[nout] <<- myG
			tOut[nout] <<- "Confidence"
			eOut[nout] <<- est2
			# the P values of the confidence are artifically too good.  LM is biased by their tight variance
			# use Mann Whitney instead
			#pOut[nout] <<- pv
			pOut[nout] <<- wilcox.test( cd ~ myGrps)$p.value
			v1Out[nout] <<- est1
			v2Out[nout] <<- est1 + est2
		}
		
		cat( "\r", myG, nout)
		return(NULL)
	})
	# now set the true length
	cat( "\nDone.\n")
	length(gOut) <- length(tOut) <- length(eOut) <- length(pOut) <- length(v1Out) <- length(v2Out) <- nout
	pOut[ is.nan(pOut)] <- 1
	pOut[ is.na(pOut)] <- 1

	# package up the results
	out1 <- data.frame( "GENE_ID"=gOut, "PRODUCT"=gene2Product(gOut), "Difference_Metric"=tOut, 
			"P_VALUE"=pOut, "LM_Estimate"=round(eOut,digits=4), "Value1"= round(v1Out,digits=2), 
			"Value2"=round(v2Out,digits=2), stringsAsFactors=F)
	colnames(out1)[ 6:7] <- grpLevels[1:2]

	# sort by P
	ord <- diffExpressRankOrder( abs(out1$LM_Estimate), out1$P_VALUE, wt.fold=wt.estimate, wt.pvalue=1)
	out1 <- out1[ ord, ]
	rownames(out1) <- 1:nrow(out1)

	# also a meta result for each gene, that shows how it fared by all 4 metrics 
	gset <- sort( unique( out1$GENE_ID))
	NG2 <- length(gset)
	edSub <- subset.data.frame( out1, Difference_Metric == "EditDistance")
	whED <- match( gset, edSub$GENE_ID)
	dlSub <- subset.data.frame( out1, Difference_Metric == "DeltaLength")
	whDL <- match( gset, dlSub$GENE_ID)
	pdSub <- subset.data.frame( out1, Difference_Metric == "PctDifferent")
	whPD <- match( gset, pdSub$GENE_ID)
	scSub <- subset.data.frame( out1, Difference_Metric == "InternalStopCodons")
	whSC <- match( gset, scSub$GENE_ID)
	igSub <- subset.data.frame( out1, Difference_Metric == "Indel_Gaps")
	whIG <- match( gset, igSub$GENE_ID)
	cdSub <- subset.data.frame( out1, Difference_Metric == "Confidence")
	whCD <- match( gset, cdSub$GENE_ID)

	wh <- matrix( c(whED,whDL,whPD,whSC,whIG,whCD), nrow=NG2, ncol=6)
	colnames(wh) <- c( "EditDist", "DeltaLen", "PctDiff", "StopCodons", "IndelGaps", "Confidence")
	# how to score/rank those that are missing...
	wh[ is.na(wh)] <- NG2 + 1
	#wh[ is.na(wh)] <- NA
	avg <- apply( wh, 1, logmean)
	out2 <- data.frame( "GENE_ID"=gset, "PRODUCT"=gene2Product(gset), "AVG_RANK"=avg, wh, stringsAsFactors=F)
	ord <- order( out2$AVG_RANK)
	out2 <- out2[ ord, ]
	rownames(out2) <- 1:nrow(out2)
	out2 <- addNameIDterms( out2)

	finalAns <- list( "difference.details"=tbl, "model.results"=out1, "meta.ranks"=out2)

	# allow the tool to write results and plots to a subfolder
	if ( ! is.null( folder)) {
		out.path <- file.path( results.path, "ConsensusProteins", paste( "BAM.Protein.Compare", folder, sep="_"))
		if ( ! file.exists( out.path)) dir.create( out.path, recursive=T)
		cat( "\nWriting BAM protein results files to: ", out.path)
		f1 <- file.path( out.path, paste( folder, "DifferenceDetails.txt", sep="."))
		write.table( tbl, f1, sep="\t", quote=F, row.names=F)
		f2 <- file.path( out.path, paste( folder, "ModelResults.txt", sep="."))
		write.table( out1, f2, sep="\t", quote=F, row.names=F)
		f3 <- file.path( out.path, paste( folder, "MetaRanks.txt", sep="."))
		write.table( out2, f3, sep="\t", quote=F, row.names=F)
		# also create a HTML of that final meta ranks answer, with plots
		f3 <- file.path( out.path, paste( folder, "MetaRanks.html", sep="."))
		localLinkPath <- "pngPlots"
		globalLinkPath <- file.path( out.path, localLinkPath)
		if ( ! file.exists( globalLinkPath)) dir.create( globalLinkPath, recursive=T)
		table2html( out2, f3, title=paste( "Top Genes in 'BAM Protein Group Compare': ", "&nbsp; &nbsp;", folder),
			maxRows=Ngenes, linkPaths=localLinkPath)
		# lastly, make those plots
		nShow <- min( nrow(out2), Ngenes)
		genesToPlot <- out2$GENE_ID[ 1:nShow]
		cat( "\nMaking BAM protein comparison plots..\n")
		for ( g in genesToPlot) {
			plotBAMproteinDifferenceOneGene( g, tbl=tbl, show.labels=T)
			plotfile <- file.path( globalLinkPath, paste( g, "png", sep="."))
			dev.print( png, plotfile, width=1200)
			cat( "\r", g)
		}
		cat( "\nDone.\n")
	}

	return( finalAns)
}


`plotBAMproteinDifferenceOneGene` <- function( geneID, tbl, show.labels=T) {

	# given the large table of BAM protein difference details
	if ( names(tbl)[1] == "difference.details") {
		tbl <- tbl[[1]]
	}

	if ( ! all( c( "GENE_ID", "Group", "EditDistance", "DeltaLength", "Pct_Different", 
			"InternalStopCodons", "Indel_Gaps", "Confidence") %in% colnames(tbl))) {
		cat( "\nIncorrect column names in BAM protein difference details object...")
		return( NULL)
	}

	sml <- subset.data.frame( tbl, GENE_ID == geneID)
	if ( !nrow(sml)) {
		cat( "\nNo details data found for gene: ", geneID)
		return(NULL)
	}
	prod <- sml$PRODUCT[1]
	grps <- factor( sml$Group)
	grpPos <- as.integer(grps) * 2
	labels <- sml$SampleID
	NGrp <- nlevels(grps)
	colSet <- rainbow( NGrp, end=0.6)

	# OK, make some plots to show what we see for this one gene
	checkX11()
	savMAI <- par( "mai")
	on.exit( par( "mai"=savMAI))
	par( mai=c( 0.4, 0.6, 0.6, 0.2))
	par( mfrow=c( 2, 3))
	on.exit( par( mfrow=c(1,1)))

	vCol <- colSet[ as.numeric( grps)]
	xx <- jitter( as.numeric( grps), amount=0.1)

	v <- as.numeric( sml$EditDistance)
	yLim <- range(v,0,5) + diff(range(v,0,5)) * c( -0.05, 0.05)
	boxplot( v ~ grps, col=colSet, ylab="Edit Distance", main=paste( "Edit Distance:  ", geneID, "\n", prod),
			ylim=yLim, pars=list( boxwex=0.4, font.axis=2, font.lab=2))
	lines( c(0,3), c(0,0), lty=2, lwd=2, col='gray50'); text( 1.5, 0, "Reference", cex=0.95, pos=3, adj=0.2)
	points( xx, v, pch=21, bg=vCol, cex=2)
	if (show.labels) {
		show <- which( v > max(v)*0.1)
		if (length(show)) text( xx[show], v[show], labels[show], pos=grpPos[show], cex=0.9)
	}

	v <- as.numeric( sml$DeltaLength)
	yLim <- range(v,0,5) + diff(range(v,0,5)) * c( -0.05, 0.05)
	boxplot( v ~ grps, col=colSet, ylab="Delta Length", main=paste( "Delta Length:  ", geneID, "\n", prod),
			ylim=yLim, pars=list( boxwex=0.4, font.axis=2, font.lab=2))
	lines( c(0,3), c(0,0), lty=2, lwd=2, col='gray50'); text( 1.5, 0, "Reference", cex=0.95, pos=3, adj=0.2)
	points( xx, v, pch=21, bg=vCol, cex=2)
	if (show.labels) {
		show <- which( abs(v) > max(abs(v))*0.1)
		if (length(show)) text( xx[show], v[show], labels[show], pos=grpPos[show], cex=0.9)
	}

	v <- as.numeric( sml$Pct_Different)
	yLim <- range(v,0,5) + diff(range(v,0,5)) * c( -0.05, 0.05)
	boxplot( v ~ grps, col=colSet, ylab="Pct Different", main=paste( "Percent Different:  ", geneID, "\n", prod),
			ylim=yLim, pars=list( boxwex=0.4, font.axis=2, font.lab=2))
	lines( c(0,3), c(0,0), lty=2, lwd=2, col='gray50'); text( 1.5, 0, "Reference", cex=0.95, pos=3, adj=0.2)
	points( xx, v, pch=21, bg=vCol, cex=2)
	if (show.labels) {
		show <- which( v > max(v)*0.1)
		if (length(show)) text( xx[show], v[show], labels[show], pos=grpPos[show], cex=0.9)
	}

	v <- as.numeric( sml$InternalStopCodons)
	yLim <- range(v,0,5) + diff(range(v,0,5)) * c( -0.05, 0.05)
	boxplot( v ~ grps, col=colSet, ylab="Internal Stop Codons", main=paste( "Internal Stop Codons:  ", geneID, "\n", prod),
			ylim=yLim, pars=list( boxwex=0.4, font.axis=2, font.lab=2))
	lines( c(0,3), c(0,0), lty=2, lwd=2, col='gray50'); text( 1.5, 0, "Reference", cex=0.95, pos=3, adj=0.2)
	points( xx, v, pch=21, bg=vCol, cex=2)
	if (show.labels) {
		show <- which( v > max(v)*0.1)
		if (length(show)) text( xx[show], v[show], labels[show], pos=grpPos[show], cex=0.9)
	}

	v <- as.numeric( sml$Indel_Gaps)
	yLim <- range(v,0,5) + diff(range(v,0,5)) * c( -0.05, 0.05)
	boxplot( v ~ grps, col=colSet, ylab="Indel Gaps", main=paste( "Indel Gaps:  ", geneID, "\n", prod),
			ylim=yLim, pars=list( boxwex=0.4, font.axis=2, font.lab=2))
	lines( c(0,3), c(0,0), lty=2, lwd=2, col='gray50'); text( 1.5, 0, "Reference", cex=0.95, pos=1, adj=0.2)
	points( xx, v, pch=21, bg=vCol, cex=2)
	if (show.labels) {
		show <- which( v > max(v)*0.1)
		if (length(show)) text( xx[show], v[show], labels[show], pos=grpPos[show], cex=0.9)
	}

	v <- as.numeric( sml$Confidence)
	yLim <- range(v,50,100) + diff(range(v,50,100)) * c( -0.05, 0.05)
	boxplot( v ~ grps, col=colSet, ylab="Confidence Score", main=paste( "Confidence Score:  ", geneID, "\n", prod),
			ylim=yLim, pars=list( boxwex=0.4, font.axis=2, font.lab=2))
	lines( c(0,3), c(100,100), lty=2, lwd=2, col='gray50'); text( 1.5, 100, "Reference", cex=0.95, pos=1, adj=0.2)
	points( xx, v, pch=21, bg=vCol, cex=2)
	if (show.labels) {
		show <- which( v < max(v)*0.9)
		if (length(show)) text( xx[show], v[show], labels[show], pos=grpPos[show], cex=0.9)
	}

	return( nrow(sml))
}


`pipe.BAMproteins.GatherOneGene` <- function( sampleIDset, geneID, groupText=NULL, optionsFile="Options.txt", 
						speciesID=getCurrentSpecies(), results.path=NULL, verbose=F) {

	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	if (speciesID != getCurrentSpecies()) setCurrentSpecies( speciesID)

	#geneID <- alias2Gene( geneID[1])
	shortID <- shortGeneName( geneID[1], keep=1)

	seqs <- descs <- vector()

	# find this one gene in the FASTA file sets
	# first the reference genome
	if (verbose) cat( "\nGathering..\n")
	refFile <- paste( speciesID, "ReferenceProteins.fasta", sep=".")
	if ( ! file.exists( refFile)) {
		cat( "\nMaking FASTA of Reference Proteins..\n")
		genomicFastaFile <- getOptionValue( optT, "genomicFastaFile", notfound="Pf_genomicDNA.fasta", verbose=F)
		allGenes <- unique( getCurrentCdsMap()$GENE_ID)
		fa <- gene2Fasta( allGenes, genomicFastaFile, mode="aa", verbose=T)
		writeFasta( fa, refFile, line.width=100)
		cat( "\nDone.\n")
	}
	if ( file.exists( refFile)) {
		fa <- loadFasta( refFile, short=T, verbose=F)
		wh <- match( geneID, fa$desc, nomatch=0)
		if (wh) {
			descs <- c( descs, paste( "Reference", sep="_"))
			seqs <- c( seqs, fa$seq[wh])
			if (verbose) cat( "\rReference", length(seqs))
		}
	}
	# then in each sample 
	if ( ! is.null(groupText)) groupText <- rep( groupText, length.out=length(sampleIDset))
	for ( i in 1:length(sampleIDset)) {
		sid <- sampleIDset[i]
		f <- file.path( results.path, "ConsensusProteins", sid, paste( "All.BAM.Proteins", sid, "fasta", sep="."))
		if ( file.exists( f)) {
			fa <- loadFasta( f, short=T, verbose=F)
			wh <- match( geneID, fa$desc, nomatch=0)
			if (wh) {
				descs <- c( descs, if ( is.null(groupText)) sid else paste( sid, groupText[i], sep="_"))
				seqs <- c( seqs, fa$seq[wh])
				if (verbose) cat( "\r", sid, length(seqs))
			}
		}
	}

	# package up the results
	if ( ! length(seqs)) {
		cat( "\nFound zero proteins for gene: ", geneID)
		return( NULL)
	} else {
		if (verbose) cat( "\nFound", length(seqs), "proteins for gene: ", geneID)
	}
	descs <- paste( descs, shortID, sep="_")
	faOut <- as.Fasta( descs, seqs)
	faFile <- paste( geneID, "BAM.Proteins.fasta", sep=".")
	outPath <- file.path( results.path, "ConsensusProteins", "BAM.Proteins.By.Gene")
	faFile <- file.path( outPath, faFile)
	if ( ! file.exists( outPath)) dir.create( outPath, recursive=T)
	writeFasta( faOut, faFile, line.width=100)
	if ( length(descs) > 1) {
		if (verbose) cat( "\nAligning..")
		alnFile <- file.path( outPath, paste( geneID, "BAM.Proteins.aln", sep="."))
		#aln <- mafft( faFile, alnFile, mode='local', mafftArgs=" --anysymbol ", verbose=verbose)
		aln <- mafft( faFile, alnFile, mode='local', mafftArgs="", verbose=verbose)
		writeALN( aln, alnFile, line.width=100)
	}

	# done.
	return( invisible(faOut))
}


`pipe.BAMproteins.FindDifferencesOneGene` <- function( sampleIDset, geneID, groupSet, optionsFile="Options.txt", 
						results.path=NULL, verbose=F) {

	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=verbose)
	}

	# find this one gene's FASTA and ALN under the BAM proteins folder
	bamProteinPath <- file.path( results.path, "ConsensusProteins", "BAM.Proteins.By.Gene")
	faFile <- file.path( bamProteinPath, paste( geneID, "BAM.Proteins.fasta", sep="."))
	if ( ! file.exists( faFile)) {
		cat( "\nExpected FASTA file not found: ", faFile)
		cat( "\nPerhaps run 'pipe.BAMproteins.GatherOneGene()' first..")
		return(NULL)
	}
	fa <- loadFasta( faFile, short=F, verbose=verbose)
	alnFile <- file.path( bamProteinPath, paste( geneID, "BAM.Proteins.aln", sep="."))
	if ( ! file.exists( alnFile)) {
		cat( "\nExpected ALN file not found: ", alnFile)
		cat( "\nPerhaps run 'pipe.BAMproteins.GatherOneGene()' first..")
		return(NULL)
	}
	aln <- readALN( alnFile, verbose=verbose)
	alnM <- aln$alignment

	# get the alignment rows that match the samples we want, noting that the reference is in row #1
	# since ALN can crop sample names, be more careful finding the sample IDs
	N <- length( sampleIDset)
	if ( length( groupSet) != N) stop( "'groupSet' must be same length as 'sampleIDset'")
	grpFac <- factor( groupSet)
	if ( nlevels(grpFac) != 2) stop( "'groupSet' must contain exactly 2 factor levels")
	id2alnPtr <- rep.int( 0, N)
	for ( i in 1:N) {
		sid <- sampleIDset[i]
		wh <- pmatch( sid, rownames(alnM), nomatch=0)
		if (wh == 0) {
			sid <- substr( sampleIDset[i], 1, 15)
			wh <- pmatch( sid, rownames(alnM), nomatch=0)
		}
		id2alnPtr[i] <- wh
	}
	if ( any( id2alnPtr == 0)) {
		keep <- which( id2alnPtr > 0)
		missing <- which( id2alnPtr == 0)
		cat( "\nSome samples not found in BAM protein FASTA files: ", length(missing))
		cat( "\nMissing = ", sampleIDset[missing])
		sampleIDset <- sampleIDset[keep]
		groupSet <- groupSet[keep]
		grpFac <- factor( groupSet)
		if ( nlevels(grpFac) != 2) stop( "'groupSet' must contain exactly 2 factor levels")
		id2alnPtr <- id2alnPtr[keep]
	}
	refM <- alnM[ 1, ]
	alnM <- alnM[ id2alnPtr, ]

	# set up to do AA comparisons at every location along the protein, using the reference for AA counting 
	is1 <- which( as.numeric(grpFac) == 1)
	is2 <- which( as.numeric(grpFac) == 2)
	if ( length(is1) < 2 || length(is2) < 2) {
		cat( "\nNot enough samples in some groups. Unable to compare:\n")
		print( table( groupSet))
		return(NULL)
	}
	NAA <- ncol( alnM)
	refAApos <- 0
	outPos <- outPval <- outCons1 <- outCons2 <- outStr1 <- outStr2 <- vector()
	
	PASTE <- base::paste
	SORT <- sort.default
	TABLE <- base::table
	WHICH.MAX <- base::which.max

	for ( i in 1:NAA) {
		refAA <- refM[ i]
		if ( refAA != "-") refAApos <- refAApos + 1
		outPos[i] <- refAApos
		aaV <- alnM[ , i]
		allAA <- sort.default( unique.default( aaV))
		if ( length(allAA) > 1) {
			cnts1 <- TABLE( factor( aaV[ is1], levels=allAA))
			cnts2 <- TABLE( factor( aaV[ is2], levels=allAA))
			cntsM <- matrix( c(cnts1,cnts2), nrow=length(allAA), ncol=2)
			outPval[i] <- suppressWarnings( prop.test( cntsM, correct=F))$p.value
			outCons1[i] <- names(cnts1)[WHICH.MAX(cnts1)]
			outCons2[i] <- names(cnts2)[WHICH.MAX(cnts2)]
			cnts1 <- SORT( cnts1, decreasing=T)
			cnts2 <- SORT( cnts2, decreasing=T)
			outStr1[i] <- PASTE( names(cnts1), as.numeric(cnts1), sep=":", collapse="; ")
			outStr2[i] <- PASTE( names(cnts2), as.numeric(cnts2), sep=":", collapse="; ")
		} else {
			outPval[i] <- 1
			outCons1[i] <- outCons2[i] <- allAA
			outStr1[i] <- PASTE( allAA,  length(is1), sep=":")[1]
			outStr2[i] <- PASTE( allAA,  length(is2), sep=":")[1]
		}
	}

	# package up the results
	out <- data.frame( "ALN.POS"=1:NAA, "REF.POS"=outPos, "REF.AA"=refM, "Consensus.1"=outCons1, 
			"Consensus.2"=outCons2, "Distribution.1"=outStr1, "Distribution.2"=outStr2,
			"P.Value"=outPval, stringsAsFactors=F)

	# put the group names in explicitly
	colnames(out)[4:5] <- paste( "Consensus", levels(grpFac), sep="_")
	colnames(out)[6:7] <- paste( "Distribution", levels(grpFac), sep="_")

	# done.
	return( out)
}



`pipe.BAMproteins.FindGeneDifferences` <- function( sampleIDset, groupSet, geneIDset=NULL, optionsFile="Options.txt", 
					results.path=NULL, speciesID=getCurrentSpecies(), 
					min.pvalue=0.05, verbose=F) {

	if ( speciesID != getCurrentSpecies()) setCurrentSpecies(speciesID)
	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=verbose)
	}

	# find all gene FASTA and ALN under the BAM proteins folder
	bamProteinPath <- file.path( results.path, "ConsensusProteins", "BAM.Proteins.By.Gene")
	if ( ! file.exists( bamProteinPath)) dir.create( bamProteinPath, recursive=T)

	# get the list of genes to query
	if ( is.null( geneIDset)) {
		cmap <- getCurrentCdsMap()
		geneIDset <- unique( cmap$GENE_ID)
	}

	cat( "\nVisiting", length(geneIDset), "proteins..\n")

	# do it as a local function so we can parallelize
	gatherOneBamProteinGene <- function( g) {

		# visit each gene, using existing data if it's there
		faFile <- file.path( bamProteinPath, paste( g, "BAM.Proteins.fasta", sep="."))
		alnFile <- file.path( bamProteinPath, paste( g, "BAM.Proteins.aln", sep="."))
		needBuild <- TRUE
		if ( all( file.exists( c(faFile,alnFile)))) {
			# make sure all the sample we want are present
			fa <- loadFasta( faFile, short=T, verbose=verbose)
			expectDesc <- paste( sampleIDset, g, sep="_")
			if ( all( expectDesc %in% fa$desc)) needBuild <- FALSE
		}
		if ( needBuild) {
			ans1 <- pipe.BAMproteins.GatherOneGene( sampleIDset, geneID=g, optionsFile=optionsFile, 
						results.path=results.path, speciesID=speciesID, verbose=F)
			if ( is.null(ans1)) return(NULL)
		}

		# do that compare
		ans2 <- pipe.BAMproteins.FindDifferencesOneGene( sampleIDset, geneID=g, groupSet=groupSet, optionsFile=optionsFile, 
						results.path=results.path, verbose=F)
		if ( is.null(ans2)) return(NULL)

		# only keep significant rows
		ans2 <- subset( ans2, P.Value <= min.pvalue)
		if ( ! nrow(ans2)) return(NULL)
		cat( "\r", g, nrow(ans2))
		return(ans2)
	}

	# accumulate results from all genes
	bigAns <- multicore.lapply( geneIDset, gatherOneBamProteinGene)
	bigOut <- data.frame()
	for ( k in 1:length(geneIDset)) {
		smlAns <- bigAns[[k]]
		if ( is.null( smlAns)) next
		if ( ! nrow( smlAns)) next
		g <- geneIDset[k]
		smlOut <- data.frame( "GENE_ID"=g, "PRODUCT"=gene2Product(g), smlAns, stringsAsFactors=F)
		bigOut <- rbind( bigOut, smlOut)
	}
	cat( "\nDone.  Total significant protein AA difference sites: ", nrow(bigOut))

	# lastly, sort into P value ordering
	if ( nrow( bigOut)) {
		ord <- order( bigOut$P.Value, bigOut$GENE_ID, bigOut$REF.POS)
		bigOut <- bigOut[ ord, ]
		rownames(bigOut) <- 1:nrow( bigOut)
	}
	return( bigOut)
}

