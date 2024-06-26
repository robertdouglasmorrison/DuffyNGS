# pipe.HLA.ConsensusProteins.R -- functions to examine HLA locus genes, using the
#				Consensus Proteins Pileups (CPP) tool

# globally define the universe of HLA genes we may track, as full current Human annotation names and short names
ALL_HLA_GeneIDs <- c( "HLA-A:GI3105:06:29942470", "HLA-B:GI3106:06:31353868", "HLA-C:GI3107:06:31268749", 
			"HLA-E:GI3133:06:30489406", "HLA-DRA:GI3122:06:32439842", "HLA-DRB1:GI3123:06:32578769", 
			"HLA-DQA1:GI3117:06:32637396", "HLA-DQB1:GI3119:06:32659464", "HLA-DPA1:GI3113:06:33064569", 
			"HLA-DPB1:GI3115:06:33075926")
ALL_HLA_GeneNames <- c( "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-DRA", "HLA-DRB1", 
			"HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1")


`pipe.HLA.ConsensusProteins` <- function( sampleID=NULL, HLAgenes=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt",
				results.path=NULL, IMGT.HLA.path="~/IMGT_HLA", max.pileup.depth=80, pct.aligned.depth=0.9,
				maxNoHits.pileup=0, maxNoHits.setup=0, min.minor.pct=15, 
				doPileups=FALSE, doExtractions=doPileups, intronMaskFasta=NULL, verbose=TRUE) {

	# path for all results
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	}
	if ( length(sampleID) > 1) {
		cat( "\nWarning: HLA consensus requires a single 'sampleID'.  Using first..\n")
		sampleID <- sampleID[1]
	}
	HLAresults.path <- file.path( results.path, "HLA.ProteinCalls", sampleID)
	if ( ! file.exists( HLAresults.path)) dir.create( HLAresults.path, recursive=T)
	consensusProteins.path <- file.path( results.path, "ConsensusProteins", sampleID)
	if ( ! file.exists( consensusProteins.path)) dir.create( consensusProteins.path, recursive=T)

	# force human as the current species
	setCurrentSpecies( "Hs_grc")

	# start from the list of Hs_grc HLA genes to harvest
	HLAgeneIDs <- ALL_HLA_GeneIDs
	HLAgeneNames <- ALL_HLA_GeneNames

	# allow being given a subset of HLA genes
	if ( ! is.null( HLAgenes)) {
		hitsID <- which( HLAgeneIDs %in% HLAgenes)
		hitsName <- which( HLAgeneNames %in% HLAgenes)
		hits <- sort( unique( c(hitsID, hitsName)))
		if ( ! length(hits)) {
			cat( "\nWarning: HLA genes must be from: ", HLAgeneNames, HLAgeneIDs)
			return(NULL)
		}
		HLAgeneIDs <- HLAgeneIDs[ hits]
		HLAgeneNames <- HLAgeneNames[ hits]
	}
	N_HLA <- length( HLAgeneIDs)
	
	# if given intron masking info, require the tool to only operate on a single HLA gene
	if ( ! is.null(intronMaskFasta)) {
		if ( N_HLA > 1) {
			cat( "\nCaution: using an Intron Mask currently requires processing only one HLA gene at a time.")
			cat( "\n  Use the 'HLAgenes' argument to down-select to a single HLA gene.")
			return(NULL)
		}
	}

	# the HLA genes are messy, so preload the reference AA sequence as a guide
	if (verbose) cat( "\nGathering HLA reference protein sequences..")
	genomicFastaFile <- getOptionValue( optionsFile, "genomicFastaFile", notfound="Hs_genomicDNA.fasta", verbose=verbose)
	HLArefAA <- vector()
	for ( ig in 1:N_HLA) {
		hlaFA <- gene2Fasta( HLAgeneIDs[ig], genomicFastaFile, mode="aa")
		HLArefAA[ig] <- hlaFA$seq[1]
	}
	
	# we will do each HLA locus all the way through
	outLocus <- outName <- outDist <- outSeq <- vector()
	for ( i in 1:N_HLA) {
		thisGene <- HLAgeneIDs[i]
		thisName <- HLAgeneNames[i]
		thisRefAA <- HLArefAA[i]

		# step 1: call the Consensus Pileups tool
		consensusFile <- paste( sampleID, thisName, "ConsensusProteinSummary.txt", sep=".")
		consensusFile <- file.path( consensusProteins.path, consensusFile)

		# the HLA name may have used "_" instead of "-" in the past.  Try to catch and allow
		thisNameIn <- thisName
		if ( ! file.exists( consensusFile)) {
			thisNameIn <- sub( "\\-", "_", thisNameIn)
			consensusFile2 <- paste( sampleID, thisNameIn, "ConsensusProteinSummary.txt", sep=".")
			consensusFile2 <- file.path( consensusProteins.path, consensusFile2)
			if ( file.exists( consensusFile2)) {
				consensusFile <- consensusFile2
				cat( "\nUsing older results for: ", thisNameIn)
			} else {
				thisNameIn <- thisName
			}
		}

		didNewPileup <- FALSE
		if ( ! file.exists( consensusFile) || doPileups) {
			cat( "\n\nCalling 'Consensus Protein Pileups' tool..  ", sampleID, " ", thisName)
			pipe.ConsensusProteinPileups( sampleID, thisGene, thisNameIn, results.path=results.path,
						max.depth=max.pileup.depth, pct.aligned.depth=pct.aligned.depth, chunkSize.pileup=50000, 
						maxNoHits.pileup=maxNoHits.pileup, maxNoHits.setup=maxNoHits.setup,
						showFrameShiftPeptides=F, referenceAA=thisRefAA, intronMaskFasta=intronMaskFasta)
			didNewPileup <- TRUE
		}
		# if the file still not found, must be some error, skip it
		if ( ! file.exists( consensusFile)) {
			cat( "\nFailed to create consensus protein result for: ", thisName, "  Skipping..")
			next
		}

		# step 2:  Extract up to 2 proteins from this one result
		#  use the existing results, unless not found or if we did the pileup step
		extractedFile <- paste( sampleID, thisName, "FinalExtractedAA.fasta", sep=".")
		extractedFile <- file.path( consensusProteins.path, extractedFile)
		if ( ! file.exists(extractedFile) || didNewPileup || doExtractions) {
			cat( "\n\nExtracting 'Consensus Protein Pileups' sequences..  ", sampleID, " ", thisName)
			ans <- pipe.ConsensusProteinExtraction( sampleID, thisGene, thisNameIn, results.path=results.path,
						min.minor.pct=min.minor.pct, max.proteins=2, verbose=FALSE, intronMaskFasta=intronMaskFasta)
		}
			
		# step 3:  make the HLA type calls for these
		ans3 <- pipe.HLA.Calls( sampleID, thisNameIn, results.path=results.path, IMGT.HLA.path=IMGT.HLA.path, verbose=verbose)

		# accumulate results
		outLocus <- c( outLocus, ans3$Locus)
		outName <- c( outName, ans3$IMGT_Name)
		outDist <- c( outDist, ans3$EditDist)
		outSeq <- c( outSeq, ans3$Sequence)
	}

	out <- data.frame( "SampleID"=sampleID, "Locus"=outLocus, "IMGT_Name"=outName, "EditDist"=outDist, "Sequence"=outSeq, stringsAsFactors=F)
	
	# since we have no idea how many HLA genes just got processed, and that each gene got written
	# to it's own file, do not write out a 'all genes' File.  Instead, call the function to merge all
	pipe.HLA.MergeGeneCalls( sampleID, results.path=results.path)
	
	return(out)
}


`pipe.HLA.MergeGeneCalls` <- function( sampleID=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt",
				results.path=NULL) {

	# path for all results
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	}
	if ( length(sampleID) > 1) {
		cat( "\nWarning: HLA merging requires a single 'sampleID'.  Using first..\n")
		sampleID <- sampleID[1]
	}
	HLAresults.path <- file.path( results.path, "HLA.ProteinCalls", sampleID)
	
	# find all the single gene results
	fset <- dir( HLAresults.path, pattern="Allele.Calls.csv$", full=T)
	if ( ! length(fset)) {
		cat( "\nWarning: found no HLA call files for sample: ", sampleID)
		return(NULL)
	}
	out <- data.frame()
	for ( f in fset) { 
		tmp <- read.csv( f, as.is=T)
		out <- rbind( out, tmp)
	}
	outfile <- paste( sampleID, "Merged.HLA.Calls.csv", sep=".")
	outfile <- file.path( HLAresults.path, outfile)
	write.csv( out, outfile, row.names=F)
	return(out)
}


`pipe.HLA.Calls` <- function( sampleID=NULL, HLAgene=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt",
				results.path=NULL, IMGT.HLA.path="~/IMGT_HLA", verbose=TRUE) {

	# set up to do fast sequence pattern comparison
	require( Biostrings)
	require( pwalign)
	data(BLOSUM62)

	# path for all results
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	}
	if ( length(sampleID) > 1) {
		cat( "\nWarning: HLA consensus requires a single 'sampleID'.  Using first..\n")
		sampleID <- sampleID[1]
	}
	HLAresults.path <- file.path( results.path, "HLA.ProteinCalls", sampleID)
	if ( ! file.exists( HLAresults.path)) dir.create( HLAresults.path, recursive=T)
	consensusProteins.path <- file.path( results.path, "ConsensusProteins", sampleID)

	# do the gather of FASTA and calling of HLA locus for one gene
	outLocus <- outName <- outDist <- outSeq <- NA
	
	# build the expected filenames we want/need
	proteinFile <- file.path( consensusProteins.path, paste( sampleID, HLAgene, "FinalExtractedAA.fasta", sep="."))
	if ( ! file.exists(proteinFile)) {
		cat( "\nError:  final consensus HLA protein file not found.  Tried: ", proteinFile)
		return(NULL)
	}
	hlaFA <- loadFasta( proteinFile, verbose=F)
	proteins <- hlaFA$seq
	if ( is.null( proteins) || !length(proteins)) return(NULL)
	# the proteins may have gaps, stops, etc.  And never let more than top 2 alleles through
	if ( length(proteins) > 2) proteins <- proteins[ 1:2]
	proteins <- gsub( "*", "", proteins, fixed=T)
	proteins <- gsub( "-", "", proteins, fixed=T)
	proteins <- gsub( "X", "", proteins, fixed=T)
	# small chance of getting back just one sequence
	if ( is.na(proteins[2])) proteins[2] <- proteins[1]
	# bail out if we got nothing
	if ( nchar( proteins[1]) < 10) return(NULL)

	# ready to make the HLA type calls for these
	referenceAAfile <- paste( HLAgene, "AA.fasta", sep=".")
	referenceAAfile <- file.path( IMGT.HLA.path, referenceAAfile)
	if ( ! file.exists( referenceAAfile)) {
		cat( "\nError:  failed to find Reference AA FASTA file: ", referenceAAfile)
		return(NULL)
	}
	refAA <- loadFasta( referenceAAfile, short=T, verbose=F)
	imgtIDs <- refAA$desc
	imgtSeqs <- refAA$seq
	# see which is closest match to what we have for our protein call
	pa1 <- pairwiseAlignment( imgtSeqs, proteins[1], type="local", scoreOnly=T, substitutionMatrix=BLOSUM62)
	best1 <- which.max( pa1)
	d1 <- adist( proteins[1], imgtSeqs[best1])[1]
	nam1 <- paste( sampleID, "|", imgtIDs[best1], " EditDist=", d1, sep="")
	pa2 <- pairwiseAlignment( imgtSeqs, proteins[2], type="local", scoreOnly=T, substitutionMatrix=BLOSUM62)
	best2 <- which.max( pa2)
	d2 <- adist( proteins[2], imgtSeqs[best2])[1]
	nam2 <- paste( sampleID, "|", imgtIDs[best2], " EditDist=", d2, sep="")

	# Done:  Write the results
	outFA <- as.Fasta( c( nam1, nam2), proteins)
	outfile <- paste( sampleID, HLAgene, "AA.fasta", sep=".")
	outfile <- file.path( HLAresults.path, outfile)
	writeFasta( outFA, outfile, line=100)

	outLocus <- rep.int( HLAgene, 2)
	outName <- c( imgtIDs[best1], imgtIDs[best2])
	outDist <- c( d1, d2)
	outSeq <- c( proteins[1], proteins[2])

	out <- data.frame( "SampleID"=sampleID, "Locus"=outLocus, "IMGT_Name"=outName, "EditDist"=outDist, "Sequence"=outSeq, stringsAsFactors=F)
	outfile <- paste( sampleID, HLAgene, "Allele.Calls.csv", sep=".")
	outfile <- file.path( HLAresults.path, outfile)
	write.csv( out, outfile, row.names=F)

	return(out)
}


`cropHLAsuffix` <- function( hlaNames, max.suffix=2) {

	# some HLA name specifiers can have multiple suffixes of specifity.  Allow cropping off the far right ends
	suffix.pattern <- ":"
	out <- hlaNames
	N <- length(hlaNames)

	hits <- gregexpr( suffix.pattern, hlaNames)
	for ( i in 1:N) {
		thisAns <- hits[[i]]
		if ( length( thisAns) < max.suffix) next
		lastWanted <- thisAns[ max.suffix]
		out[i] <- substr( out[i], 1, lastWanted - 1)
	}
	out
}


`pipe.HLA.CallsOverview` <- function( sampleIDset, annotationFile="Annotation.txt", optionsFile="Options.txt",
				results.path=NULL, max.suffix=2) {

	# path for all results
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	}
	HLAresults.path <- file.path( results.path, "HLA.ProteinCalls")

	# visit every sample, and gather up all the HLA calls
	N <- length( sampleIDset)
	hlaTbl <- data.frame()
	for ( i in 1:N) {
		sid <- sampleIDset[i]
		my.file <- file.path( HLAresults.path, sid, paste( sid, "Merged.HLA.Calls.csv", sep="."))
		if ( ! file.exists( my.file)) {
			cat( "\nHLA Results file not found: ", my.file, "  Skipping..")
			next
		}
		smlTbl <- read.csv( my.file, as.is=T)
		if ( ! nrow(smlTbl)) next
		hlaTbl <- rbind( hlaTbl, smlTbl)
	}

	# Ok, now summmarize each HLA gene by what alleles were seen.
	geneFac <- factor( hlaTbl$Locus)
	hlaGenes <- levels(geneFac)
	NG <- length( hlaGenes)

	mOut <- matrix( "", nrow=N, ncol=NG)
	colnames(mOut) <- hlaGenes
	rownames(mOut) <- sampleIDset

	for ( i in 1:N) {
		mySID <- sampleIDset[i]
		smlTbl <- subset( hlaTbl, SampleID == mySID)
		if ( ! nrow(smlTbl)) next
		for (j in 1:NG) {
			myGID <- hlaGenes[j]
			tiny <- subset( smlTbl, Locus == myGID)
			if ( ! nrow(tiny)) next
			imgtNames <- tiny$IMGT_Name
			# remove the HLA prefix?...
			imgtNames <- sub( "^HLA\\-", "", imgtNames)
			imgtNames <- cropHLAsuffix( imgtNames, max.suffix=max.suffix)
			mOut[ i, j] <- paste( unique( imgtNames), collapse="/")
		}
	}
		
	out <- data.frame( "SampleID"=sampleIDset, mOut, stringsAsFactors=F)
	rownames(out) <- 1:nrow(out)
	out
}


`pipe.HLA.CombineSubjectReplicates` <- function( sampleIDset, annotationFile="Annotation.txt", optionsFile="Options.txt",
				results.path=NULL, max.suffix=2) {

	# path for all results
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	}
	HLAresults.path <- file.path( results.path, "HLA.ProteinCalls")

	# visit every replicate for this one subject
	N <- length( sampleIDset)
	hlaTbl <- data.frame()
	for ( i in 1:N) {
		sid <- sampleIDset[i]
		my.file <- file.path( HLAresults.path, sid, paste( sid, "Merged.HLA.Calls.csv", sep="."))
		if ( ! file.exists( my.file)) {
			cat( "\nHLA Results file not found: ", my.file, "  Skipping..")
			next
		}
		smlTbl <- read.csv( my.file, as.is=T)
		if ( ! nrow(smlTbl)) next
		hlaTbl <- rbind( hlaTbl, smlTbl)
	}

	# Ok, now summmarize each HLA gene by what alleles were most often seen.
	geneFac <- factor( hlaTbl$Locus)
	namesOut <- cntsOut <- fullOut <- vector()
	nOut <- 0
	tapply( hlaTbl$IMGT_Name, geneFac, function(x) {
		# given all the IMGT allele names for one gene from one subject
		
		# the alleles can have any number of ":xx:xx" suffix specifiers, that make calling common hits 
		# very difficult.  Allow trimming off the far right ends to make them more constant
		myAlleles <- cropHLAsuffix( x, max.suffix=max.suffix)
		
		# now count them, where we expect 2 calls for each gene (one from each parent)
		alleleCnts <- sort( table( myAlleles), decreasing=T)
		bestName <- names(alleleCnts)[1]
		bestCnt <- as.integer( alleleCnts[1])
		secondName <- names(alleleCnts)[2]
		secondCnt <- as.integer( alleleCnts[2])
		
		# catch rare case of getting same allele from both parents
		if ( is.na(secondCnt) || secondCnt < (N*0.29) || bestCnt > (N*1.35)) {
			secondName <- bestName
			secondCnt <- 0
		}
		
		# make the combined call, after some text cleanup
		nOut <<- nOut + 1
		namesOut[nOut] <<- paste( sub( "^HLA\\-", "", bestName), sub( "^.+\\*", "*", secondName), sep="/")
		cntsOut[nOut] <<- paste( bestCnt, secondCnt, sep=" / ")
		fullOut[nOut] <<- paste( sub(  "^.+\\*", "*", names(alleleCnts)), "(", as.integer(alleleCnts), ")", 
					sep="", collapse="; ")
	})
	
	out <- data.frame( "Locus"=levels(geneFac), "Allele.Calls"=namesOut, "Allele.Counts"=cntsOut, 
			"All.Allele.Frequencies"=fullOut, stringsAsFactors=F)
	return(out)
}


`pipe.HLA.Allele.Proportions.ByGroup` <- function( sampleIDset, groupColumn="Group", HLAgenes=NULL, annotationFile="Annotation.txt", 
						optionsFile="Options.txt", results.path=NULL, max.suffix=1, col=c(2,4), 
						nFDR=1000, label=NULL) {

	# path for all results
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	}
	HLAresults.path <- file.path( results.path, "HLA.ProteinCalls")

	# verify the grouping column exists
	annT <- readAnnotationTable( annotationFile)
	if ( ! (groupColumn %in% colnames(annT))) {
		cat( "\nWarning: wanted grouping column not found in Annotation file: ", groupColumn)
		return(NULL)
	}
	
	# visit every sample, and gather up all the HLA calls
	N <- length( sampleIDset)
	hlaTbl <- data.frame()
	for ( i in 1:N) {
		sid <- sampleIDset[i]
		grp <- annT[[ groupColumn]][ match( sid, annT$SampleID)]
		my.file <- file.path( HLAresults.path, sid, paste( sid, "Merged.HLA.Calls.csv", sep="."))
		if ( ! file.exists( my.file)) {
			cat( "\nHLA Results file not found: ", my.file, "  Skipping..")
			next
		}
		smlTbl <- read.csv( my.file, as.is=T)
		if ( ! nrow(smlTbl)) next
		smlTbl$Group <- grp
		hlaTbl <- rbind( hlaTbl, smlTbl)
	}
	# shorten the allele calls to wanted specificity
	hlaTbl$Short_IMGT_Name <- cropHLAsuffix( hlaTbl$IMGT_Name, max.suffix=max.suffix)
	if ( ! is.null( HLAgenes)) {
		hlaTbl <- subset( hlaTbl, Locus %in% HLAgenes)
	}
	
	# set up to know breakdowns by group
	grpFac <- factor( hlaTbl$Group)
	grpNames <- levels(grpFac)
	nGrp <- nlevels(grpFac)
	hlaNames <- sort( unique( hlaTbl$Locus))
	nHLA <- length( hlaNames)
	colUse <- rep( col, length.out=nGrp)
	
	# do each HLA gene separately
	checkX11()
	padFactor <- 1.5
	for (hlagene in hlaNames) {
		smlTbl <- subset( hlaTbl, Locus == hlagene)
		nAlleles <- length( unique( smlTbl$Short_IMGT_Name))
		cntsM <- tapply( 1:nrow(smlTbl), list(factor(smlTbl$Group,levels=grpNames),factor(smlTbl$Short_IMGT_Name)), FUN=length)
		cntsM[ is.na(cntsM)] <- 0
		colnames(cntsM) <- sub( "^HLA\\-", "", colnames(cntsM))
		# turn counts to percentages, within each group
		pctsM <- cntsM
		for ( k in 1:nGrp) pctsM[ k, ] <- round( cntsM[ k, ] * 100 / sum(cntsM[k,]), digits=2)
		
		mainText <- paste( "HLA Allele Proportions:   ", hlagene)
		if ( ! is.null(label)) mainText <- paste( mainText, label, sep="\n")
		barAns <- barplot( pctsM, beside=T, col=colUse, ylim=c( 0, max( max(pctsM,na.rm=T)*1.2, 50)), 
				xlim=c(0.4,(nrow(cntsM)*ncol(cntsM)*padFactor)+1),
				main=mainText, ylab="Percentage of allele calls", 
				xlab=NA, las=3, font.axis=2, font.lab=2, cex.axis=1.1, cex.lab=1.1, yaxt="n")
		axis( side=2, seq( 0, 100, by=20))
		
		# if asked for, now do FDR simulations that permute the grouping calls
		if ( nFDR) {
			tmpTbl <- smlTbl
			trueDiff <- apply( pctsM, 2, function(x) diff( range(x)))
			diffsByGrp <- matrix( NA, nrow=nFDR, ncol=nAlleles)
			for ( ifdr in 1:nFDR) {
				tmpTbl$Group <- sample( smlTbl$Group)
				cntsM <- tapply( 1:nrow(tmpTbl), list(factor(tmpTbl$Group,levels=grpNames),factor(tmpTbl$Short_IMGT_Name)), FUN=length)
				cntsM[ is.na(cntsM)] <- 0
				colnames(cntsM) <- sub( "^HLA\\-", "", colnames(cntsM))
				pcts2 <- cntsM
				for ( k in 1:nGrp) pcts2[ k, ] <- round( cntsM[ k, ] * 100 / sum(cntsM[k,]), digits=2)
				thisDiffs <<- apply( pcts2, 2, function(x) diff( range(x)))
				diffsByGrp[ ifdr, ] <- thisDiffs
			}
			# now see how often random was at least this different as the real data
			for ( j in 1:nAlleles) {
				nRandBetter <- sum( diffsByGrp[ , j] >= trueDiff[j])
				fdr <- round( nRandBetter / nFDR, digits=3)
				if ( fdr <= 0.2) {
					ptxt <- paste( "P=", fdr, sep="")
					ytxt <- max( pctsM[ , j])
					xtxt <- mean( barAns[ ,j])
					if ( nAlleles <= 12) {
						text( xtxt, ytxt, ptxt, cex=0.85, pos=3)
					} else {
						text( xtxt-(nAlleles/20), ytxt+0.25, ptxt, cex=0.85, srt=90, pos=4, offset=1)
					}
				}
			}
		}
		
		# show number of samples per group
		grpIDcnts <- tapply( smlTbl$SampleID, factor(smlTbl$Group), function(x) length( unique(x)))
		legendText <- paste( rownames(cntsM), " (N=", grpIDcnts, ")", sep="")
		legend( 'topright', legendText, fill=colUse, bg='white', cex=1.1)
		dev.flush();  Sys.sleep(1)
		plotFile <- paste( "HLA.Allele.Proportions_", hlagene, "_By.", groupColumn, sep="")
		printPlot( plotFile)
	}
}


`pipe.HLA.AminoAcid.Proportions.ByGroup` <- function( sampleIDset, groupColumn="Group", HLAgenes=NULL, annotationFile="Annotation.txt", 
						optionsFile="Options.txt", results.path=NULL, label=NULL) {

	# path for all results
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	}
	HLAresults.path <- file.path( results.path, "HLA.ProteinCalls")

	# verify the grouping column exists
	annT <- readAnnotationTable( annotationFile)
	if ( ! (groupColumn %in% colnames(annT))) {
		cat( "\nWarning: wanted grouping column not found in Annotation file: ", groupColumn)
		return(NULL)
	}
	
	# visit every sample, and gather up all the HLA calls & actual AA sequences
	N <- length( sampleIDset)
	hlaTbl <- data.frame()
	for ( i in 1:N) {
		sid <- sampleIDset[i]
		grp <- annT[[ groupColumn]][ match( sid, annT$SampleID)]
		my.file <- file.path( HLAresults.path, sid, paste( sid, "Merged.HLA.Calls.csv", sep="."))
		if ( ! file.exists( my.file)) {
			cat( "\nHLA Results file not found: ", my.file, "  Skipping..")
			next
		}
		smlTbl <- read.csv( my.file, as.is=T)
		if ( ! nrow(smlTbl)) next
		smlTbl$Group <- grp
		hlaTbl <- rbind( hlaTbl, smlTbl)
	}
	if ( ! is.null( HLAgenes)) {
		hlaTbl <- subset( hlaTbl, Locus %in% HLAgenes)
	}
	
	# set up to know breakdowns by group
	grpFac <- factor( hlaTbl$Group)
	grpNames <- levels(grpFac)
	nGrp <- nlevels(grpFac)
	hlaNames <- sort( unique( hlaTbl$Locus))
	nHLA <- length( hlaNames)
	
	#  local function to assess one HLA group
	`HLA.AA.Differences` <- function( hla="HLA-A") {

		# grab the AA seqs and do the MSA
		cat( "\nDoing MSA for: ", hla)
		tmpFasta <- paste( hla, "AllSamples.AA.fasta", sep=".")
		tmpALN <- paste( hla,"AllSamples.AA.aln", sep=".")
		use <- which( hlaTbl$Locus == hla)
		# combine the sample ID and the IMGT allele name into the descriptor
		desc <- paste( hlaTbl$SampleID[use], sub( "^HLA\\-", "", hlaTbl$IMGT_Name[use]), sep="_")
		tmpFA <- as.Fasta( desc, hlaTbl$Sequence[use])
		writeFasta( tmpFA, tmpFasta, line=100)
		aln <- mafft( tmpFasta, tmpALN, verbose=F)
		writeALN( aln, tmpALN, line=100)

		aaM <- aln$alignment
		nAA <- ncol(aaM)
		nSEQ <- nrow(aaM)
		if ( ! nAA || ! nSEQ) return(NULL)
		
		# now we can look at the AA calls versus the group calls
		grp <- hlaTbl$Group[use]
		aaLevels <- sort( unique( as.character( aaM)))
		nAAlevels <- length(aaLevels)
		grpFac <- factor( grp)
		grpLevels <- levels(grpFac)
		nGrps <- length(grpLevels)

		# visit every AA and assess the proprotion differences
		outAA <- outPos <- outPval <- vector( length=nAA)
		outPctStrs <- matrix( "", nrow=nAA, ncol=nGrps)
		colnames(outPctStrs) <- grpLevels
		
		for ( i in 1:nAA) {
			thisVec <- aaM[ , i]
			cntsM <- tapply( thisVec, list(grpFac,factor(thisVec,levels=aaLevels)), length)
			cntsM[ is.na(cntsM)] <- 0
			# drop AA not seen by anyone
			bigCnt <- apply( cntsM, 2, max)
			cntsM <- cntsM[ , bigCnt > 0, drop=F]
			SAV2 <<- cntsM
			pv <- 1
			if ( ncol(cntsM) > 1) {
				# we got 2+ different AA detected, so we can ask if th groups are different
				if ( nrow(cntsM) == 2) {
					test <- suppressWarnings( prop.test( t(cntsM)))
					pv <- test$p.value
				} else {
					# loop over all the pairs of groups
					pvVec <- vector()
					for (j in 1:(nrow(cntsM)-1)) {
						for (k in (j+1):nrow(cntsM)) {
							tmpM <- cntsM[ c(j,k), ]
							bigCnt <- apply( tmpM, 2, max)
							tmpM <- tmpM[ , bigCnt > 0, drop=F]
							if ( ncol(tmpM) < 2) {
								pvVec <- c( pvVec, 1)
								next
							}
							test <- suppressWarnings( prop.test( t( tmpM)))
							pvVec <- c( pvVec, test$p.value)
						}
					}
					pv <- min( pvVec)
					if ( pv < 0.05) SAV5 <<- pvVec
				}
			}
			outPos[i] <- i
			outAA[i] <- names( sort( table(thisVec), decreasing=T))[1]
			outPval[i] <- pv
			# let's show the top K AA calls in each group, as percentages
			MAX_K <- 4
			aaPctCalls <- apply( cntsM, 1, function(x) {
					pcts <- sort( round( x * 100 / sum(x), digits=0), decreasing=T)
					nKeep <- min( MAX_K, sum( pcts >= 2))
					if ( length(pcts) > nKeep) pcts <- pcts[ 1:nKeep]
					txtStr <- paste( names(pcts), ":", as.numeric(pcts), "%", sep="", collapse="; ")
					return( txtStr)
				})
			outPctStrs[ i, ] <- aaPctCalls
			if ( i %% 50 == 0) cat( "\r", hla, i, pv, aaPctCalls)
		}
		out <- data.frame( "Locus"=hla, "Position"=outPos, "Consensus.AA"=outAA, "Pvalue"=round(outPval,digits=5),
				outPctStrs, stringsAsFactors=F)
		rownames(out) <- 1:nrow(out)
		return( out)
	}

	# now call it for each HLA locus
	bigOut <- data.frame()
	for (hla in hlaNames) {
		sml <- HLA.AA.Differences( hla)
		if ( is.null(sml)) {
			cat( "\nError:  no rows of AA results for HLA group: ", hla)
			next
		}
		bigOut <- rbind( bigOut, sml)
	}
	
	# render the data as a Manhattan plot
	colSet <- rainbow( nHLA, end=0.65)
	ptCol <- colSet[ match( bigOut$Locus, hlaNames)]
	y <- -log10(bigOut$Pvalue)
	bigY <- max( c( 3, y), na.rm=T) * 1.05
	smlY <- bigY * -0.15
	mainText <- paste( "HLA Locus Amino Acid Differences by: ", groupColumn, "\n(", paste( grpNames, collapse=" .vs. "), ")")
	if ( ! is.null(label)) mainText <- paste( mainText, label, sep="\n")
	plot( 1:nrow(bigOut), y, main=mainText, xlab="Amino Acid location with each HLA Locus", ylab="-Log10( P.value)", 
		ylim=c(smlY, bigY), pch=19, col=ptCol, cex=0.8, yaxt="n")
	axis( side=2, at=pretty( c(0,bigY)))
	pv05 <- -log10( 0.05)
	lines( c(-500,nrow(bigOut)*2), rep.int(pv05,2), lty=3, lwd=1, col='grey50')
	text( nrow(bigOut)/2, pv05, "'P = 0.05' threshold", pos=3, col='grey50', cex=0.85)
	starts <- match( hlaNames, bigOut$Locus)
	stops <- c( (starts-1)[2:nHLA], nrow(bigOut))
	rect( starts, smlY, stops, smlY*0.3, border='black', col=colSet, lwd=2)
	showNames <- hlaNames
	tooLong <- which( nchar(showNames) > 5)
	showNames[tooLong] <- sub( "^HLA\\-", "", showNames[tooLong])
	text( (starts+stops)/2, smlY*0.65, showNames, font=2, cex=0.85)
	
	# done
	return( bigOut)
}
