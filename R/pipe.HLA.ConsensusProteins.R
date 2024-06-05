# pipe.HLA.ConsensusProteins.R -- wrapper to the Consensus Proteins tool for HLA loci

`pipe.HLA.ConsensusProteins` <- function( sampleID=NULL, HLAgenes=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt",
				results.path=NULL, IMGT.HLA.path="~/IMGT_HLA", max.pileup.depth=80, 
				maxNoHits.pileup=0, maxNoHits.setup=0,
				min.minor.pct=15, doPileups=FALSE, intronMaskFasta=NULL, verbose=TRUE) {

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

	# list of Hs_grc HLA genes to harvest
	HLAgeneIDs <- c( "HLA-A:GI3105:06:29942470", "HLA-B:GI3106:06:31353868", 
			"HLA-C:GI3107:06:31268749", "HLA-E:GI3133:06:30489406",
			"HLA-DRA:GI3122:06:32439842", "HLA-DRB1:GI3123:06:32578769", 
			"HLA-DQA1:GI3117:06:32637396", "HLA-DQB1:GI3119:06:32659464", 
			"HLA-DPA1:GI3113:06:33064569", "HLA-DPB1:GI3115:06:33075926")
	HLAgeneNames <- c( "HLA-A", "HLA-B", "HLA-C", "HLA-E",
			"HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", 
			"HLA-DPA1", "HLA-DPB1")

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

	require( Biostrings)
	data(BLOSUM62)

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

		if ( ! file.exists( consensusFile) || doPileups) {
			cat( "\n\nCalling 'Consensus Protein Pileups' tool..  ", sampleID, " ", thisName)
			pipe.ConsensusProteinPileups( sampleID, thisGene, thisNameIn, results.path=results.path,
						max.depth=max.pileup.depth, chunkSize.pileup=50000, 
						maxNoHits.pileup=maxNoHits.pileup, maxNoHits.setup=maxNoHits.setup,
						showFrameShiftPeptides=F, referenceAA=thisRefAA, intronMaskFasta=intronMaskFasta)
		}
		# if the file still not found, must be some error, skip it
		if ( ! file.exists( consensusFile)) {
			cat( "\nFailed to create consensus protein result for: ", thisName, "  Skipping..")
			next
		}

		# step 2:  Extract up to 2 proteins from this one result
		#proteins <- top2proteins( consensusFile, min.heterozygousPct=min.heterozygousPct)
		cat( "\n\nExtracting 'Consensus Protein Pileups' sequences..  ", sampleID, " ", thisName)
		ans <- pipe.ConsensusProteinExtraction( sampleID, thisGene, thisNameIn, results.path=results.path,
						min.minor.pct=min.minor.pct, max.proteins=2, verbose=FALSE, intronMaskFasta=intronMaskFasta)
		proteins <- ans$AA.Fasta$seq
		if ( is.null( proteins) || !length(proteins)) next
		# the proteins may have gaps, stops, etc
		proteins <- gsub( "*", "", proteins, fixed=T)
		proteins <- gsub( "-", "", proteins, fixed=T)
		proteins <- gsub( "X", "", proteins, fixed=T)
		# small chance of getting back just one sequence
		if ( is.na(proteins[2])) proteins[2] <- proteins[1]
		# bail out if we got nothing
		if ( nchar( proteins[1]) < 10) next

		# step 3:  make the HLA type calls for these
		referenceAAfile <- paste( thisName, "AA.fasta", sep=".")
		referenceAAfile <- file.path( IMGT.HLA.path, referenceAAfile)
		if ( ! file.exists( referenceAAfile)) {
			cat( "\nError:  failed to find Reference AA FASTA file: ", referenceAAfile)
			next
		}
		refAA <- loadFasta( referenceAAfile, short=T, verbose=F)
		# the IMGT names we want are the second field
		#imgtTerms <- strsplit( refAA$desc, split=" +")
		#imgtIDs <- sapply( imgtTerms, FUN=`[`, 2)
		imgtIDs <- refAA$desc
		imgtSeqs <- refAA$seq
		# see which is closest
		#d1 <- adist( proteins[1], imgtSeqs)
		#best1 <- which.min( d1)
		pa1 <- pairwiseAlignment( imgtSeqs, proteins[1], type="local", scoreOnly=T, substitutionMatrix=BLOSUM62)
		best1 <- which.max( pa1)
		d1 <- adist( proteins[1], imgtSeqs[best1])[1]
		nam1 <- paste( sampleID, "|", imgtIDs[best1], " EditDist=", d1, sep="")

		#d2 <- adist( proteins[2], imgtSeqs)
		#best2 <- which.min( d2)
		pa2 <- pairwiseAlignment( imgtSeqs, proteins[2], type="local", scoreOnly=T, substitutionMatrix=BLOSUM62)
		best2 <- which.max( pa2)
		d2 <- adist( proteins[2], imgtSeqs[best2])[1]
		nam2 <- paste( sampleID, "|", imgtIDs[best2], " EditDist=", d2, sep="")

		# step 4:  Write the resuls
		outFA <- as.Fasta( c( nam1, nam2), proteins)
		outfile <- paste( sampleID, thisName, "AA.fasta", sep=".")
		outfile <- file.path( HLAresults.path, outfile)
		writeFasta( outFA, outfile, line=100)

		outLocus <- c( outLocus, thisName, thisName)
		outName <- c( outName, imgtIDs[best1], imgtIDs[best2])
		outDist <- c( outDist, d1, d2)
		outSeq <- c( outSeq, proteins[1], proteins[2])
	}

	out <- data.frame( "SampleID"=sampleID, "Locus"=outLocus, "IMGT_Name"=outName, "EditDist"=outDist, "Sequence"=outSeq, stringsAsFactors=F)
	outfile <- paste( sampleID, "HLA.Calls.csv", sep=".")
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
		my.file <- file.path( HLAresults.path, sid, paste( sid, "HLA.Calls.csv", sep="."))
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
		my.file <- file.path( HLAresults.path, sid, paste( sid, "HLA.Calls.csv", sep="."))
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
			"All.Alleles.Frequencies"=fullOut, stringsAsFactors=F)
	return(out)
}
