# pipe.HLA.ConsensusProteins.R -- wrapper to the Consensus Proteins tool for HLA loci

`pipe.HLA.ConsensusProteins` <- function( sampleID=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt",
				results.path=NULL, IMGTreferenceAA.path="~/seq2HLA/references_Dec2017", 
				min.heterozygousPct=5, doPileups=FALSE, verbose=TRUE) {

	# path for all results
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	}
	HLAresults.path <- file.path( results.path, "HLA.ProteinCalls", sampleID)
	if ( ! file.exists( HLAresults.path)) dir.create( HLAresults.path, recursive=T)
	consensusProteins.path <- file.path( results.path, "ConsensusProteins", sampleID)

	# force human as the current species
	setCurrentSpecies( "Hs_grc")

	# list of Hs_grc HLA genes to harvest
	HLAgeneIDs <- c( "HLA-A:GI3105:06:29942470", "HLA-B:GI3106:06:31353868", "HLA-C:GI3107:06:31268749",
			"HLA-DRA:GI3122:06:32439842", "HLA-DRB1:GI3123:06:32578769", 
			"HLA-DQA1:GI3117:06:32637396", "HLA-DQB1:GI3119:06:32659464", 
			"HLA-DPA1:GI3113:06:33064569", "HLA-DPB1:GI3115:06:33075926")
			#"HLA-DPA2:GI646702:06:33091482", "HLA-DPB2:GI3116:06:33112516")
	HLAgeneNames <- c( "HLA_A", "HLA_B", "HLA_C", "HLA_DRA", "HLA_DRB1", "HLA_DQA1", "HLA_DQB1", 
			"HLA_DPA1", "HLA_DPB1")   # "HLA_DPA2", "HLA_DPB2")

	# we will do each HLA locus all the way through
	outLocus <- outName <- outDist <- vector()

	N_HLA <- length( HLAgeneIDs)
	for ( i in 1:N_HLA) {
		thisGene <- HLAgeneIDs[i]
		thisName <- HLAgeneNames[i]

		# step 1: call the Consensus Pileups tool
		consensusFile <- paste( sampleID, thisName, "ConsensusProteinSummary.txt", sep=".")
		consensusFile <- file.path( consensusProteins.path, consensusFile)
		if ( ! file.exists( consensusFile) || doPileups) {
			cat( "\n\nCalling 'Consensus Protein Pileups' tool..")
			pipe.ConsensusProteinPileups( sampleID, thisGene, thisName, results.path=results.path,
						max.depth=80, chunkSize.pileup=50000, maxNoHits.pileup=100000)
		}
		# if the file still not found, must be some error, skip it
		if ( ! file.exists( consensusFile)) {
			cat( "\nFailed to create consensus protein result for: ", thisName, "  Skipping..")
			next
		}

		# step 2:  Extract up to 2 proteins from this one result
		proteins <- top2proteins( consensusFile, min.heterozygousPct=min.heterozygousPct)
		if ( ! length(proteins)) next

		# step 3:  make the HLA type calls for these
		referenceAAfile <- paste( sub( "HLA_", "", thisName), "aa.fasta", sep="_")
		referenceAAfile <- file.path( IMGTreferenceAA.path, referenceAAfile)
		if ( ! file.exists( referenceAAfile)) {
			cat( "\nError:  failed to find Reference AA FASTA file: ", referenceAAfile)
			next
		}
		refAA <- loadFasta( referenceAAfile, short=F, verbose=F)
		# the IMGT names we want are the second field
		imgtTerms <- strsplit( refAA$desc, split=" +")
		imgtIDs <- sapply( imgtTerms, FUN=`[`, 2)
		imgtSeqs <- refAA$seq
		# see which is closest
		d1 <- adist( proteins[1], imgtSeqs)
		best1 <- which.min( d1)
		nam1 <- paste( sampleID, " ", imgtIDs[best1], " EditDist=", d1[best1], sep="")
		d2 <- adist( proteins[2], imgtSeqs)
		best2 <- which.min( d2)
		nam2 <- paste( sampleID, " ", imgtIDs[best2], " EditDist=", d2[best2], sep="")

		# step 4:  Write the resuls
		outFA <- as.Fasta( c( nam1, nam2), proteins)
		outfile <- paste( sampleID, thisName, "AA.fasta", sep=".")
		outfile <- file.path( HLAresults.path, outfile)
		writeFasta( outFA, outfile, line=100)

		outLocus <- c( outLocus, thisName, thisName)
		outName <- c( outName, imgtIDs[best1], imgtIDs[best2])
		outDist <- c( outDist, d1[best1], d2[best2])
	}

	out <- data.frame( "Locus"=outLocus, "IMGT_Name"=outName, "EditDist"=outDist, stringsAsFactors=F)
	return(out)
}


`top2proteins` <- function( f, min.heterozygousPct=5) {

	tbl <- read.delim( f, as.is=T)

	# the most likely one protein is already called
	bestProtein <- tbl$ConsensusAA

	# to know the second protein, look to the table of details
	pctDetails <- tbl$Percentages
	pctTerms <- strsplit( pctDetails, split="; ", fixed=T)
	nTerms <- sapply( pctTerms, length)
	has2plus <- which( nTerms > 1)

	# start with the best, and then substitute where the 2nd call is deep enough
	secondProtein <- bestProtein
	for ( j in has2plus) {
		thisTerm <- pctTerms[[j]][2]
		thisAA <- sub( ":.+", "", thisTerm)
		thisPct <- as.numeric(sub( ".+:", "", thisTerm))
		if ( ! is.na( thisPct) && thisPct >= min.heterozygousPct) {
			secondProtein[j] <- thisAA
		}
	}

	out <- c(paste(bestProtein,collapse=""), paste(secondProtein,collapse=""))
	return(out)
}
