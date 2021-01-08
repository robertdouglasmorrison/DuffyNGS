# pipe.HLA.ConsensusProteins.R -- wrapper to the Consensus Proteins tool for HLA loci

`pipe.HLA.ConsensusProteins` <- function( sampleID=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt",
				results.path=NULL, IMGT.HLA.path="~/IMGT_HLA",
				min.minor.pct=5, doPileups=FALSE, verbose=TRUE) {

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
	HLAgeneIDs <- c( "HLA-A:GI3105:06:29942470", "HLA-B:GI3106:06:31353868", 
			"HLA-C:GI3107:06:31268749", "HLA-E:GI3133:06:30489406",
			"HLA-DRA:GI3122:06:32439842", "HLA-DRB1:GI3123:06:32578769", 
			"HLA-DQA1:GI3117:06:32637396", "HLA-DQB1:GI3119:06:32659464", 
			"HLA-DPA1:GI3113:06:33064569", "HLA-DPB1:GI3115:06:33075926")
	HLAgeneNames <- c( "HLA-A", "HLA-B", "HLA-C", "HLA-E",
			"HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", 
			"HLA-DPA1", "HLA-DPB1")

	require( Biostrings)
	data(BLOSUM62)

	# we will do each HLA locus all the way through
	outLocus <- outName <- outDist <- outSeq <- vector()

	N_HLA <- length( HLAgeneIDs)
	for ( i in 1:N_HLA) {
		thisGene <- HLAgeneIDs[i]
		thisName <- HLAgeneNames[i]

		# step 1: call the Consensus Pileups tool
		consensusFile <- paste( sampleID, thisName, "ConsensusProteinSummary.txt", sep=".")
		consensusFile <- file.path( consensusProteins.path, consensusFile)
		if ( ! file.exists( consensusFile) || doPileups) {
			cat( "\n\nCalling 'Consensus Protein Pileups' tool..  ", sampleID, " ", thisName)
			pipe.ConsensusProteinPileups( sampleID, thisGene, thisName, results.path=results.path,
						max.depth=80, chunkSize.pileup=50000, maxNoHits.pileup=0, maxNoHits.setup=0,
						showFrameShiftPeptides=F)
		}
		# if the file still not found, must be some error, skip it
		if ( ! file.exists( consensusFile)) {
			cat( "\nFailed to create consensus protein result for: ", thisName, "  Skipping..")
			next
		}

		# step 2:  Extract up to 2 proteins from this one result
		#proteins <- top2proteins( consensusFile, min.heterozygousPct=min.heterozygousPct)
		cat( "\n\nExtracting 'Consensus Protein Pileups' sequences..  ", sampleID, " ", thisName)
		ans <- pipe.ConsensusProteinExtraction( sampleID, thisGene, thisName, results.path=results.path,
						min.minor.pct=min.minor.pct, max.proteins=2, verbose=FALSE)
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

