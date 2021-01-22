# pipe.ConsensusProteinExtraction.R -- tools to extract the final protein calls from CPP
#				allowing for 2+ proteins of different proportions to be resolved by deconvolution


`pipe.ConsensusProteinExtraction` <- function( sampleID, geneID="PF3D7_1200600", geneName="Var2csa", optionsFile="Options.txt",
						results.path=NULL, max.proteins=NULL, min.minor.pct=6, min.mutation.pct=1.0, min.minor.read.count=2,
						verbose=TRUE ) {
						
	# min.minor.pct - at any one amino acid, how frequent must a minor AA call be to not be treated as just noise
	# min.mutation.pct -  for the full length protein, what fraction of AA must be mutated to call it a separate protein call


	# make sure we were given valid parameters..
	if ( length(sampleID) > 1) {
		sampleID <- sampleID[1]
		cat( "\nWarning:  only one sampleID allowed")
	}
	gmap <- getCurrentGeneMap()
	if ( ! (geneID %in% gmap$GENE_ID)) {
		cat( "\nNot a known geneID: ", geneID)
		return(NULL)
	}
	if ( is.null(results.path)) results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	peptide.path <- file.path( results.path, "ConsensusProteins", sampleID)
	if ( ! file.exists( peptide.path)) {
		cat( "\nNo 'Consensus Protein' results found for sample: ", sampleID)
		return(NULL)
	}

	# build the filename of the expected results, and get the consensus protein call details 
	summaryFile <- file.path( peptide.path, paste( sampleID, geneName, "ConsensusProteinSummary.txt", sep="."))						
	if ( ! file.exists( summaryFile)) {
		cat( "\nProtein summary file not found: ", summaryFile)
		return(NULL)
	}
	protDF <- read.delim( summaryFile, as.is=T)
	NAA <- nrow( protDF)

	# there is a small chance that this file is old, from befoer the 'Counts' field was added.
	# check and perhaps remake it
	if ( ! "Counts" %in% colnames(protDF)) {
		sumFile <- sub( "ProteinSummary.txt$", "Answer.rda", summaryFile)
		tmpAns <- proteinConstructPileupSummary( sumFile, sampleID=sampleID, geneName=geneName, doPlot=F)
		protDF <- read.delim( summaryFile, as.is=T)
		NAA <- nrow( protDF)
	}
	
	# get the consensus protein(s) themselves, and extract the initial guess about proportions from their names
	topProteins <- cpp.ExtractTopProteins( summaryFile, min.minor.pct=min.minor.pct, min.mutation.pct=min.mutation.pct,
						drop.gaps=FALSE)
	NPROT <- length( topProteins)
	protNames <- names( topProteins)

	# allow the user to put a upper limit on how many proteins we consider
	if ( ! is.null( max.proteins)) {
		topProteins <- topProteins[ 1:min(NPROT,max.proteins)]
		NPROT <- length( topProteins)
		protNames <- names( topProteins)
	}
	protProportions <- rep.int( NA, NPROT)
	protProportions[1] <- 100
	if (NPROT > 1) {
		for ( i in 2:NPROT) protProportions[i] <- as.numeric( sub( "(^MinorVariant.+Pct=)([0-9]+)(%$)", "\\2", protNames[i]))
		protProportions[1] <- 100 - sum( protProportions[2:NPROT], na.rm=T)
	}
	if (verbose) {
		cat( "\nInitial Proportions: \n")
		print( paste( protNames, protProportions, sep=" = "))
	}
	# as defined, all these constructs should be the same length.  Make sure
	if ( any( nchar(topProteins) != NAA)) {
		stop( "\nSize error:  consensus protein lengths not all equal..")
	}
	# turn all the protein calls into one AA matrix
	protV <- strsplit( topProteins, split="")
	aaM <- matrix( unlist(protV), nrow=NAA, ncol=NPROT)
	colnames(aaM) <- protNames
	
	# prep all the protein calls and all the summary percentages & counts
	protV <- strsplit( topProteins, split="")
	pcts <- protDF$Percentages
	cnts <- protDF$Counts
	pctTerms <- strsplit( pcts, split="; ")
	cntTerms <- strsplit( cnts, split="; ")
	pctAAs <- lapply( pctTerms, function(x) sub( ":[0-9]+", "", x))
	pctVals <- lapply( pctTerms, function(x) as.numeric( sub( "[-*A-Z]:", "", x)))
	cntAAs <- lapply( cntTerms, function(x) sub( ":[0-9]+", "", x))
	cntVals <- lapply( cntTerms, function(x) as.numeric( sub( "[-*A-Z]:", "", x)))
	
	# we can tell who is conserved or not by how many AA were called
	nAAterms <- sapply( pctAAs, length)
	useForScoring <- which( nAAterms > 1)
	
	# we will ignore very small read counts as being not seen enough to be considered a true minor variant
	for ( k in useForScoring) {
		myCnts <- cntVals[[k]]
		tooLow <- which( myCnts < min.minor.read.count)
		if ( length(tooLow)) {
			# remove these from both percents and counts
			cntVals[[k]] <- cntVals[[k]][ -tooLow]
			cntAAs[[k]] <- cntAAs[[k]][ -tooLow]
			pctVals[[k]] <- pctVals[[k]][ -tooLow]
			pctAAs[[k]] <- pctAAs[[k]][ -tooLow]
			# re-evaluate the percentages after this trim of low count minor variants
			pctVals[[k]] <- round( cntVals[[k]] * 100 / sum( cntVals[[k]], na.rm=T))
		}
	}
	
	# now we need to reassess who is conserved or not by how many AA were called
	nAAterms <- sapply( pctAAs, length)
	useForScoring <- which( nAAterms > 1)
	if (verbose) cat( "\nN AA with minor variants: ", length(useForScoring))
	
	# use one universal set of all possible AA for doing our compares and distance scoring
	ALL.AA <- sort( unique( c( unlist( pctAAs), unlist(cntAAs))))
	
	# prebuild the static distributions of the called AAs
	aaPctTables <- vector( mode="list", length=NAA)
	for ( k in useForScoring) {
		myAA <- rep( pctAAs[[k]], times=pctVals[[k]])
		myTbl <- table( factor(myAA, levels=ALL.AA))
		aaPctTables[[k]] <- myTbl
	}
	
	
	# local function that scores how well the consensus proteins match/explain the call details
	# puts some results into global storage for speed / later use...
	aaDist <- rep.int( 0, NAA)
	aaM.out <- aaM
	
	consensus.AA.Distance <- function( aaM, protPcts) {
	
		# given one matrix of all the protein AA calls for every location, and the percentage 
		# that each protein contributes to the total, see how this current consensus and percentage
		# fits to the called data from the CPP tool
		protPcts <- as.integer( protPcts)
			
		# visit each AA that has any minor forms, and see how distant is current from what was called
		for ( k in useForScoring) {
			myAA <- rep( aaM[ k, ], times=protPcts)
			myTbl <- table( factor(myAA, levels=ALL.AA))
			calledTbl <- aaPctTables[[k]]
			myDistTbl <- myTbl - calledTbl
			myDist <- sum( abs( myDistTbl), na.rm=T) 
			# store these where the caller can use them later
			aaDist[k] <<- myDist 
		}
		return( sum(aaDist) / nrow(aaM))
	}
	
	
	optimize.Proportions <- function( aaM, protPcts, verbose=T) {
	
		# given one matrix of all the protein AA calls for every location, and the percentage 
		# that each protein contributes to the total, optimize those percentages to give the best
		# fit to the called data
		
		# ready to do the optimization.   Evaluate at the initial estimate
		startingDist <- consensus.AA.Distance( aaM, protProportions)
		if (verbose) cat( "\nInitial Distance: ", startingDist)

		bestDist <- startingDist
		bestProportion <- protProportions
		iterCount <- 0
		while ( TRUE ) {
			# step over all alternative pairs of proteins
			thisProportion <- lastProportion <- bestProportion
			iterCount <- iterCount + 1
			bestI <- bestJ <- 0
			for (i in 1:NPROT) {
				for ( j in 1:NPROT) {
					if ( i == j) next
					# tweak the proportions between all 2 choices of protein
					thisProportion <- lastProportion
					thisProportion[i] <- thisProportion[i] + 1
					thisProportion[j] <- thisProportion[j] - 1
					# rescore the distance
					thisDist <- consensus.AA.Distance( aaM, thisProportion)
					if ( thisDist < bestDist) {
						bestDist <- thisDist
						bestProportion <- thisProportion
						bestI <- i
						bestJ <- j
					}
				}
			}
			# we did all possible steps, was any better?
			if ( bestI == 0 || bestJ == 0) break
			if (verbose) cat( "\nIter: ", iterCount, "  Dist: ", bestDist, "  Proportions: ", bestProportion)
		}
		return( list( "Distance"=bestDist, "Proportions"=bestProportion))
	}


	optimize.ResidueCalls <- function( aaM, protPcts, verbose=T) {
	
		# given one matrix of all the protein AA calls for every location, and the percentage 
		# that each protein contributes to the total, optimize the AA residues between the proteins
		# to give the best fit to the called data
		
		# the current distance for each minor site already sits in global storage...

		aaM.out <- aaM
		dist.out <- aaDist
		
		# visit each AA that has any minor forms, and see how distant is affected by permuting the residues
		for ( k in useForScoring) {
		
			# start with the set of AA the proteins were on function entry
			myAAorig <- aaM[ k, ]
			aaSet <- sort( unique( myAAorig))
			nChoice <- length( aaSet)
			nNeed <- ncol(aaM)
			
			# step over all possible combinations of AA
			bestDist <- origDist <- aaDist[k]
			setPtr <- rep.int( 1, nNeed)
			bestAAtry <- NULL
			while (TRUE) {
				# use the set pointers to make this possible choice of AA calls
				myAAtry <- aaSet[ setPtr]
				myAA <- rep( myAAtry, times=protPcts)
				myTbl <- table( factor(myAA, levels=ALL.AA))
				calledTbl <- aaPctTables[[k]]
				myDistTbl <- myTbl - calledTbl
				thisDist <- sum( abs( myDistTbl), na.rm=T) 
				if ( thisDist < bestDist) {
					bestDist <- thisDist
					bestAAtry <- myAAtry
				}
				# now iterate by incrementing the pointers to eventually try all possible combinations
				setPtr[1] <- setPtr[1] + 1
				for (ich in 1:(nNeed-1)) {
					if (setPtr[ich] > nChoice) {
						setPtr[ich] <- 1
						setPtr[ich+1] <- setPtr[ich+1] + 1
					}
				}
				# we are done when the very last pointer is too big
				if ( setPtr[nNeed] > nChoice) break
			}
			# we did all possible steps, was any better?
			if ( is.null(bestAAtry)) next
			if (verbose && any( myAAorig != bestAAtry)) cat( "\nAA Swap: ", k, "  OldDist:", origDist, "  OldAA:", myAAorig, 
									"|  NewDist:", bestDist, "  NewAA:", bestAAtry)
			
			# update the proteins
			aaM.out[ k, ] <- bestAAtry
			dist.out[k] <- bestDist
		}
		
		newDist <- sum( dist.out) / nrow(aaM)
		return( list( "Distance"=newDist, "AA_Matrix"=aaM.out))
	}
	
	
	topDistSites <- function( aaM, aaDist, cutoffDist, nExtraAA=2) {
	
		# gather the top residual distance sites as extra info
		keep <- which( aaDist > cutoffDist)
		badDist <- aaDist[ keep]
		badMotif <- aaM[ keep, 1]
		for ( i in 1:length(keep)) {
			k <- keep[i]
			badMotif[i] <- paste( aaM[(max(1,k-nExtraAA)):(min(nrow(aaM),k+nExtraAA)), 1], collapse="")
		}
		sml <- data.frame( "AA.Position"=keep, "AA.Motif"=badMotif, "Residual.Distance"=badDist, stringsAsFactors=F)
		ord <- order( sml$Residual.Distance, decreasing=T)
		sml <- sml[ ord, ]
		rownames(sml) <- 1:nrow(sml)
		return(sml)
	}

	
	# ready to do the work...
	# Step 1:  optimize the starting proportions, given what the extractor guessed
	myProportions <- protProportions
	ans1 <- optimize.Proportions( aaM, myProportions, verbose=verbose)
	myProportions <- ans1$Proportions
	myDist <- ans1$Distance
	
	if ( NPROT > 1) {
		# Step 2:  visit each minor mutation site, and see if swapping residues improves the distance scoring
		ans2 <- optimize.ResidueCalls( aaM, myProportions, verbose=verbose)
		aaM.out <- ans2$AA_Matrix
	
		# Step 3:  re-optimize the  proportions, given the improved residue calls
		ans3 <- optimize.Proportions( aaM.out, myProportions, verbose=verbose)
		myProportions <- ans3$Proportions
		myDist <- ans3$Distance
	}
	
	# make the final answer
	finalProportions <- myProportions
	finalDist <- myDist
	finalAAseqs <- apply( aaM.out, MARGIN=2, FUN=function(x) paste( x, collapse=""))
	finalNames <- c( "Consensus")
	if (NPROT > 1) finalNames <- c( finalNames, paste( "Variant", 2:NPROT, sep=""))
	finalSuffix <- paste( round( finalProportions), "%", sep="")
	finalNames <- paste( finalNames, finalSuffix, sep="_")
	
	if (verbose) cat( "\nFinal Extraction: \n", finalNames, "\nFinal Distance: ", finalDist)

	outSeqs <- as.Fasta( finalNames, finalAAseqs)
	finalFile <- file.path( peptide.path, paste( sampleID, geneName, "FinalExtractedAA.fasta", sep="."))
	writeFasta( outSeqs, finalFile, line=100)

	alnFile <- file.path( peptide.path, paste( sampleID, geneName, "FinalExtractedAA.aln", sep="."))
	if ( NPROT > 1) {
		aln <- mafft( finalFile, alnFile)
		writeALN( aln, alnFile, line=100)	
	} else {
		# make sure any earlier result gets removed
		file.delete( alnFile)
	}
	
	ans4 <- topDistSites( aaM.out, aaDist, finalDist*2)
	residSiteFile <- file.path( peptide.path, paste( sampleID, geneName, "FinalExtracted.ResidualSites.csv", sep="."))
	write.table( ans4, residSiteFile, sep=",", quote=T, row.names=F)	
	
	out <- list( "AA.Fasta"=outSeqs, "Residual"=finalDist, "Problem.Sites"=ans4)
	return( invisible(out))
}
	


`cpp.ExtractTopProteins` <- function( proteinSummaryFile, min.minor.pct=5.0, min.mutation.pct=1.0, 
					drop.gaps=TRUE) {

	# min.minor.pct - at any one amino acid, how frequent must a minor AA call be to not be treated as just noise
	# min.mutation.pct -  for the full length protein, what fraction of AA must be mutated to call it a separate protein call

	tbl <- read.delim( proteinSummaryFile, as.is=T)

	# the most likely one protein is already called
	bestProtein <- tbl$ConsensusAA
	NAA <- length( bestProtein)
	
	# always return the best
	out <- paste( bestProtein, collapse="")
	names(out) <- "BestConsensus"

	# to know the minor proteins, look to the table of details
	pctDetails <- tbl$Percentages
	pctTerms <- strsplit( pctDetails, split="; ", fixed=T)
	nTerms <- sapply( pctTerms, length)
	has2plus <- which( nTerms > 1)

	# start with the best, and then substitute where the minor call is deep enough
	minorPctsFound <- rep.int( NA, 5)
	for ( minor in 2:5) {
		nMutat <- 0
		myPcts <- vector()
		minorProtein <- bestProtein
		for ( j in has2plus) {
			if ( nTerms[j] < minor) next
			thisTerm <- pctTerms[[j]][minor]
			thisAA <- sub( ":.+", "", thisTerm)
			thisPct <- as.numeric(sub( ".+:", "", thisTerm))
			if ( ! is.na( thisPct) && thisPct >= min.minor.pct) {
				minorProtein[j] <- thisAA
				nMutat <- nMutat + 1
				myPcts[ nMutat] <- thisPct
			}
		}
		
		# did we get enought mutations to call this a new protein variant?
		if ( nMutat < (NAA * (min.mutation.pct/100))) next
		
		minorProtein <- paste( minorProtein, collapse="")
		out <- c( out, minorProtein)
		avgPct <- round( mean( myPcts))
		minorPctsFound[ minor] <- avgPct
	}

	# there is a chance of gaps being in the final set, drop them
	if ( drop.gaps) {
		out <- gsub( "-", "", out, fixed=T)
	}
	
	if ( length(out) > 1) {
		names(out) <- c( "BestConsensus", paste( "MinorVariant", 2:length(out), "_Pct=", minorPctsFound[2:length(out)], "%", sep=""))
	}
	return(out)
}



`pipe.RealignExtractedProteins` <- function( sampleID, geneID="PF3D7_1200600", geneName="Var2csa", optionsFile="Options.txt",
						results.path=NULL) {
						
	# after hand curation of the extracted .ALN file against the pileup image, remake the final FASTA and ALN files

	if ( is.null(results.path)) results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	peptide.path <- file.path( results.path, "ConsensusProteins", sampleID)
	if ( ! file.exists( peptide.path)) {
		cat( "\nNo 'Consensus Protein' results found for sample: ", sampleID)
		return(NULL)
	}

	finalFile <- file.path( peptide.path, paste( sampleID, geneName, "FinalExtractedAA.fasta", sep="."))
	alnFile <- file.path( peptide.path, paste( sampleID, geneName, "FinalExtractedAA.aln", sep="."))
	if ( ! file.exists( alnFile)) {
		cat( "\nError:  Can not find required Final Extracted ALN file: ", alnFile)
		return(NULL)
	}
	
	# convert the ALN that was hand modified
	ALNtoFasta( alnFile, outfile=finalFile, gap.character="", line.width=100, reorder=NULL)
	
	# now redo the MSA alignment again
	aln <- mafft( finalFile, alnFile)
	writeALN( aln, alnFile, line=100)	
}
