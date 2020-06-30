# pipe.ConsensusProteinExtraction.R -- tools to extract the final protein calls from CPP
#				allowing for 2+ proteins of different proportions to be resolved by deconvolution


`pipe.ConsensusProteinExtraction` <- function( sampleID, geneID, geneName=geneID, optionsFile="Options.txt",
						results.path=NULL, min.minor.pct=5.0, min.mutation.pct=1.0 ) {
						
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

	# buiild the filename of the expected results, and get the protein call details and sequences
	summaryFile <- file.path( peptide.path, paste( sampleID, geneName, "ConsensusProteinSummary.txt", sep="."))						
	if ( ! file.exists( summaryFile)) {
		cat( "\nProtein summary file not found: ", summaryFile)
		return(NULL)
	}
	protDF <- read.delim( summaryFile, as.is=T)
	NAA <- nrow( protDF)
	
	topProteins <- cpp.ExtractTopProteins( summaryFile, min.minor.pct=min.minor.pct, min.mutation.pct=min.mutation.pct,
											drop.gaps=FALSE)
	NPROT <- length( topProteins)
	protNames <- names( topProteins)
	
	# get
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
