# spliceDetailSummary.R


`spliceDetailSummary` <- function( genePattTable) {

	# given a table of splice counts, where the names are one or more gene::splice IDs separated by '|'

	# extract the individual gene names, etc
	nams <- names( genePattTable)
	cnts <- as.vector( genePattTable)

	# multi hits need to be partitioned back to the separete genes
	needSplit <- grep( "|", nams, fixed=TRUE)
	moreNams <- moreCnts <- vector()
	nmore <- 0

	if ( length(needSplit)) {
		splitNames <- strsplit( nams[needSplit], split="|", fixed=TRUE)
		sapply( 1:length(needSplit), function(i)  {
				myNewNames <- splitNames[[i]]
				myOldCount <- cnts[ needSplit[i]]
				N <- length(myNewNames)
				moreNams[ (nmore+1):(nmore+N)] <<- myNewNames
				moreCnts[ (nmore+1):(nmore+N)] <<- myOldCount / N
				nmore <<- nmore + N
			})
		# we have the genes from multi hits, so add those to the ones that were not multis
		nams <- c( nams[ -needSplit], moreNams)
		cnts <- c( cnts[ -needSplit], moreCnts)
	}

	# now combine the counts for genes that appear more than once
	tapply( 1:length(nams), factor(nams), function(x) {
			if ( length(x) < 2) return()
			cnts[ x[1]] <<- sum( cnts[ x])
			cnts[ x[ 2:length(x)]] <<- 0
			return()
		})
	keep <- which( cnts > 0)
	nams <- nams[ keep]
	cnts <- cnts[ keep]

	# now we need to figure out what species these all are
	terms <- strsplit( nams, split="::", fixed=T)
	genes <- sapply( terms, function(x) x[1])
	splices <- sapply( terms, function(x) x[2])
	species <- rep.int( "", length(genes))

	curSpecies <- getCurrentSpecies()
	for (spec in getCurrentTargetSpecies()) {
		setCurrentSpecies( spec)
		gmap <- getCurrentGeneMap()
		where <- match( genes, gmap$GENE_ID, nomatch=0)
		species[ where > 0] <- spec
		where <- match( genes, gmap$NAME, nomatch=0)
		species[ where > 0] <- spec
	}
	setCurrentSpecies( curSpecies)

	#OK, make a table of all these
	out <- data.frame( "GENE_ID"=genes, "SPLICE_ID"=splices, "N_READS"=round( cnts), "SPECIES_ID"=species,
			stringsAsFactors=FALSE)

	keep <- which( out$N_READS >= 1)
	out <- out[ keep, ]
	ord <- order( out$N_READS, decreasing=TRUE)
	out <- out[ ord, ]
	rownames(out) <- 1:nrow(out)

	return( out)
}

