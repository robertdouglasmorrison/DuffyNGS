# calcAlignSummary.R

# feel up the distribution of alignments to the various species and genes

`calcAlignSummary` <- function( mode=c("setup", "addData", "report"), 
			chunk, geneIDs, speciesIDs, filename, rawReadCount=NULL,
			alignPhase=c("Genomic","RiboClear", "Splicing"), verbose=!interactive()) {

	mode <- match.arg( mode)
	alignPhase <- match.arg( alignPhase)


	#local functions to get just uniques from a multi-read set
	`forceUnique` <- function(x) { if ( length(x) == 1) x else sort.int( unique.default(x)) } 
	`testUnique` <- function(x) { if ( length(x) == 1) TRUE else FALSE } 
	

	# start building the result
	if ( mode == "setup") {

		outText <- vector()
		outText <- c( paste( "\n\nSummary: \t\t\t\tAlignment Step:  ", alignPhase), 
				paste( "\nFilename:            \t\t", filename, "\n"))

		readStatsOutText <<- outText
		readStatsNreads <<- readStatsNaligns <<- 0
		readStatsNmultiS <<- 0
		readStatsGenePattTable <<- NULL
		readStatsSpeciesTable <<- NULL
		readStatsRawReadCount <<- rawReadCount
		readStatsAlignPhase <<- alignPhase
		return()
	}

	if (mode == "addData") {

		if ( ! is.bamChunk(chunk)) {
			cat( "\nWarning: BAM content does not look unsorted. Unable to assess BAM alignment summary")
			return(NULL)
		}
		N <- size( chunk)
		cat( "  summarize..")
		fac <- factor.align( chunk)
		nReads <- nlevels( fac)
		PASTE <- base::paste

		# let's try to do this faster
		geneIDpatterns <- nSpecies <- speciesIDpatterns <- vector( length=nReads)
		nNow <- 0
		if ( readStatsAlignPhase != "Splicing") geneIDs <- shortGeneName(geneIDs, keep=1)
		base::tapply( 1:N, INDEX=fac, function(x) {
				x1 <- x[1]
				nNow <<- nNow + 1
				if ( length(x) > 1) {
					geneIDpatterns[nNow] <<- PASTE( sort.int(unique.default(geneIDs[x])), collapse="|")
					mySpecies <- sort.int(unique.default(speciesIDs[x]))
					nSpecies[nNow] <<- ns <- length( mySpecies)
					if ( ns == 1) speciesIDpatterns[nNow] <<- mySpecies
				} else {
					geneIDpatterns[nNow] <<- geneIDs[x1]
					speciesIDpatterns[nNow] <<- speciesIDs[x1]
					nSpecies[nNow] <<- 1
				}
			})
		genePattTable <- base::table( geneIDpatterns)
		speciesTable <- base::table( speciesIDpatterns[ nSpecies == 1])
		nMultiSpecies <- sum( nSpecies > 1)

		# old way...
		#if ( readStatsAlignPhase == "Splicing") {
		#	geneIDsets <- base::tapply( geneIDs, INDEX=fac, forceUnique, simplify=FALSE)
		#} else {
		#	geneIDsets <- base::tapply( shortGeneName(geneIDs, keep=1), INDEX=fac, forceUnique, simplify=FALSE)
		#}
		#geneIDpatterns <- sapply( geneIDsets, base::paste, collapse="|")
		#genePattTable <- base::table( geneIDpatterns)
		#speciesIDsets <- tapply( speciesIDs, INDEX=fac, forceUnique, simplify=FALSE)
		#uniqueSpecies <- which( sapply( speciesIDsets, FUN=testUnique, simplify=TRUE))
		#speciesTable <- base::table( base::unlist( speciesIDsets[ uniqueSpecies]))
		#nMultiSpecies <- nReads - length( uniqueSpecies)

		# push these results back to storage
		readStatsNreads <<- readStatsNreads + nReads
		readStatsNaligns <<- readStatsNaligns + N
		readStatsNmultiS <<- readStatsNmultiS + nMultiSpecies
		if ( is.null( readStatsGenePattTable)) {
			readStatsGenePattTable <<- genePattTable
		} else {
			readStatsGenePattTable <<- mergeTables( readStatsGenePattTable, genePattTable)
		}
		if ( is.null( readStatsSpeciesTable)) {
			readStatsSpeciesTable <<- speciesTable
		} else {
			readStatsSpeciesTable <<- mergeTables( readStatsSpeciesTable, speciesTable)
		}
		return()
	}

	# reporting...
	outText <- readStatsOutText
	nReads <- readStatsNreads
	nAligns <- readStatsNaligns
	nMultiSpecies <- readStatsNmultiS
	genePattTable <- readStatsGenePattTable
	speciesTable <- readStatsSpeciesTable
	if ( is.null( genePattTable)) {
		cat( "Warning: no gene data collected from BAM file.")
		return( NULL)
	}

	outText <- base::append( outText, base::paste( "\nN_Reads:             \t", 
				formatC(nReads, format="d", width=12, big.mark=",")))
	outText <- base::append( outText, base::paste( "\nN_Alignments:        \t", 
				formatC(nAligns, format="d", width=12, big.mark=","),"\n"))

	# tell percents of the original file if that count is given...
	if ( ! is.null( rawReadCount)) readStatsRawReadCount <<- rawReadCount
	nRawReads <- if ( is.null( readStatsRawReadCount)) nReads else readStatsRawReadCount

	# gather how many hit various species
	if ( ! is.null( speciesTable)) {
	    for ( i in 1:length( speciesTable)) {
		outText <- base::append( outText, base::paste( "\nReads Hitting Only:   ", 
				format( names(speciesTable)[i], width=10, justify="right"), " \t", 
				formatC( speciesTable[i], width=10, big.mark=",", format="d"), 
				"\t", as.percent( speciesTable[i], big.value=nReads),
				"\t", as.percent( speciesTable[i], big.value=nRawReads)))
	    }
	}
	outText <- base::append( outText, base::paste( "\nMore Than One Species:            \t", 
				formatC( nMultiSpecies, width=10, big.mark=",", format="d"), "\t", 
				as.percent( nMultiSpecies, big.value=nReads),
				"\t", as.percent( nMultiSpecies, big.value=nRawReads),"\n"))

	# show the summary, with a special format for ribo cleared reads
	nAllGlobin <- 0
	if ( readStatsAlignPhase == "RiboClear") {

	    cat( "\nSummarize by group..")

	    if ( ! is.null(speciesTable) && ! is.null( genePattTable)) {
		# for all gene level counts, we can count up how often each gene pattern occurs and do each once
		gpLevels <- names( genePattTable)
		NgpLevels <- length( genePattTable)
		gpCnts <- genePattTable

		# changing from GREP method to set overlap.  Some species have many riboclear genes, that
		# makes grep pattern string a huge logical or.
		gpLevelGenes <- strsplit( gpLevels, split="|", fixed=T)
		# we may need to mask characters that have special meaning to 'grep'
		#gpLevels <- gsub( "(", ":", gpLevels, fixed=T)
		#gpLevels <- gsub( ")", ":", gpLevels, fixed=T)
		#gpLevels <- gsub( "[", ":", gpLevels, fixed=T)
		#gpLevels <- gsub( "]", ":", gpLevels, fixed=T)

	    	for( i in 1:length( speciesTable)) {
			thisSpecies <- names( speciesTable)[i]
			setCurrentSpecies( thisSpecies)
			rrnaMap <- getCurrentRrnaMap()
			if ( ! ( "GROUP" %in% colnames( rrnaMap))) next
			# only report on the things we are trying to clear
			if ( "CLEAR" %in% colnames( rrnaMap)) {
				keep <- which( as.logical( rrnaMap$CLEAR))
				if ( length( keep) > 0) {
					rrnaMap <- rrnaMap[ keep, ]
				}
			}
			rrna <- subset.data.frame( rrnaMap, subset=(GROUP != ""), select=c(GENE_ID,GROUP))
			grps <- factor( rrna$GROUP)
			ptrs <- tapply( 1:nrow(rrna), INDEX=grps, FUN=NULL)
			rrnaGrpGenes <- tapply( shortGeneName(rrna$GENE_ID,keep=1), grps, c, simplify=FALSE)
			for( ig in 1:nlevels(grps)) {
				thisGrp <- levels(grps)[ig]
				cat( " ", thisGrp)
				# changing from GREP method to set overlap.  Some species have many riboclear genes, that
				# makes grep pattern string a huge logical or.
				#theseGenes <- base::paste( shortGeneName( rrna$GENE_ID[ ptrs == ig], keep=1), collapse="|")
				#theseGenes <- gsub( "(", ":", theseGenes, fixed=T)
				#theseGenes <- gsub( ")", ":", theseGenes, fixed=T)
				#theseGenes <- gsub( "[", ":", theseGenes, fixed=T)
				#theseGenes <- gsub( "]", ":", theseGenes, fixed=T)
				#hits <- which( sapply( gpLevels, function(x) return( length( grep( theseGenes, x)) > 0)))
				thisGeneList <- rrnaGrpGenes[[ig]]
				hits <- which( sapply( gpLevelGenes, function(x) return( sum( x %in% thisGeneList) > 0)))
				nhits <- sum( gpCnts[hits])
				outText <- base::append( outText, base::paste( "\nCleared: ", format(thisSpecies, width=8), " ",
						format( thisGrp, width=10), "\t", 
						formatC( nhits, big.mark=",", width=10, format="d"), "\t",
						as.percent( nhits, big.value=nReads),
						"\t", as.percent( nhits, big.value=nRawReads)))
				if( thisGrp == "Globin") nAllGlobin <- nAllGlobin + nhits
			}
	    	}
		out <- list( "textSummary"=outText, "nReads"=nReads, "SpeciesTable"=speciesTable, 
				"nMultiSpecies"=nMultiSpecies, "nGlobin"=nAllGlobin)
	    } else {
		out <- list( "textSummary"=outText, "nReads"=nReads, "SpeciesTable"=NULL, 
				"nMultiSpecies"=nMultiSpecies, "nGlobin"=nAllGlobin)
	    }

	} else if ( readStatsAlignPhase == "Genomic") {

		# only show the top N of these
		N <- min( 20, length( genePattTable))
		ord <- order( genePattTable, decreasing=T)
		genePattTable <- genePattTable[ ord[ 1:N]]
		maxCharShow <- 28
		if (N) for( ig in 1:N) {
			thisGrp <- names(genePattTable)[ig]
			if ( nchar(thisGrp) > maxCharShow) thisGrp <- paste( substr( thisGrp, 1, maxCharShow), "...", sep="")
			nhits <- genePattTable[ig]
			outText <- base::append( outText, base::paste( "\nTop Hits:\t", 
					format( thisGrp, width=maxCharShow+3), "\t", 
					formatC( nhits, big.mark=",", width=10, format="d"), "\t",
					as.percent( nhits, big.value=nReads),
					"\t", as.percent( nhits, big.value=nRawReads)))
		}
		nhits <- sum( genePattTable)
		outText <- base::append( outText, base::paste( "\n\nTop", N, "Genes Combined:\t\t\t\t\t",
				formatC( nhits, big.mark=",", width=10, format="d"), "\t",
				as.percent( nhits, big.value=nReads),
				"\t", as.percent( nhits, big.value=nRawReads), "\n"))

		out <- list( "textSummary"=outText, "nReads"=nReads, "SpeciesTable"=speciesTable, 
				"nMultiSpecies"=nMultiSpecies, "nGlobin"=nAllGlobin)

	} else {   # last case is 'Splicing'

		fullGenePattTable <- genePattTable
		if ( ! is.null( genePattTable)) {
			# only show the top N splicies
			N <- min( 20, length( genePattTable))
			ord <- order( genePattTable, decreasing=T)
			genePattTable <- genePattTable[ ord[ 1:N]]
			maxCharShow <- 32
			if (N) for( ig in 1:N) {
				thisGrp <- names(genePattTable)[ig]
				if ( nchar(thisGrp) > maxCharShow) thisGrp <- paste( substr( thisGrp, 1, maxCharShow), "...", sep="")
				nhits <- genePattTable[ig]
				outText <- base::append( outText, base::paste( "\nTop Hits:\t", 
						format( thisGrp, width=maxCharShow+3), "\t", 
						formatC( nhits, big.mark=",", width=10, format="d"), "\t",
						as.percent( nhits, big.value=nReads),
						"\t", as.percent( nhits, big.value=nRawReads)))
			}
			nhits <- sum( genePattTable)
			outText <- base::append( outText, base::paste( "\n\nTop", N, "Splices Combined:\t\t\t\t\t",
					formatC( nhits, big.mark=",", width=10, format="d"), "\t",
					as.percent( nhits, big.value=nReads),
					"\t", as.percent( nhits, big.value=nRawReads), "\n"))

			detailSummary <- spliceDetailSummary( fullGenePattTable)

			out <- list( "textSummary"=outText, "nReads"=nReads, "SpeciesTable"=speciesTable, 
					"nMultiSpecies"=nMultiSpecies, "nGlobin"=nAllGlobin, "spliceDetailSummary"=detailSummary )
		} else {
			out <- list( "textSummary"=outText, "nReads"=nReads, "SpeciesTable"=NULL, 
					"nMultiSpecies"=nMultiSpecies, "nGlobin"=nAllGlobin, "spliceDetailSummary"=NULL )
		}
	}

	# clean up
	try( rm( readStatsNreads, readStatsNaligns, readStatsNmultiS, readStatsSpeciesTable, 
			readStatsOutText, readStatsRawReadCount, readStatsGenePattTable, 
			readStatsAlignPhase, envir=.GlobalEnv))

	return( out)
}

