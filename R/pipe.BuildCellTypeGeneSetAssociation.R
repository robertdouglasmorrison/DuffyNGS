# pipe.BuildCellTypeGeneSetAssociation.R -- turn a CellType TargetMatrix and DE results into a mapping
#				from gene sets to cell types, thus a "cell type for each gene set" resource
#
#		Note that Gene Set to Cell Type mappings are now based primarily on
#		differential expression results from the K-wise DE comparison

`pipe.BuildCellTypeGeneSetAssociation` <- function( reference=NULL, folderName="", optionsFile="Options.txt",
						results.path=NULL, min.fold=0.1, min.expression=1, max.pvalue=0.05) {

	speciesID <- getCurrentSpecies()
	prefix <- getCurrentSpeciesFilePrefix()
	
	# use the target matrix of gene expression from the same reference as the starting data
	if ( is.null( reference)) stop( "'reference' argument of a cell type target matrix must be specified")
	CellTypeSetup( reference=reference, optionsFile=optionsFile)
	targetM <- getCellTypeMatrix()
	cellNames <- colnames(targetM)
	geneNames <- rownames(targetM)
	NG <- nrow(targetM)
	NC <- ncol(targetM)

	# use the DE results folder path, to extract all the DE details we want for every gene
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound="results", verbose=F)
	}
	de.path <- file.path( results.path, "MetaResults", paste( prefix, folderName, sep="."), "MetaGeneSets")
	if ( ! file.exists( de.path)) {
		cat( "\nError:  unable to find the MetaGeneSets folder: ", de.path)
		return(NULL)
	}
	deFiles <- dir( de.path, pattern=paste( prefix, "MetaGeneSets.UP.txt$", sep="."), full.name=T)
	if ( ! length( deFiles)) {
		cat( "\nError:  unable to find any 'MetaGeneSet.UP.txt' DE results files")
		return(NULL)
	}
	deNames <- sub( paste( ".", prefix, ".MetaGeneSets.UP.txt", sep=""), "", basename(deFiles))
	if ( !all( cellNames %in% deNames)) {
		cat( "\nError:  failed to find all Cell Type names in DE meta gene set results: ")
		cat( "\n  Missing: ", setdiff( cellNames, deNames))
		return(NULL)
	}


	# local function that does all the work...
	extractGeneSetCellRanks <- function( defile, min.fold=0.1, max.pvalue=0.05) {

		# Collect gene sets that are up-regulated, using thresholds for both DE fold and significance
		# extract the cell type  from the filename
		myCellType <- sub( ".MetaGeneSets.UP.txt", "", basename(defile))
		myCellType <- sub( paste( ".", prefix, sep=""), "", myCellType)
	
		# skip any that do not belong to the target
		if ( ! (myCellType %in% cellNames)) {
			cat( "\nWarning:  skipping DE results not in target matrix: ", myCellType)
			return(NULL)
		}

		tbl <- read.delim( defile, as.is=T)
		myGeneSets <- cleanGeneSetModuleNames( tbl$PathName, wrapParen=F)
		fold <- tbl$LOG2FOLD
		pval <- tbl$AVG_PVALUE
		ngenes <- tbl$N_Genes
		# no expression level in these files

		# we only want those gene sets that are UP and significant enough
		keep <- which( fold >= min.fold & pval <= max.pvalue)

		out <- data.frame( "GeneSetName"=myGeneSets[keep], "CellType"=myCellType, 
				"N_Genes"=ngenes[keep], "Fold"=round( fold[keep], digits=3), 
				"Pvalue"=pval[keep], "Rank"=1:length(keep), 
				stringsAsFactors=FALSE)

		# never let any duplicates stay
		dups <- which( duplicated( out$GeneSetName))
		if ( length(dups)) out <- out[ -dups, ]

		cat( "\n", myCellType, nrow(out))
		
		return( out)
	}

	# OK, run the gene set collector
	cat( "\nExtracting Gene Set cell types from DE folder: ", de.path)
	ans <- data.frame()
	for ( i in 1:length( deFiles)) {
		sml <- extractGeneSetCellRanks( deFiles[i], min.fold=min.fold, max.pvalue=max.pvalue)
		if( ! is.null(sml) && nrow(sml)) ans <- rbind( ans, sml)
	}
	
	SAV_ANS <<- ans

	# get the giant set of all gene sets as a flat list
	gsObj <- gatherGeneSets( defaultGeneSets())
	allgenesets <- gsObj[[3]][[1]]
	geneSetNames <- cleanGeneSetModuleNames( names( allgenesets), wrapParen=F)

	# now we can factor on the gene set name, and find which cell type(s) is best for each
	cat( "\nFactor for best from each..")
	genesetFac <- factor( ans$GeneSetName)
	ansPct <- rep.int( 0, nrow(ans))

	# with K dimensions, so make threshold be 2x that to count as having a cell that is best
	MIN_PCT <- 100 / NC * 2
	MIN_PCT_ANY <- 100 / NC
	MAX_GENE_SETS <- round( length( allgenesets) / NC)


	bestAnsRows <- tapply( 1:nrow(ans), genesetFac, function(x) {

			myGS <- ans$GeneSetName[x[1]]
			myFolds <- ans$Fold[x]
			myPvals <- ans$Pvalue[x]
			myCells <- ans$CellType[x]
			myRanks <- ans$Rank[x]
			
			# to know the relative expression, we need to build that subset from global data
			whereGS <- match( myGS, geneSetNames)
			myGenes <- allgenesets[[whereGS]]
			whereG <- match( myGenes, rownames(targetM), nomatch=0)

			# if not enough genes in the target matrix, it is by definition very low expression
			if ( sum( whereG > 0) > 0) {
				smlM <- targetM[ whereG, , drop=F]
				avgExpress <- apply( smlM, 2, sqrtmean, na.rm=T)
				if ( max(avgExpress) < min.expression) return(0)
				allPcts <- round( avgExpress * 100 / sum(avgExpress,na.rm=T), digits=1)
			} else {
				allPcts <- rep.int( round(100 / ncol(targetM), digits=1), ncol(targetM))
			}
			
			# save the pcts that go with these cell types, and order by expression (not fold)
			use <- match( myCells, colnames(targetM), nomatch=NA)
			myPcts <- allPcts[use]
			# gracefully catch missing cell types
			if ( all( is.na( myPcts))) return(0)
			
			ansPct[ x] <<- myPcts
			myOrder <- order( myPcts, decreasing=T)
			#myOrder <- order( myFolds, decreasing=T)

			# allow a gene set to record more than one cell type?
			# and if not at least 2x the background, skip it completely
			# Let's be a bit looser here, such that if not over the threshold perhaps return the one highest anyway
			if ( max(myPcts,na.rm=T) < MIN_PCT) {
				if ( max(myPcts,na.rm=T) > MIN_PCT_ANY) {
					return( x[ myOrder[1]])
				}
				return( 0)
			}

			# if any above the minimum, send back all that are
			keep <- which( myPcts >= MIN_PCT)
			myOrder <- order( myPcts[keep], decreasing=T)
			#myOrder <- order( myFolds[keep], decreasing=T)
			return( x[ keep[ myOrder]])
		})
		
	# store this Pct by cell type data into the global version
	ans$PctExpression <- round( ansPct)

	# might have a list now, so revert to a vector
	bestAnsRows <- unlist( bestAnsRows)

	# extract the subset of the big answer, that is just the best cell type(s) of each
	geneSetCellTypes <- ans[ bestAnsRows, ]
	ord <- order( toupper( geneSetCellTypes$GeneSetName))
	geneSetCellTypes <- geneSetCellTypes[ ord, ]

	# don't let certain junk get through...
	drops <- which( geneSetCellTypes$GeneSetName == "")
	if ( length(drops)) geneSetCellTypes <- geneSetCellTypes[ -drops, ]
	rownames(geneSetCellTypes) <- 1:nrow(geneSetCellTypes)
	
	# lastly, don't let any cell type have too many gene sets
	keepRows <- tapply( 1:nrow(geneSetCellTypes), factor( geneSetCellTypes$CellType), function(x) {
				# given the rows for one cell type, decide which are the best UP ones
				if ( length(x) <= MAX_GENE_SETS) return(x)
				
				ord <- order( geneSetCellTypes$Rank[x], -geneSetCellTypes$PctExpression[x])
				return( x[ord[ 1:MAX_GENE_SETS]])
			})
	geneSetCellTypes <- geneSetCellTypes[ sort( unlist(keepRows)), ]
	rownames(geneSetCellTypes) <- 1:nrow(geneSetCellTypes)

	# give it the final name
	cat( "\nSaving gene set association results..")
	finalFile <- paste( reference, "GeneSetAssociation.txt", sep=".")
	write.table( geneSetCellTypes, finalFile, sep="\t", quote=F, row.names=F)
	finalFile <- paste( reference, "GeneSetAssociation.rda", sep=".")
	save( geneSetCellTypes, file=finalFile)

	cat( "\nDone.  \nCopy this new 'GeneSetAssociation.rda' file to your DuffyTools/data/ subfolder and then remake the R package.\n")
	cat( "\nNote:  to fully incorporate these newly created 'Top N' cell type gene sets, then next also rerun the 'pipe.MetaGeneSets()'",
		"steps and then rerun this 'pipe.BuildCellTypeGeneSetAssociation()' tool.  This will add the new 'Top N' gene sets to the",
		"final gene set association database.")

	return( invisible( geneSetCellTypes))
}
