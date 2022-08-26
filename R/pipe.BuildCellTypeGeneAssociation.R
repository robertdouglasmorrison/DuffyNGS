# pipe.BuildCellTypeGeneAssociation.R -- turn a CellType TargetMatrix into gene associations
#				thus a "cell type for each GeneID" resource
#
#		Note that Gene to Cell Type mappings are now based on Differential Expression results,
#		no longer on differential expression results.  They are usually the same but need not be.

`pipe.BuildCellTypeGeneAssociation` <- function( de.folder, referenceName=NULL, optionsFile="Options.txt",
						results.path=NULL, min.fold=0.1, min.expression=1.0) {

	speciesID <- getCurrentSpecies()
	prefix <- getCurrentSpeciesFilePrefix()
	
	# use the target matrix of gene expression from the same reference as the starting data
	verifyCellTypeSetup( referenceName=referenceName, optionsFile=optionsFile)
	targetM <- getCellTypeMatrix()
	cellNames <- colnames(targetM)
	geneNames <- rownames(targetM)
	NG <- nrow(targetM)
	NC <- ncol(targetM)

	# use the DE results folder path, to extract all the DE details we want for every gene
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound="results", verbose=F)
	}
	de.path <- file.path( results.path, "MetaResults", paste( prefix, de.folder, sep="."))
	if ( ! file.exists( de.path)) {
		cat( "\nError:  unable to find the MetaResults folder: ", de.path)
		return(NULL)
	}
	deFiles <- dir( de.path, pattern=paste( prefix, "Meta.UP.txt$", sep="."), full.name=T)
	if ( ! length( deFiles)) {
		cat( "\nError:  unable to find any 'Meta.UP.txt' DE results files")
		return(NULL)
	}
	deNames <- sub( paste( prefix, "Meta.UP.txt", sep="."), "", basename(deFiles))
	if ( !all( cellNames %in% deNames)) {
		cat( "\nError:  failed to find all Cell Type names in DE results: ")
		cat( "\n  Missing: ", setdiff( cellNames, deNames))
		return(NULL)
	}


	# local function that does all the extraction work per cell type file...
	extractGeneCellRanks <- function( defile, min.fold=0.1, min.expression=1.0) {

		# Collect genes that are up-regulated, using thresholds for both DE fold and expression level
		# extract the cell type terms from the filename
		myCellType <- sub( paste( ".", prefix, "Meta.UP.txt", sep="."), "", basename(defile))

		tbl <- read.delim( defile, as.is=T)
		myGenes <- shortGeneName( tbl$GENE_ID, keep=1)
		fold <- tbl$LOG2FOLD
		expres <- tbl[[ myCellType]]

		# we only want those real genes that are UP and expressed high enough
		keep <- which( fold >= min.fold)
		keep2 <- which( expres >= min.expression)
		nonGenes <- grep( "(ng)", tbl$GENE_ID, fixed=T)
		keep <- sort( setdiff( intersect(keep,keep2), nonGenes))

		out <- data.frame( "GENE_ID"=myGenes[keep], "CellType"=myCellType, 
				"Expression"=expres[keep], "Fold"=round( fold[keep], digits=3), 
				"Pvalue"=tbl$AVG_PVALUE[keep], "Rank"=1:length(keep), 
				stringsAsFactors=FALSE)
				
		# never let any duplicates stay
		dups <- which( duplicated( out$GENE_ID))
		if ( length(dups)) out <- out[ -dups, ]

		cat( "\n", i, myCellType, nrow(out))
		return( out)
	}


	getManualGeneFixes <- function() {

		cat( "\nManually fixing hemoglobins, etc...")
		# since we filter out hemoglobin, it may not even be seen here
		hgbNames <- c( "HBA1", "HBA2", "HBB", "HBBP1", "HBD", "HBE1", "HBG1",  
					"HBG2", "HBM", "HBQ1", "HBZ", "HBZP1")
		if ( speciesID %in% MAMMAL_SPECIES && speciesID != "Hs_grc") {
			hgbNames <- ortholog( hgbNames, from="Hs_grc", to=speciesID)
			hgbNames <- setdiff( hgbNames, "")
			if ( ! length(hgbNames)) return( data.frame())
		}
		
		smlDF <- data.frame( "GENE_ID"=hgbNames, "CellType"="RBC", 
					"Expression"=100, "Fold"=10, 
					"Pvalue"=1e-10, "Rank"=1, 
					stringsAsFactors=FALSE)
		return( smlDF)
	}
	# end of local functions..


	# OK, run the gene collector
	cat( "\nExtracting Gene cell types..\n")
	ans <- data.frame()
	for ( i in 1:length( deFiles)) {
		sml <- extractGeneCellRanks( deFiles[i], min.fold=min.fold, min.expression=min.expression)
		if( nrow(sml)) ans <- rbind( ans, sml)
	}

	# apply manual currated fixes, like for HGB
	if (species %in% MAMMAL_SPECIES) {
		sml <- getManualGeneFixes()
		if( nrow(sml)) ans <- rbind( ans, sml)
	}
	
	# only keep those genes that are in the target matrix
	keep <- which( ans$GENE_ID %in% geneNames)
	ans <- ans[ keep, ]

	# now we can factor on the gene name, and find which cell type(s) is best for each gene
	cat( "\nFactor for best cell types for each gene..")
	geneFac <- factor( ans$GENE_ID)
	ansPct <- rep.int( NA, nrow(ans))

	# use the global gene expression data where possible
	pctM <- matrix( 0, nrow(targetM), ncol(targetM))
	rownames(pctM) <- rownames(targetM)
	colnames(pctM) <- colnames(targetM)
	for ( i in 1:nrow(targetM)) { x <- targetM[ i, ]; pctM[ i, ] <- round( x * 100 / sum(x), digits=1)}

	# define some minimum metrics about if a gene can be called cell associated
	MIN_PCT <- 100 / NC * 2
	MIN_PCT_ANY <- 100 / (NC-1)

	bestAnsRows <- tapply( 1:nrow(ans), geneFac, function(x) {

			myExpress <- round( ans$Expression[x], digits=1)
			myFolds <- ans$Fold[x]
			myPvals <- ans$Pvalue[x]
			myRanks <- ans$Rank[x]
			myCells <- ans$CellType[x]
			myPcts <- round( myExpress * 100 / sum(myExpress,na.rm=T), digits=1)
			
			# whenever possible, use the global expression and Percentage data to know these key details
			myGene <- ans$GENE_ID[x[1]]
			whereGene <- match( myGene, geneNames, nomatch=0)
			if ( whereGene > 0) {
				whereCell <- match( myCells, cellNames)
				myExpress <- targetM[ whereGene, whereCell]
				myPcts <- pctM[ whereGene, whereCell]
			}
			
			# save those pcts, and order by expression, more than by fold change alone
			ansPct[ x] <<- myPcts
			myOrder <- order( myExpress, myFolds, decreasing=T)

			# allow a gene to record more than one cell type?
			# and if not at least 2x the expected background, skip it completely
			# Let's be a bit looser here, such that if not over the threshold perhaps return the one highest anyway
			if ( max(myPcts) < MIN_PCT) {
				if ( max(myPcts) > MIN_PCT_ANY) {
					return( x[ myOrder[1]])
				}
				return( 0)
			}

			# if any above the minimum, send back all that are
			keep <- which( myPcts >= MIN_PCT)
			myOrder <- order( myExpress[keep], myFolds[keep], decreasing=T)
			return( x[ keep[ myOrder]])
		})
		
	# store this Pct by cell type data into the global version
	ans$PctExpression <- round( ansPct)

	# might have a list now, so revert to a vector
	bestAnsRows <- unlist( bestAnsRows)

	# extract the subset of the big answer, that is just the best cell type(s) of each
	geneCellTypes <- ans[ bestAnsRows, ]
	ord <- order( toupper( geneCellTypes$GENE_ID))
	geneCellTypes <- geneCellTypes[ ord, ]

	# don't let certain junk get through...
	drops <- which( geneCellTypes$GENE_ID == "")
	if ( length(drops)) geneCellTypes <- geneCellTypes[ -drops, ]
	rownames(geneCellTypes) <- 1:nrow(geneCellTypes)

	# give it the final object name
	# .CSV ruins Human gene names!!
	cat( "\nSaving gene association results..")
	finalFile <- paste( prefix, referenceName, "GeneAssociation.txt", sep=".")
	write.table( geneCellTypes, finalFile, sep="\t", quote=F, row.names=F)
	save( geneCellTypes, file="Hs.CoreGeneCellTypes.rda")
	finalFile <- paste( prefix, referenceName, "GeneAssociation.rda", sep=".")
	save( geneCellTypes, file=finalFile)


	# also make a "allGeneSets" gene set version that keeps all genes
	# Refine this idea a bit, to put an upper limit on how many genes can be in each cell type
	# and a minimum fold change to be in the gene set
	max.genes <- nrow( subset( getCurrentGeneMap(), REAL_G == TRUE)) / NC
	
	allGeneSets <- tapply( 1:nrow(geneCellTypes), factor( geneCellTypes$CellType), function(x) {

			#order by fold change vs other cell types
			myFold <- geneCellTypes$Fold[x]
			myExpress <- geneCellTypes$Expression[x]
			myGenes <- geneCellTypes$GENE_ID[x]
			keep <- which( myFold >= min.fold & myExpress >= min.expression)
			myFold <- myFold[keep]
			myExpress <- myExpress[keep]
			myGenes <- myGenes[keep]
			ord <- order( myFold, myExpress, decreasing=T)
			finalGenes <- myGenes[ord]
			if ( length(finalGenes) > max.genes) finalGenes <- finalGenes[ 1:max.genes]
			return( unique( finalGenes))
		})
	save( allGeneSets, file=paste( prefix, referenceName, "GeneSets.rda", sep="."))


	# lastly, make orthologged versions...
	saveGeneCellTypes <- geneCellTypes
	
	orthoSpecies <- orthoPrefix <- vector()
	if ( speciesID %in% MAMMAL_SPECIES) {
		orthoSpecies <- MAMMAL_SPECIES
		orthoPrefix <- MAMMAL_PREFIXES
	} else if ( speciesID %in% PARASITE_SPECIES) {
		orthoSpecies <- PARASITE_SPECIES
		orthoPrefix <- PARASITE_PREFIXES
	} else if ( speciesID %in% BACTERIA_SPECIES) {
		orthoSpecies <- BACTERIA_SPECIES
		orthoPrefix <- BACTERIA_PREFIXES
	}

	if ( length( orthoSpecies) > 0) {
		cat( "\nOrthologging..")
		for (j in 1:length(orthoSpecies)) {
			destSpecies <- orthoSpecies[j]
			if ( destSpecies == speciesID) next
			geneCellTypes <- saveGeneCellTypes
			destGenes <- ortholog( saveGeneCellTypes$GENE_ID, from=speciesID, to=destSpecies)
			keep <- which( destGenes != "")
			if ( length(keep) < 10) { 
				cat( "\n ", destSpecies, "  Too few ortholog genes..  Skip.")
				next
			}
			geneCellTypes <- geneCellTypes[ keep, ]
			geneCellTypes$GENE_ID <- destGenes[ keep]
			
			# orthologging could induce duplicates
			key <- paste( geneCellTypes$GENE_ID, geneCellTypes$CellType, sep=":")
			drops <- which( duplicated( key))
			if ( length(drops)) geneCellTypes <- geneCellTypes[ -drops, ]
			# revert to alphabetical ordering
			geneCellTypes <- geneCellTypes[ order( geneCellTypes$GENE_ID), ]
			rownames(geneCellTypes) <- 1:nrow(geneCellTypes)

			finalFile <- paste( orthoPrefix[j], referenceName, "GeneAssociation.rda", sep=".")
			cat( "\n ", destSpecies, "  N_Genes =", length( unique( geneCellTypes$GENE_ID)))
			save( geneCellTypes, file=finalFile)

		}
	}

	# when all done, restore the human answer
	geneCellTypes <- saveGeneCellTypes
	return( invisible( geneCellTypes))
}
