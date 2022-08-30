# pipe.BuildCellTypeTargetmatrx.R -- turn a set of samples into a CellType TargetMatrix 
#
#		This is intended to create the gene expression dataset used by all the Cell Type tools
#		This is called the 'target matrix', a table of gene expression that defines an
#		K-dimensional surface.


`pipe.BuildCellTypeTargetMatrix` <- function( sids, reference, annotationFile="Annotation.txt", optionsFile="Options.txt",
						groupColumn="Group", colorColumn="Color", results.path=NULL,
						expressionUnits="TPM_M", min.expression=1.0, finalColumnOrder=NULL) {
	
	annT <- readAnnotationTable( annotationFile)
	if ( ! all( sids %in% annT$SampleID)) stop( "All sample IDs not found in annotation file")
	speciesID <- getCurrentSpecies()
	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound="results", verbose=F)
	}
	
	cat( "\nGathering gene expression..")
	files <- file.path( results.path, "transcript", paste( sids, prefix, "Transcript.txt", sep="."))
	expM <- expressionFileSetToMatrix( files, fids=sids, intensityColumn=expressionUnits, verbose=T)
	rownames( expM) <- shortGeneName( rownames(expM), keep=1)
	bigM <- expM

	cat( "\nCombining replicates..")
	grps <- annT[[ groupColumn]][ match( sids, annT$SampleID)]
	cFac <- factor( grps)
	NC <- nlevels(cFac)
	NG <- nrow(bigM)
	newM <- matrix( NA, nrow=NG, ncol=NC)
	colnames(newM) <- levels(cFac)
	rownames(newM) <- rownames( bigM)
	for ( i in 1:nrow(newM)) newM[ i, ] <- tapply( bigM[i, ], cFac, sqrtmean)
	targetM <- newM
	rm( expM, bigM, newM)

	cat( "\nDropping very low intensity genes..")
	maxV <- apply( targetM, 1, max, na.rm=T)
	whoTooLow <- which( maxV < min.expression)
	targetM <- targetM[ -whoTooLow, ]
	cat( "   N_Drop for too low expression (", expressionUnits, "<", min.expression, ") =", length( whoTooLow))
	
	# since we are using the short gene names as the rownames, make sure we only keep one copy of each
	cat( "\nRemoving duplicate named genes..")
	dups <- which( duplicated( rownames( targetM)))
	if ( length( dups)) {
		cat( "   N_Duplicate row names =", length( dups), "\n")
		targetM <- targetM[ -dups, ]
	}
	
	# if we are using TPM units, sum of all genes should be 1 million
	# averaging above may have distorted that.  Reimpose that ideal.
	if (grepl( "TPM", expressionUnits)) {
		cat( "\nRescaling TPM to 1M total reads")
		for ( k in 1:NC) {
			v <- targetM[ , k]
			sf <- 1000000 / sum(v)
			targetM[ , k] <- v * sf
		}
	}

	# rearrange into the order wanted
	if ( is.null( finalColumnOrder)) finalColumnOrder <- 1:NC
	targetM <- targetM[ , finalColumnOrder]
	cat( "\nFinal Target Matrix Columns: ", NC, "\n  ", colnames(targetM))

	# we can round the precision a bit, as it's all RPKM or TPM
	targetM <- round( targetM, digits=3)

	# define the colors for each cell type
	targetColors <- annT[[ colorColumn]][ match( colnames(targetM), annT[[ groupColumn]])]
	names(targetColors) <- colnames(targetM)

	cat( "\nSaving files..")
	finalFilename <- paste( prefix, reference, "TargetMatrix.rda", sep=".")
	save( targetM, targetColors, file=finalFilename)
	write.table( targetM, sub( "rda$", "txt", finalFilename), sep="\t", quote=F, row.names=T)

	# make a few images and files
	checkX11()
	cat( "\nPlotting cluster images..")
	plot( expressionCluster(targetM), which=2, main=paste( "Cell Type Expression data: ", reference))
	printPlot( sub( "TargetMatrix.rda$", "ExpresssionClustering", finalFilename))
	plot( expressionCluster( log2(targetM+1)), which=2, main=paste( "Cell Type Log2 Expression data: ", reference))
	printPlot( sub( "TargetMatrix.rda$", "Log2.ExpresssionClustering", finalFilename))

	# make the correlation table too
	cat( "\nEvaluating Correlations..")
	ccm <- matrix( 0, NC, NC)
	colnames(ccm) <- rownames(ccm) <- colnames(targetM)
	for ( i in 1:NC) for( j in 1:NC) ccm[i,j] <- cor( targetM[,i], targetM[,j])
	ccm <- round( ccm, digits=3)
	write.csv( ccm, sub( "TargetMatrix.rda$", "CorrelationMatrix.csv", finalFilename), row.names=T)
	for ( i in 1:NC) for( j in 1:NC) ccm[i,j] <- cor( log2(targetM[,i]+1), log2(targetM[,j]+1))
	ccm <- round( ccm, digits=3)
	write.csv( ccm, sub( "TargetMatrix.rda$", "Log2.CorrelationMatrix.csv", finalFilename), row.names=T)

	# also save orthologged versions for our other systems
	saveGenes <- rownames(targetM)
	saveM <- targetM
	
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
			destGenes <- ortholog( saveGenes, from=speciesID, to=destSpecies)
			keep <- which( destGenes != "")
			if ( length(keep) < 10) { 
				cat( "\n ", destSpecies, "  Too few ortholog genes..  Skip.")
				next
			}
			targetM <- saveM[ keep, ]
			rownames(targetM) <- destGenes[keep]
			# orthologging could induce duplicates
			targetM <- targetM[ ! duplicated( rownames(targetM)), ]
			# revert to alphabetical ordering
			targetM <- targetM[ order( rownames(targetM)), ]
			# since we are using TPM units, sum of all genes should be 1 million.  Reimpose that ideal.
			if (grepl( "TPM", expressionUnits)) {
				for ( k in 1:NC) {
					v <- targetM[ , k]
					sf <- 1000000 / sum(v)
					targetM[ , k] <- v * sf
				}
			}
			cat( "\n ", destSpecies, "  N_Genes =", nrow(targetM))
			save( targetM, targetColors, file=paste( orthoPrefix[j], reference, "TargetMatrix.rda", sep="."))
		}
	}

	# all done
	targetM <- saveM
	setCurrentSpecies( speciesID)
	
	cat( "\nDone.  \nCopy these new 'TargetMatrix.rda' files to your DuffyTools/data/ subfolder and then remake the R package.\n")
	return( invisible( targetM))
}


