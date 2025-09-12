# pipe.BaselineAdjustTranscripts.R

# scale a set of transcriptomes to have identical time zero expression levels

`pipe.BaselineAdjustTranscripts` <- function( sampleIDs, subjectIDs, timeIDs, baselineID=timeIDs[1], 
				adjusted.path="adjustedResults/transcript",
				annotationFile="Annotation.txt", optionsFile="Options.txt", 
				speciesID=getCurrentSpecies(), results.path=NULL, 
				max.scale.factor=4, min.expression.for.scaling=1, expression.units=NULL,
				verbose=TRUE, debug.genes=NULL) {

	if ( speciesID != getCurrentSpecies()) setCurrentSpecies( speciesID)
	speciesPrefix <- getCurrentSpeciesFilePrefix()

	if ( is.null( results.path)) {
		optT <- readOptionsTable( optionsFile)
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	if ( ! file.exists( adjusted.path)) dir.create( adjusted.path, recursive=TRUE, showWarnings=FALSE)

	# get the expression units to use for scaling
	if ( is.null( expression.units)) expression.units <- getExpressionUnitsColumn( optionsFile)

	# get the set of existing transcripts
	old.path <- file.path( results.path, "transcript")
	fileSet <- file.path( old.path, paste( sampleIDs, speciesPrefix, "Transcript.txt", sep="."))

	if ( !all( file.exists( fileSet))) {
		cat( "\nError:  Some required transciptome files not found: \n")
		print( fileSet[ ! file.exists(fileSet)])
		return( NULL)
	}


	# grab all the original transcriptomes into a matrix to aid calculations
	cat( "\nGathering original transcriptomes..")
	expressM <- expressionFileSetToMatrix( fileSet, sampleIDs, intensityColumn=expression.units,
					missingGenes="na", keepIntergenic=F, verbose=T)
	
	# deduce the median Day 0 expression level for each gene
	whoTime0 <- which( timeIDs == baselineID)
	avgTime0 <- apply( expressM[ , whoTime0], MARGIN=1, sqrtmean, na.rm=T)
	avgGeneNames <- rownames( expressM)
	
	# now read every transcriptome, adjust all expression columns and write it out
	ExpressColumns <- c("RPKM_M", "SIGMA_M", "RPKM_U", "SIGMA_U", "READS_M", "READS_U", "TPM_M", "TPM_U")

	# the scaling is based on time zero, and applied to all other of the same subjectID
	allUsedScaleFactors <- vector()
	curSubject <- ""

	for ( i in 1:length(fileSet)) {
		thisFile <- fileSet[i]
		thisSubj <- subjectIDs[i]
		newFile <- sub( old.path, adjusted.path, thisFile, fixed=T)

		cat( "\n", i, "ReadFile:", sampleIDs[i])
		tbl <- read.delim( thisFile, as.is=T)
		myColumns <- match( ExpressColumns, colnames(tbl))
		if ( any( is.na( myColumns))) {
			cat( "\nUnable to find required columns in transcriptome file: \n")
			print( ExpressColumns[ is.na(myColumns)])
			cat( "\nSkipping..")
			next
		}
	
		# we will clip the uncertainies at their unscaled lower bound
		minSigma <- min( tbl$SIGMA_M, na.rm=T)
	
		# grab the time 0 expression values
		genePtr <- match( tbl$GENE_ID, avgGeneNames, nomatch=0)
		if ( curSubject != thisSubj) {
			cat( "  grabTime0..")
			whoMy0 <- which( subjectIDs == thisSubj & timeIDs == baselineID)[1]
			if ( is.na( whoMy0)) {
				cat( "\nNo Time=0 data for subject: ", thisSubj)
				cat( "\nUnable to baseline adjust..   Using unadjusted..")
				write.table( tbl, newFile, sep="\t", quote=F, row.names=F)
				next
			}
			myTime0 <- expressM[ , whoMy0]
			curSubject <- thisSubj
		}
	
		# calculate the scale factor for every gene, but only on the time 0 samples
		cat( "  calcScaling..")
		scaleFac <- avgTime0 / myTime0
		# allow scale clipping
		if ( ! is.null( max.scale.factor)) {
			if ( max.scale.factor < 1.0) stop( "'max.scale.factor' must be greater than 1.0")
			isBIG <- which( scaleFac > max.scale.factor)
			if ( length(isBIG)) scaleFac[isBIG] <- max.scale.factor
			isSMALL <- which( scaleFac < (1.0/max.scale.factor))
			if ( length(isSMALL)) scaleFac[isSMALL] <- (1.0/max.scale.factor)
		}
		# prevent divide by zero and extreme scaling
		scaleFac[ myTime0 <= min.expression.for.scaling] <- 1
		allUsedScaleFactors <- c( allUsedScaleFactors, scaleFac)
	
		# make the new adjusted transcriptome
		cat( "  applyScaling..")
		newtbl <- tbl
		for ( j in myColumns) {
			v <- tbl[[j]]
			vnew <- v
			vnew[ genePtr > 0] <- v[ genePtr > 0] * scaleFac[genePtr]
			newtbl[[j]] <- vnew
		}
		newtbl$SIGMA_M <- ifelse( newtbl$SIGMA_M < minSigma, minSigma, newtbl$SIGMA_M)
		newtbl$SIGMA_U <- ifelse( newtbl$SIGMA_U < minSigma, minSigma, newtbl$SIGMA_U)

		# perhaps re-force the TPM units to sum to a million?
		# this seems to "undo" the effect of the adjustment.  Turn back off for now...
		#if ( "TPM_M" %in% colnames(newtbl)) {
			#newtbl$TPM_M <- tpm( newtbl$READS_M, newtbl$N_EXON_BASES)
			#newtbl$TPM_U <- tpm( newtbl$READS_U, newtbl$N_EXON_BASES)
		#}
		if ( all( c("TPM_M","RANK_M") %in% colnames(newtbl))) {
        	newtbl$RANK_M <- rank( newtbl$TPM_M, ties.method="min")
        	newtbl$RANK_U <- rank( newtbl$TPM_U, ties.method="min")
		}
		
		# inspect details about some genes?
		if ( ! is.null(debug.genes)) {
			for (dbgG in debug.genes) {
				whereTbl <- match( dbgG, shortGeneName(tbl$GENE_ID,keep=1))
				whereScale <- match( dbgG, shortGeneName(avgGeneNames,keep=1))
				cat( "\nDebug: ", dbgG, "(whT,whS,oldV,my0,avg0,scale,newV)", whereTbl, whereScale, tbl[[expression.units]][whereTbl],
					"|", myTime0[whereScale], avgTime0[whereScale], scaleFac[whereScale], "|", newtbl[[expression.units]][whereTbl])
			}
		}

		# put the new transcripome back into expression order
		exOrd <- order( newtbl$RPKM_M, decreasing=T)
		newtbl <- newtbl[ exOrd, ]
	
		cat( "  writeNewFile..")
		write.table( newtbl, newFile, sep="\t", quote=F, row.names=F)
	}
	cat( "\nDone.\n")

	# pass back the set of scaling used, invisibly
	out <- allUsedScaleFactors
	return( invisible(out))
}
