# pipe.ScaleFactorAdjustTranscripts.R

# scale a set of transcriptomes given an explicit set of gene scaling factors
# most typically used to compensate for some types of global DNA/RNA amplification methods

`pipe.ScaleFactorAdjustTranscripts` <- function( sampleIDset,  scaleFactorFile, geneColumn="GENE_ID", 
				scaleFactorColumn="SCALE_FACTOR", sep="\t", adjusted.path="adjustedResults/transcript",
				annotationFile="Annotation.txt", optionsFile="Options.txt", speciesID=getCurrentSpecies(), results.path=NULL, 
				min.scale.factor=0.001, max.scale.factor=100, 
				expression.units=NULL, verbose=TRUE, debug.genes=NULL, exemptSampleIDset=NULL) {

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
	fileSet <- file.path( old.path, paste( sampleIDset, speciesPrefix, "Transcript.txt", sep="."))
	if ( !all( file.exists( fileSet))) {
		cat( "\nError:  Some required transciptome files not found: \n")
		print( fileSet[ ! file.exists(fileSet)])
		return( NULL)
	}

	# read in the scale factor data
	if ( ! file.exists( scaleFactorFile)) {
		cat( "\nError:  Scale Factor file not found: ", scaleFactorFile)
		return( NULL)
	}
	sf <- read.delim( scaleFactorFile, as.is=T, sep=sep)
	if ( ! all( c( geneColumn, scaleFactorColumn) %in% colnames(sf))) {
		cat( "\nError:  Some needed columns not found.  Looked for: ", geneColumn, scaleFactorColumn)
		cat( "\n        Found: ", colnames(sf))
		return( NULL)
	}
	scaleGenes <- as.character( sf[[ geneColumn]])
	scaleFactors <- as.numeric( sf[[ scaleFactorColumn]])
	scaleFactors[ is.na(scaleFactors)] <- 1
	tooBig <- which( scaleFactors > max.scale.factor)
	if ( length(tooBig)) {
		cat( "\nWarn: some scale factors larger than", max.scale.factor, "  N =", length(tooBig))
		scaleFactors[ tooBig] <- max.scale.factor
	}
	tooSmall <- which( scaleFactors < min.scale.factor)
	if ( length(tooSmall)) {
		cat( "\nWarn: some scale factors smaller than", min.scale.factor, "  N =", length(tooSmall))
		scaleFactors[ tooSmall] <- min.scale.factor
	}
	NSG <- length(scaleGenes)
		
	# now read every transcriptome, adjust all expression columns and write it out
	ExpressColumns <- c("RPKM_M", "SIGMA_M", "RPKM_U", "SIGMA_U", "READS_M", "READS_U", "TPM_M", "TPM_U")

	# the scaling is given explicitly, and applies to all given samples
	for ( i in 1:length(fileSet)) {
		thisFile <- fileSet[i]
		thisSID <- sampleIDset[i]
		newFile <- sub( old.path, adjusted.path, thisFile, fixed=T)

		cat( "\n", i, "ReadFile:", thisSID)
		tbl <- read.delim( thisFile, as.is=T)
		
		# if this sample is exempt from scaling, just write it across to new locations now and skip over it
		if ( ! is.null( exemptSampleIDset) && thisSID %in% exemptSampleIDset) {
			cat( "  writeExemptFile..")
			write.table( tbl, newFile, sep="\t", quote=F, row.names=F)
			next
		}
		
		myColumns <- match( ExpressColumns, colnames(tbl))
		if ( any( is.na( myColumns))) {
			cat( "\nSome expression columns not in transcriptome file: ", ExpressColumns[ is.na(myColumns)])
			myColumns <- myColumns[ ! is.na(myColumns)]
		}
	
		# we will clip the uncertainies at their unscaled lower bound
		minSigma <- min( tbl$SIGMA_M, na.rm=T)
		
		# find where every gene is 
		myGenes <- tbl$GENE_ID
		N <- nrow(tbl)
		myScaleFacs <- rep.int( 1, N)
		where <- match( myGenes, scaleGenes, nomatch=0)
		if ( sum( where > 0) < min(N,NSG)/2) {
			where <- match( shortGeneName(myGenes,keep=1), shortGeneName(scaleGenes,keep=1), nomatch=0)
			if ( sum( where > 0) < min(N,NSG)/2) {
				cat( "\nError:  Not enough Gene IDs match between scale factor file and transcriptome...")
				stop()
			}
		}
		myScaleFacs[ where > 0] <- scaleFactors[ where]
		
		# make the new adjusted transcriptome
		cat( "  applyScaling..")
		newtbl <- tbl
		for ( j in myColumns) {
			v <- tbl[[j]]
			vnew <- as.numeric(v)
			vnew <- vnew * myScaleFacs
			newtbl[[j]] <- vnew
		}
		newtbl$SIGMA_M <- ifelse( newtbl$SIGMA_M < minSigma, minSigma, newtbl$SIGMA_M)
		newtbl$SIGMA_U <- ifelse( newtbl$SIGMA_U < minSigma, minSigma, newtbl$SIGMA_U)

		# force the TPM units to sum to a million ???
		#if ( "TPM_M" %in% colnames(newtbl)) {
			#newtbl$TPM_M <- tpm( newtbl$READS_M, newtbl$N_EXON_BASES)
			#newtbl$TPM_U <- tpm( newtbl$READS_U, newtbl$N_EXON_BASES)
		#}
	
		# put the new transcripome back into expression order
		exOrd <- order( newtbl$RPKM_M, decreasing=T)
		newtbl <- newtbl[ exOrd, ]
	
		cat( "  writeNewFile..")
		write.table( newtbl, newFile, sep="\t", quote=F, row.names=F)
	}
	cat( "\nDone.\n")

	# pass back the set of samples we scaled
	out <- sampleIDset
	return( invisible(out))
}
