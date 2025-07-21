# genericDiffExpression.R - alternative method for making "Ratio" files from non-RNA-seq type data...

# make the DE using the transcript files.

`pipe.GenericDiffExpression` <- function( sampleIDset, groupSet=NULL, annotationFile="Annotation.txt", 
		optionsFile="Options.txt", results.path=NULL, speciesID=NULL, 
		intensityColumn="RPKM_M", sep="\t", minRPKM=NULL, missingOnly=FALSE, verbose=TRUE) {

	# pre-catch the case of MissingOnly=TRUE and exactly 2 samples that are already done
	# so we can silently return in a hurry
	# when called in multicore mode, we get a lot of console status for nothing..
	if ( missingOnly && length(sampleIDset) == 2) {
		if ( is.null( results.path)) {
			results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
		}
		if ( ! is.null(speciesID) && getCurrentSpecies() != speciesID) setCurrentSpecies( speciesID)
		speciesPrefix <- getCurrentSpeciesFilePrefix()
		sampI <- sampleIDset[1]
		sampJ <- sampleIDset[2]
		thisPair1 <- paste( sampI, "v", sampJ, sep=".")
		thisPair2 <- paste( sampJ, "v", sampI, sep=".")
		file1 <- paste( thisPair1, speciesPrefix, "Ratio.txt", sep=".")
		file2 <- paste( thisPair2, speciesPrefix, "Ratio.txt", sep=".")
		bothFiles <- file.path( results.path, "ratios", c(file1,file2))
		if ( all( file.exists( bothFiles))) {
			cat( "\nRatio files already made.  Skip: ", sampI, "vs", sampJ)
			return()
		}
	}

	
	# the samples could come as a list of vectors of samples...
	if ( is.list( sampleIDset)) {
		listOfSets <- TRUE
		flatIDset <- unlist( sampleIDset)
		if ( ! is.null( groupSet)) {
			cat( "\n'groupSet' only valid with a vector of samples, not a list of vectors")
			groupSet <- NULL
		}
	} else {
		listOfSets <- FALSE
		flatIDset <- sampleIDset
	}

	if (verbose) {
		cat( verboseOutputDivider)
		if ( is.list( sampleIDset)) {
			cat( "\nStarting pipe 'DiffExpression' on Sample Set:\n")
			print( sampleIDset)
		} else {
			cat( "\nStarting pipe 'DiffExpression' on Sample Set:     ", sampleIDset, "\n")
		}
	}

	startTotalT <- proc.time()
	grandTotalReads <- 0

	Nsamples <- length( flatIDset)
	if ( Nsamples < 2) {
		cat( "\nDifferential Expression requires at least 2 sampleIDs.  Skipping...")
		return("Failure")
	}
	if ( ! is.null( groupSet)) {
		Ngroups <- length( groupSet)
		if ( Ngroups != Nsamples) {
			cat( "\nGroup names must be same length as SampleIDs: ", Ngroups, Nsamples)
			cat( "\nGroups: ", groupSet)
			return("Failure")
		}
	}

	optT <- readOptionsTable( optionsFile)
	target <- getAndSetTarget( optionsFile)

	if ( is.null(speciesID)) {
		allSpecies <- getCurrentTargetSpecies()
	} else {
		allSpecies <- speciesID
	}

	if ( is.null( results.path)) {
		resultsPath <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	} else {
		resultsPath <- results.path
	}

	if ( is.null( minRPKM)) {
		minRPKM <- as.numeric( getOptionValue( optT, "DE.minimumRPKM", notfound="1"))
	}
	wt.folds <- as.numeric( getOptionValue( optT, "DE.weightFoldChange", notfound="1"))
	wt.pvalues <- as.numeric( getOptionValue( optT, "DE.weightPvalue", notfound="1"))

	
	findSampleInSets <- function( x) {
		for ( k in 1:length(sampleIDset)) {
			if ( any( sampleIDset[[k]] == x)) return(k)
		}
		return(0)
	}

	# do each species all the way through
	for( speciesID in allSpecies) {
	
	cat( "\n\nCalculating DiffExpress for Species:  ", speciesID,"\n")
	startT <- proc.time()
	gc()
	setCurrentSpecies( speciesID)
	speciesPrefix <- getCurrentSpeciesFilePrefix()

	# allow for a species specific minimum
	thisMinRPKM <- as.numeric( getOptionValue( optT, "DE.minimumRPKM", speciesID=speciesID, 
			notfound=as.character( minRPKM), verbose=T))

	gmap <- getCurrentGeneMap()
	pathIn <- file.path( resultsPath, "transcript")
	pathOut <- file.path( resultsPath, "ratios")
	if ( ! file.exists( pathOut)) dir.create( pathOut, recursive=TRUE, showWarnings=FALSE)

	# do all pairs of samples
	for ( i in 1:Nsamples) {

		sampI <- flatIDset[i]
		if ( listOfSets) {
			grpI <- findSampleInSets( sampI)
		} else {
			grpI <- 1
		}
		if ( ! is.null(groupSet)) grpNameI <- groupSet[i]

		# get the transcript for this sample
		fileIn1 <- paste( sampI, speciesPrefix, "Transcript.txt", sep=".")
		fileIn1 <- file.path( pathIn, fileIn1)
		if ( ! file.exists( fileIn1)) {
			cat( "\nSkipping DE for this sampleID/species.   file not found: ", fileIn1)
			break
		}

		for ( j in 1:Nsamples) {
			if ( i == j) next

			# only do if in the same subset
			sampJ <- flatIDset[j]
			if ( listOfSets) {
				grpJ <- findSampleInSets( sampJ)
			} else {
				grpJ <- 1
			}
			if ( grpI != grpJ) next

			# if given group membership per sample, don't do DE for 2 samples in the same group
			if ( ! is.null(groupSet)) {
				grpNameJ <- groupSet[j]
				if ( grpNameI == grpNameJ) next
			}

			# build the name for the DE file
			thisPair <- paste( sampI, "v", sampJ, sep=".")
			fileOutDE <- paste( thisPair, speciesPrefix, "Ratio.txt", sep=".")
			fileOutDE <- file.path( pathOut, fileOutDE)

			if ( missingOnly) {
				if ( file.exists( fileOutDE)) {
					cat( "\nRatios file already made.  Skip: ", basename(fileOutDE))
					next
				}
			}

			# now get this sample's expression data
			fileIn2 <- paste( sampJ, speciesPrefix, "Transcript.txt", sep=".")
			fileIn2 <- file.path( pathIn, fileIn2)
			if ( ! file.exists( fileIn2)) {
				cat( "\nSkipping DE for this sampleID/species.   file not found: ", fileIn2)
				break
			}

			thisDE <- calcGenericDiffExpressFromTranscripts( fileIn1, fileIn2, minRPKM=thisMinRPKM, 
					wt.folds=wt.folds, wt.pvalues=wt.pvalues,
					intensityColumn=intensityColumn, sep=sep, 
					fileout=fileOutDE, verbose=verbose)
		}
	}

	if (verbose) {
		cat( "\n\nFinished 'DiffExpression' for Species: ", speciesID)
	}
	gc()

	}  # end of allSpecies...


	if (verbose) {
		cat( verboseOutputDivider)
		if (listOfSets) {
			cat( "\n\nFinished pipe 'DiffExpression' on Sample Set:\n")
			print( sampleIDset)
		} else {
			cat( "\n\nFinished pipe 'DiffExpression' on Sample Set:     ", flatIDset, "\n")
		}
	}

	return( paste(flatIDset, collapse="|"))
}


`calcGenericDiffExpressFromTranscripts` <- function( file1, file2, minRPKM=2, clipFold=10.0, 
				wt.folds=1.0, wt.pvalues=1.0, intensityColumn="RPKM_M", sep="\t",
				fileout=NA, verbose=TRUE) {

	if (verbose) {
		cat( "\n\nCalculating DE for files:\n\t", file1, "\n\t", file2, "\n")
	}

	# turn the 2 transcripts into matching order...
	trans1 <- read.delim( file1, as.is=T, sep=sep)
	trans2 <- read.delim( file2, as.is=T, sep=sep)

	allG <- sort( intersect( trans1$GENE_ID, trans2$GENE_ID))
	wh1 <- match( allG, trans1$GENE_ID)
	wh2 <- match( allG, trans2$GENE_ID)
	trans1 <- trans1[ wh1, ]
	trans2 <- trans2[ wh2, ]
	NG <- length( allG)

	# grab those two sets of expression values
	val1 <- trans1[[ intensityColumn]]
	val2 <- trans2[[ intensityColumn]]

	# allocate strorage for returned data
	geneSet <- gProdSet <- nBaseSet <- vector( length=NG)
	PvalueSet <- vector( length=NG)
	FoldSet <- rpkmSetA <- rpkmSetB <- vector( length=NG)

	# visit every gene
	ngenes <- nDE <- 0

	for( ig in 1:NG) {

		ans <- calcGenericDE( val1[ig], val2[ig], minRPKM=minRPKM, clipFold=clipFold)

		# load the vectors and do the DE calc
		ngenes <- ngenes + 1
		geneSet[ ngenes] <- gene <- allG[ig]
		gProdSet[ ngenes] <- trans1$PRODUCT[ig]
		PvalueSet[ ngenes] <- ans$Pvalue
		FoldSet[ ngenes] <- ans$Fold
		rpkmSetA[ ngenes] <- ans$RPKM1
		rpkmSetB[ ngenes] <- ans$RPKM2
		if ( ig %% 1000 == 0) cat( "\r", ig, "\t", gene, "\t", ans$rpkmFold)
	}

	PvalueSet[ PvalueSet > 1] <- 1
	PIvalueSet <- piValue( FoldSet, PvalueSet)

	out <- data.frame( geneSet, gProdSet, PvalueSet, FoldSet, PIvalueSet, rpkmSetA, rpkmSetB, 
				stringsAsFactors=FALSE)
	colnames( out) <- c("GENE_ID", "PRODUCT", "PVALUE_M", "LOG2FOLD_M", "PIVALUE_M", "RPKM_1_M", "RPKM_2_M")

	# now sort based on Pvalue and fold..
	ord <- diffExpressRankOrder( out$LOG2FOLD_M, out$PVALUE_M, wt.folds, wt.pvalues)
	out <- out[ ord, ]
	rownames( out) <- 1:nrow(out)
	
	if ( ! is.na( fileout)) {
		if ( getCurrentSpecies() %in% MAMMAL_SPECIES) out <- addHumanIDterms( out)
		#if ( getCurrentSpecies() %in% ORIGID_PARASITE_SPECIES) out <- addOrigIDterms( out)
		if ( getCurrentSpecies() %in% BACTERIA_SPECIES) out <- addNameIDterms( out)

		write.table( out, file=fileout, sep="\t", quote=FALSE, row.names=FALSE)
		cat( "\nWrote Differential Expression file:  \t", fileout)
	}

	if (verbose) {
		cat( "\nN_Gene regions processed:         \t", nrow(out))
	}

	return( out)
}


`calcGenericDE` <- function( v1, v2, minRPKM=1, clipFold=10.0) {

	# evaluate the differential expression between two samples for one gene...
	# extract the values from the 2 transcript subsets...
	rawA <- v1
	rawB <- v2

	# we don;t have any sense of sigma, and assume the values given are like RPKM, alreadly normalized

	# make a Pvalue... create a distribution around each one point
	otherA <- rnorm( 5, mean=rawA, sd=if(rawA > 1) sqrt(rawA) else 1)
	otherB <- rnorm( 5, mean=rawB, sd=if(rawB > 1) sqrt(rawB) else 1)
	Pvalue <- t.test( c(rawA,otherA), c(rawB,otherB))$p.value

	# RPKM stuff
	Fold <- log2( (rawA+minRPKM) / (rawB+minRPKM))

	out <- list( "Pvalue"=Pvalue, "Fold"=Fold, "RPKM1"=rawA, "RPKM2"=rawB)
	return( out)
}

