# pipe.DiffExpression.R

# make the DE using the transcript files.

`pipe.DiffExpression` <- function( sampleIDset, groupSet=NULL, annotationFile="Annotation.txt", 
		optionsFile="Options.txt", results.path=NULL, speciesID=NULL, 
		altGeneMap=NULL, altGeneMapLabel=NULL, 
		minRPKM=NULL, missingOnly=FALSE, verbose=TRUE) {

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

	# during 'altGeneMap' runs, skip over species that are not covered...
	doingAltGeneMap <- FALSE
	if ( ! is.null( altGeneMap)) {
		doingAltGeneMap <- TRUE
		altSpeciesSet <- unique.default( getSpeciesFromSeqID( altGeneMap$SEQ_ID))
	}

	# do each species all the way through
	for( speciesID in allSpecies) {
	
	if ( doingAltGeneMap && (!( speciesID %in% altSpeciesSet))) next

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
	if ( doingAltGeneMap) {
		pathIn <- file.path( resultsPath, "transcript", altGeneMapLabel)
		pathOut <- file.path( resultsPath, "ratios", altGeneMapLabel)
		if ( ! file.exists( pathOut)) dir.create( pathOut, recursive=TRUE, showWarnings=FALSE)
	}

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
		if ( doingAltGeneMap) {
			fileIn1 <- paste( sampI, speciesPrefix, altGeneMapLabel, "Transcript.txt", sep=".")
		}
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
			if ( doingAltGeneMap) {
				fileOutDE <- paste( thisPair, speciesPrefix, altGeneMapLabel, "Ratio.txt", sep=".")
			}
			fileOutDE <- file.path( pathOut, fileOutDE)

			if ( missingOnly) {
				if ( file.exists( fileOutDE)) {
					cat( "\nRatios file already made.  Skip: ", basename(fileOutDE))
					next
				}
			}

			# now get this sample's expression data
			fileIn2 <- paste( sampJ, speciesPrefix, "Transcript.txt", sep=".")
			if ( doingAltGeneMap) {
				fileIn2 <- paste( sampJ, speciesPrefix, altGeneMapLabel, "Transcript.txt", sep=".")
			}
			fileIn2 <- file.path( pathIn, fileIn2)
			if ( ! file.exists( fileIn2)) {
				cat( "\nSkipping DE for this sampleID/species.   file not found: ", fileIn2)
				break
			}

			thisDE <- calcDiffExpressFromTranscripts( fileIn1, fileIn2, minRPKM=thisMinRPKM, 
					wt.folds=wt.folds, wt.pvalues=wt.pvalues,
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

