# buildAllAnnotationFiles.R

# try to build every file after an annotation map set gets make...

`buildAllAnnotationFiles` <- function( speciesIDSet, genomicFastaFileSet, fastaPatternSet="\\.fasta$",
			optionsFileSet="Options.txt", outPath=".", 
			steps=1:4, riboTailSize=36, exonFragmentSize=36, 
			maxExonSkipSet=2, intronSplicesTooSet=FALSE, 
			detectableReadSize=32, chunkSize=1000000, 
			otherSpeciesIDSet="Hs_grc", otherGenomicIndexSet="Hs.genomic_idx",
			verbose=TRUE, debug=FALSE) {


	# let's try doing combo maps in place...
	bigComboFastaSet <- riboComboFastaSet <- spliceComboFastaSet <- vector()
	riboComboMap <- data.frame()
	comboGenes <- 0

	# several args may now be vectors...
	N_SpeciesID <- length(speciesIDSet)
	doingComboIndex <- ( N_SpeciesID > 1)

	# make sure all Sets are the right size
	if ( length( fastaPatternSet) < N_SpeciesID) fastaPatternSet <- rep( fastaPatternSet, length.out=N_SpeciesID)
	if ( length( maxExonSkipSet) < N_SpeciesID) maxExonSkipSet <- rep( maxExonSkipSet, length.out=N_SpeciesID)
	if ( length( intronSplicesTooSet) < N_SpeciesID) 
				intronSplicesTooSet <- rep( intronSplicesTooSet, length.out=N_SpeciesID)

	if ( doingComboIndex && length( optionsFileSet) != (N_SpeciesID + 1)) {
		stop( "For buildind 'Combo Indexes', argument 'optionsFileSet' must be exactly 1 longer than 'speciesIDset'")
	}

	
	# most steps are done for each species first
	for( iSpecies in 1:N_SpeciesID) {


	speciesID <- speciesIDSet[ iSpecies]
	cat( "\n\nBuilding All Index files for:  ", speciesID, "\n")
	genomicFastaFile <- genomicFastaFileSet[ iSpecies]
	optionsFile <- optionsFileSet[ iSpecies]
	otherSpeciesID <- otherSpeciesIDSet[ iSpecies]
	if ( is.list( otherSpeciesID)) otherSpeciesID <- otherSpeciesIDSet[[ iSpecies]]
	otherGenomicIndex <- otherGenomicIndexSet[ iSpecies]
	if ( is.list(otherGenomicIndexSet)) otherGenomicIndex <- otherGenomicIndexSet[[ iSpecies]]
	maxExonSkip <- maxExonSkipSet[ iSpecies]
	intronSplicesToo <- intronSplicesTooSet[ iSpecies]

	
	# load that species MapSet data...
	setCurrentSpecies( speciesID=speciesID)
	prefix <- outfilePrefix <- getCurrentSpeciesFilePrefix()

	# prepare to use Bowtie:  most all steps use Bowtie2,  detectability needs Bowtie1
	if ( 4 %in% steps) bowtiePar.defaults( optionsFile)
	bowtie2Par.defaults( optionsFile)

	cat( "\n\nGetting names for new Indexes and Maps..")
	optT <- readOptionsTable( optionsFile)
	riboIndexFile <- getOptionValue( optT, "RiboIndex", notfound=paste( prefix, "ribo_idx", sep="."))
	riboMapFile <- getOptionValue( optT, "RiboMap", notfound=paste( prefix, "riboMap.txt", sep="."))
	genomicIndexFile <- getOptionValue( optT, "GenomicIndex", notfound=paste( prefix, "genomic_idx", sep="."))
	spliceIndexFile <- getOptionValue( optT, "SpliceIndex", notfound=paste( prefix, "splice_idx", sep="."))
	spliceMapPrefix <- getOptionValue( optT, "SpliceMapPrefix", notfound="spliceMap")
	bowtie1IndexPath <- getOptionValue( optT, "bowtieIndex.path", notfound=outPath)
	bowtie2IndexPath <- getOptionValue( optT, "bowtie2Index.path", notfound=outPath)

	startTime <- proc.time()

	if ( 1 %in% steps) {

		cat( "\n\nPart 1:  Genomic Bowtie Index...")

		outfile <- file.path( bowtie2IndexPath, genomicIndexFile)
		cat( "\n\nIndex files being written to:    ", outfile)

		# the 'fasta' file can be a directory...
		myFastaFile <- genomicFastaFile
		cat( "\nUsing Fasta file name/directory: ", myFastaFile)
		fInfo <- file.info( myFastaFile)
		if ( all( is.na( fInfo))) stop( "Fasta File/Folder not found...")
		if ( fInfo$isdir) {
			cat( "\nFinding fasta files in folder:    ", myFastaFile)
			pattern <- fastaPatternSet[ iSpecies]
			cat( "\nPattern string:                ", pattern)
			allFastas <- dir( path=myFastaFile, pattern=pattern, full.names=TRUE)
			cat( "\nN_files found: ", length( allFastas))
			bigComboFastaSet <- base::append( bigComboFastaSet, allFastas)
			# bowtie wants comma sep list...
			myFastaFile <- paste( allFastas, collapse=",")
		} else {
			bigComboFastaSet <- base::append( bigComboFastaSet, genomicFastaFile)
		}

		if ( ! doingComboIndex) {
			cmdline <- buildBowtie2BuildCommandLine( inputFastaFile=myFastaFile, 
					outputIndexFile=outfile, optionsFile=optionsFile, 
					verbose=FALSE, debug=debug)
			callBowtie2Build( cmdline, wait=TRUE, verbose=verbose)
		} else {
			cat( "\nSkipping 'non-Combo' genomic index step for:  ", speciesID)
		}
	} else {
		cat( "\n\nPart 1:  Skipping build of Genomic Bowtie Index...")
	}


	nRRNA <- nrow( getCurrentRrnaMap())
	if ( 2 %in% steps) {
	    if ( nRRNA > 0 && sum( getCurrentRrnaMap()$CLEAR) > 0) {

		cat( "\n\nPart 2a:  Ribo Clearing files...")
		ans <- buildRiboClearingFiles( speciesID=speciesID, genomicFastaFile=genomicFastaFile, 
			outPath=outPath, mapFile=riboMapFile, tailSize=riboTailSize, verbose=T)
		riboFastaFile <- ans$FastaFile
		if ( ! file.exists( riboFastaFile)) stop( "RiboClear Fasta file not created... ")
		riboComboFastaSet <- base::append( riboComboFastaSet, riboFastaFile)
		riboComboMap <- rbind( riboComboMap, read.delim( ans$MapFile, as.is=TRUE))

		cat( "\n\nPart 2b:  Ribo Clearing Bowtie Index...")
		outfile <- file.path( bowtie2IndexPath, riboIndexFile)
		cmdline <- buildBowtie2BuildCommandLine( inputFastaFile=riboFastaFile, outputIndexFile=outfile, 
				optionsFile=optionsFile, verbose=FALSE, debug=debug)
		callBowtie2Build( cmdline, wait=TRUE, verbose=verbose)

		# move a copy of the map file...
		cmd <- paste( "cp ", ans$MapFile, file.path( bowtie2IndexPath, riboMapFile))
		cat( "\nMoving copy of Ribo Map...\n", cmd, "\n")
		system( cmd)

	    } else {
		cat( "\n\nNo Genes flagged for clearing in the rRNA Map...  No RiboClearing files made.\n")
	    }
	} else {
		cat( "\n\nPart 2:  Skipping build of Ribo Clearing files...")
	}


	emap <- getCurrentExonMap()
	nExons <- nrow( emap)
	nGenes <- length( unique.default( emap$GENE_ID))
	comboGenes <- comboGenes + nGenes

	if ( 3 %in% steps) {
	    if ( nExons > nGenes ) {

		cat( "\n\nPart 3a:  Splice Junction files...")
		spliceFastaFile <- paste( prefix, spliceMapPrefix, "fasta", sep=".")

		#if ( !doingComboIndex || !file.exists( spliceFastaFile)) {
			ans <- buildSpliceJunctionFiles( speciesID=speciesID, genomicFastaFile=genomicFastaFile, 
				outPath=outPath, spliceMapPrefix=spliceMapPrefix, 
				fragmentSize=exonFragmentSize, maxExonSkip=maxExonSkip, 
				intronSplicesToo=intronSplicesToo, verbose=T)
			spliceFastaFile <- ans$FastaFile
			spliceMapFile <- ans$MapFile

			cat( "\n\nPart 3b:  Compressing Splice Maps...")
			ans2 <- compressSpliceMaps( outPath=bowtie2IndexPath, spliceMaps=spliceMapFile,
					spliceMapPrefix=spliceMapPrefix)
		#} else {
		#	cat( "\nSkipping 'non-Combo' splice map step for:  ", speciesID)
		#	spliceFastaFile <- paste( prefix, spliceMapPrefix, "fasta", sep=".")
		#	spliceMapFile <- paste( prefix, spliceMapPrefix, "txt", sep=".")
		#}

		# make sure we have files we need
		if ( ! file.exists( spliceFastaFile)) stop( "Splice Junctions Fasta file not created... ")
		if ( ! file.exists( spliceMapFile)) stop( "Splice Junctions Map file not created... ")
		spliceComboFastaSet <- base::append( spliceComboFastaSet, spliceFastaFile)


		if ( ! doingComboIndex) {
			cat( "\n\nPart 3c:  Splice Junction Bowtie Index...")
			outfile <- file.path( bowtie2IndexPath, spliceIndexFile)
			cmdline <- buildBowtie2BuildCommandLine( inputFastaFile=spliceFastaFile, outputIndexFile=outfile, 
					optionsFile=optionsFile, verbose=verbose, debug=debug)
			callBowtie2Build( cmdline, wait=TRUE, verbose=verbose)
		} else {
			cat( "\nSkipping 'non-Combo' splice index step for:  ", speciesID)
		}

	    } else {
		cat( "\n\nNo Exon Junctions in the Exon Map...  No Splice Junction files made.\n")
	    }
	} else {
		cat( "\n\nPart 3:  Skipping build of Splice Junction files...")
	}


	if ( 4 %in% steps) {
	    if ( ! doingComboIndex) {
		
		# verify the index we need is there
		idxFiles <- dir( bowtie1IndexPath, pattern=paste( "^",genomicIndexFile,sep=""))
		if ( ! length( idxFiles)) {
			cat( "\n\nNo Bowtie1 index found for:  ", genomicIndexFile)
			cat( "\nBuilding now..")
			# the 'fasta' file can be a directory...
			myFastaFile <- genomicFastaFile
			cat( "\nUsing Fasta file name/directory: ", myFastaFile)
			outfile <- file.path( bowtie1IndexPath, genomicIndexFile)
			cmdline <- buildBowtieBuildCommandLine( inputFastaFile=myFastaFile, 
					outputIndexFile=outfile, optionsFile=optionsFile, 
					verbose=FALSE, debug=debug)
			callBowtieBuild( cmdline, wait=TRUE, verbose=verbose)
		} else {
			cat( "\nBowtie1 index found and ready..")
		}
		# also verify the index we need for all the 'other' genomes are there
		nOtherIndex <- length( otherSpeciesIDSet)
		if (nOtherIndex) for ( k in 1:nOtherIndex) {
			otherIndex <- otherGenomicIndexSet[k]
			otherSpecies <- otherSpeciesIDSet[k]
			if ( is.list( otherSpeciesIDSet)) otherSpecies <- otherSpeciesIDSet[[k]]
			if ( is.list(otherGenomicIndexSet)) otherIndex <- otherGenomicIndexSet[[k]]
			idxFiles <- dir( bowtie1IndexPath, pattern=paste( "^",otherIndex,sep=""))
			if ( ! length( idxFiles)) {
				cat( "\n\nNo Bowtie1 index found for:  ", otherIndex)
				cat( "\nTrying to Build now..")
				# the 'fasta' file can be a directory...
				myOtherFastaPath <- dirname(genomicFastaFile)
				myOtherFastaFile <- paste( otherSpecies, "genomicDNA.fasta", sep="_")
				myOtherFastaFile <- file.path( myOtherFastaPath, myOtherFastaFile)
				cat( "\nUsing Fasta file name/directory: ", myOtherFastaFile)
				outfile <- file.path( bowtie1IndexPath, otherIndex)
				cmdline <- buildBowtieBuildCommandLine( inputFastaFile=myOtherFastaFile, 
						outputIndexFile=outfile, optionsFile=optionsFile, 
						verbose=FALSE, debug=debug)
				callBowtieBuild( cmdline, wait=TRUE, verbose=verbose)
			} else {
				cat( "\nOther Bowtie1 index found and ready..", otherSpecies)
			}
		}

		cat( "\n\nPart 4a:  Detectability files...")
		ans <- buildDetectabilityFiles( speciesID=speciesID, genomicFastaFile=genomicFastaFile, 
				selfGenomicIndex=genomicIndexFile, otherGenomicIndex=otherGenomicIndex,
				chunkSize=chunkSize, verbose=verbose)

		cat( "\n\nPart 4b:  Convert Detectability to WiggleBin files...")
		selfFile <- ans$SelfUniqueFile
		otherFile <- ans$OtherDetectableFile
		wb <- loadDetectabilityToWB( filein=selfFile, WB=NULL, speciesID=speciesID, 
				type="selfUnique", otherSpeciesID="", readSize=detectableReadSize)
		if ( ! is.null( otherFile)) {
			# this guy can now be a vector of 'others'
			for( k in 1:length(otherFile)) {
				wb <- loadDetectabilityToWB( filein=otherFile[k], WB=NULL, speciesID=speciesID, 
						type="otherDetectable", otherSpeciesID=otherSpeciesID[k], 
						readSize=detectableReadSize)
			}
		}
		file.delete( selfFile)
		file.delete( otherFile)
	    } else {
		cat( "\nSkipping detectability step for 'Combo' indexes.")
	    }
	}


	cat( "\n\nAll 'Annotation Files' for speciesID: ", speciesID, " successfully made.\n")

	stopTime <- proc.time()
	print( elapsedProcTime( startTime, stopTime, N=nGenes, what="Gene"))

	} # end of combo species 'by each' loop


	# do the combo pass ??
	if ( ! doingComboIndex) return()


	optionsFile <- optionsFileSet[ N_SpeciesID + 1]
	prefix <- outfilePrefix <- getOptionValue( optionsFile, "targetID", notfound="HsPf")

	# prepare to use Bowtie
	bowtie2Par.defaults( optionsFile)
	optT <- readOptionsTable( optionsFile)
	riboIndexFile <- getOptionValue( optT, "RiboIndex", notfound=paste( prefix, "ribo_idx", sep="."))
	riboMapFile <- getOptionValue( optT, "RiboMap", notfound=paste( prefix, "riboMap.txt", sep="."))
	genomicIndexFile <- getOptionValue( optT, "GenomicIndex", notfound=paste( prefix, "genomic_idx", sep="."))
	spliceIndexFile <- getOptionValue( optT, "SpliceIndex", notfound=paste( prefix, "splice_idx", sep="."))
	spliceMapPrefix <- getOptionValue( optT, "SpliceMapPrefix", notfound="spliceMap")
	bowtie2IndexPath <- getOptionValue( optT, "bowtie2Index.path", notfound=outPath)

	startTime <- proc.time()
	cat("\n")

	if ( 1 %in% steps) {

		cat( "\n\nPart 1:  'Combo' Genomic Bowtie Index...")

		outfile <- file.path( bowtie2IndexPath, genomicIndexFile)
		cat( "\n\nIndex files being written to: ", outfile)

		# the 'fasta' files are now a vector
		myFastaFile <- paste( bigComboFastaSet, collapse=",")
		cmdline <- buildBowtie2BuildCommandLine( inputFastaFile=myFastaFile, 
				outputIndexFile=outfile, optionsFile=optionsFile, 
				verbose=verbose, debug=debug)
		callBowtie2Build( cmdline, wait=TRUE, verbose=verbose)
	} else {
		cat( "\n\nPart 1:  Skipping build of 'Combo' Genomic Bowtie Index...")
	}

	nRRNA <- nrow( riboComboMap)
	if ( 2 %in% steps) {
	    if ( nRRNA > 0) {

		cat( "\n\nPart 2a:  'Combo' Ribo Clearing files...")
		write.table( riboComboMap, file=file.path( ".", riboMapFile), sep="\t", quote=F)

		cat( "\n\nPart 2b:  'Combo' Ribo Clearing Bowtie Index...")
		outfile <- file.path( bowtie2IndexPath, riboIndexFile)
		riboFastaFile <- paste( riboComboFastaSet, collapse=",")
		cmdline <- buildBowtie2BuildCommandLine( inputFastaFile=riboFastaFile, outputIndexFile=outfile, 
				optionsFile=optionsFile, verbose=verbose, debug=debug)
		callBowtie2Build( cmdline, wait=TRUE, verbose=verbose)

		# move a copy of the map file...
		write.table( riboComboMap, file=file.path( bowtie2IndexPath, riboMapFile), sep="\t", quote=F)

	    } else {
		cat( "\n\nNo Genes in the rRNA Map...  No 'Combo' RiboClearing files made.\n")
	    }
	} else {
		cat( "\n\nPart 2:  Skipping build of 'Combo' Ribo Clearing files...")
	}


	if ( 3 %in% steps) {

		cat( "\n\nPart 3c:  'Combo' Splice Junction Bowtie Index...")
		if ( length( spliceComboFastaSet) > 1) {
			outfile <- file.path( bowtie2IndexPath, spliceIndexFile)
			spliceFastaFile <- paste( spliceComboFastaSet, collapse=",")
			cmdline <- buildBowtie2BuildCommandLine( inputFastaFile=spliceFastaFile, outputIndexFile=outfile, 
					optionsFile=optionsFile, verbose=verbose, debug=debug)
			callBowtie2Build( cmdline, wait=TRUE, verbose=verbose)
		} else {
			cat( "\nNot enough splice files for a 'Combo' index...  Skipping.")
		}

	} else {
		cat( "\n\nPart 3:  Skipping build of 'Combo' Splice Junction files...")
	}


	cat( "\n\nAll 'Annotation Files' for combo species:  (", speciesIDSet, ")  successfully made.\n\n")

	stopTime <- proc.time()
	print( elapsedProcTime( startTime, stopTime, N=comboGenes, what="Gene"))

	return()
}
