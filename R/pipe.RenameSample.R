# pipe.RenameSample.R -- change a 'SampleID' in ALL results files


`pipe.RenameSample` <- function( sampleID,  newSampleID, results.path=NULL,
				annotationFile="Annotation.txt", optionsFile="Options.txt",
				verbose=TRUE) {

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nRenaming files for sample:     ", sampleID)
		cat( "\nNew 'sampleID' will be:        ", newSampleID, "\n")
	}
	gc()

	# character substitution can fail if old is a direct subset of the new
	if ( regexpr( sampleID, newSampleID, fixed=T) > 0) {
		cat( "\nWarning:\n'sampleID' can not be an exact substring of 'newSampleID'")
		cat( "\nCauses recursive character substitution in filenames.")
		stop()
	}
	
	# make folders for all results...
	if ( is.null( results.path)) {
		resultsPath <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	}
	cat( "\n")
	subfolders <- c( "align", "fastq", "html", "ratios", "riboClear", "splicing", "summary", 
			"transcript", "wig", "AlignStats", "USR", "VariantCalls", 
			"VelvetContigs", "VelvetPeptides", "CR", "HLA.typing", "ConsensusProteins",
			"SieveAnalysis", "Viral.BAMS", "LineageCalls")
	dataType <- getAnnotationValue( annotationFile, key=sampleID, columnArg="DataType", notfound="RNA-seq")
	if ( dataType == "ChIP-seq") subfolders <- c( subfolders, "ChIPpeaks")

	ndone <- 0
	for ( folder in subfolders) {
		thisFolder <- file.path( resultsPath, folder)

		# find all files that start with this ID
		files <- dir( thisFolder, pattern=paste( "^", sampleID, sep=""), full.name=TRUE)
		for ( f in files) {
			if ( file.info(f)$isdir) next
			newfile <- gsub( sampleID, newSampleID, f, fixed=TRUE)
			if ( ! file.exists( newfile)) {
				file.rename( f, newfile)
				cat( "\nFolder: ", basename(thisFolder), "\tFile: ", basename(f))
				ndone <- ndone + 1
				renameFileContents( newfile, sampleID, newSampleID)
			}
		}

		# repeat this pass, watching for directories of subfolders for this sample
		files <- dir( thisFolder, pattern=paste( "^", sampleID, sep=""), full.name=TRUE)
		for ( f in files) {
		    if ( ! file.info(f)$isdir) next
		    newfolder <- gsub( sampleID, newSampleID, f)
		    file.rename( f, newfolder)
		    cat( "\nRename Folder: ", f)
		    ndone <- ndone + 1
		    subfiles <- dir( newfolder, pattern=paste( "^", sampleID, sep=""), full.name=TRUE)
		    for ( subf in subfiles) {
			newfile <- gsub( sampleID, newSampleID, subf, fixed=TRUE)
			if ( ! file.exists( newfile)) {
				file.rename( subf, newfile)
				cat( "\nFolder: ", basename(thisFolder), "\tFile: ", basename(subf))
				ndone <- ndone + 1
				renameFileContents( newfile, sampleID, newSampleID)
			}
		    }
		}

		# some folders can have SampleIDs not at the front of the name
		if ( basename( thisFolder) %in% c( "ratios")) {
			files <- dir( thisFolder, pattern=sampleID, full.name=TRUE)
			for ( f in files) {
				if ( file.info(f)$isdir) next
				newfile <- gsub( sampleID, newSampleID, f, fixed=TRUE)
				if ( ! file.exists( newfile)) {
					file.rename( f, newfile)
					cat( "\nFolder: ", basename(thisFolder), "\tFile: ", basename(f))
					ndone <- ndone + 1
				}
			}
		}
	}

	# all subfolders done...

	# top level files ??
	topfiles <- c( paste( sampleID, "log.txt", sep="."), 
			paste( ".FileLineCounts", sampleID, "txt", sep="."))
	for ( f in topfiles) {
		if ( ! file.exists(f)) next
		if ( file.info(f)$isdir) next
		newfile <- gsub( sampleID, newSampleID, f, fixed=TRUE)
		if ( ! file.exists( newfile)) {
			file.rename( f, newfile)
			cat( "\n\tFile: ", basename(f))
			ndone <- ndone + 1
			renameFileContents( newfile, sampleID, newSampleID)
		}
	}

	cat( "\nDone.\nN_Files Renamed: ", ndone,"\n")

	return()
}


renameFileContents <- function( newfile, sampleID, newSampleID) {

	# only a few file types can should get their content adjusted
	doText <- FALSE
	if (regexpr( "\\.html$", newfile) > 0) doText <- TRUE
	if (regexpr( "Summary.txt$", newfile) > 0) doText <- TRUE
	if (regexpr( "FileLineCounts.+txt$", newfile) > 0) doText <- TRUE
	if (regexpr( "log.txt$", newfile) > 0) doText <- TRUE

	if ( doText) {
		cat( "\nRenaming contents in file: ", basename(newfile))
		txt <- readLines( newfile)
		txtout <- gsub( sampleID, newSampleID, txt, fixed=T)
		writeLines( txtout, con=newfile)
		return()
	}

	doWIG <- FALSE
	if (regexpr( "\\.WIG.rda$", newfile) > 0) doWIG <- TRUE

	if (doWIG) {
		who <- load( newfile)
		if( who != "wiggles") return()
		myInfo <- wiggles$Info
		if (is.null( myInfo)) return()
		mySubFolder <- wiggles$SubWigFolder
		if (is.null( mySubFolder)) return()
		cat( "\nRenaming contents in WIG file: ", basename(newfile))

		myInfo$FileName <- gsub( sampleID, newSampleID, myInfo$FileName)
		myInfo$SampleID <- newSampleID
		mySubFolder <- gsub( sampleID, newSampleID, mySubFolder)
		wiggles$Info <- myInfo
		wiggles$SubWigFolder <- mySubFolder

		save( wiggles, file=newfile)
		return()
	}

}
