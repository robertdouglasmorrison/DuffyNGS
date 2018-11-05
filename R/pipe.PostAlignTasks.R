# pipe.PostAlignTasks.R

# do any and all post alignment tasks for a sample

`pipe.PostAlignTasks` <- function( sampleID,  annotationFile="Annotation.txt", optionsFile="Options.txt",
			dataType=c("RNA-seq", "DNA-seq", "ChIP-seq", "RIP-seq"), verbose=TRUE) {

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nStarting 'PostAlignTasks' on Sample:     ", sampleID, "\n\n")
	}
	gc()
	dataType <- match.arg( dataType)


	# convert/append any BAM file details...

	# currently nothing to do at this point

	# assess No Hits
	doNoHits <- getOptionTrue( optionsFile, "investigate.NoHits", notfound=TRUE)
	if ( doNoHits) {
		pipe.NoHits_Investigate( sampleID, annotationFile=annotationFile,
				optionsFile=optionsFile)
		pipe.CR_Investigate( sampleID, annotationFile=annotationFile,
				optionsFile=optionsFile)
	}
	
	cat( "\n\nFinished 'PostAlignTasks' on Sample:     ", sampleID, "\n\n")
	gc()
	return()
}

