# readTrimming - standardize the lookup of read base trimming requested for a sample


getReadTrimming <- function( sampleID="", annotationFile="Annotation.txt", optionsFile="Options.txt", verbose=T) {

	# trimming info is in the options file
	optT <- readOptionsTable( optionsFile)
	trim5 <- as.integer( getOptionValue( optT, "trim5", notfound="0", verbose=verbose))
	trim3 <- as.integer( getOptionValue( optT, "trim3", notfound="0", verbose=verbose))

	# there may also be trim arguments in the annotation file, and these would override
	if ( sampleID != "") {
		originalID <- originalSamplePairID(sampleID, annotationFile=annotationFile) 
		annT <- readAnnotationTable( annotationFile)
		if ( "Trim5" %in% colnames(annT)) {
			trim5.a <- getAnnotationValue( annT, key=originalID, columnArg="Trim5", notfound="0", verbose=F)
			trim5.a <- as.integer( trim5.a)
			if ( ! is.na( trim5.a)) {
				trim5 <- trim5.a
				if (verbose) cat( "\nOverride: Using 'Trim5' from Annotation file: ", trim5)
			}
		}
		if ( "Trim3" %in% colnames(annT)) {
			trim3.a <- getAnnotationValue( annT, key=originalID, columnArg="Trim3", notfound="0", verbose=F)
			trim3.a <- as.integer( trim3.a)
			if ( ! is.na( trim3.a)) {
				trim3 <- trim3.a
				if (verbose) cat( "\nOverride: Using 'Trim3' from Annotation file: ", trim3)
			}
		}
	}

	return( list( "trim5"=trim5, "trim3"=trim3))
}

