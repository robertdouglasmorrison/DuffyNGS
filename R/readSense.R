# readSense.R - standardize the string for the sense of RNA-seq reads


getReadSense <- function( sampleID, annotationFile="Annotation.txt") {

	# get the string from the annotation file
	originalID <- originalSamplePairID(sampleID, annotationFile) 
	readSense <- getAnnotationValue( annotationFile, key=originalID, 
					columnArg="ReadSense", notfound="antisense", verbose=F)

	# second part of paired reads should land on the opposite strand
	paired <- getAnnotationTrue( annotationFile, key=originalID, 
					columnArg="PairedEnd", verbose=F)
	if ( paired && (regexpr( "_2$", sampleID) > 0)) {
		readSense <- as.readAntisense( readSense)
	} else {
		readSense <- as.readSense( readSense)
	}
	return( readSense)
}


as.readSense <- function( sense) {

	choices <- c( "sense", "antisense")

	out <- base::tolower( base::gsub( "-| ","",sense))
	if ( ! (out %in% choices)) {
		warning( paste( "Invalid 'ReadSense':  ", sense, "\tChoices:   sense, antisense "))
		out <- "sense"
	}
	return( out)
}


as.readAntisense <- function( sense) {

	choices <- c( "sense", "antisense")
	given <- as.readSense( sense)
	return( choices[ 3 - match( given, choices)])
}

