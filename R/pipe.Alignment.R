# pipe.Alignment.R

# do the A to Z alignment pipeline on a sample

`pipe.Alignment` <- function( sampleID=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt", ...) {

	datatype <- getAnnotationValue( annotationFile, key=sampleID, columnArg="DataType", notfound="RNA-seq")

	if ( datatype == "RNA-seq") {
		ans <- pipe.RNAalignment( sampleID, annotationFile, optionsFile, ...)
	} else if ( datatype == "DNA-seq") {
		ans <- pipe.DNAalignment( sampleID, annotationFile, optionsFile, ...)
	} else if ( datatype == "ChIP-seq") {
		ans <- pipe.ChIPalignment( sampleID, annotationFile, optionsFile, ...)
	} else if ( datatype == "RIP-seq") {
		ans <- pipe.RIPalignment( sampleID, annotationFile, optionsFile, ...)
	} else {
		cat( "\nInvalid DuffyNGS data type: ", datatype)
		cat( "\nValid types:   RNA-seq,  DNA-seq,  ChIP-seq,  RIP-seq \n")
		stop()
	}

	return( ans)
}

