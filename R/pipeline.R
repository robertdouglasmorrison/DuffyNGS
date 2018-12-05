# pipeline.R

# do the A to Z alignment pipeline on a sample

`pipeline` <- function( sampleID=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt", ...) {

	checkSampleName( sampleID)
	datatype <- getAnnotationValue( annotationFile, key=sampleID, columnArg="DataType", notfound="RNA-seq")

	if ( datatype == "RNA-seq") {
		ans <- pipe.RNAseq( sampleID, annotationFile, optionsFile, ...)
	} else if ( datatype == "DNA-seq") {
		ans <- pipe.DNAseq( sampleID, annotationFile, optionsFile, ...)
	} else if ( datatype == "ChIP-seq") {
		ans <- pipe.ChIPseq( sampleID, annotationFile, optionsFile, ...)
	} else if ( datatype == "RIP-seq") {
		ans <- pipe.RIPseq( sampleID, annotationFile, optionsFile, ...)
	} else {
		cat( "\nInvalid DuffyNGS data type: ", datatype)
		cat( "\nValid types:   RNA-seq,  DNA-seq,  ChIP-seq,  RIP-seq \n")
		stop()
	}

	return( ans)
}

