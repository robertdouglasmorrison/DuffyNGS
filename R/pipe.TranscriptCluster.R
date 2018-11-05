# pipe.TranscriptCluster.R

# run a set of samples transcripts  through a clustering tool

`pipe.TranscriptCluster` <- function( sampleIDset=NULL, annotationFile="Annotation.txt", 
		optionsFile="Options.txt", results.path=NULL, speciesID=getCurrentSpecies(),
		intensityColumn="RPKM_M", useLog=TRUE, label="") {


	require( cluster)


	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		resultsPath <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	} else {
		resultsPath <- results.path
	}

	# other annotation fact we need...
	annT <- readAnnotationTable( annotationFile)
	if ( is.null( sampleIDset)) {
		sampleIDset <- annT$SampleID
	}

	setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()

	# gather the transcriptome files
	transcriptFiles <- paste( sampleIDset, prefix, "Transcript.txt", sep=".")
	transcriptFiles <- file.path( resultsPath, "transcript", transcriptFiles)
	allTranscripts <- expressionFileSetToMatrix( fnames=transcriptFiles, fids=sampleIDset, 
					intensityColumn=intensityColumn)

	clusterAns <- expressionCluster( allTranscripts, useLog=useLog)

	mainText <- paste( "Transcript Clustering:  ", label)
	plot( clusterAns, which.plot=2, main=mainText, font=2)

	return( invisible( list( "expressionMatrix"=allTranscripts, "cluster"=clusterAns)))
}

