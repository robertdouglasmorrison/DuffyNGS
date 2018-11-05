# pipe.PreAlignTasks.R

# do any and all post alignment tasks for a sample

`pipe.PreAlignTasks` <- function( sampleID,  annotationFile="Annotation.txt", optionsFile="Options.txt",
			dataType=c("RNA-seq", "ChIP-seq", "DNA-seq", "RIP-seq"), verbose=TRUE) {

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nStarting 'PreAlignTasks' on Sample:     ", sampleID, "\n")
	}
	gc()
	dataType <- match.arg( dataType)
	
	# make folders for all results...
	resultsPath <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	cat( "\n")
	subfolders <- c( "align", "fastq", "html", "ratios", "riboClear", "splicing", "summary", "transcript", "wig")
	if ( dataType == "DNA-seq") subfolders <- c( "align", "fastq", "html", "summary", "wig")
	if ( dataType == "ChIP-seq") subfolders <- c( "align", "fastq", "html", "summary", "transcript", "wig", "ChIPpeaks")
	if ( dataType == "RIP-seq") subfolders <- c( "align", "fastq", "html", "riboClear", "splicing", "summary", 
							"transcript", "wig", "RIPpeaks")

	for ( f in subfolders) {
		thisFolder <- file.path( resultsPath, f)
		if ( ! file.exists(thisFolder)) {
			dir.create( thisFolder, recursive=TRUE, showWarnings=TRUE)
			cat( "\nCreating new results subfolder: ", thisFolder)
		}
	}

	# pre-delete the files we expect to make, to help catch any pipeline crashes
	pairedEnd <- getAnnotationTrue( annotationFile, key=sampleID, "PairedEnd", notfound=FALSE)
	strandSpecific <- getAnnotationTrue( annotationFile, key=sampleID, "StrandSpecific", notfound=FALSE)
	doPairs <- (pairedEnd && strandSpecific)

	# part 1: get all the species and prefixes we may need
	mySpecies <- getCurrentTargetSpecies()
	Nspecies <- length( mySpecies)
	allPrefix <- getAllSpeciesFilePrefixes()
	where <- base::match( mySpecies, names( allPrefix))
	myPrefixes <- allPrefix[ where]

	# build all the file names we may create
	sampleIDs <- sampleID
	if (doPairs) {
		sampleIDs <- c( sampleID, getSamplePairIDs( sampleID, annotationFile))
	}

	for (sampleID in sampleIDs) {
	    allfiles <- list( 
			"align"=c( paste( sampleID, "genomic", c("bam","bam.tmpBufferFile",
					"sorted.bam","sorted.bam.bai"), sep=".")),
			"fastq"=c( paste( sampleID, c("noHits","not.genomic"), "fastq", sep="."),
				paste( sampleID, c("noHits","not.genomic"), "fastq.gz", sep=".")),
			"html"=NULL, 
			"ratios"=c( paste( sampleID, "v.*", myPrefixes, "Ratio.txt", sep=".")),  
			"splicing"=c( paste( sampleID, "splice", 
					rep(myPrefixes,times=Nspecies), "Details.txt", sep="."),
				 paste( sampleID, "splice", 
				 	rep(myPrefixes,times=Nspecies), "Summary.txt", sep="."),
				paste( sampleID, "splice.bam", sep="."),
				paste( sampleID, "splice.converted.bam", sep="."),
				paste( sampleID, "splice.converted.sorted.bam", sep="."),
				paste( sampleID, "splice.converted.sorted.bam.bai", sep=".")),
			"transcript"=c( paste( sampleID, myPrefixes, "Transcript.txt", sep=".")),  
			"riboClear"=c( paste( sampleID, "ribo.bam", sep="."),
				paste( sampleID, "ribo.converted.bam", sep=".")),
			"summary"=c( paste( sampleID, "pipeline.Summary.txt", sep="."), 
				paste( sampleID, "genomic.Summary.txt", sep="."),
				paste( sampleID, "noHits.Summary.txt", sep="."),
				paste( sampleID, "ribo.Summary.txt", sep="."),
				paste( sampleID, "splice.Summary.txt", sep="."),
				paste( sampleID, myPrefixes, "Metrics.txt", sep=".")),
			"wig"=c( paste( sampleID, myPrefixes, "Plus.wig", sep="."),
				paste( sampleID, myPrefixes, "Minus.wig", sep="."),
				paste( sampleID, myPrefixes, "WIG.rda", sep=".")),
			"."=c( paste( sampleID, "not.genomic.fq.gz", sep="."),
				paste( sampleID, "not.ribo.fq.gz", sep="."),
				paste( sampleID, "not.splice.fq.gz", sep="."),
				paste( sampleID, "convertGenomicBAM.log.txt", sep="."),
				paste( sampleID, "convertRiboBAM.log.txt", sep="."),
				paste( sampleID, "convertSpliceBAM.log.txt", sep="."))
			)

	    for( i in 1:length( allfiles)) {
		thisset <- allfiles[[ i]]
		thispath <- names( allfiles)[i]
		for( f in thisset) {
			onefile <- file.path( resultsPath, thispath, f)
			if ( file.exists( onefile)) {
				file.remove( onefile)
				cat( "\nRemoving old result file: ", onefile)
			}
			# also catch any compressed results...
			onefile <- paste( onefile, "gz", sep=".")
			if ( file.exists( onefile)) {
				file.remove( onefile)
				cat( "\nRemoving old result file: ", onefile)
			}
		}
	    }

	    # also any entire subfolders...
	    wigfolder <- file.path( resultsPath, "wig", sampleID)
	    if ( file.exists( wigfolder)) {
		cat( "\nRemoving old WIG sub-folder: ", wigfolder)
		system( paste( "rm -rf ", wigfolder))
	    }
	    vcffolder <- file.path( resultsPath, "VariantCalls", sampleID)
	    if ( file.exists( vcffolder)) {
		cat( "\nRemoving old VariantCalls sub-folder: ", vcffolder)
		system( paste( "rm -rf ", vcffolder))
	    }
	} # end of all pairs

	cat( "\n")

	return()
}

