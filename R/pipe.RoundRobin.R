# pipe.RoundRobin.R

# run a N-way round robin comparison of the differential expression between a set of samples

`pipe.RoundRobin` <- function( sampleIDset, speciesID=getCurrentSpecies(), annotationFile="Annotation.txt",
				optionsFile="Options.txt", useMultiHits=TRUE, results.path=NULL,  folderName="", 
				groupColumn="Group", colorColumn="Color", method=c("Matrix", "RatioFiles"),
				altGeneMap=NULL, altGeneMapLabel=NULL, targetID=NULL, Ngenes=100, 
				geneColumnHTML=if (speciesID %in% MAMMAL_SPECIES) "NAME" else "GENE_ID", 
				keepIntergenics=FALSE, verbose=TRUE, doDE=TRUE, PLOT.FUN=NULL, 
				forceMulticore=FALSE, ...)
{

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nStarting pipe 'RoundRobin' on Sample Set: \n")
		print(sampleIDset)
		cat("\n\nUsing results from Species:  ", speciesID,"\n")
	}

	if ( is.null( targetID)) targetID <- getOptionValue( optionsFile, "targetID", notfound="HsPf", verbose=F)
	setCurrentTarget( targetID)

	annT <- readAnnotationTable( annotationFile)
	if ( ! (groupColumn %in% colnames(annT))) {
		cat( "\n\nError: Sample grouping column: '", groupColumn, "' not found in annotation file.", sep="")
		stop()
	}

	# the new 'Matrix' method needs a few option values
	method <- match.arg( method)
	if ( method == "Matrix") {
		optT <- readOptionsTable( optionsFile)
		# get the weights for ranking the results
		wt.fold <- as.numeric( getOptionValue( optT, "DE.weightFoldChange", notfound="1"))
		wt.pvalue <- as.numeric( getOptionValue( optT, "DE.weightPvalue", notfound="1"))
		# allow for a species specific minimum
		minRPKM <- as.numeric( getOptionValue( optT, "DE.minimumRPKM", speciesID=speciesID, 
				notfound="1", verbose=T))

		# do the Round Robin the new way from expression matrix data
		roundRobinCompare.MatrixMethod( sampleIDset, speciesID, annotationFile=annotationFile, optionsFile=optionsFile, 
				useMultiHits=useMultiHits, results.path=results.path, folderName=folderName, 
				groupColumn=groupColumn, colorColumn=colorColumn, Ngenes=Ngenes,
				altGeneMap=altGeneMap, altGeneMapLabel=altGeneMapLabel, 
				geneColumnHTML=geneColumnHTML, keepIntergenics=keepIntergenics, 
				wt.fold=wt.fold, wt.pvalue=wt.pvalue, minRPKM=minRPKM,
				verbose=verbose, doDE=doDE, PLOT.FUN=PLOT.FUN, forceMulticore=forceMulticore, ...)

	} else {

		# do the Round Robin the old way, from premade Ratio files
		roundRobinCompare.RatioMethod( sampleIDset, speciesID, annotationFile=annotationFile, optionsFile=optionsFile, 
				useMultiHits=useMultiHits, results.path=results.path, folderName=folderName, 
				groupColumn=groupColumn, colorColumn=colorColumn, Ngenes=Ngenes,
				altGeneMap=altGeneMap, altGeneMapLabel=altGeneMapLabel, 
				geneColumnHTML=geneColumnHTML, keepIntergenics=keepIntergenics, 
				verbose=verbose, doDE=doDE, PLOT.FUN=PLOT.FUN, ...)
	}
	geneTableToHTMLandPlotsCleanup()

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\n\nFinished pipe 'RoundRobin' on Sample Set:     ", unlist(sampleIDset), 
			"\nSpecies: ", speciesID,"\n")
	}

	return()
}

