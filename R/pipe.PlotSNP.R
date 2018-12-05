# pipe.PlotSNP.R

`pipe.PlotSNP` <- function( sampleIDs, seqID, position, geneID=NULL, groups=sampleIDs,
				annotationFile="Annotation.txt", optionsFile="Options.txt", 
				results.path=NULL, tailWidth=NULL, 
				plotFormat=c("","png","pdf"), plotFileName=NULL, plot.path="SNP_Plots", 
				label="", SNPtablePath="~/SNPs/", mf=NULL, start=NULL, stop=NULL, 
				verbose=TRUE, ...) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	curSpecies <- getCurrentSpecies()
	geneMap <- getCurrentGeneMap()
	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	fastaFile <- getOptionValue( optT, "genomicFastaFile", verbose=F)

	if ( ! is.null(SNPtablePath)) SNP_curSNPpath <<- SNPtablePath

	if ( !is.null(start) && !is.null(stop)) {
		start <- as.numeric( start)
		stop <- as.numeric( stop)
		dLeft <- abs( position - start)
		dRight <- abs( stop - position)
		tailWidth <- max( dLeft, dRight)
	}

	plotFormat <- match.arg( plotFormat)
	if (plotFormat != "") {
		if ( ! file.exists( plot.path)) dir.create( plot.path, recursive=T, showWarnings=F)
	}
	checkX11( bg="white", width=10, height=8)

	bamfiles <- file.path( results.path, "align", paste( sampleIDs, "genomic.sorted.bam", sep="."))
	vcffiles <- file.path( results.path, "VariantCalls", sampleIDs, 
				paste( sampleIDs, seqID, "VCF.txt", sep="."))

	# make sure we know all these sampleIDs
	annT <- readAnnotationTable( annotationFile)
	doSamples <- intersect( sampleIDs, annT$SampleID)
	if ( length( doSamples) < 1) {
		cat( "\nNo Samples found that match annotation: ", annT$SampleID, "\n")
		return()
	}
	if ( length( doSamples) < length( sampleIDs)) {
		cat( "\nSome Samples not found: ", setdiff( sampleIDs, doSamples), "\n")
	}
	doMultiSample <- ( length( doSamples) > 1)
	doMultiPosition <- ( length( position) > 1)

	# allow auto guess of how wide a tail to use
	if ( is.null( tailWidth)) {
		nplots <- max( length(doSamples), length(positions))
		tailWidth <- max( 3, round( 100 / nPlots))
	}

	# pre-get the gene map and geneID if we can
	if ( ! doMultiPosition) {
		if ( is.null( geneID)) {
			gmap <- subset.data.frame( geneMap, SEQ_ID == seqID & position >= POSITION & position <= END)
			if ( nrow( gmap) > 1) {
				bestOne <- which.min( gmap$N_EXON_BASES)
				gmap <- gmap[ bestOne, ]
			}
			geneID <- gmap$GENE_ID[1]
		} else if ( is.null( seqID)) {
			gmap <- subset.data.frame( geneMap, GENE_ID == geneID | NAME == geneID)
			geneID <- gmap$GENE_ID[1]
			seqID <- gmap$SEQ_ID[1]
		} else {
			gmap <- subset.data.frame( geneMap, GENE_ID == geneID | NAME == geneID)
		}
	} else {
		gmap <- geneID <- NULL
	}

	checkX11( bg="white", width=12, height=9)

	if (doMultiSample) {
		multiSample.plotSNP( position=position, seqID=seqID, sampleSet=sampleIDs, bamfileSet=bamfiles, 
				vcffileSet=vcffiles, fastaFile=fastaFile, groupSet=groups, 
				tailWidth=tailWidth, geneID=geneID, gmap=gmap, mf=mf, verbose=verbose, ...)
	} else if (doMultiPosition) {
		multiPosition.plotSNP( positionSet=position, seqIDset=seqID, sample=sampleIDs[1], bamfile=bamfiles[1], 
				vcffile=vcffiles[1], fastaFile=fastaFile, labelSet=label, 
				tailWidth=tailWidth, geneID=geneID, gmap=gmap, mf=mf, verbose=verbose, ...)
	} else {
		sid <- sampleIDs[1]
		bamfile <- bamfiles[1]
		vcffile <- vcffiles[1]
		plotSNP( position=position, seqID=seqID, sampleID=sid, bamfile=bamfile, 
				vcffile=vcffile, fastaFile=fastaFile,
				tailWidth=tailWidth, mode="single", geneID=geneID, gmap=gmap, 
				label=label, verbose=verbose, ...)
	}
	if ( plotFormat != "") {
		if ( is.null( plotFileName)) {
			geneID <- shortGeneName( gmap$GENE_ID[1], keep=1)
			plotfile <- paste( geneID, position, plotFormat, sep=".")
		} else {
			plotfile <- plotFileName
		}
		plotfile <- file.path( plot.path, plotfile)
		if ( plotFormat == "png") dev.print( png, plotfile, width=1000, "bg"="white") #, type="Xlib")
		if ( plotFormat == "pdf") dev.print( pdf, plotfile, width=10, "bg"="white")
	}
	return()
}


`pipe.PlotGeneSNPs` <- function( sampleIDs, geneID, tailWidth=100, groups=sampleIDs,
				annotationFile="Annotation.txt", optionsFile="Options.txt", results.path=NULL, 
				plotFormat=c("","png","pdf"), plot.path="SNP_Plots", label="", SNPtablePath="~/SNPs/", 
				mf=NULL, start=NULL, stop=NULL, verbose=TRUE, ...) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	curSpecies <- getCurrentSpecies()
	geneMap <- getCurrentGeneMap()
	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	fastaFile <- getOptionValue( optT, "genomicFastaFile", verbose=F)

	if ( ! is.null(SNPtablePath)) SNP_curSNPpath <<- SNPtablePath

	# map from a gene name to the seq, position that the SNP tool wants
	gmap <- subset.data.frame( geneMap, GENE_ID == geneID)
	if ( nrow( gmap) < 1) {
		gmap <- subset.data.frame( geneMap, NAME == geneID)
		if ( nrow( gmap) < 1) {
			cat( "\nNo Gene found that matches: ", geneID)
			return()
		}
		geneID <- gmap$GENE_ID[1]
	}
	seqID <- gmap$SEQ_ID[1]
	midpt <- round( (gmap$POSITION[1] + gmap$END[1]) / 2)
	geneHalfWidth <- round( abs(gmap$END[1] - gmap$POSITION[1] + 1) / 2)

	if ( !is.null(start) && !is.null(stop)) {
		start <- as.numeric( start)
		stop <- as.numeric( stop)
		midpt <- round( (start+stop) / 2)
		geneHalfWidth <- round( (stop-start+1) / 2)
	}

	plotFormat <- match.arg( plotFormat)
	if (plotFormat != "") {
		if ( ! file.exists( plot.path)) dir.create( plot.path, recursive=T, showWarnings=F)
	}
	checkX11( bg="white", width=12, height=9)

	# get the data we need
	bamfiles <- file.path( results.path, "align", paste( sampleIDs, "genomic.sorted.bam", sep="."))
	vcffiles <- file.path( results.path, "VariantCalls", sampleIDs, 
				paste( sampleIDs, seqID, "VCF.txt", sep="."))

	# make sure we know all these sampleIDs
	annT <- readAnnotationTable( annotationFile)
	doSamples <- intersect( sampleIDs, annT$SampleID)
	if ( length( doSamples) < 1) {
		cat( "\nNo Samples found that match annotation: ", annT$SampleID, "\n")
		return()
	}
	if ( length( doSamples) < length( sampleIDs)) {
		cat( "\nSome Samples not found: ", setdiff( sampleIDs, doSamples), "\n")
	}
	doMulti <- ( length( doSamples) > 1)

	myTailWidth <- tailWidth + geneHalfWidth
	if (doMulti) {
		multiSample.plotSNP( position=midpt, seqID=seqID, sampleSet=sampleIDs, bamfileSet=bamfiles, 
				vcffileSet=vcffiles, fastaFile=fastaFile, groupSet=groups, 
				tailWidth=myTailWidth, geneID=geneID, gmap=gmap, show.legends="none", mf=mf, verbose=verbose, ...)
	} else {
		sid <- sampleIDs[1]
		bamfile <- bamfiles[1]
		plotSNP( position=midpt, seqID=seqID, sampleID=sid, bamfile=bamfiles, 
				vcffile=vcffiles, fastaFile=fastaFile,
				tailWidth=myTailWidth, mode="single", geneID=geneID, gmap=gmap, label=label, 
				show.legends="gene", verbose=verbose, ...)
	}
	if ( plotFormat != "") {
		geneID <- shortGeneName( gmap$GENE_ID[1], keep=1)
		plotfile <- paste( geneID, "FullLength", plotFormat, sep=".")
		plotfile <- file.path( plot.path, plotfile)
		if ( plotFormat == "png") dev.print( png, plotfile, width=1000, "bg"="white") #, type="Xlib")
		if ( plotFormat == "pdf") dev.print( pdf, plotfile, width=10, "bg"="white")
	}
	return()
}

