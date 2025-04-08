# pipe.LimmaVoom.R

# run the LimmaVoom tool to assess differential expression

`pipe.LimmaVoom` <- function( sampleIDset, speciesID=getCurrentSpecies(), annotationFile="Annotation.txt",
				optionsFile="Options.txt", useMultiHits=TRUE, results.path=NULL, 
				groupColumn="Group", colorColumn="Color", folderName="", 
				altGeneMap=NULL, altGeneMapLabel=NULL, targetID=NULL,
				Ngenes=100, geneColumnHTML=if (speciesID %in% MAMMAL_SPECIES) "NAME" else "GENE_ID", 
				keepIntergenics=FALSE, verbose=TRUE, label="", 
				doDE=TRUE, PLOT.FUN=NULL, adjust.lowReadCounts=TRUE, ...)
{

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nStarting pipe 'LimmaVoom' on Sample Set: \n")
		print(sampleIDset)
		cat("\n", label, "\n\nUsing results from Species:  ", speciesID,"\n")
	}

	# setting the target can reset the speciesID, so do that explicitly
	optT <- readOptionsTable( optionsFile)
	wantedSpeciesID <- speciesID
	if ( is.null( targetID)) targetID <- getOptionValue( optT, "targetID", notfound="HsPf", verbose=F)
	setCurrentTarget( targetID)
	setCurrentSpecies( wantedSpeciesID)
	
	# sanity check on the inputs...
	if ( length( unlist(sampleIDset)) < 2) stop( "LimmaVoom requires at least 2 sampleIDs")
	if ( base::nchar( folderName) < 1) stop( "LimmaVoom needs an explicit 'folderName' argument...")

	annT <- readAnnotationTable( annotationFile)
	if ( ! (groupColumn %in% colnames(annT))) {
		cat( "\n\nError: Sample grouping column: '", groupColumn, "' not found in annotation file.", sep="")
		stop()
	}
	flatSamples <- unlist( sampleIDset)
	where <- match( flatSamples, annT$SampleID, nomatch=0)
	if ( any( where == 0)) stop( "Some named SampleIDs not in Annotation File")
	LV_samples <- unique( flatSamples)
	where <- match( LV_samples, annT$SampleID)
	LV_groups <- checkGroupNames( annT[[ groupColumn]][where])
	LV_colors <- annT[[ colorColumn]][where]
	LV_species <- speciesID
	LV_prefix <- getCurrentSpeciesFilePrefix()
	unique_LV_groups <- base::sort( unique.default( LV_groups))
	N_LV_groups <- length( unique_LV_groups)
	if ( N_LV_groups < 2) {
		stop( paste( "\nDifferential Expression needs at least 2 groups.  Check Annotation column:", groupColumn))
	} else {
		cat( "\nSample counts by group:\n")
		print( table( LV_groups))
	}

	LV_altGeneMapLabel <- altGeneMapLabel
	HTML_geneColumn <- geneColumnHTML

	# set up directories to read from or write to...
	if ( is.null( results.path)) {
		resultsPath <- getOptionValue( optT, "results.path", notfound=".")
	} else {
		resultsPath <- results.path
	}
	LV_path <- file.path( resultsPath, "LimmaVoom", paste( LV_prefix, folderName, sep="."))
	if ( ! file.exists( LV_path)) dir.create( LV_path, recursive=T, showWarnings=F)

	transFIDs <- LV_samples
	# we could be given a non-standard geneMap to use...
	if ( is.null( altGeneMap)) {

		# regular...
		gmap <- getCurrentGeneMap()
		LV_altGeneMapLabel <- NULL
		ratiosPath <- file.path( resultsPath, "ratios")
		transPath <- file.path( resultsPath, "transcript")
		transFileSet <- paste( LV_samples, LV_prefix, "Transcript.txt", sep=".")

	} else {
		# alternate...
		gmap <- altGeneMap
		if ( is.character( altGeneMap)) {
			gmap <- read.delim( file=altGeneMap, as.is=TRUE)
		}
		if ( ! all( c("GENE_ID", "SEQ_ID") %in% colnames(gmap))) 
			stop( paste( "Invalid alternate geneMap",
				"does not have required GENE_ID, SEQ_ID columns."))
		if ( is.null( altGeneMapLabel) || base::nchar( altGeneMapLabel) < 1) 
			stop( "Missing required file name text term 'altGeneMapLabel' ")

		cat( "\n\nDoing Alternate Gene Map LimmaVoom for:  ", altGeneMapLabel, "\n")
		LV_altGeneMapLabel <- altGeneMapLabel
		ratiosPath <- file.path( resultsPath, "ratios", altGeneMapLabel)
		transPath <- file.path( resultsPath, "transcript", altGeneMapLabel)
		transFileSet <- paste( LV_samples, LV_prefix, altGeneMapLabel, "Transcript.txt", sep=".")
	}
	transFileSet <- file.path( transPath, transFileSet)

	# get the weights for ranking the results
	wt.folds <- as.numeric( getOptionValue( optT, "DE.weightFoldChange", notfound="1"))
	wt.pvalues <- as.numeric( getOptionValue( optT, "DE.weightPvalue", notfound="1"))

	intensityColumn <- if (useMultiHits) "READS_M" else "READS_U"

	# ready to do the LimmaVoom...
	genesToPlot <- vector()
	dev.type <- getPlotDeviceType( optT)

	if (doDE) {
	    for ( targetgroup in sort( unique( LV_groups))) {

		cat( "\n\nDoing LimmaVoom on:    ", targetgroup, "\n")

		out <- LimmaVoom.DiffExpress( transFileSet, transFIDs, groupSet=LV_groups, targetGroup=targetgroup, 
				geneColumn="GENE_ID", intensityColumn=intensityColumn,
				keepIntergenics=keepIntergenics, 
				average.FUN=sqrtmean,
				missingGenes="fill", wt.folds=wt.folds, wt.pvalues=wt.pvalues, wt.dists=1,
				adjust.lowReadCounts=adjust.lowReadCounts)

		outfile <- paste( targetgroup, LV_prefix, "LV.Ratio.txt", sep=".")
		if ( ! is.null( altGeneMap)) {
			outfile <- paste( targetgroup, LV_prefix, altGeneMapLabel, "LV.Ratio.txt", sep=".")
		}
		outfile <- file.path( LV_path, outfile)

		# add Human ID terms if needed
		extraCols <- 0
		if ( LV_prefix %in% MAMMAL_PREFIXES) {
			out <- addHumanIDterms( out)
			extraCols <- 2
		}
		if ( LV_prefix %in% ORIGID_PARASITE_PREFIXES) {
			out <- addOrigIDterms( out)
			extraCols <- 1
		}
		if ( LV_prefix %in% BACTERIA_PREFIXES) {
			out <- addNameIDterms( out)
			extraCols <- 1
		}

		write.table( out, outfile, sep="\t", quote=F, row.names=F)
		cat( "\nWrote LimmaVoom Gene Data:  ", outfile)

		# HTML too...
		htmlFile1 <- sub( "Ratio.txt$", "UP.html", basename(outfile))
		htmlFile2 <- sub( "Ratio.txt$", "DOWN.html", basename(outfile))
		htmlPath <- LV_path

		# simplify the names?
		fullGname <- out$GENE_ID
		if ( HTML_geneColumn != "GENE_ID") {
			where <- base::match( fullGname, gmap$GENE_ID, nomatch=0)
			newGname <- fullGname
			newGname[ where > 0] <-  gmap[ , HTML_geneColumn][ where]
			out$GENE_ID <- newGname
		}

		# special mods for altGeneMap...
		Nshow <- Ngenes
		title1 <- paste( "LimmaVoom: &nbsp;  Genes most up-regulated in group: &nbsp;  ", targetgroup)
		title2 <- paste( "LimmaVoom: &nbsp;  Genes most down-regulated in group: &nbsp;  ", targetgroup)
		if ( ! is.null( LV_altGeneMapLabel)) {
			title1 <- paste( "LimmaVoom: &nbsp;  ", LV_altGeneMapLabel, 
					"  most up-regulated in group: &nbsp;  ", targetgroup)
			title2 <- paste( "LimmaVoom: &nbsp;  ", LV_altGeneMapLabel, 
					"  most down-regulated in group: &nbsp;  ", targetgroup)
			# for plots of varGenes we need to fudge a few items...
			out <- addAltGeneIdentifierColumn( out)
			extraCols <- extraCols + 1
			fullGname <- out$GENE_ID
			HTML_geneColumn <- "GENE_ID"
			#if ( regexpr( "vargene", tolower(LV_altGeneMapLabel)) > 0) {
			#	out <- cbind( "DOMAIN_ID"=fullGname, out)
			#	out$GENE_ID <- sub( "::.*", "", fullGname)
			#	extraCols <- extraCols + 1
			#	fullGname <- out$GENE_ID
			#	HTML_geneColumn <- "GENE_ID"
			#	if ( "ORIG_ID" %in% colnames(out)) out$ORIG_ID <- gene2OrigID( out$GENE_ID)
			#}
		}

		nColShow <- ncol(out)
		# only keep those that are UP
		out1 <- out[ out$LOG2FOLD > 0, ]
		if ( nrow(out1) > 0) {
			# clean up formats...
			out1$PRODUCT <- gsub( "   ", " &nbsp; ", out1$PRODUCT)
			out1$LOG2FOLD <- formatC( out1$LOG2FOLD, format="f", digits=3, flag="+")
			out1$PVALUE <- formatC( out1$PVALUE, format="e", digits=2)
			out1$FDR <- formatC( out1$FDR, format="e", digits=2)
			out1$PIVALUE <- formatC( out1$PIVALUE, format="f", digits=3, flag="+")
			for ( k in (7+extraCols):nColShow) {
				out1[[k]] <- formatC( out1[[k]], format="f", digits=2, big.mark=",")
			}
			colnames(out1)[c(3:6 + extraCols)] <- c( "Log2 Fold", "P Value", "FDR", "PI Value")
			colnames(out1)[(7+extraCols):nColShow] <- gsub( "_", " ", colnames(out1)[(7+extraCols):nColShow])
			colnames(out1)[(7+extraCols):nColShow] <- gsub( "\\.", " ", colnames(out1)[(7+extraCols):nColShow])
			# write it
			geneTableToHTMLandPlots( geneDF=out1[ , 1:nColShow], LV_samples, LV_colors, N=Nshow, 
				title=title1, 
				htmlFile=htmlFile1, html.path=htmlPath, results.path=resultsPath, 
				makePlots=FALSE)
			genesToPlot <- base::union( genesToPlot, unique.default( fullGname[1:Nshow]))
		}

		# for the DOWN table, flip it and use the DOWN Pvalues and adjust the ranks
		out2 <- out[ rev( 1:nrow(out)), ]
		# only keep those that are DOWN
		out2 <- out2[ out2$LOG2FOLD < 0, ]
		if ( nrow(out2) > 0) {
			# clean up formats...
			out2$PRODUCT <- gsub( "   ", " &nbsp; ", out2$PRODUCT)
			out2$LOG2FOLD <- formatC( out2$LOG2FOLD, format="f", digits=3, flag="+")
			out2$PVALUE <- formatC( out2$PVALUE, format="e", digits=2)
			out2$FDR <- formatC( out2$FDR, format="e", digits=2)
			out2$PIVALUE <- formatC( out2$PIVALUE, format="f", digits=3, flag="+")
			for ( k in (6+extraCols):nColShow) {
				out2[[k]] <- formatC( out2[[k]], format="f", digits=2, big.mark=",")
			}
			colnames(out2)[c(3:6 + extraCols)] <- c( "Log2 Fold", "P Value", "FDR", "PI Value")
			colnames(out2)[(7+extraCols):nColShow] <- gsub( "_", " ", colnames(out2)[(7+extraCols):nColShow])
			colnames(out2)[(7+extraCols):nColShow] <- gsub( "\\.", " ", colnames(out2)[(7+extraCols):nColShow])
			# write it
			geneTableToHTMLandPlots( geneDF=out2[ , 1:nColShow], LV_samples, LV_colors, N=Nshow, 
				title=title2, 
				htmlFile=htmlFile2, html.path=htmlPath, results.path=resultsPath, 
				makePlots=FALSE)
			genesToPlot <- base::union( genesToPlot, unique.default( rev(fullGname)[1:Nshow]))
		}
	    }
	} else {
		cat( "\nSkipping DE...  Gathering genes for plots..")
		htmlPath <- LV_path
		for (targetgroup in sort( unique( LV_groups))) {
			outfile <- paste( targetgroup, LV_prefix, "LV.Ratio.txt", sep=".")
			if ( ! is.null( altGeneMap)) {
				outfile <- paste( targetgroup, LV_prefix, altGeneMapLabel, "LV.Ratio.txt", sep=".")
			}
			outfile <- file.path( LV_path, outfile)
			if ( ! file.exists( outfile)) {
				cat( "\nLimmaVoom Ratios file not found..  Skipping.")
				next
			}
			tmp <- read.delim( outfile, as.is=T)
			genesToPlot <- c( genesToPlot, tmp$GENE_ID[1:Ngenes], rev(tmp$GENE_ID)[1:Ngenes])
		}
		genesToPlot <- unique.default( genesToPlot)
		# if we were just re-plotting and found no genes, quit now
		if ( length( genesToPlot) == 0) return()
	}

	# after all tables and results, make those gene plots
	genesToPlot <- sort( genesToPlot)
	if ( ! is.null(altGeneMap)) genesToPlot <- sort( unique( sub( "::.+", "", genesToPlot)))

	# put in chromosomal order?
	whereGmap <- match( genesToPlot, gmap$GENE_ID, nomatch=NA)
	genesToPlot <- genesToPlot[ order( whereGmap)]

	if ( is.null(PLOT.FUN) || is.function(PLOT.FUN)) {
		geneTableToHTMLandPlots( geneDF=NULL, LV_samples, LV_colors, N=Ngenes, htmlFile=htmlFile, 
				html.path=htmlPath, results.path=resultsPath, makePlots=TRUE, 
				genesToPlot=genesToPlot, label=folderName, 
				geneNameColumn=HTML_geneColumn, PLOT.FUN=PLOT.FUN, ...)
	}

	tm <- expressionFileSetToMatrix( fnames=transFileSet, fids=transFIDs, intensityColumn=intensityColumn)
	gnames <- rownames(tm)

	if ( is.null( LV_altGeneMapLabel)) {
		gprods <- gene2Product( gnames)
		outfile <- file.path( LV_path, paste( "All", LV_prefix, "GeneData.txt",sep="."))
		cat( "\nWriting 'all transcripts' gene data:  ", outfile, "\n")
	} else {
		gprods <- gmap$PRODUCT[ base::match( gnames, gmap$GENE_ID)]
		outfile <- file.path( LV_path, paste( "All.", LV_prefix, ".", LV_altGeneMapLabel, 
					"Data.txt", sep=""))
		cat( "\nWriting 'all transcripts' ", LV_altGeneMapLabel, " data:  ", outfile, "\n")
	}
		
	outTM <- data.frame( "GENE_ID"=gnames, "PRODUCT"=gprods, tm, stringsAsFactors=FALSE)
	rownames(outTM) <- 1:nrow(outTM)
	write.table( outTM, file=outfile, sep="\t", quote=F, row.names=F)
	
	# make some cluster images...
	if ( ncol(tm) > 2 && (is.null(PLOT.FUN) || is.function(PLOT.FUN))) {

	    # there are two good cluster tools...  let's do both
	    require( cluster)
	    func <- list( diana, agnes)
	    funcName <- c( "Divide", "Aggregate")
	    subtitle <- c( "Divisive hierarchical clustering (DIANA)", 
	    			"Agglomerative hierarchical clustering (AGNES)")

	    for ( i in 1:2) {

		if ( is.null( LV_altGeneMapLabel)) {
			pltText <- paste( "Transcriptome Clustering:   ", folderName,
					"\nSpecies: ", speciesID, "    Expression Units: ", intensityColumn)
			plotFile <- file.path( LV_path, paste( LV_prefix,"Cluster",funcName[i], dev.type, sep="."))
		} else {
			pltText <- paste( "Transcriptome Clustering:   ", folderName, 
					"\nSpecies: ", speciesID, "    Expression Units: ", intensityColumn,
					"    using geneMap:  ", LV_altGeneMapLabel)
			plotFile <- file.path( LV_path, paste( LV_prefix, LV_altGeneMapLabel, 
					"Cluster", funcName[i], dev.type, sep="."))
		}

		clusterAns <- expressionCluster( tm, useLog=TRUE, FUN=func[[i]])
		plot( clusterAns, which=2, main=pltText, sub=subtitle[i], font=2)

		cat( "\nMaking cluster plot: ", plotFile)
		openPlot( plotFile, bg='white')
		plot( clusterAns, which=2, main=pltText, sub=subtitle[i], font=2)
		dev.off()
	    }

	    # PCA plot too...
		pltText <- paste( "Transcriptome PCA:   ", folderName,
				"\nSpecies: ", speciesID, "    Expression Units: ", intensityColumn)
		plotFile <- file.path( LV_path, paste( LV_prefix,"PCA", dev.type, sep="."))
		matrix.PCAplot( tm, main=pltText, col=LV_colors)
		openPlot( plotFile, bg='white')
		matrix.PCAplot( tm, main=pltText, col=LV_colors)
		dev.off()

	} else {
		if ( ncol(tm) < 3) cat( "\nToo few samples to cluster...")
	}
	
	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\n\nFinished pipe 'LimmaVoom' on Sample Set:     ", unlist(sampleIDset), 
			"\nSpecies: ", speciesID,"\n", label, "\n")
	}
	return()
}

