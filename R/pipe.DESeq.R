# pipe.DESeq.R

# run the DESeq tool to assess differential expression

`pipe.DESeq` <- function( sampleIDset, speciesID=getCurrentSpecies(), annotationFile="Annotation.txt",
				optionsFile="Options.txt", useMultiHits=TRUE, results.path=NULL, 
				groupColumn="Group", colorColumn="Color", folderName="", 
				altGeneMap=NULL, altGeneMapLabel=NULL, targetID=NULL,
				Ngenes=100, geneColumnHTML=if (speciesID %in% MAMMAL_SPECIES) "NAME" else "GENE_ID", 
				keepIntergenics=FALSE, verbose=!interactive(), label="", 
				doDE=TRUE, PLOT.FUN=NULL, ...)
{

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nStarting pipe 'DESeq' on Sample Set: \n")
		print(sampleIDset)
		cat("\n", label, "\n\nUsing results from Species:  ", speciesID,"\n")
	}

	# set up for this species...
	optT <- readOptionsTable( optionsFile)
	if ( is.null( targetID)) targetID <- getOptionValue( optT, "targetID", notfound="HsPf", verbose=F)
	setCurrentTarget( targetID)
	setCurrentSpecies( speciesID)
	
	# sanity check on the inputs...
	if ( length( unlist(sampleIDset)) < 2) stop( "DESeq requires at least 2 sampleIDs")
	if ( base::nchar( folderName) < 1) stop( "DESeq needs an explicit 'folderName' argument...")

	annT <- readAnnotationTable( annotationFile)
	if ( ! (groupColumn %in% colnames(annT))) {
		cat( "\n\nError: Sample grouping column: '", groupColumn, "' not found in annotation file.", sep="")
		stop()
	}
	flatSamples <- unlist( sampleIDset)
	where <- match( flatSamples, annT$SampleID, nomatch=0)
	if ( any( where == 0)) stop( "Some named SampleIDs not in Annotation File")
	RP_samples <- unique( flatSamples)
	where <- match( RP_samples, annT$SampleID)
	RP_groups <- annT[[ groupColumn]][where]
	RP_colors <- annT[[ colorColumn]][where]
	RP_species <- speciesID
	RP_prefix <- getCurrentSpeciesFilePrefix()
	unique_RP_groups <- base::sort( unique.default( RP_groups))
	N_RP_groups <- length( unique_RP_groups)
	if ( N_RP_groups < 2) {
		stop( paste( "\nDifferential Expression needs at least 2 groups.  Check Annotation column:", groupColumn))
	} else {
		cat( "\nSample counts by group:\n")
		print( table( RP_groups))
	}

	RP_altGeneMapLabel <- altGeneMapLabel
	HTML_geneColumn <- geneColumnHTML

	# set up directories to read from or write to...
	if ( is.null( results.path)) {
		resultsPath <- getOptionValue( optT, "results.path", notfound=".")
	} else {
		resultsPath <- results.path
	}
	RP_path <- file.path( resultsPath, "DESeq", paste( RP_prefix, folderName, sep="."))
	if ( ! file.exists( RP_path)) dir.create( RP_path, recursive=T, showWarnings=F)

	transFIDs <- RP_samples
	# we could be given a non-standard geneMap to use...
	if ( is.null( altGeneMap)) {

		# regular...
		gmap <- getCurrentGeneMap()
		RP_altGeneMapLabel <- NULL
		ratiosPath <- file.path( resultsPath, "ratios")
		transPath <- file.path( resultsPath, "transcript")
		transFileSet <- paste( RP_samples, RP_prefix, "Transcript.txt", sep=".")

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

		cat( "\n\nDoing Alternate Gene Map DESeq for:  ", altGeneMapLabel, "\n")
		RP_altGeneMapLabel <- altGeneMapLabel
		ratiosPath <- file.path( resultsPath, "ratios", altGeneMapLabel)
		transPath <- file.path( resultsPath, "transcript", altGeneMapLabel)
		transFileSet <- paste( RP_samples, RP_prefix, altGeneMapLabel, "Transcript.txt", sep=".")
	}
	transFileSet <- file.path( transPath, transFileSet)

	# get the weights for ranking the results
	wt.folds <- as.numeric( getOptionValue( optT, "DE.weightFoldChange", notfound="1"))
	wt.pvalues <- as.numeric( getOptionValue( optT, "DE.weightPvalue", notfound="1"))

	# allow for a species specific minimum
	minRPKM <- as.numeric( getOptionValue( optT, "DE.minimumRPKM", speciesID=speciesID, 
			notfound="1", verbose=T))

	intensityColumn <- if (useMultiHits) "READS_M" else "READS_U"

	# ready to do the DESeq...
	genesToPlot <- vector()

	if (doDE) {
	    for ( targetgroup in sort( unique( RP_groups))) {

		cat( "\n\nDoing DESeq on:    ", targetgroup, "\n")

		out <- DESeq.DiffExpress( transFileSet, transFIDs, groupSet=RP_groups, targetGroup=targetgroup, 
				geneColumn="GENE_ID", intensityColumn=intensityColumn,
				minimumRPKM=minRPKM, keepIntergenics=keepIntergenics, 
				average.FUN=sqrtmean,
				missingGenes="fill", wt.folds=wt.folds, wt.pvalues=wt.pvalues,
				wt.dists=1, ...)

		outfile <- paste( targetgroup, RP_prefix, "DESeq.Ratio.txt", sep=".")
		if ( ! is.null( altGeneMap)) {
			outfile <- paste( targetgroup, RP_prefix, altGeneMapLabel, "DESeq.Ratio.txt", sep=".")
		}
		outfile <- file.path( RP_path, outfile)

		# add Human ID terms if needed
		extraCols <- 0
		if ( RP_prefix %in% MAMMAL_PREFIXES) {
			out <- addHumanIDterms( out)
			extraCols <- 2
		}
		if ( RP_prefix %in% ORIGID_PARASITE_PREFIXES) {
			out <- addOrigIDterms( out)
			extraCols <- 1
		}
		if ( RP_prefix %in% BACTERIA_PREFIXES) {
			out <- addNameIDterms( out)
			extraCols <- 1
		}

		write.table( out, outfile, sep="\t", quote=F, row.names=F)
		cat( "\nWrote DESeq Gene Data:  ", outfile)

		# HTML too...
		htmlFile1 <- sub( "Ratio.txt$", "UP.html", basename(outfile))
		htmlFile2 <- sub( "Ratio.txt$", "DOWN.html", basename(outfile))
		htmlPath <- RP_path

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
		title1 <- paste( "DESeq: &nbsp;  Genes most up-regulated in group: &nbsp;  ", targetgroup)
		title2 <- paste( "DESeq: &nbsp;  Genes most down-regulated in group: &nbsp;  ", targetgroup)
		if ( ! is.null( RP_altGeneMapLabel)) {
			title1 <- paste( "DESeq: &nbsp;  ", RP_altGeneMapLabel, 
					"  most up-regulated in group: &nbsp;  ", targetgroup)
			title2 <- paste( "DESeq: &nbsp;  ", RP_altGeneMapLabel, 
					"  most down-regulated in group: &nbsp;  ", targetgroup)
			# for plots of varGenes we need to fudge a few items...
			out <- addAltGeneIdentifierColumn( out)
			extraCols <- extraCols + 1
			fullGname <- out$GENE_ID
			HTML_geneColumn <- "GENE_ID"
			#if ( regexpr( "vargene", tolower(RP_altGeneMapLabel)) > 0) {
			#	out <- cbind( "DOMAIN_ID"=fullGname, out)
			#	out$GENE_ID <- sub( "::.*", "", fullGname)
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
			for ( k in (5+extraCols):nColShow) {
				out1[[k]] <- formatC( out1[[k]], format="f", digits=2, big.mark=",")
			}
			colnames(out1)[c(3:4 + extraCols)] <- c( "Log2 Fold", "P Value")
			colnames(out1)[(5+extraCols):nColShow] <- gsub( "_", " ", colnames(out1)[(5+extraCols):nColShow])
			colnames(out1)[(5+extraCols):nColShow] <- gsub( "\\.", " ", colnames(out1)[(5+extraCols):nColShow])
			# write it
			geneTableToHTMLandPlots( geneDF=out1[ , 1:nColShow], RP_samples, RP_colors, N=Nshow, 
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
			for ( k in (5+extraCols):nColShow) {
				out2[[k]] <- formatC( out2[[k]], format="f", digits=2, big.mark=",")
			}
			colnames(out2)[c(3:4 + extraCols)] <- c( "Log2 Fold", "P Value")
			colnames(out2)[(5+extraCols):nColShow] <- gsub( "_", " ", colnames(out2)[(5+extraCols):nColShow])
			colnames(out2)[(5+extraCols):nColShow] <- gsub( "\\.", " ", colnames(out2)[(5+extraCols):nColShow])
			# write it
			geneTableToHTMLandPlots( geneDF=out2[ , 1:nColShow], RP_samples, RP_colors, N=Nshow, 
				title=title2, 
				htmlFile=htmlFile2, html.path=htmlPath, results.path=resultsPath, 
				makePlots=FALSE)
			genesToPlot <- base::union( genesToPlot, unique.default( rev(fullGname)[1:Nshow]))
		}
	    }
	} else {
		cat( "\nSkipping DE...  Gathering genes for plots..")
		htmlPath <- RP_path
		for (targetgroup in sort( unique( RP_groups))) {
			outfile <- paste( targetgroup, RP_prefix, "DESeq.Ratio.txt", sep=".")
			if ( ! is.null( altGeneMap)) {
				outfile <- paste( targetgroup, RP_prefix, altGeneMapLabel, "DESeq.Ratio.txt", sep=".")
			}
			outfile <- file.path( RP_path, outfile)
			tmp <- read.delim( outfile, as.is=T)
			genesToPlot <- c( genesToPlot, tmp$GENE_ID[1:Ngenes], rev(tmp$GENE_ID)[1:Ngenes])
		}
		genesToPlot <- unique.default( genesToPlot)
	}

	# after all tables and results, make those gene plots
	genesToPlot <- sort( genesToPlot)
	if ( ! is.null(altGeneMap)) genesToPlot <- sort( unique( sub( "::.+", "", genesToPlot)))
	if ( is.null(PLOT.FUN) || is.function( PLOT.FUN)) {
		geneTableToHTMLandPlots( geneDF=NULL, RP_samples, RP_colors, N=Ngenes, htmlFile=htmlFile, 
				html.path=htmlPath, results.path=resultsPath, makePlots=TRUE, 
				genesToPlot=genesToPlot, label=folderName, 
				geneNameColumn=HTML_geneColumn, PLOT.FUN=PLOT.FUN, ...)
	}

	tm <- expressionFileSetToMatrix( fnames=transFileSet, fids=transFIDs, intensityColumn=intensityColumn)
	gnames <- rownames(tm)

	if ( is.null( RP_altGeneMapLabel)) {
		gprods <- gene2Product( gnames)
		outfile <- file.path( RP_path, paste( "All", RP_prefix, "GeneData.txt",sep="."))
		cat( "\nWriting 'all transcripts' gene data:  ", outfile, "\n")
	} else {
		gprods <- gmap$PRODUCT[ base::match( gnames, gmap$GENE_ID)]
		outfile <- file.path( RP_path, paste( "All.", RP_prefix, ".", RP_altGeneMapLabel, 
					"Data.txt", sep=""))
		cat( "\nWriting 'all transcripts' ", RP_altGeneMapLabel, " data:  ", outfile, "\n")
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

		if ( is.null( RP_altGeneMapLabel)) {
			pltText <- paste( "Transcriptome Clustering:   ", folderName,
					"\nTranscriptomes for species:   ", speciesID)
			pngFile <- file.path( RP_path, paste( RP_prefix,"Cluster",funcName[i],"png",sep="."))
		} else {
			pltText <- paste( "Transcriptome Clustering:   ", folderName, 
					"\nTranscriptomes for species:   ", speciesID,
					"    using geneMap:  ", RP_altGeneMapLabel)
			pngFile <- file.path( RP_path, paste( RP_prefix, RP_altGeneMapLabel, 
					"Cluster", funcName[i], "png", sep="."))
		}

		clusterAns <- expressionCluster( tm, useLog=TRUE, FUN=func[[i]])
		plot( clusterAns, which=2, main=pltText, sub=subtitle[i], font=2)

		cat( "\nMaking cluster plot: ", pngFile)
		png( filename=pngFile, width=1000, height=700, bg="white")
		plot( clusterAns, which=2, main=pltText, sub=subtitle[i], font=2)
		dev.off()
	    }

	    # PCA plot too...
		pltText <- paste( "Transcriptome PCA:   ", folderName,
				"\nTranscriptomes for species:   ", speciesID)
		pngFile <- file.path( RP_path, paste( RP_prefix,"PCA.png",sep="."))
		matrix.PCAplot( tm, main=pltText, col=RP_colors)
		png( filename=pngFile, width=800, height=800, bg="white")
		matrix.PCAplot( tm, main=pltText, col=RP_colors)
		dev.off()

	} else {
		if ( ncol(tm) < 3) cat( "\nNot able to cluster fewer than 3 transcripts...")
	}
	
	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\n\nFinished pipe 'DESeq' on Sample Set:     ", unlist(sampleIDset), 
			"\nSpecies: ", speciesID,"\n", label, "\n")
	}
	return()
}

