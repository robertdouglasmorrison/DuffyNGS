# pipe.RankProduct.R

# run the Rank Product tool to assess differential expression

`pipe.RankProduct` <- function( sampleIDset, speciesID=getCurrentSpecies(), annotationFile="Annotation.txt",
				optionsFile="Options.txt", 
				useMultiHits=TRUE, results.path=NULL,  
				groupColumn="Group", colorColumn="Color", folderName="", 
				altGeneMap=NULL, altGeneMapLabel=NULL, targetID=NULL,
				Ngenes=100, geneColumnHTML=if (speciesID %in% MAMMAL_SPECIES) "NAME" else "GENE_ID", 
				keepIntergenics=FALSE, verbose=!interactive(), 
				nSimulations=100, label="", doDE=TRUE, PLOT.FUN=NULL, forceMulticore=FALSE, ...)
{

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nStarting pipe 'RankProduct' on Sample Set: \n")
		print(sampleIDset)
		cat("\n", label, "\n\nUsing results from Species:  ", speciesID,"\n")
	}

	# set up for this species...
	optT <- readOptionsTable( optionsFile)
	if ( is.null( targetID)) targetID <- getOptionValue( optT, "targetID", notfound="HsPf", verbose=F)
	setCurrentTarget( targetID)
	setCurrentSpecies( speciesID)
	
	# sanity check on the inputs...
	if ( length( unlist(sampleIDset)) < 2) stop( "RankProduct requires at least 2 sampleIDs")
	if ( base::nchar( folderName) < 1) stop( "RankProduct needs an explicit 'folderName' argument...")

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
	RP_groups <- checkGroupNames( annT[[ groupColumn]][where])
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
	RP_path <- file.path( resultsPath, "RankProduct", paste( RP_prefix, folderName, sep="."))
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

		cat( "\n\nDoing Alternate Gene Map RankProduct for:  ", altGeneMapLabel, "\n")
		RP_altGeneMapLabel <- altGeneMapLabel
		ratiosPath <- file.path( resultsPath, "ratios", altGeneMapLabel)
		transPath <- file.path( resultsPath, "transcript", altGeneMapLabel)
		transFileSet <- paste( RP_samples, RP_prefix, altGeneMapLabel, "Transcript.txt", sep=".")
	}
	transFileSet <- file.path( transPath, transFileSet)

	# see if we are doing all as one big set, or are they by subset...
	poolSet <- rep( 1, times=length( transFileSet))
	if ( is.list( sampleIDset)) {
		for ( k in 1:length(sampleIDset)) {
			thisset <- sampleIDset[[k]]
			poolSet[ RP_samples %in% thisset] <- k
		}
	}

	# get the weights for ranking the results
	wt.folds <- as.numeric( getOptionValue( optT, "DE.weightFoldChange", notfound="1"))
	wt.pvalues <- as.numeric( getOptionValue( optT, "DE.weightPvalue", notfound="1"))

	# allow for a species specific minimum
	minRPKM <- as.numeric( getOptionValue( optT, "DE.minimumRPKM", speciesID=speciesID, 
			notfound="1", verbose=T))

	#intensityColumn <- if (useMultiHits) "RPKM_M" else "RPKM_U"
	intensityColumn <- getExpressionUnitsColumn( optionsFile, useMultiHits=useMultiHits)

	# ready to do the Rank Product...
	genesToPlot <- vector()
	htmlPath <- RP_path

	if (doDE) {

	    cat( "\n\nLoading Transcriptomes..")
	    m <- expressionFileSetToMatrix( fnames=transFileSet, fids=transFIDs, geneColumn="GENE_ID",
			intensityColumn=intensityColumn, missingGenes="fill", keepIntergenics=keepIntergenics)
	    cat( "  Done.\n")


	    # use a local function so we can call in parallel if needed
	    processOneRankProductGroup <- function( targetgroup) {

		cat( "\n\nDoing Rank Product on:  ", targetgroup)
		myGenesToPlot <- vector()

		out <- rankProductDiffExpress( transFileSet, fids=transFIDs, RP_groups, targetgroup, 
				m=m, geneColumn="GENE_ID", intensityColumn=intensityColumn,
				offset=minRPKM, keepIntergenics=keepIntergenics, 
				poolSet=poolSet, average.FUN=sqrtmean, nSimulations=nSimulations, 
				missingGenes="fill")

		outfile <- paste( targetgroup, RP_prefix, "RP.Ratio.txt", sep=".")
		if ( ! is.null( altGeneMap)) {
			outfile <- paste( targetgroup, RP_prefix, altGeneMapLabel, "RP.Ratio.txt", sep=".")
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
		cat( "\nWrote RankProduct Gene Data:  ", outfile)

		# HTML too...
		htmlFile1 <- sub( "Ratio.txt$", "UP.html", basename(outfile))
		htmlFile2 <- sub( "Ratio.txt$", "DOWN.html", basename(outfile))

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
		title1 <- paste( "Rank Product: &nbsp;  Genes most up-regulated in group: &nbsp;  ", targetgroup)
		title2 <- paste( "Rank Product: &nbsp;  Genes most down-regulated in group: &nbsp;  ", targetgroup)
		if ( ! is.null( RP_altGeneMapLabel)) {
			title1 <- paste( "Rank Product: &nbsp;  ", RP_altGeneMapLabel, 
					"  most up-regulated in group: &nbsp;  ", targetgroup)
			title2 <- paste( "Rank Product: &nbsp;  ", RP_altGeneMapLabel, 
					"  most down-regulated in group: &nbsp;  ", targetgroup)
			# for plots of alternate gene maps, we need to fusge a few items...
			out <- addAltGeneIdentifierColumn( out)
			extraCols <- extraCols + 1
			fullGname <- out$GENE_ID
			HTML_geneColumn <- "GENE_ID"
		}

		out1 <- out[ out$LOG2FOLD > 0, ]
		if ( nrow(out1) > 0) {
			# clean up formats...
			# only keep those that are UP, and just the 'UP' columns
			downColumns <- grep( "DOWN", colnames(out1))
			if (length(downColumns)) out1 <- out1[ , -downColumns]
			nColShow <- ncol(out1)

			out1$PRODUCT <- gsub( "   ", " &nbsp; ", out1$PRODUCT)
			out1$LOG2FOLD <- formatC( out1$LOG2FOLD, format="f", digits=3, flag="+")
			out1$RP_VALUE <- formatC( out1$RP_VALUE, format="e", digits=2)
			out1$AVG_RANK <- formatC( out1$AVG_RANK, format="f", digits=1, big.mark=",")
			out1$E_VALUE <- formatC( out1$E_VALUE, format="f", digits=4)
			out1$FP_RATE <- formatC( out1$FP_RATE, format="f", digits=4)
			for ( k in (8+extraCols):nColShow) {
				out1[[k]] <- formatC( out1[[k]], format="f", digits=2, big.mark=",")
			}
			colnames(out1)[c(3:7 + extraCols)] <- c( "Log2 Fold", "RP Value", "Avg Rank", "E Value", "FDR")
			colnames(out1)[(8+extraCols):nColShow] <- gsub( "_", " ", colnames(out1)[(8+extraCols):nColShow])
			colnames(out1)[(8+extraCols):nColShow] <- gsub( "\\.", " ", colnames(out1)[(8+extraCols):nColShow])
			# write it
			geneTableToHTMLandPlots( geneDF=out1[ , 1:nColShow], RP_samples, RP_colors, N=Nshow, 
				title=title1, 
				htmlFile=htmlFile1, html.path=htmlPath, results.path=resultsPath, 
				makePlots=FALSE)
			myGenesToPlot <- base::union( myGenesToPlot, unique.default( fullGname[1:Nshow]))
		}

		# for the DOWN table, we now have separate explicit DOWN columns
		# move them into place
		out2 <- out
		out2$RP_VALUE <- out2$RP_VALUE_DOWN
		out2$AVG_RANK <- out2$AVG_RANK_DOWN
		out2$E_VALUE <- out2$E_VALUE_DOWN
		out2$FP_VALUE <- out2$FP_VALUE_DOWN
		# only keep those that are DOWN
		out2 <- out2[ out2$LOG2FOLD < 0, ]
		if ( nrow(out2) > 0) {
			# clean up formats...
			ord <- order( out2$RP_VALUE)
			out2 <- out2[ ord, ]
			rownames(out2) <- 1:nrow(out2)
			# since we move these 'down' columns, remove them
			downColumns <- grep( "DOWN", colnames(out2))
			if (length(downColumns)) out2 <- out2[ , -downColumns]
			nColShow <- ncol(out2)

			out2$PRODUCT <- gsub( "   ", " &nbsp; ", out2$PRODUCT)
			out2$LOG2FOLD <- formatC( out2$LOG2FOLD, format="f", digits=3, flag="+")
			out2$RP_VALUE <- formatC( out2$RP_VALUE, format="e", digits=2)
			out2$AVG_RANK <- formatC( out2$AVG_RANK, format="f", digits=1, big.mark=",")
			out2$E_VALUE <- formatC( out2$E_VALUE, format="f", digits=4)
			out2$FP_RATE <- formatC( out2$FP_RATE, format="f", digits=4)
			for ( k in (8+extraCols):nColShow) {
				out2[[k]] <- formatC( out2[[k]], format="f", digits=2, big.mark=",")
			}
			colnames(out2)[c(3:7 + extraCols)] <- c( "Log2 Fold", "RP Value", "Avg Rank", "E Value", "FDR")
			colnames(out2)[(8+extraCols):nColShow] <- gsub( "_", " ", colnames(out2)[(8+extraCols):nColShow])
			colnames(out2)[(8+extraCols):nColShow] <- gsub( "\\.", " ", colnames(out2)[(8+extraCols):nColShow])
			# write it
			geneTableToHTMLandPlots( geneDF=out2[ , 1:nColShow], RP_samples, RP_colors, N=Nshow, 
				title=title2, 
				htmlFile=htmlFile2, html.path=htmlPath, results.path=resultsPath, 
				makePlots=FALSE)
			myGenesToPlot <- base::union( myGenesToPlot, unique.default( rev(fullGname)[1:Nshow]))
		}

		# all done.  Send back the genes worth plotting
		return( myGenesToPlot)
	    } 

	
	    # either do all the groups in a for loop, or parallel
	    allGroups <- sort( unique( RP_groups))
	    if ( forceMulticore) {
		mcAns <- multicore.lapply( allGroups, processOneRankProductGroup)
		genesToPlot <- unlist( mcAns)
	    } else {
		for ( targetgroup in allGroups) {
			genesOneGroup <- processOneRankProductGroup( targetgroup)
			genesToPlot <- c( genesToPlot, genesOneGroup)
		}
	    }
	    genesToPlot <- unique.default( genesToPlot)

	} else {
		cat( "\nSkipping DE...  Gathering genes for plots..")
		for (targetgroup in sort( unique( RP_groups))) {
			outfile <- paste( targetgroup, RP_prefix, "RP.Ratio.txt", sep=".")
			if ( ! is.null( altGeneMap)) {
				outfile <- paste( targetgroup, RP_prefix, altGeneMapLabel, "RP.Ratio.txt", sep=".")
			}
			outfile <- file.path( RP_path, outfile)
			tmp <- read.delim( outfile, as.is=T)
			genesToPlot <- c( genesToPlot, tmp$GENE_ID[1:Ngenes], rev(tmp$GENE_ID)[1:Ngenes])
		}
		genesToPlot <- unique.default( genesToPlot)
	}

	# make the table of all genen expression
	cat( "\nGathering Expression data for clustering..")
	tm <- expressionFileSetToMatrix( fnames=transFileSet, fids=transFIDs, intensityColumn=intensityColumn, 
					keepIntergenics=keepIntergenics)
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
		if (ncol(tm) < 3) cat( "\nToo few samples to cluster...")
	}
	
	# after all tables and results, make those gene plots
	genesToPlot <- sort( genesToPlot)
	if ( ! is.null(altGeneMap)) genesToPlot <- sort( unique( sub( "::.+", "", genesToPlot)))

	# put in chromosomal order?
	whereGmap <- match( genesToPlot, gmap$GENE_ID, nomatch=NA)
	genesToPlot <- genesToPlot[ order( whereGmap)]

	if ( is.null(PLOT.FUN) || is.function( PLOT.FUN)) {
		geneTableToHTMLandPlots( geneDF=NULL, RP_samples, RP_colors, N=Ngenes, htmlFile=htmlFile, 
				html.path=htmlPath, results.path=resultsPath, makePlots=TRUE, 
				genesToPlot=genesToPlot, label=folderName, 
				geneNameColumn=HTML_geneColumn, PLOT.FUN=PLOT.FUN, ...)
	}

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\n\nFinished pipe 'RankProduct' on Sample Set:     ", unlist(sampleIDset), 
			"\nSpecies: ", speciesID,"\n", label, "\n")
	}
	return()
}

