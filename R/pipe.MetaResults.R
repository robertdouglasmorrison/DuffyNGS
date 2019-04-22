# pipe.MetaResults.R

# run the 5 different DE tools and combine their results

`pipe.MetaResults` <- function( sampleIDset, speciesID=getCurrentSpecies(), annotationFile="Annotation.txt",
				optionsFile="Options.txt", 
				useMultiHits=TRUE, results.path=NULL,  folderName="", 
				groupColumn="Group", colorColumn="Color", average.FUN=sqrtmean, 
				tools=c("RoundRobin", "RankProduct", "SAM", "EdgeR", "DESeq"), 
				altGeneMap=NULL, altGeneMapLabel=NULL, targetID=NULL,
				Ngenes=100, geneColumnHTML=if (speciesID %in% MAMMAL_SPECIES) "NAME" else "GENE_ID", 
				keepIntergenics=FALSE, verbose=TRUE, doDE=TRUE, makePlots=doDE, copyPlots=makePlots,
				nFDRsimulations=0,
				addCellTypes=(speciesID %in% MAMMAL_SPECIES), 
				addLifeCycle=(speciesID %in% PARASITE_SPECIES), PLOT.FUN=NULL, ...)
{

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nStarting pipe 'MetaResults' on Sample Set: \n")
		print(sampleIDset)
		cat("\n\nUsing results from Species:  ", speciesID,"\n")
	}

	# set up for this species...
	annT <- readAnnotationTable( annotationFile)
	if ( ! (groupColumn %in% colnames(annT))) {
		cat( "\n\nError: Sample grouping column: '", groupColumn, "' not found in annotation file.", sep="")
		stop()
	}
	groupIDs <- checkGroupNames( annT[[ groupColumn]])

	optT <- readOptionsTable( optionsFile)
	if ( is.null( targetID)) targetID <- getOptionValue( optT, "targetID", notfound="HsPf", verbose=F)
	setCurrentTarget( targetID)
	setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".")
	}
	
	# first call all 5 DE tools

	# plotting no longer works in multicore mode...
	specialPlotMode <- FALSE
	save.PLOT.FUN <- PLOT.FUN
	save.multicore <- multicore.currentCoreCount()
	if ( doDE && multicore.currentCoreCount() > 1 && (is.null(PLOT.FUN) || is.function(PLOT.FUN))) {
		cat( "\nWarning:  no gene plotting in multicore mode..")
		PLOT.FUN <- NA
		specialPlotMode <- TRUE
	}

	# put the slowest tools first, so they all end asap
	toolFuncList <- vector( mode="list")
	if ( "RoundRobin" %in% tools) toolFuncList <- c( toolFuncList, pipe.RoundRobin)
	if ( "RankProduct" %in% tools) toolFuncList <- c( toolFuncList, pipe.RankProduct)
	if ( "DESeq" %in% tools) toolFuncList <- c( toolFuncList, pipe.DESeq)
	if ( "EdgeR" %in% tools) toolFuncList <- c( toolFuncList, pipe.EdgeR)
	if ( "SAM" %in% tools) toolFuncList <- c( toolFuncList, pipe.SAM)

	if ( doDE) {

		saveMClapplyAns <<- multicore.lapply( toolFuncList, FUN=function(x) x( sampleIDset, speciesID, 
					annotationFile, optionsFile, useMultiHits=useMultiHits, results.path=results.path, 
					groupColumn=groupColumn, colorColumn=colorColumn,
					folderName=folderName, altGeneMap=altGeneMap, altGeneMapLabel=altGeneMapLabel,
					Ngenes=Ngenes, geneColumnHTML=geneColumnHTML, keepIntergenics=keepIntergenics,
					targetID=targetID, verbose=verbose, PLOT.FUN=PLOT.FUN, 
					doDE=doDE, ...))
	} # end of 'doDE'

	if ( (specialPlotMode || makePlots) && exists( "toolFuncList")) {
		cat( "\n\nNow recalling DE tools just to make plots..")
		multicore.setup(1)
		lapply( toolFuncList, FUN=function(x) x( sampleIDset, speciesID, 
					annotationFile, optionsFile, useMultiHits=useMultiHits, results.path=results.path, 
					groupColumn=groupColumn, colorColumn=colorColumn,
					folderName=folderName, altGeneMap=altGeneMap, altGeneMapLabel=altGeneMapLabel,
					Ngenes=Ngenes, geneColumnHTML=geneColumnHTML, keepIntergenics=keepIntergenics,
					targetID=targetID, verbose=verbose, PLOT.FUN=save.PLOT.FUN, 
					doDE=FALSE, ...))
		PLOT.FUN <- save.PLOT.FUN
		multicore.setup( save.multicore)
	}

	# now we can do meta analysis
	if (verbose) cat( "\nAll DE steps done.\n\nStarting 'MetaResults' phase...\n")

	wt.fold <- 1; #as.numeric( getOptionValue( optT, "DE.weightFoldChange", notfound="1"))
	wt.pval <- 1; #as.numeric( getOptionValue( optT, "DE.weightPvalue", notfound="1"))
	wt.rank <- 1; #as.numeric( getOptionValue( optT, "DE.weightRank", notfound="1"))

	metaPath <- file.path( results.path, "MetaResults", paste( prefix, folderName, sep="."))
	if ( ! file.exists( metaPath)) dir.create( metaPath, recursive=TRUE)
	pngPath <- "./pngPlots"

	flatSamples <- sort( unique( unlist( sampleIDset)))
	annT2 <- subset( annT, SampleID %in% flatSamples)
	myGrps <- sort( unique( annT2[[ groupColumn]]))
	for ( grp in myGrps) {

		cat( "\n\nDoing MetaResults on:    ", grp)
		for (direction in c( "UP", "DOWN")) {

		out <- metaResults( targetGroup=grp, results.path=results.path, speciesID=speciesID,
				geneColumn="GENE_ID", subfolderName=folderName, tools=tools,
				altGeneMapLabel=altGeneMapLabel,
				rank.average.FUN=average.FUN, value.average.FUN=mean, 
				keepIntergenics=keepIntergenics, topFolder=NULL, 
				other.DE.files=NULL, missingGenes="na", nFDRsimulations=nFDRsimulations,
				direction=direction)
		if ( nrow(out) < 1) next

		
		# lastly, we may have a better ordering by using all 3 features
		myFold <- out$LOG2FOLD
		if ( direction == "DOWN") myFold <- -myFold
		myPval <- out$AVG_PVALUE
		myRank <- out$AVG_RANK
		ord <- diffExpressMetaResultOrder( myFold, myPval, myRank, wt.fold=wt.fold, wt.pvalue=wt.pval, wt.rank=wt.rank)
		out <- out[ ord, ]
		rownames(out) <- 1:nrow(out)

		# round to sensible digits of resolution
		out$LOG2FOLD <- round( out$LOG2FOLD, digits=4)
		out$AVG_RANK <- round( out$AVG_RANK, digits=2)

		if (addCellTypes) {
			cellType <- gene2CellType( out$GENE_ID)
			out <- cbind( out[,1:2], "CellType"=cellType, out[,3:ncol(out)], stringsAsFactors=F)
		} else if (addLifeCycle) {
			lifeCycle <- gene2LifeCycle( out$GENE_ID)
			out <- cbind( out[,1:2], "LifeCycle"=lifeCycle, out[,3:ncol(out)], stringsAsFactors=F)
		}

		# make the text file version
		outDF <- out
		if ( speciesID %in% MAMMAL_SPECIES) outDF <- addHumanIDterms( outDF)
		if ( speciesID %in% ORIGID_PARASITE_SPECIES) outDF <- addOrigIDterms( outDF)
		if ( speciesID %in% BACTERIA_SPECIES) outDF <- addNameIDterms( outDF)

		fileout <- paste( grp, prefix, "Meta", direction, "txt", sep=".")
		if ( ! is.null( altGeneMapLabel)) {
			fileout <- paste( grp, prefix, altGeneMapLabel, "Meta", direction, "txt", sep=".")
		}
		fileout <- file.path( metaPath, fileout)
		write.table( outDF, fileout, sep="\t", quote=F, row.names=F)
		rm( outDF)

		if (addCellTypes) {
			fout <- paste( grp, prefix, "Meta", direction, "GeneCellTypeEnrichment.csv", sep=".")
			fout <- file.path( metaPath, fout)
			if ( direction == "UP") {
				who <- which( out$AVG_PVALUE < 0.05 & out$LOG2FOLD > 0.20)
			} else {
				who <- which( out$AVG_PVALUE < 0.05 & out$LOG2FOLD < -0.20)
			}
			n5pct <- max( 5, round( nrow(out) * 0.05))
			who <- sort( union( who, 1:n5pct))
			enrich <- cellTypeEnrichment( out$CellType[who], mode="genes", speciesID=speciesID,
						upOnly=F, minEnrich=1, maxPval=1, correct=T, 
						geneUniverse=out$GENE_ID, verbose=F)
			write.table( enrich, fout, sep=",", quote=T, row.names=F)
		} else if (addLifeCycle) {
			fout <- paste( grp, prefix, "Meta", direction, "GeneLifeCycleEnrichment.csv", sep=".")
			fout <- file.path( metaPath, fout)
			if ( direction == "UP") {
				who <- which( out$AVG_PVALUE < 0.05 & out$LOG2FOLD > 0.20)
			} else {
				who <- which( out$AVG_PVALUE < 0.05 & out$LOG2FOLD < -0.20)
			}
			n5pct <- max( 5, round( nrow(out) * 0.05))
			who <- sort( union( who, 1:n5pct))
			enrich <- lifeCycleEnrichment( out$LifeCycle[who], mode="genes", speciesID=speciesID,
						upOnly=F, minEnrich=1, maxPval=1, correct=T, 
						geneUniverse=out$GENE_ID, verbose=F)
			write.table( enrich, fout, sep=",", quote=T, row.names=F)
		}
			

		# simplify the names?
		fullGname <- out$GENE_ID
		if ( geneColumnHTML != "GENE_ID") {
			gmap <- getCurrentGeneMap()
			wh <- match( fullGname, gmap$GENE_ID, nomatch=0)
			fullGname[ wh > 0] <- gmap[[ geneColumnHTML]][wh]
			out$GENE_ID <- fullGname
		}

		#if ( speciesID %in% c("Hs_grc", "Mmu_grc", "MacMu")) out <- addHumanIDterms( out)
		if ( speciesID %in% ORIGID_PARASITE_SPECIES) out <- addOrigIDterms( out)
		if ( speciesID %in% BACTERIA_SPECIES) out <- addNameIDterms( out)

		# special mods for altGeneMap...
		Nshow <- Ngenes
		otherGrps <- setdiff( myGrps, grp)
		if ( length(otherGrps) > 1) otherGrps <- paste( "{", paste(sort(otherGrps),collapse=" + "), "}", sep="")
		title <- paste( "MetaResults: &nbsp;  Genes most", direction, "regulated in group: &nbsp;  ", grp,
				" &nbsp; vs &nbsp; ", otherGrps)
		if ( ! is.null( altGeneMapLabel)) {
			title <- paste( "MetaResults: &nbsp;  ", altGeneMapLabel, 
					"  most", direction, "regulated in group: &nbsp;  ", grp,
					" &nbsp; vs &nbsp; ", otherGrps)
			# for plots of alt gene maps, we need to fudge a few items...
			out <- addAltGeneIdentifierColumn( out)
			fullGname <- out$GENE_ID
			geneColumnHTML <- "GENE_ID"
			#if ( regexpr( "vargene|vardomain", tolower(altGeneMapLabel)) > 0) {
			#	out <- cbind( "DOMAIN_ID"=fullGname, out)
			#	out$GENE_ID <- sub( "::.*", "", fullGname)
			#	fullGname <- out$GENE_ID
			#	geneColumnHTML <- "GENE_ID"
			#	if ( "ORIG_ID" %in% colnames(out)) out$ORIG_ID <- gene2OrigID( out$GENE_ID)
			#}
		}

		if ( direction == "UP") {
			out <- out[ out$LOG2FOLD > 0, ]
		} else {
			out <- out[ out$LOG2FOLD < 0, ]
		}
		N <- min( nrow( out), Ngenes)
		htmlout <- sub( ".txt$", ".html", fileout)
		metaResultsToHTML( out, htmlout, title, maxRows=N, linkColumnName="GENE_ID", linkPaths=pngPath)
		rm( out)
	}}  # for direction and group

	# copy all the gene plots to this new results location
	if ( copyPlots) {
		if (verbose) cat( "\nCopying plots.. ")
		for ( folder in c( "RoundRobin", "RankProduct", "SAM", "EdgeR", "DESeq")) {
			pathFrom <- file.path( results.path, folder, paste( prefix, folderName, sep="."), 
						"pngPlots")
			if ( ! file.exists( pathFrom)) next
	
			cmdLine <- paste( " cp -r ", pathFrom, "  ", metaPath)
			system( cmdLine)

			# also grab and copy any Cluster and PCA plots
			pathFrom <- file.path( results.path, folder, paste( prefix, folderName, sep="."))
			fset <- dir( pathFrom, pattern="(Cluster|PCA).+png$", full.name=T)
			if ( length(fset)) {
				fout <- file.path( metaPath, basename(fset))
				file.copy( fset, fout)
			}
			if (verbose) cat( "  ", folder)
		}
		if (verbose) cat( "  Done.\n")
	}

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\n\nFinished pipe 'MetaResults' on Sample Set:     ", unlist(sampleIDset), 
			"\nSpecies: ", speciesID,"\n")
	}
	return()
}

