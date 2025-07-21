# roundRobinCompare.R

# run a N-way round robin comparison of the differential expression between a set of samples
# where each sample is a member of some 'group'.   Final results are group vs group gene rankings.


`roundRobinCompare.RatioMethod` <- function( sampleIDset, speciesID=getCurrentSpecies(), annotationFile="Annotation.txt",
		optionsFile="Options.txt", useMultiHits=TRUE, keepIntergenics=FALSE, 
		results.path=NULL, folderName="", groupColumn="Group", colorColumn="Color",
		altGeneMap=NULL, altGeneMapLabel=NULL, geneColumnHTML="GENE_ID", 
		average.FUN=sqrtmean, Ngenes=100, useLog=FALSE, 
		verbose=!interactive(), doDE=TRUE, PLOT.FUN=NULL, ...) {

	# set up for this species...  Target setup done up in the calling routine!
	setCurrentSpecies( speciesID)
	
	# sanity check on the inputs...
	if ( length( unlist(sampleIDset)) < 2) stop( "RoundRobinCompare requires at least 2 sampleIDs")
	if ( base::nchar( folderName) < 1) stop( "Round Robin needs an explicit 'folderName' argument...")

	# setup persistent storage for RR
	annT <- readAnnotationTable( annotationFile)
	flatSamples <- unlist( sampleIDset)
	where <- match( flatSamples, annT$SampleID, nomatch=0)
	if ( any( where == 0)) stop( "Some named SampleIDs not in Annotation File")
	flatGroups <- annT[[ groupColumn]][ where]
	flatColors <- annT[[ colorColumn]][ where]
	RR_samples <- unique( flatSamples)
	where <- match( RR_samples, annT$SampleID)
	RR_groups <- checkGroupNames( annT[[ groupColumn]][where])
	RR_colors <- annT[[ colorColumn]][where]
	RR_species <- speciesID
	RR_prefix <- getCurrentSpeciesFilePrefix()
	unique_RR_groups <- base::sort( unique.default( RR_groups))
	unique_RR_colors <- RR_colors[ match( unique_RR_groups, RR_groups)]
	N_RR_groups <- length( unique_RR_groups)
	if ( N_RR_groups < 2) {
		stop( paste( "\nDifferential Expression needs at least 2 groups.  Check Annotation column:", groupColumn))
	} else {
		cat( "\nSample counts by group:\n")
		print( table( RR_groups))
	}

	# set up storage for all these RR groups' data
	RR_List <- vector( mode="list", length=N_RR_groups)
	RR_GroupCounts <- rep( 0, times=N_RR_groups)
	names( RR_List) <- unique_RR_groups
	
	# making average transcripts for each group too...
	RR_TransList <- vector( mode="list", length=N_RR_groups)
	names( RR_TransList) <- unique_RR_groups

	RR_altGeneMapLabel <- altGeneMapLabel
	HTML_geneColumn <- geneColumnHTML

	# set up directories to read from or write to...
	if ( is.null( results.path)) {
		resultsPath <- getOptionValue( optionsFile, "results.path", notfound=".")
	} else {
		resultsPath <- results.path
	}
	RR_path <- file.path( resultsPath, "RoundRobin", paste( RR_prefix, folderName, sep="."))
	if ( ! file.exists( RR_path)) dir.create( RR_path, recursive=T, showWarnings=F)

	# we could be given a non-standard geneMap to use...
	if ( is.null( altGeneMap)) {

		# regular...
		gmap <- getCurrentGeneMap()
		RR_altGeneMapLabel <- NULL
		ratiosPath <- file.path( resultsPath, "ratios")
		transPath <- file.path( resultsPath, "transcript")

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

		cat( "\n\nDoing Alternate Gene Map RoundRobinCompare for:  ", altGeneMapLabel, "\n")
		RR_altGeneMapLabel <- altGeneMapLabel
		ratiosPath <- file.path( resultsPath, "ratios", altGeneMapLabel)
		transPath <- file.path( resultsPath, "transcript", altGeneMapLabel)
	}

	# get the weights for ranking the results
	wt.folds <- as.numeric( getOptionValue( optionsFile, "DE.weightFoldChange", notfound="1"))
	wt.pvalues <- as.numeric( getOptionValue( optionsFile, "DE.weightPvalue", notfound="1"))

	# also, build one giant matrix of all transcriptomes as another part of results...
	transFileSet <- transFIDs <- vector()
	#intensityColumn <- if (useMultiHits) "RPKM_M" else "RPKM_U"
	intensityColumn <- getExpressionUnitsColumn( optionsFile, useMultiHits=useMultiHits)

	dev.type <- getPlotDeviceType( optionsFile)


	# local functions to do the Round Robin...

	`empty.RRlistData` <- function() {
		
		# set up storage that will grow for each pair added
		g <- gprod <- vector( "mode"="character")
		x <- v1 <- v2 <- ranks <- vector( "mode"="numeric")
		pval <- pvalUp <- pvalDown <- vector( "mode"="numeric")
		return( list( "GENE_ID"=g, "PRODUCT"=gprod, "LOG2FOLD"=x, "PVALUE"=pval, 
				"RANK"=ranks, "VALUE_1"=v1, "VALUE_2"=v2, "P_UP"=pvalUp, 
				"P_DOWN"=pvalDown))
	}


	`add.RRlistData` <- function( i, g, gprod, x, pval, ranks, v1, v2, pvalUp, pvalDown) {

		# add this data to the growing list
		Nold <- length( RR_List[[i]]$GENE_ID)
		Nnew <- length( g)
		now <- (Nold+1) : (Nold+Nnew)
		RR_List[[i]]$GENE_ID[ now] <<- g
		RR_List[[i]]$PRODUCT[ now] <<- gprod
		RR_List[[i]]$LOG2FOLD[ now] <<- x
		RR_List[[i]]$PVALUE[ now] <<- pval
		RR_List[[i]]$RANK[ now] <<- ranks
		RR_List[[i]]$VALUE_1[ now] <<- v1
		RR_List[[i]]$VALUE_2[ now] <<- v2
		RR_List[[i]]$P_UP[ now] <<- pvalUp
		RR_List[[i]]$P_DOWN[ now] <<- pvalDown
		return()
	}


	`finalize.RRlistData` <- function(i) {

		mylist <- RR_List[[i]]
		out <- data.frame( "GENE_ID"=mylist$GENE_ID, "PRODUCT"=mylist$PRODUCT, "LOG2FOLD"=mylist$LOG2FOLD, 
				"PVALUE"=mylist$PVALUE, "RANK"=mylist$RANK, "VALUE_1"=mylist$VALUE_1, 
				"VALUE_2"=mylist$VALUE_2, "P_UP"=mylist$P_UP, "P_DOWN"=mylist$P_DOWN, 
				stringsAsFactors=FALSE)
		RR_List[[i]] <<- out
		return()
	}


	`roundRobinAddData` <- function( thisDE, i, useMultiHits=TRUE) {

		# we may have a gene set with the fake intergenics... remove them
		if ( ! keepIntergenics) {
			isNG <- grep( "(ng)", thisDE$GENE_ID, fixed=TRUE)
			if ( length( isNG) > 0) {
				thisDE <- thisDE[ -isNG, ]
			}
			nonGenes <- subset( getCurrentGeneMap(), REAL_G == FALSE)$GENE_ID
			drops <- which( thisDE$GENE_ID %in% nonGenes)
			if ( length(drops) > 0) {
				thisDE <- thisDE[ -drops, ]
			}
		}

		# extract the wanted terms from this Diff Expression dataset
		# these Ratio objects still call them 'RPKM' for now...
		g <- thisDE$GENE_ID
		gProd <- thisDE$PRODUCT
		x <- thisDE$LOG2FOLD_U
		pval <- thisDE$PVALUE_U
		v1 <- thisDE$RPKM_1_U
		v2 <- thisDE$RPKM_2_U

		if ( useMultiHits) {
			x <- thisDE$LOG2FOLD_M
			pval <- thisDE$PVALUE_M
			v1 <- thisDE$RPKM_1_M
			v2 <- thisDE$RPKM_2_M
		}

		# clip very small Pvalues
		SMALL_PVALUE <- 1e-20
		GIANT_PVALUE <- 1
		pval[ pval < SMALL_PVALUE] <- SMALL_PVALUE

		# set the rank order of each gene in this DE object, by P-value
		ord <- diffExpressRankOrder( x, pval, wt.folds, wt.pvalues)
		ranks <- vector( length=length(x))
		ranks[ ord] <- 1:length(x)

		# turn all the down regulated genes into terrible P-values
		pvalUp <- pvalDown <- pval
		pvalUp[ x <= 0] <- 1.0 / pvalUp[ x <= 0]
		pvalDown[ x >= 0] <- 1.0 / pvalDown[ x >= 0]

		add.RRlistData( i, g, gProd, x, pval, ranks, v1, v2, pvalUp, pvalDown)
		return()
	}


	`roundRobinAddTranscriptData` <- function( thisT, i, useMultiHits=TRUE) {

		# we may have a gene set with the fake intergenics... remove them
		if ( ! keepIntergenics) {
			isNG <- grep( "(ng)", thisT$GENE_ID, fixed=TRUE)
			if ( length( isNG) > 0) {
				thisT <- thisT[ -isNG, ]
			}
			nonGenes <- subset( getCurrentGeneMap(), REAL_G == FALSE)$GENE_ID
			drops <- which( thisT$GENE_ID %in% nonGenes)
			if ( length(drops) > 0) {
				thisT <- thisT[ -drops, ]
			}
		}

		# extract the wanted terms from this transcript dataset
		g <- thisT$GENE_ID
		gProd <- thisT$PRODUCT
		#x <- thisT$RPKM_U
		x <- thisT[ , intensityColumn]

		#if ( useMultiHits) {
		#	x <- thisT$RPKM_M
		#}

		# set the rank order of each gene in this transcript, by RPKM
		ord <- base::order( x, decreasing=TRUE)
		ranks <- vector( length=length(g))
		ranks[ ord] <- 1:length(g)
		sml <- data.frame( g, gProd, x, ranks, stringsAsFactors=FALSE)
		colnames( sml) <- c( "GENE_ID", "PRODUCT", intensityColumn, "RANK")

		# now add this set to it's group
		if ( is.null( RR_TransList[[i]])) {
			RR_TransList[[i]] <<- sml
		} else {
			RR_TransList[[i]] <<- rbind( RR_TransList[[i]], sml)
		}
		return()
	}



	`roundRobinResults` <- function( Ngenes=100, tailWidth=2000, PLOT.FUN=NULL, ...) {

		# now we have all round robin sets in one place, build the final consensus
		genesToPlot <- vector()

		cat( "\n\nExtracting RR Differential Expression Ratios...")

		do.oneRRresult <- function(k) {

			mydf <- RR_List[[ k]]
			grpName <- unique_RR_groups[k]
			otherGroups <- setdiff( unique_RR_groups, grpName)
			if ( length(otherGroups) > 1) otherGroups <- paste("Not",grpName,sep=".")

			if (verbose) cat( "\nExtracting: ", grpName)

			# factor by geneID, to get all entries for each gene
			gfac <- factor( mydf$GENE_ID)
			N <- nlevels( gfac)
			gout <- gProd <- vector( "mode"="character", length=N)
			fout <- pout <- rout <- rpkm1 <- rpkm2 <- vector( "mode"="numeric", length=N)
			poutUp <- poutDown <- vector( "mode"="numeric", length=N)
			nout <- 0
			# build the consensus average for each
			tapply( 1:nrow(mydf), gfac, function(who) {
				i <- (nout+1)
				fout[i] <<- mean.default( mydf$LOG2FOLD[ who])
				rout[i] <<- average.FUN( mydf$RANK[ who])
				pout[i] <<- p.combine( mydf$PVALUE[ who])
				poutUp[i] <<- p.combine( mydf$P_UP[ who])
				poutDown[i] <<- p.combine( mydf$P_DOWN[ who])
				rpkm1[i] <<- average.FUN( mydf$VALUE_1[ who])
				rpkm2[i] <<- average.FUN( mydf$VALUE_2[ who])
				gout[i] <<- mydf$GENE_ID[ who[1]]
				gProd[i] <<- mydf$PRODUCT[ who[1]]
				nout <<- i
			})

			# final order by average of fold, Pvalue,  ... leave out the ranks...
			# use the final fold change to decide which P-value to carry forward
			pout[ fout > 0] <- poutUp[ fout > 0]
			pout[ fout < 0] <- poutDown[ fout < 0]
			ord <- diffExpressRankRankOrder( fout, pout, rout, wt.folds, wt.pvalues, wt.ranks=1)

			# after the final rankings, then clip the terrible P-values
			pout <- ifelse( pout > 1, 1, pout)
			piout <- piValue( fout, pout)

			# round to sensible digits of resolution
			fout <- round( fout, digits=4)
			rout <- round( rout, digits=2)
			piout <- round( piout, digits=3)
			rpkm1 <- round( rpkm1, digits=2)
			rpkm2 <- round( rpkm2, digits=2)

			out <- data.frame( gout, gProd, fout, pout, rout, piout, rpkm1, rpkm2, poutUp, poutDown, 
					stringsAsFactors=FALSE)
			colnames( out) <- c( "GENE_ID", "PRODUCT", "LOG2FOLD", "PVALUE", "RANK", "PIVALUE",
					"VALUE_1", "VALUE_2", "P_UP", "P_DOWN")
			out <- out[ ord, ]
			rownames( out) <- 1:nrow(out)
			nColShow <- 8

			# write it out
			outfile <- paste( grpName, RR_prefix, "RR.Ratio.txt", sep=".")
			if ( !is.null( RR_altGeneMapLabel)) outfile <- paste( grpName, RR_prefix, 
						RR_altGeneMapLabel, "RR.Ratio.txt", sep=".")
			outfile <- file.path( RR_path, outfile)

			# add Human ID terms if needed
			extraCols <- 0
			if ( RR_prefix %in% MAMMAL_PREFIXES) {
				out <- addHumanIDterms( out)
				extraCols <- 2
			}
			#if ( RR_prefix %in% ORIGID_PARASITE_PREFIXES) {
			#	out <- addOrigIDterms( out)
			#	extraCols <- 1
			#}
			if ( RR_prefix %in% BACTERIA_PREFIXES) {
				out <- addNameIDterms( out)
				extraCols <- 1
			}
			nColShow <- nColShow + extraCols

			outText <- out
			colnames(outText)[ colnames(outText) == "VALUE_1"] <- grpName
			colnames(outText)[ colnames(outText) == "VALUE_2"] <- otherGroups
			write.table( outText, file=outfile, sep="\t", quote=FALSE, row.names=F)
			if (verbose) cat( "\nWrote RoundRobin Gene Data:  ", outfile)

			# HTML too...
			genesToPlot <- vector()
			htmlFile1 <- sub( "Ratio.txt$", "UP.html", basename(outfile))
			htmlFile2 <- sub( "Ratio.txt$", "DOWN.html", basename(outfile))
			htmlPath <- RR_path

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
			title1 <- paste( "Round Robin:   Genes most up-regulated in group:   ", grpName)
			title2 <- paste( "Round Robin:   Genes most down-regulated in group:   ", grpName)
			if ( ! is.null( RR_altGeneMapLabel)) {
				title1 <- paste( "Round Robin:   ", RR_altGeneMapLabel, 
						"  most up-regulated in group:   ", grpName)
				title2 <- paste( "Round Robin:   ", RR_altGeneMapLabel, 
						"  most down-regulated in group:   ", grpName)
				# for plots of varGenes we need to fudge a few items...
				out <- addAltGeneIdentifierColumn( out)
				extraCols <- extraCols + 1
				fullGname <- out$GENE_ID
				HTML_geneColumn <- "GENE_ID"
				#if ( regexpr( "vargene", tolower(RR_altGeneMapLabel)) > 0) {
				#	out <- cbind( "DOMAIN_ID"=fullGname, out)
				#	out$GENE_ID <- sub( "::.*", "", fullGname)
				#	extraCols <- extraCols + 1
				#	fullGname <- out$GENE_ID
				#	HTML_geneColumn <- "GENE_ID"
				#	if ( "ORIG_ID" %in% colnames(out)) out$ORIG_ID <- gene2OrigID( out$GENE_ID)
				#}
			}

			# for the UP table, use that Pvalue and hide the 2 directional Pvalue columns
			out1 <- out
			out1$PVALUE <- out1$P_UP
			# only keep those that are UP
			out1 <- out1[ out1$LOG2FOLD > 0, ]
			if ( nrow(out1) > 0) {
				# clean up formats...
				out1$PRODUCT <- gsub( "   ", " &nbsp; ", out1$PRODUCT)
				out1$LOG2FOLD <- formatC( out1$LOG2FOLD, format="f", digits=3, flag="+")
				out1$PVALUE <- formatC( out1$PVALUE, format="e", digits=2)
				out1$RANK <- formatC( out1$RANK, format="f", digits=2)
				out1$PIVALUE <- formatC( out1$PIVALUE, format="f", digits=3, flag="+")
				out1$VALUE_1 <- formatC( out1$VALUE_1, format="f", digits=2)
				out1$VALUE_2 <- formatC( out1$VALUE_2, format="f", digits=2)
				colnames(out1)[3:6 + extraCols] <- c( "Log2 Fold", "Avg Pvalue", "Avg Rank", "Avg PIvalue")
				colnames(out1)[ colnames(out1) == "VALUE_1"] <- groupName
				colnames(out1)[ colnames(out1) == "VALUE_2"] <- sub( "_|\\."," ",otherGroups)
				# write it
				geneTableToHTMLandPlots( geneDF=out1[ , 1:nColShow], RR_samples, RR_colors, N=Nshow, 
					title=title1, 
					htmlFile=htmlFile1, html.path=htmlPath, results.path=resultsPath, 
					makePlots=FALSE)
				genesToPlot <- base::union( genesToPlot, unique.default( fullGname[1:Nshow]))
			}

			# for the DOWN table, flip it and use the DOWN Pvalues and adjust the ranks
			out2 <- out[ rev( 1:nrow(out)), ]
			out2$PVALUE <- out2$P_DOWN
			# only keep those that are DOWN
			out2 <- out2[ out2$LOG2FOLD < 0, ]
			if ( nrow(out2) > 0) {
				# clean up formats...
				out2$PRODUCT <- gsub( "   ", " &nbsp; ", out2$PRODUCT)
				out2$LOG2FOLD <- formatC( out2$LOG2FOLD, format="f", digits=3, flag="+")
				out2$PVALUE <- formatC( out2$PVALUE, format="e", digits=2)
				out2$RANK <- formatC( out2$RANK, format="f", digits=2)
				out2$PIVALUE <- formatC( out2$PIVALUE, format="f", digits=3, flag="+")
				out2$VALUE_1 <- formatC( out2$VALUE_1, format="f", digits=2)
				out2$VALUE_2 <- formatC( out2$VALUE_2, format="f", digits=2)
				colnames(out2)[3:6 + extraCols] <- c( "Log2 Fold", "Avg Pvalue", "Avg Rank", "Avg PIvalue")
				colnames(out2)[ colnames(out2) == "VALUE_1"] <- groupName
				colnames(out2)[ colnames(out2) == "VALUE_2"] <- sub( "_|\\."," ",otherGroups)
				# write it
				geneTableToHTMLandPlots( geneDF=out2[ , 1:nColShow], RR_samples, RR_colors, N=Nshow, 
					title=title2, 
					htmlFile=htmlFile2, html.path=htmlPath, results.path=resultsPath, 
					makePlots=FALSE)
				genesToPlot <- base::union( genesToPlot, unique.default( rev(fullGname)[1:Nshow]))
			}

			# send back the genes to plot
			return( genesToPlot)
		}

		ans <- multicore.lapply( 1:length( RR_List), FUN=do.oneRRresult)
		genesToPlot <- sort( unique( unlist( ans)))
		if ( ! is.null(altGeneMap)) genesToPlot <- sort( unique( sub( "::.+", "", genesToPlot)))

		# put in chromosomal order?
		whereGmap <- match( genesToPlot, gmap$GENE_ID, nomatch=NA)
		genesToPlot <- genesToPlot[ order( whereGmap)]

		# after all tables and results, make those gene plots
		htmlPath <- RR_path
		if ( is.null(PLOT.FUN) || is.function(PLOT.FUN)) {
			geneTableToHTMLandPlots( geneDF=NULL, RR_samples, RR_colors, N=Ngenes, htmlFile=htmlFile, 
					html.path=htmlPath, results.path=resultsPath, makePlots=TRUE, 
					genesToPlot=genesToPlot, label=folderName, tailWidth=tailWidth,
					geneNameColumn=HTML_geneColumn, useLog=useLog, PLOT.FUN=PLOT.FUN, ...)
		}

		# all done, nothing to pass back...
		return()
	}


	`roundRobinTranscriptResults` <- function() {

		# now we have all round robin sets in one place, build the final consensus

		# for a consensus transcriptome, this is only useful if there was more than one
		# dataset per group, so watch for that...
		cat( "\n\nExtracting RR Average Transcripts...")

		do.oneRRtranscript <- function(k) {

			mydf <- RR_TransList[[ k]]
			grpName <- unique_RR_groups[k]

			# skip if nothing to combine...
			#if ( RR_GroupCounts[k] < 2) next

			if (verbose) cat( "\nExtracting: ", grpName)

			# factor by geneID, to get all entries for each gene
			gfac <- factor( mydf$GENE_ID)
			N <- nlevels( gfac)
			gout <- gProd <- vector( mode="character", length=N)
			tout <- rout <- vector( mode="numeric", length=N)
			nout <- 0

			# build the consensus average for each
			tapply( 1:nrow(mydf), gfac, function(who) {
				i <- nout + 1
				tout[i] <<- average.FUN( mydf[ , intensityColumn][ who])
				rout[i] <<- average.FUN( mydf$RANK[ who])
				gout[i] <<- mydf$GENE_ID[ who[1]]
				gProd[i] <<- mydf$PRODUCT[ who[1]]
				nout <<- i
			})

			# final order by RPKM
			ord <- base::order( tout, decreasing=TRUE)
			out <- data.frame( gout, gProd, tout, rout, stringsAsFactors=FALSE)
			colnames( out) <- c( "GENE_ID", "PRODUCT", intensityColumn, "RANK")
			out <- out[ ord, ]
			rownames( out) <- 1:nrow(out)

			# write it out
			outfile <- paste( grpName, RR_prefix, "RR.Transcript.txt", sep=".")
			if ( !is.null( RR_altGeneMapLabel)) outfile <- paste( grpName, RR_prefix, 
						RR_altGeneMapLabel, "RR.Transcript.txt", sep=".")
			outfile <- file.path( RR_path, outfile)

			# add Human ID terms if needed
			extraCols <- 0
			if ( RR_prefix %in% MAMMAL_PREFIXES) {
				out <- addHumanIDterms( out)
				extraCols <- 2
			}
			#if ( RR_prefix %in% ORIGID_PARASITE_PREFIXES) {
			#	out <- addOrigIDterms( out)
			#	extraCols <- 1
			#}
			if ( RR_prefix %in% BACTERIA_PREFIXES) {
				out <- addNameIDterms( out)
				extraCols <- 1
			}

			write.table( out, file=outfile, sep="\t", quote=FALSE, row.names=F)
			if (verbose) cat( "\nWrote RoundRobin Gene Data:  ", outfile)
			
			# add this to the final matrix of transcripts for clustering...
			# (but only if its a true average of more than 1 transcriptome)
			if ( RR_GroupCounts[k] > 1) {
				return( list( "file"=outfile, "fid"=grpName))
			} else {
				return( list( "file"="", "fid"=""))
			}
		}

		ans <- multicore.lapply( 1:length( RR_TransList), FUN=do.oneRRtranscript)
		transFileSet <<- c( transFileSet, sapply(ans, function(x) x$file))
		transFIDs <<- c( transFIDs, sapply(ans, function(x) x$fid))
		transFileSet <<- setdiff( unique( transFileSet), "")
		transFIDs <<- setdiff( unique( transFIDs), "")

		# if we want a clickable HTML file, do that part here...
		# not yet done...

		# now we can make one big matrix of all transcripts, and do a cluster map...
		cat( "\nGathering all transcriptomes..")
		tm <- expressionFileSetToMatrix( fnames=transFileSet, fids=transFIDs, intensityColumn=intensityColumn)
		gnames <- rownames(tm)
	
		if ( is.null( RR_altGeneMapLabel)) {
			gprods <- gene2Product( gnames)
			outfile <- file.path( RR_path, paste( "All", RR_prefix, "GeneData.txt",sep="."))
			cat( "\nWriting 'all transcripts' gene data:  ", outfile, "\n")
		} else {
			gprods <- gmap$PRODUCT[ base::match( gnames, gmap$GENE_ID)]
			outfile <- file.path( RR_path, paste( "All.", RR_prefix, ".", RR_altGeneMapLabel, 
						"Data.txt", sep=""))
			cat( "\nWriting 'all transcripts' ", RR_altGeneMapLabel, " data:  ", outfile, "\n")
		}
			
		outTM <- data.frame( "GENE_ID"=gnames, "PRODUCT"=gprods, tm, stringsAsFactors=FALSE)
		rownames(outTM) <- 1:nrow(outTM)
		write.table( outTM, file=outfile, sep="\t", quote=F, row.names=F)
		
		roundRobinExtraPlots( tm=tm)

		return()
	}


	`roundRobinExtraPlots` <- function( tm=NULL) {

		if ( is.null( tm)) {
			transFileSet1 <- paste( RR_samples, RR_prefix, "Transcript.txt", sep=".")
			if ( ! is.null( RR_altGeneMapLabel)) {
				transFileSet1 <- paste( RR_samples, RR_prefix, RR_altGeneMapLabel, "Transcript.txt", sep=".")
			}
			transFileSet1 <- file.path( transPath, transFileSet1)
			transFIDs1 <- RR_samples
			transFileSet2 <- paste( unique_RR_groups, RR_prefix, "RR.Transcript.txt", sep=".")
			if ( ! is.null( RR_altGeneMapLabel)) {
				transFileSet2 <- paste( unique_RR_groups, RR_prefix, RR_altGeneMapLabel, "RR.Transcript.txt", sep=".")
			}
			transFileSet2 <- file.path( RR_path, transFileSet2)
			transFIDs2 <- unique_RR_groups
			transFileSet <- c( transFileSet1, transFileSet2)
			transFIDs <- c( transFIDs1, transFIDs2)
	
			# now we can make one big matrix of all transcripts, and do a cluster map...
			cat( "\nRe-gathering all transcriptomes..")
			tm <- expressionFileSetToMatrix( fnames=transFileSet, fids=transFIDs, intensityColumn=intensityColumn)
		}
	
		# make some cluster images...
		if ( ncol(tm) > 2 && (is.null(PLOT.FUN) || is.function(PLOT.FUN))) {
	
	    		# there are two good cluster tools...  let's do both
	    		require( cluster)
	    		func <- list( diana, agnes)
	    		funcName <- c( "Divide", "Aggregate")
	    		subtitle <- c( "Divisive hierarchical clustering (DIANA)", 
	    					"Agglomerative hierarchical clustering (AGNES)")
		
	    		for ( i in 1:2) {
				if ( is.null( RR_altGeneMapLabel)) {
					pltText <- paste( "Transcriptome Clustering:   ", folderName,
							"\nTranscriptomes for species:   ", speciesID)
					plotFile <- file.path( RR_path, paste( RR_prefix,"Cluster",funcName[i], dev.type, sep="."))
				} else {
					pltText <- paste( "Transcriptome Clustering:   ", folderName, 
							"\nTranscriptomes for species:   ", speciesID,
							"    using geneMap:  ", RR_altGeneMapLabel)
					plotFile <- file.path( RR_path, paste( RR_prefix, RR_altGeneMapLabel, 
							"Cluster", funcName[i], dev.type, sep="."))
				}
		
				clusterAns <- expressionCluster( tm, useLog=TRUE, FUN=func[[i]])
				plot( clusterAns, which=2, main=pltText, sub=subtitle[i], font=2)
		
				cat( "\nMaking cluster plot: ", plotFile)
				openPlot( plotFile, bg="white")
				plot( clusterAns, which=2, main=pltText, sub=subtitle[i], font=2)
				dev.off()
	    		}

	    		# PCA plot too...
			# knowing the colors to show can be tricky, as the number of "group average" columns can vary
			pcaColors <- rep.int( 1, ncol(tm))
			wh <- match( colnames(tm), RR_samples, nomatch=0)
			pcaColors[ wh > 0] <-RR_colors[wh]
			wh <- match( colnames(tm), annT[[ groupColumn]], nomatch=0)
			pcaColors[ wh > 0] <- annT[[colorColumn]][wh]
			pltText <- paste( "Transcriptome PCA:   ", folderName,
					"\nSpecies: ", speciesID, "    Expression Units: ", intensityColumn)
			plotFile <- file.path( RR_path, paste( RR_prefix,"PCA", dev.type, sep="."))
			matrix.PCAplot( tm, main=pltText, col=pcaColors)
			openPlot( filename=plotFile, bg="white")
			matrix.PCAplot( tm, main=pltText, col=pcaColors)
			dev.off()
		} else {
			if (ncol(tm) < 3) cat( "\nNot able to make cluster plots...")
		}
		return()
	}
	# end of local functions...


	# the 'sampleIDset' used to be a vector, but now it MAY be a list of vectors
	if ( is.list( sampleIDset)) {
		mySampleSets <- lapply( sampleIDset, function(x) match( x, flatSamples))
	} else {
		mySampleSets <- vector( mode="list", length=1)
		mySampleSets[[1]] <- match( sampleIDset, flatSamples)
	}


	# if not doing DE, just re-draw
	if ( ! doDE) {
		cat( "\nSkipping DE...  Gathering genes for plots..")
		grps <- sort( unique( flatGroups))
		genesToPlot <- vector()
		for( grp in grps) {
			outfile <- paste( grp, RR_prefix, "RR.Ratio.txt", sep=".")
			if ( ! is.null( altGeneMap)) {
				outfile <- paste( grp, RR_prefix, altGeneMapLabel, "RR.Ratio.txt", sep=".")
			}
			outfile <- file.path( RR_path, outfile)
			tmp <- read.delim( outfile, as.is=T)
			genesToPlot <- c( genesToPlot, tmp$GENE_ID[1:Ngenes], rev(tmp$GENE_ID)[1:Ngenes])
		}
		genesToPlot <- unique.default( genesToPlot)
		if ( ! is.null(altGeneMap)) genesToPlot <- sort( unique( sub( "::.+", "", genesToPlot)))

		# put in chromosomal order?
		whereGmap <- match( genesToPlot, gmap$GENE_ID, nomatch=NA)
		genesToPlot <- genesToPlot[ order( whereGmap)]

		roundRobinExtraPlots()

		# after all tables and results, make those gene plots
		htmlPath <- RR_path
		geneTableToHTMLandPlots( geneDF=NULL, RR_samples, RR_colors, N=Ngenes, htmlFile=htmlFile, 
					html.path=htmlPath, results.path=resultsPath, makePlots=TRUE, 
					genesToPlot=genesToPlot, label=folderName, 
					geneNameColumn=HTML_geneColumn, useLog=useLog, PLOT.FUN=PLOT.FUN, ...)
		return()
	}

	# initialize the storage that will grow for each pair of samples
	for ( i in 1:N_RR_groups) RR_List[[i]] <- empty.RRlistData()

	# do all pairs of samples
	if (verbose) cat("\n")
	for ( k in 1:length( mySampleSets)) {
	    thisSet <- mySampleSets[[k]]
	    for ( i in thisSet) {
		iGrp <- flatGroups[i]
		iListPtr <- base::match( iGrp, unique_RR_groups)

		for ( j in thisSet) {
			if ( i == j) next

			# skip if they are from the same group
			jGrp <- flatGroups[j]
			if ( iGrp == jGrp) next

			thisPair <- paste( flatSamples[i], "v", flatSamples[j], sep=".")
			if ( is.null( RR_altGeneMapLabel)) {
				fileInDE <- paste( thisPair, RR_prefix, "Ratio.txt", sep=".")
			} else {
				fileInDE <- paste( thisPair, RR_prefix, RR_altGeneMapLabel, "Ratio.txt", 
						sep=".")
			}
			fileInDE <- file.path( ratiosPath, fileInDE)
			if ( ! file.exists( fileInDE)) {
				cat( "\nDiffExpression result file not found: ", fileInDE,
					"\nSkipping...\n")
				next
			}
			 
			if (verbose) cat( "\rAdding:  ", thisPair, " to group: ", unique_RR_groups[ iListPtr])
			thisDE <- read.delim( fileInDE, as.is=TRUE)
			roundRobinAddData( thisDE, iListPtr, useMultiHits=useMultiHits)
		}

		# also make an 'average transcript' for each group...
		if ( is.null( RR_altGeneMapLabel)) {
			fileInT <- paste( flatSamples[i], RR_prefix, "Transcript.txt", sep=".")
		} else {
			fileInT <- paste( flatSamples[i], RR_prefix, RR_altGeneMapLabel, 
					"Transcript.txt", sep=".")
		}
		fileInT <- file.path( transPath, fileInT)
		if ( ! file.exists( fileInT)) {
			cat( "\nTranscript result file not found: ", fileInT, "\nSkipping...\n")
			next
		}
		if (verbose) cat( "\rAdding:  ", flatSamples[i], " to group: ", unique_RR_groups[ iListPtr])

		# add this transcriptome to the list for gene matrix and clustering..
		transFileSet <- base::append( transFileSet, fileInT)
		transFIDs <- base::append( transFIDs, flatSamples[i])

		# OK, read that transcriptome and add it to the growing dataset
		thisT <- read.delim( fileInT, as.is=TRUE)
		roundRobinAddTranscriptData( thisT, iListPtr, useMultiHits=useMultiHits)
		RR_GroupCounts[iListPtr] <- RR_GroupCounts[iListPtr] + 1
	}
	} 	# end 'for K in list of sets

	# finalize the storage that was growing for each pair of samples
	for (i in 1:N_RR_groups) finalize.RRlistData(i)
	if ( any( (nSize <- sapply( RR_List, nrow)) == 0)) {
		cat( "\n\nSome Round Robin groups have no data:  ", names(RR_List)[nSize == 0], "\n")
		stop( "Round Robin quitting..")
	}

	# all datasets have been added into to round robin pool, now do the tally to make results.
	roundRobinTranscriptResults()
	roundRobinResults( Ngenes=Ngenes, PLOT.FUN=PLOT.FUN, ...)

	cat( "\nRound Robin done.\n")
	return()
}

