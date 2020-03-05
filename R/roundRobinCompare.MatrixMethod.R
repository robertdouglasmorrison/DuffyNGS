# roundRobinCompare.R

# run a N-way round robin comparison of the differential expression between a set of samples
# where each sample is a member of some 'group'.   Final results are group vs group gene rankings.


`roundRobinCompare.MatrixMethod` <- function( sampleIDset, speciesID=getCurrentSpecies(), annotationFile="Annotation.txt",
		optionsFile="Options.txt", useMultiHits=TRUE, keepIntergenics=FALSE, 
		results.path=NULL, folderName="", groupColumn="Group", colorColumn="Color",
		altGeneMap=NULL, altGeneMapLabel=NULL, geneColumnHTML="GENE_ID", 
		average.FUN=sqrtmean, Ngenes=100, useLog=FALSE, 
		wt.fold=1, wt.pvalue=2, minRPKM=1,
		verbose=!interactive(), doDE=TRUE, PLOT.FUN=NULL, ...) {

	# set up for this species...  Target setup done up in the calling routine!
	setCurrentSpecies( speciesID)
	
	# sanity check on the inputs...
	if ( length( unlist(sampleIDset)) < 2) stop( "RoundRobinCompare requires at least 2 sampleIDs")
	if ( base::nchar( folderName) < 1) stop( "Round Robin needs an explicit 'folderName' argument...")

	# setup persistent storage for RR
	annT <- readAnnotationTable( annotationFile)

	# we used to allow lists of samples...  Not any more...
	flatSamples <- unlist( sampleIDset)
	RR_samples <- base::sort( unique( flatSamples))
	NS <- length( RR_samples)
	where <- match( RR_samples, annT$SampleID, nomatch=0)
	if ( any( where == 0)) stop( "Some named SampleIDs not in Annotation File")
	RR_groups <- checkGroupNames( annT[[ groupColumn]][where])
	RR_colors <- annT[[ colorColumn]][where]
	RR_species <- speciesID
	RR_prefix <- getCurrentSpeciesFilePrefix()
	unique_RR_groups <- base::sort( unique.default( RR_groups))
	unique_RR_colors <- RR_colors[ match( unique_RR_groups, RR_groups)]
	N_RR_groups <- length( unique_RR_groups)
	if ( N_RR_groups < 2) {
		stop( paste( "\nDifferential Expression needs at least 2 groups.  Check Annotation column: ", groupColumn))
	} else {
		cat( "\nSample counts by group:\n")
		print( table( RR_groups))
	}

	RR_altGeneMapLabel <- altGeneMapLabel
	HTML_geneColumn <- geneColumnHTML

	# set up directories to read from or write to...
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".")
	}
	RR_path <- file.path( results.path, "RoundRobin", paste( RR_prefix, folderName, sep="."))
	if ( ! file.exists( RR_path)) dir.create( RR_path, recursive=T, showWarnings=F)

	# we could be given a non-standard geneMap to use...
	if ( is.null( altGeneMap)) {
		# regular...
		gmap <- getCurrentGeneMap()
		RR_altGeneMapLabel <- NULL
		trans.path <- file.path( results.path, "transcript")
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
		trans.path <- file.path( results.path, "transcript", altGeneMapLabel)
	}

	# get the weights for ranking the results
	wt.folds <- as.numeric( getOptionValue( optionsFile, "DE.weightFoldChange", notfound="1"))
	wt.pvalues <- as.numeric( getOptionValue( optionsFile, "DE.weightPvalue", notfound="1"))

	#intensityColumn <- if (useMultiHits) "RPKM_M" else "RPKM_U"
	intensityColumn <- getExpressionUnitsColumn( optionsFile, useMultiHits=useMultiHits)

	# get the set of transcriptome filenames we will work from
	transFiles <- file.path( trans.path, paste( RR_samples, RR_prefix, "Transcript.txt", sep="."))
	if ( ! all( exists <- file.exists( transFiles))) {
		cat( "\nError.  Missing some transcriptome files: ", sum( !exists))
		cat( "\nCan't find:  ", basename( transFiles[ ! exists]))
		return(NULL)
	}
	
	# read in one transcriptome as a reference, for extracting the other details we will need later
	refTrans <- read.delim( transFiles[1], as.is=T)

	# load matrices of all the values we need to do Round Robin in place, without using the old "Ratio" files
	if (useMultiHits) {
		readsColumn <- "READS_M"
		sigmaColumn <- "SIGMA_M"
	} else {
		readsColumn <- "READS_U"
		sigmaColumn <- "SIGMA_U"
	}

	# we only need these matrices if we are doing the DE step
	if (doDE) {
	cat( "\nLoading Gene Expression: (units=", intensityColumn, ")..", sep="")
		intenM <- expressionFileSetToMatrix( transFiles, RR_samples, intensityColumn=intensityColumn, 
					keepIntergenics=keepIntergenics, verbose=F)
		cat( "\nLoading Read Counts..")
		readsM <- expressionFileSetToMatrix( transFiles, RR_samples, intensityColumn=readsColumn, 
					keepIntergenics=keepIntergenics, verbose=F)
		cat( "\nLoading Read Sigmas..")
		sigmaM <- expressionFileSetToMatrix( transFiles, RR_samples, intensityColumn=sigmaColumn, 
					keepIntergenics=keepIntergenics, verbose=F)
		# get the other facts we will need, in matrix row order
		whereRef <- match( rownames(intenM), refTrans$GENE_ID)
		refGeneIDs <- refTrans$GENE_ID[ whereRef]
		refProds <- refTrans$PRODUCT[ whereRef]
		refExonBases <- refTrans$N_EXON_BASES[ whereRef]
	}

	# and grow the list of genes that were DE enough to get plotted for the HTML results
	genesToPlot <- vector()



	# local functions to do the Round Robin...

	`empty.RR.data` <- function() {
		
		# set up storage that will grow for each pair added
		# do smarter allocation and re-allocation
		chunkSize <- nrow(intenM) * 10
		g <- vector( "mode"="character")
		x <- v1 <- v2 <- ranks <- vector( "mode"="numeric", length=chunkSize)
		pval <- pvalUp <- pvalDown <- vector( "mode"="numeric", length=chunkSize)
		return( list( "N_USED"=0, "CUR_SIZE"=chunkSize, "GENE_ID"=g, "LOG2FOLD"=x, "PVALUE"=pval, 
				"RANK"=ranks, "VALUE_1"=v1, "VALUE_2"=v2, "P_UP"=pvalUp, "P_DOWN"=pvalDown))
	}


	`add.RR.data` <- function( RR_List, g, x, pval, ranks, v1, v2, pvalUp, pvalDown) {

		# add this new data from one 2-sample ratio DF to the growing list
		# step 0:  unpack the list
		bigG <- RR_List$GENE_ID
		bigX <- RR_List$LOG2FOLD
		bigP <- RR_List$PVALUE
		bigR <- RR_List$RANK
		bigV1 <- RR_List$VALUE_1
		bigV2 <- RR_List$VALUE_2
		bigPU <- RR_List$P_UP
		bigPD <- RR_List$P_DOWN

		# step 1:  see if we need to extend
		Nold <- RR_List$N_USED
		curSize <- RR_List$CUR_SIZE
		Nnew <- length( g)
		if ( (Nold + Nnew) > curSize) {
			extra <- length(g) * 10
			newSize <- curSize + extra
			length(bigG) <- length(bigX) <- length(bigP) <- length(bigR) <- newSize
			length(bigV1) <- length(bigV2) <- length(bigPU) <- length(bigPD) <- newSize
		}

		# step 2:  stuff the data in
		now <- (Nold+1) : (Nold+Nnew)
		bigG[ now] <- g
		bigX[ now] <- x
		bigP[ now] <- pval
		bigR[ now] <- ranks
		bigV1[ now] <- v1
		bigV2[ now] <- v2
		bigPU[ now] <- pvalUp
		bigPD[ now] <- pvalDown
		newSize <- length(bigX)
		newUsed <- max(now)

		# step 3:  make the result
		outList <- list( "N_USED"=newUsed, "CUR_SIZE"=newSize, "GENE_ID"=bigG, "LOG2FOLD"=bigX, "PVALUE"=bigP, 
				"RANK"=bigR, "VALUE_1"=bigV1, "VALUE_2"=bigV2, "P_UP"=bigPU, "P_DOWN"=bigPD)
		return( outList)
	}


	`finalize.RR.data` <- function( RR_List) {

		# turn this fully grown list into a data frame

		# step 1:   trim all vectors to their true used size
		nUsed <- RR_List$N_USED
		bigG <- RR_List$GENE_ID
		bigX <- RR_List$LOG2FOLD
		bigP <- RR_List$PVALUE
		bigR <- RR_List$RANK
		bigV1 <- RR_List$VALUE_1
		bigV2 <- RR_List$VALUE_2
		bigPU <- RR_List$P_UP
		bigPD <- RR_List$P_DOWN
		length(bigG) <- length(bigX) <- length(bigP) <- length(bigR) <- nUsed
		length(bigV1) <- length(bigV2) <- length(bigPU) <- length(bigPD) <- nUsed

		# step 2:  make that data frame
		outList <- data.frame( "GENE_ID"=bigG, "LOG2FOLD"=bigX, "PVALUE"=bigP, "RANK"=bigR, "VALUE_1"=bigV1, 
					"VALUE_2"=bigV2, "P_UP"=bigPU, "P_DOWN"=bigPD, stringsAsFactors=FALSE)
		return( outList)
	}


	`makeOneDE.dataFrame` <- function( i, j, intenM, readsM, sigmaM, wt.fold, wt.pvalue, minRPKM) {

		# given two samples from 2 different RR Groups, make the tiny data frame of Differential Expresssion
		# like the 'Ratio' files.
		x1 <- intenM[ ,i]
		x2 <- intenM[ ,j]
		r1 <- readsM[ ,i]
		r2 <- readsM[ ,j]
		s1 <- sigmaM[ ,i]
		s2 <- sigmaM[ ,j]

		# the fold change is easy
		f <- log2( (x1 + minRPKM) / (x2 + minRPKM))

		# P values use the Rosetta Sigmas
		normFac1 <- 1.0
		normFac2 <- sum( r1, na.rm=T) / sum( r2, na.rm=T)
		pv <- rosettaPvalue( r1, s1, normFac1, r2, s2, normFac2, refExonBases)

		smlDF <- data.frame( "GENE_ID"=rownames(intenM), "PVALUE"=pv, "LOG2FOLD"=f, 
					"VALUE_1"=x1, "VALUE_2"=x2, stringsAsFactors=F)

		# the RR addeer function will do the ordering again, so don't bother here
		return( smlDF)
}


	`roundRobinAddData` <- function( RR_List, thisDE) {

		# extract the wanted terms from this Diff Expression dataset
		# these Ratio objects still call them 'RPKM' for now...
		g <- thisDE$GENE_ID
		x <- thisDE$LOG2FOLD
		pval <- thisDE$PVALUE
		v1 <- thisDE$VALUE_1
		v2 <- thisDE$VALUE_2

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

		outList <- add.RR.data( RR_List, g, x, pval, ranks, v1, v2, pvalUp, pvalDown)
		return( outList)
	}


	`roundRobinResults` <- function( RR_List, groupName, Ngenes=100) {

		# now we have all round robin sets in one place, build the final consensus
		cat( "\nExtracting RR DE Results for group: ", groupName)
		mydf <- RR_List

		# factor by geneID, to get all entries for each gene
		gfac <- factor( mydf$GENE_ID)
		N <- nlevels( gfac)
		gout <- vector( "mode"="character", length=N)
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
			nout <<- i
		})
		gProd <- gene2Product( gout)

		# final order by average of fold, Pvalue, and Rank
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
		outfile <- paste( groupName, RR_prefix, "RR.Ratio.txt", sep=".")
		if ( !is.null( RR_altGeneMapLabel)) outfile <- paste( groupName, RR_prefix, 
						RR_altGeneMapLabel, "RR.Ratio.txt", sep=".")
		outfile <- file.path( RR_path, outfile)

		# add Human ID terms if needed
		extraCols <- 0
		if ( RR_prefix %in% MAMMAL_PREFIXES) {
			out <- addHumanIDterms( out)
			extraCols <- 2
		}
		if ( RR_prefix %in% ORIGID_PARASITE_PREFIXES) {
			out <- addOrigIDterms( out)
			extraCols <- 1
		}
		if ( RR_prefix %in% BACTERIA_PREFIXES) {
			out <- addNameIDterms( out)
			extraCols <- 1
		}
		nColShow <- nColShow + extraCols

		write.table( out, file=outfile, sep="\t", quote=FALSE, row.names=F)
		if (verbose) cat( "\nWrote RoundRobin Gene Data:  ", outfile, "\n")

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
		title1 <- paste( "Round Robin:   Genes most up-regulated in group:   ", groupName)
		title2 <- paste( "Round Robin:   Genes most down-regulated in group:   ", groupName)
		if ( ! is.null( RR_altGeneMapLabel)) {
			title1 <- paste( "Round Robin:   ", RR_altGeneMapLabel, 
					"  most up-regulated in group:   ", groupName)
			title2 <- paste( "Round Robin:   ", RR_altGeneMapLabel, 
					"  most down-regulated in group:   ", groupName)
			# for plots of varGenes we need to fudge a few items...
			out <- addAltGeneIdentifierColumn( out)
			extraCols <- extraCols + 1
			fullGname <- out$GENE_ID
			HTML_geneColumn <- "GENE_ID"
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
			colnames(out1)[7:8 + extraCols] <- paste( c( "", "Not "), gsub("_|\\."," ",groupName), sep="")
			# write it
			geneTableToHTMLandPlots( geneDF=out1[ , 1:nColShow], RR_samples, RR_colors, N=Nshow, title=title1, 
				htmlFile=htmlFile1, html.path=htmlPath, results.path=results.path, makePlots=FALSE)
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
			colnames(out2)[7:8 + extraCols] <- paste( c( "", "Not "), gsub("_|\\."," ",groupName), sep="")
			# write it
			geneTableToHTMLandPlots( geneDF=out2[ , 1:nColShow], RR_samples, RR_colors, N=Nshow, title=title2, 
				htmlFile=htmlFile2, html.path=htmlPath, results.path=results.path, makePlots=FALSE)
			genesToPlot <- base::union( genesToPlot, unique.default( rev(fullGname)[1:Nshow]))
		}

		# send back the genes to plot
		return( genesToPlot)
	}


	`roundRobinMakePlots` <- function( genesToPlot) {

		# after all tables and results, make those gene plots
		if ( ! is.null(altGeneMap)) genesToPlot <- sort( unique( sub( "::.+", "", genesToPlot)))

		# put in chromosomal order?
		whereGmap <- match( genesToPlot, gmap$GENE_ID, nomatch=NA)
		genesToPlot <- genesToPlot[ order( whereGmap)]

		htmlPath <- RR_path
		if ( is.null(PLOT.FUN) || is.function(PLOT.FUN)) {
			geneTableToHTMLandPlots( geneDF=NULL, RR_samples, RR_colors, N=Ngenes, htmlFile=htmlFile, 
					html.path=htmlPath, results.path=results.path, makePlots=TRUE, 
					genesToPlot=genesToPlot, label=folderName, 
					geneNameColumn=HTML_geneColumn, useLog=useLog, PLOT.FUN=PLOT.FUN, ...)
		}

		# all done, nothing to pass back...
		return()
	}


	`roundRobinGroupTranscripts` <- function( tm, rrFactors, units="RPKM_M") {

		# now we can make one big matrix of all transcripts, that shows the group averages, and do a cluster map...
		if ( is.null( RR_altGeneMapLabel)) {
			outfile <- file.path( RR_path, paste( "All", RR_prefix, "GeneData.txt",sep="."))
			cat( "\nMaking 'all transcripts' gene data:  ", outfile, "\n")
		} else {
			# the alt gene map knows these composite names...
			outfile <- file.path( RR_path, paste( "All.", RR_prefix, ".", RR_altGeneMapLabel, 
						"Data.txt", sep=""))
			cat( "\nMaking 'all transcripts' ", RR_altGeneMapLabel, " data:  ", outfile, "\n")
		}

		# make those group average transcripts... and track gene average ranks too
		ng <- nrow( tm)
		ngrps <- nlevels( rrFactors)
		grpNams <- levels( rrFactors)
		grpM <- rnkM <- matrix( 0, ng, ngrps)
		colnames(grpM) <- colnames(rnkM) <- grpNams
		grpSize <- vector( length=ngrps)

		# load the ranks based on the intensities
		ranksM <- tm
		for (i in 1:ncol(tm)) {
			ord <- order( tm[ ,i], decreasing=T)
			ranksM[ ord,i] <- 1:ng
		}

		# OK, do each group
		nDone <- 0
		tapply( 1:ncol(tm), rrFactors, function(x) {
				# given the columns of TM that all belong to one group
				avgInten <- apply( tm[ , x, drop=F], MARGIN=1, sqrtmean)
				avgRank <- apply( ranksM[ , x, drop=F], MARGIN=1, sqrtmean)
				nDone <<- nDone + 1
				grpM[ , nDone] <<- avgInten
				rnkM[ , nDone] <<- avgRank
				grpSize[nDone] <- length(x)
				return()
			})

		# Write out the average transcript for each group
		for ( i in 1:ngrps) {
			smlDF <- data.frame( "GENE_ID"=refGeneIDs, "PRODUCT"=refProds, "RPKM"=round(grpM[,i],digits=4), 
						"RANK"=round(rnkM[,i],digits=2), stringsAsFactors=F)
			ord <- order( smlDF$RPKM, decreasing=T)
			smlDF <- smlDF[ ord, ]
			rownames(smlDF) <- 1:ng
			# give the expression column the correct name
			colnames(smlDF)[3] <- units
			tfile <- file.path( RR_path, paste( grpNams[i], RR_prefix, "RR.Transcript.txt", sep="."))
			write.table( smlDF, tfile, sep="\t", quote=F, row.names=F)
			if (verbose) cat( "\nWrote Group Average Transcriptome: ", tfile)
		}

		# lastly make the one final matrix of all samples and the group averages
		outTM <- data.frame( "GENE_ID"=refGeneIDs, "PRODUCT"=refProds, round(tm,digits=4), 
						stringsAsFactors=FALSE)
		# but only add group averages if there was 2+ replicates in the group
		if( sum( grpSize > 1)) {
			grpM <- grpM[ , grpSize > 1]
			outTM <- cbind( outTM, grpM, stringsAsFactors=FALSE)
		}
		rownames(outTM) <- 1:nrow(outTM)
		write.table( outTM, file=outfile, sep="\t", quote=F, row.names=F)
		return( outTM)
	}


	`roundRobinExtraPlots` <- function( tm=NULL) {

		# make some cluster images...
		if ( is.null( tm)) {
			tmfile <- file.path( RR_path, paste( "All", RR_prefix, "GeneData.txt",sep="."))
			tm <- read.delim( tmfile, as.is=T)
		}

		# we may be given a data frame of the "All.GeneData.txt" object
		if ( colnames(tm)[1] == "GENE_ID") tm <- as.matrix( tm[ ,3:ncol(tm)])

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
					pngFile <- file.path( RR_path, paste( RR_prefix,"Cluster",funcName[i],"png",sep="."))
				} else {
					pltText <- paste( "Transcriptome Clustering:   ", folderName, 
							"\nTranscriptomes for species:   ", speciesID,
							"    using geneMap:  ", RR_altGeneMapLabel)
					pngFile <- file.path( RR_path, paste( RR_prefix, RR_altGeneMapLabel, 
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
			# knowing the colors to show can be tricky, as the number of "group average" columns can vary
			pcaColors <- RR_colors
			if ( ncol(tm) > length( pcaColors)) {
				extraGroups <- colnames(tm)[ (length(pcaColors)+1):ncol(tm)]
				wh <- match( extraGroups, annT[[ groupColumn]])
				extraColors <- annT[[colorColumn]][wh]
				pcaColors <- c( pcaColors, extraColors)
			}
			pltText <- paste( "Transcriptome PCA:   ", folderName,
					"\nTranscriptomes for species:   ", speciesID)
			pngFile <- file.path( RR_path, paste( RR_prefix,"PCA.png",sep="."))
			matrix.PCAplot( tm, main=pltText, col=pcaColors)
			png( filename=pngFile, width=800, height=800, bg="white")
			matrix.PCAplot( tm, main=pltText, col=pcaColors)
			dev.off()
		} else {
			if (ncol(tm) < 3) cat( "\nToo few samples to make cluster plots...")
		}
		return()
	}
	# end of local functions...


	# if not doing DE, just re-draw
	if ( ! doDE) {
		cat( "\nSkipping DE...  Gathering genes for plots..")
		grps <- unique_RR_groups
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
		genesToPlot <- sort( unique.default( genesToPlot))
		if ( ! is.null(altGeneMap)) genesToPlot <- sort( unique( sub( "::.+", "", genesToPlot)))

		# put in chromosomal order?
		whereGmap <- match( genesToPlot, gmap$GENE_ID, nomatch=NA)
		genesToPlot <- genesToPlot[ order( whereGmap)]

		roundRobinExtraPlots()

		# after all tables and results, make those gene plots
		htmlPath <- RR_path
		geneTableToHTMLandPlots( geneDF=NULL, RR_samples, RR_colors, N=Ngenes, htmlFile=htmlFile, 
					html.path=htmlPath, results.path=results.path, makePlots=TRUE, 
					genesToPlot=genesToPlot, label=folderName, 
					geneNameColumn=HTML_geneColumn, useLog=useLog, PLOT.FUN=PLOT.FUN, ...)
		return()
	}


	# now ready to do each group, all the way thru one at a time
	if (verbose) cat("\n")

	# build the factoring for the groups
	RR_GrpFac <- factor( RR_groups)
	RR_GrpPtrs <- tapply( 1:NS, RR_GrpFac, FUN=NULL)

	# accumulate what genes will get drawn
	allGenesToPlot <- vector()

	for ( ig in 1:N_RR_groups) {

		# get the name and member samples for this group, and the set of 'not this group' samples
		thisGroup <- unique_RR_groups[ig]
		mySamples <- which( RR_GrpPtrs == ig)
		otherSamples <- setdiff( 1:NS, mySamples)
			 
		# start empty, then add all the 2-way compares
		rrList <- empty.RR.data()
		cat( "\n")
		for ( i in mySamples) {
			if (verbose) cat( "\rAdding:  ", RR_samples[i], " to group: ", thisGroup)
			for ( j in otherSamples) {
				oneDF <- makeOneDE.dataFrame( i, j, intenM, readsM, sigmaM, wt.fold, wt.pvalue, minRPKM)
				rrList <- roundRobinAddData( rrList, oneDF)
			}
		}
		rrList <- finalize.RR.data( rrList)

		# we have all we need for this group, summarize and write it
		genesOneGroup <- roundRobinResults( rrList, thisGroup, Ngenes=Ngenes)
		allGenesToPlot <- c( allGenesToPlot, genesOneGroup)

		# clean up to save space
		rm( rrList)
	}

	# all groups have been done...  now any last items and make those plots
	# first, make the "group average transcriptomes"
	tm <- roundRobinGroupTranscripts( intenM, RR_GrpFac, units=intensityColumn)

	# second, make the few extra plots about all transcripts
	roundRobinExtraPlots( tm)

	# lastly, plot all the genes that need it
	roundRobinMakePlots( unique( allGenesToPlot))

	cat( "\nRound Robin done.\n")
	return()
}

