# pipe.SpliceCoverageAnalysis.R

`pipe.SpliceCoverageAnalysis` <- function( sampleIDset, annotationFile="Annotation.txt",
				optionsFile="Options.txt", speciesID=NULL, results.path=NULL,
				groupColumn="Group", 
				tool=c("SAM","RankProduct","EdgeR", "DESeq"),
				Nsplices=1000) {

	# get needed paths, etc. from the options file
	annT <- readAnnotationTable( annotationFile)
	optT <- readOptionsTable( optionsFile)
	target <- getAndSetTarget( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	splicePath <- file.path( results.path, "splicing")

	if ( is.null( speciesID)) {
		speciesSet <- getCurrentTargetSpecies()
		cat( "\nNo 'speciesID' specified.  Choices: ", speciesSet)
		return(NULL)
	}

	NSAMPLES <- length( sampleIDset)

	setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()

	files <- paste( sampleIDset, "splice", prefix, "Summary.txt", sep=".")
	files <- file.path( splicePath, files)
	fids <- sampleIDset
	groups <- annT[[groupColumn]][ match( fids, annT$SampleID)]
	targets <- sort( unique( groups))
	NTARGETS <- length( targets)

	tool <- match.arg( tool)
	intensityColumn <- if (tool %in% c("EdgeR","DESeq")) "N_READS" else "PCT_READS"

	# get the matrix of values
	cat( "\nGethering splice coverage depths..")
	m <- expressionFileSetToMatrix( files, fids=sampleIDset, geneColumn="KEY_ID",
			intensityColumn=intensityColumn, missingGenes="fill")
	cat( "  Done.   N_Splices: ", nrow(m))

	# if big, trim down to the most interesting bits
	if ( nrow(m) > Nsplices) {
		cat( "\nTrimming to top ", Nsplices, "variant splices..")
		dif <- apply( m, 1, function(x) diff( range(x)))
		ord <-order( dif, decreasing=T)
		keep <- ord[ 1:Nsplices]
		m <- m[ keep, ]
	}
	saveM <<- m

	# call the chosen DE tool
	ansDE <- NULL
	if ( NTARGETS > 1) {
		for ( k in 1:NTARGETS) {
			cat( "\nTesting for differential splice coverage: \tTool: ", tool, "  \tGroup: ", targets[k], "\n")
			if (tool == "SAM") {
				ansDE[[k]] <- SAM.DiffExpress( files, fids, m=m, groupSet=groups, targetGroup=targets[k],
						geneColumn="KEY_ID", intensityColumn=intensityColumn, useLog=FALSE)
			} else if ( tool == "RankProduct") {
				# rank product dies on zeros...
				smallValue <- min( m[ m > 0], na.rm=T)
				ansDE[[k]] <- rankProductDiffExpress( files, m=m, groupSet=groups, targetGroup=targets[k],
						geneColumn="KEY_ID", intensityColumn=intensityColumn, 
						missingGenes="fill", offset=smallValue*2, nSimu=100)
			} else if ( tool == "EdgeR") {
				# EdgeR gets tricked by very low read counts...
				ansDE[[k]] <- EdgeR.DiffExpress( files, fids, m=m, groupSet=groups, targetGroup=targets[k],
						geneColumn="KEY_ID", intensityColumn=intensityColumn,
						minimumIntensity=1)
			} else {
				# DESeq gets tricked by very low read counts...
				ansDE[[k]] <- DESeq.DiffExpress( files, fids, m=m, groupSet=groups, targetGroup=targets[k],
						geneColumn="KEY_ID", intensityColumn=intensityColumn,
						minimumIntensity=1)
			}
			saveAns <<- ansDE
		}
		names( ansDE) <- targets
	}

	out <- vector( mode="list", length=(NTARGETS+1))
	out[[1]] <- m
	for (i in 1:NTARGETS) out[[i+1]] <- ansDE[[i]]
	names( out) <- c( "matrix", targets)

	return( out)
}


`mergePairedEndSpliceResults` <- function( sampleIDset, annotationFile="Annotation.txt",
				optionsFile="Options.txt", speciesID=NULL, results.path=NULL) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	target <- getOptionValue( optT, "targetID", notfound="HsPf", verbose=F)
	setCurrentTarget( target)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	splicePath <- file.path( results.path, "splicing")

	if ( is.null( speciesID)) {
		speciesSet <- getCurrentTargetSpecies()
		cat( "\nNo 'speciesID' specified.  Choices: ", speciesSet)
		return(NULL)
	}
	NSAMPLES <- length( sampleIDset)
	setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()

	nmerge <- 0
	for ( i in 1:NSAMPLES) {
		files <- paste( sampleIDset[i], "_", 1:2, ".splice.", prefix, ".Summary.txt", sep="")
		files <- file.path( splicePath, files)
		fid <- sampleIDset[i]
		if ( all( file.exists( files))) {
			cat( "\nMerging splices for: ", fid)
			tbl1 <- read.delim( files[1], as.is=T)
			tbl2 <- read.delim( files[2], as.is=T)
			key1 <- tbl1$KEY_ID
			key2 <- tbl2$KEY_ID
			keys <- sort( union( key1, key2))
			NSP <- length( keys)
			wh1 <- match( keys, key1, nomatch=0)
			wh2 <- match( keys, key2, nomatch=0)
			gid <- sid <- rep.int( "", NSP)
			cnt <- rep.int( 0, NSP)
			gid[ wh1 > 0] <- tbl1$GENE_ID[wh1]
			gid[ wh2 > 0] <- tbl2$GENE_ID[wh2]
			sid[ wh1 > 0] <- tbl1$SPLICE_ID[wh1]
			sid[ wh2 > 0] <- tbl2$SPLICE_ID[wh2]
			cnt[ wh1 > 0] <- cnt[ wh1 > 0] + tbl1$N_READS[wh1]
			cnt[ wh2 > 0] <- cnt[ wh2 > 0] + tbl2$N_READS[wh2]
			pct <- cnt * 100 / sum(cnt)
			strPct <- formatC( pct, format="e", digits=2)
			out <- data.frame( "GENE_ID"=gid, "SPLICE_ID"=sid, "N_READS"=cnt, "PCT_READS"=strPct,
					"KEY_ID"=keys, stringsAsFactors=FALSE)
			ord <- order( pct, decreasing=T)
			out <- out[ ord, ]
			rownames(out) <- 1:nrow(out)

			outfile <- paste( fid, "splice", prefix, "Summary.txt", sep=".")
			outfile <- file.path( splicePath, outfile)
			write.table( out, outfile, sep="\t", quote=F, row.names=F)
			nmerge <- nmerge + 1
		} 
	}
	cat( "\nTotal Samples with merged splice file pairs: ", nmerge)
}
