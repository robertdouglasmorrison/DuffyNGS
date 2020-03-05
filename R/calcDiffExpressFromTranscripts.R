# calcDiffExpressFromTranscripts.R

# calculate the differential expression between two transcripts

`calcDiffExpressFromTranscripts` <- function( file1, file2, minRPKM=2, clipFold=10.0, 
				wt.folds=1.0, wt.pvalues=1.0, 
				fileout=NA, verbose=TRUE) {


	if (verbose) {
		cat( "\n\nCalculating DE for files:\n\t", file1, "\n\t", file2, "\n")
	}

	# turn the 2 transcripts into mathing order...
	trans1 <- read.delim( file1, as.is=T)
	trans2 <- read.delim( file2, as.is=T)

	allG <- sort( intersect( trans1$GENE_ID, trans2$GENE_ID))
	wh1 <- match( allG, trans1$GENE_ID)
	wh2 <- match( allG, trans2$GENE_ID)
	trans1 <- trans1[ wh1, ]
	trans2 <- trans2[ wh2, ]
	NG <- length( allG)

	# allocate strorage for returned data
	geneSet <- gProdSet <- nBaseSet <- vector( length=NG)
	rawCntSetA <- rawCntSetB <- PvalueSet <- vector( length=NG)
	rpkmFoldSet <- rpkmCntSetA <- rpkmCntSetB <- vector( length=NG)
	combo_rawCntSetA <- combo_rawCntSetB <- combo_PvalueSet <- vector( length=NG)
	combo_rpkmFoldSet <- combo_rpkmCntSetA <- combo_rpkmCntSetB <- vector( length=NG)

	# calc the 'net total reads' for RPKM one time
	# try to re-create what the WIG data structure knows
	u1 <- sum( trans1$READS_U)
	m1 <- sum( trans1$READS_M) 
	u2 <- sum( trans2$READS_U)
	m2 <- sum( trans2$READS_M) 
	totalReads1 <- list( "Unique"=u1, "Multi"=(m1 - u1))
	totalReads2 <- list( "Unique"=u2, "Multi"=(m2 - u2))

	# visit every gene
	ngenes <- nDE <- 0

	for( ig in 1:NG) {

		ans <- calcDEoneGeneFromTranscripts( trans1[ ig, ], trans2[ ig, ],
				totalReads1=totalReads1, totalReads2=totalReads2, 
				minRPKM=minRPKM, clipFold=clipFold)

		# load the vectors and do the DE calc
		ngenes <- ngenes + 1
		geneSet[ ngenes] <- gene <- allG[ig]
		gProdSet[ ngenes] <- trans1$PRODUCT[ig]
		nBaseSet[ ngenes] <- trans1$N_EXON_BASES[ig]
		PvalueSet[ ngenes] <- ans$Pvalue
		rawCntSetA[ ngenes] <- ans$rawReads1
		rawCntSetB[ ngenes] <- ans$rawReads2
		rpkmFoldSet[ ngenes] <- ans$rpkmFold
		rpkmCntSetA[ ngenes] <- ans$rpkm1
		rpkmCntSetB[ ngenes] <- ans$rpkm2
		combo_PvalueSet[ ngenes] <- ans$Pvalue.Multi
		combo_rawCntSetA[ ngenes] <- ans$rawReads1.Multi
		combo_rawCntSetB[ ngenes] <- ans$rawReads2.Multi
		combo_rpkmFoldSet[ ngenes] <- ans$rpkmFold.Multi
		combo_rpkmCntSetA[ ngenes] <- ans$rpkm1.Multi
		combo_rpkmCntSetB[ ngenes] <- ans$rpkm2.Multi
		if ( ig %% 1000 == 0) cat( "\r", ig, "\t", gene, "\t", ans$rpkmFold)
	}

	combo_PvalueSet[ combo_PvalueSet > 1] <- 1
	PvalueSet[ PvalueSet > 1] <- 1

	# do PI Values too
	PIvalueSet <- piValue( rpkmFoldSet, PvalueSet)
	combo_PIvalueSet <- piValue( combo_rpkmFoldSet, combo_PvalueSet)

	out <- data.frame( geneSet, gProdSet, combo_PvalueSet, combo_rpkmFoldSet, combo_rpkmCntSetA, 
			combo_rpkmCntSetB, combo_rawCntSetA, combo_rawCntSetB, combo_PIvalueSet,
			PvalueSet, rpkmFoldSet, rpkmCntSetA, rpkmCntSetB, rawCntSetA, rawCntSetB, 
			PIvalueSet, nBaseSet, stringsAsFactors=FALSE)
	colnames( out) <- c("GENE_ID", "PRODUCT", "PVALUE_M", "LOG2FOLD_M", "RPKM_1_M", "RPKM_2_M", 
				"READS_1_M", "READS_2_M", "PIVALUE_M", "PVALUE_U", "LOG2FOLD_U", "RPKM_1_U", 
				"RPKM_2_U", "READS_1_U", "READS_2_U", "PIVALUE_U", "N_EXON_BASES")

	# now sort based on Pvalue and fold..
	ord <- diffExpressRankOrder( out$LOG2FOLD_M, out$PVALUE_M, wt.folds, wt.pvalues)
	out <- out[ ord, ]
	rownames( out) <- 1:nrow(out)
	
	if ( ! is.na( fileout)) {
		if ( getCurrentSpecies() %in% MAMMAL_SPECIES) out <- addHumanIDterms( out)
		if ( getCurrentSpecies() %in% ORIGID_PARASITE_SPECIES) out <- addOrigIDterms( out)
		if ( getCurrentSpecies() %in% BACTERIA_SPECIES) out <- addNameIDterms( out)

		write.table( out, file=fileout, sep="\t", quote=FALSE, row.names=FALSE)
		cat( "\nWrote Differential Expression file:  \t", fileout)
	}

	if (verbose) {
		cat( "\nN_Gene regions processed:         \t", nrow(out))
	}

	return( out)
}

