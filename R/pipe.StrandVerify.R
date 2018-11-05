# pipe.StrandVerify.R

# look at the STRAND statistics of a transcriptome to judge correct sense vs antisense settings

`pipe.StrandVerify` <- function( sampleID, annotationFile="Annotation.txt", optionsFile="Options.txt", 
				speciesID=getCurrentSpecies(), results.path=NULL, 
				altGeneMap=NULL, altGeneMapLabel=NULL,
				min.reads=100, min.strandCC=0.5) 
{

	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	strandSpecific <- getAnnotationTrue( annT, sampleID, "StrandSpecific", notfound=FALSE, verbose=F)
	if ( ! strandSpecific) {
		cat( "\nSample is not flagged for strand-specific reads.  Ignored..")
		return(NULL)
	}

	# during 'altGeneMap' runs, skip over species that are not covered...
	doingAltGeneMap <- FALSE
	if ( ! is.null( altGeneMap)) {
		doingAltGeneMap <- TRUE
	}

	setCurrentSpecies( speciesID)
	speciesPrefix <- getCurrentSpeciesFilePrefix()

	if ( is.null( altGeneMap)) {
		# regular...
		fileTrans <- paste( sampleID, speciesPrefix, "Transcript.txt", sep=".")
		fileTrans <- file.path( results.path, "transcript", fileTrans)
	} else {
		if ( is.null( altGeneMapLabel) || base::nchar( altGeneMapLabel) < 1) 
				stop( "pipe.StrandVerify: missing file path term 'altGeneMapLabel' ")
		fileTrans <- paste( sampleID, speciesPrefix, altGeneMapLabel, "Transcript.txt", sep=".")
		fileTrans <- file.path( results.path, "transcript", altGeneMapLabel, fileTrans)
	}

	if ( ! file.exists( fileTrans)) {
		cat( "\nExisting Transcriptome file not found: ", fileTrans)
		return(NULL)
	}
	tbl <- read.delim( fileTrans, as.is=T)

	return( strandVerify( tbl, min.reads=min.reads, min.strandCC=min.strandCC))
}


`strandVerify` <- function( tbl, geneColumn="GENE_ID", readColumn="READS_M", strandColumn="STRAND_M",
				min.reads=100, min.strandCC=0.5) {

	gid <- tbl[[ geneColumn]]
	readCnt <- as.numeric( tbl[[ readColumn]])
	strandValue <- as.numeric( tbl[[ strandColumn]])
	N <- nrow(tbl)

	isNG <- grep( "(ng)", gid, fixed=T)
	isG <- setdiff( 1:N, isNG)
	hasNG <- ( length(isNG) >= N * 0.01)

	# test #1:  are the NonGenes near the bottom?
	pval1 <- NA
	if (hasNG) {
		pval1 <- wilcox.test( isG, isNG, alternate="greater")$p.value
	}
	
	# test #2, are the strand call values mostly positive and near one
	hasCount <- which( readCnt >= min.reads)
	use <- intersect( isG, hasCount)
	sVals <- strandValue[ use]
	sMean <- mean( sVals)
	sSD <- sd( sVals)
	pval2 <- t.test( sVals, mu=min.strandCC, altern="greater")$p.value

	return( list( "strandMean"=sMean, "NonGene_Pvalue"=pval1, "StrandValue_Pvalue"=pval2))
}
