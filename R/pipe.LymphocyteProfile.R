# pipe.LymphocyteProfile.R -- various tools for Bcell and Tcell constructs and expression, etc.


# turn RNA-seq data into a profile of lymphocyte expression
`pipe.LymphocyteProfile` <- function( sampleID=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt",
				results.path=NULL, kmerSize=73, doVelvet=TRUE, buildHash=doVelvet, buildContigs=doVelvet,
				doIGBlast=TRUE, verbose=TRUE) {

	if ( is.null( results.path)) results.path <- getOptionValue( optionsFile, "results.path", 
				notfound=".", verbose=F)

	# path for all results
	velvet.path <- file.path( results.path, "VelvetContigs", sampleID, "Lymphocytes")
	contigFile <- file.path( velvet.path, "contigs.fa")

	didVelvet <- FALSE
	if ( doVelvet || ! file.exists(contigFile)) {

		# step 1:   gather all aligned reads that land in/near the lymphocytes loci
		velvetFiles <- "noHits"

		igFile <- file.path( results.path, "fastq", paste( sampleID, "IgHits.fastq.gz", sep="."))
		if ( ! file.exists( igFile)) {
			cat( "\nGathering B cell aligned reads..\n")
			IG <- getLymphocyteRegions( "IG")
			pipe.GatherRegionAlignments( sampleID, IG$SEQ_ID, IG$START, IG$STOP, asFASTQ=TRUE, 
						fastq.keyword="IgHits")
		}
		velvetFiles <- c( "IgHits", velvetFiles)

		trFile <- file.path( results.path, "fastq", paste( sampleID, "TrHits.fastq.gz", sep="."))
		if ( ! file.exists( trFile)) {
			cat( "\nGathering T cell aligned reads..\n")
			TCR <- getLymphocyteRegions( "TCR")
			pipe.GatherRegionAlignments( sampleID, TCR$SEQ_ID, TCR$START, TCR$STOP, asFASTQ=TRUE, 
						fastq.keyword="TrHits")
		}
		velvetFiles <- c( "TrHits", velvetFiles)
	
		# step 2:  create contigs of anything that may be lymphocyte reads
		#		only rebuild the data structures that are missing...
		pipe.VelvetContigs( sampleID, fastqSource=velvetFiles, kmerSize=kmerSize, folderName="Lymphocytes", 
					buildHash=buildHash, buildContigs=buildContigs, makePep=T)
		didVelvet <- TRUE
	}

	# step 3:  Throw those contigs at IgBlast
	blastOutFile <- file.path( velvet.path, "IgBlast.IG.Out.txt")
	if (didVelvet || doIGBlast) {
		cat( "\n\nSearching Contigs for B cell constructs..")
		callIgBlast( fastafile=contigFile, outfile=blastOutFile, db="IG", organism="human", 
					path="~/IgBlast", outfmt=7)
	}
	tbl <- readIgBlastOutput( infile=blastOutFile, m=7, min.blast.score=100)
	summaryFile <- file.path( velvet.path, "IgBlast.IG.Summary.txt")
	write.table( tbl, summaryFile, sep="\t", quote=F, row.names=F)
	cat( "\nWrote IG Summary file: ", summaryFile)
	ans <- profileIgBlastCoverage( tbl, min.blast.score=100)
	profile <- ans$profile
	pctIdent <- ans$identity
	colnames(pctIdent) <- paste( "Ident", colnames(pctIdent), sep="_")
	out <- cbind( round(profile, digits=4), round(pctIdent, digits=4))
	profileFile <- file.path( velvet.path, "IgBlast.IG.Profile.txt")
	write.table( out, profileFile, sep="\t", quote=F, row.names=T)
	cat( "\nWrote IG Profile file: ", profileFile)

	blastOutFile <- file.path( velvet.path, "IgBlast.TCR.Out.txt")
	if (didVelvet || doIGBlast) {
		cat( "\n\nSearching Contigs for T cell constructs..")
		callIgBlast( fastafile=contigFile, outfile=blastOutFile, db="TR", organism="human", 
					path="~/IgBlast", outfmt=7)
	}
	tbl <- readIgBlastOutput( infile=blastOutFile, m=7, min.blast.score=100)
	summaryFile <- file.path( velvet.path, "IgBlast.TCR.Summary.txt")
	write.table( tbl, summaryFile, sep="\t", quote=F, row.names=F)
	cat( "\nWrote TCR Summary file: ", summaryFile)
	ans <- profileIgBlastCoverage( tbl, min.blast.score=100)
	profile <- ans$profile
	pctIdent <- ans$identity
	colnames(pctIdent) <- paste( "Ident", colnames(pctIdent), sep="_")
	out <- cbind( round(profile, digits=4), round(pctIdent, digits=4))
	profileFile <- file.path( velvet.path, "IgBlast.TCR.Profile.txt")
	write.table( out, profileFile, sep="\t", quote=F, row.names=T)
	cat( "\nWrote TCR Profile file: ", profileFile)

}


`pipe.PlotLymphocyteProfile` <- function( sampleID, type=c("IG", "TCR"), min.read.pct=0.5, 
					annotationFile="Annotation.txt", optionsFile="Options.txt",
					results.path=NULL, label="", sublabel=NULL) {

	type <- match.arg( type)
	if ( is.null( results.path)) results.path <- getOptionValue( optionsFile, "results.path", 
				notfound=".", verbose=F)
	velvet.path <- file.path( results.path, "VelvetContigs", sampleID, "Lymphocytes")
	profileFile <- file.path( velvet.path, paste( "IgBlast", type, "Profile.txt", sep="."))

	# the data is stored as a cbind of the profile and then the germline percentages 
	ans <- getIgBlastProfiles( profileFile)
	plotIgBlastProfile( ans, sampleID=sampleID, label=label, sublabel=sublabel, min.read.pct=min.read.pct)
}


`pipe.PlotLymphocyteDiffProfile` <- function( sampleID1, sampleID2, type=c("IG", "TCR"), min.read.pct=0.5, 
					annotationFile="Annotation.txt", optionsFile="Options.txt",
					results.path=NULL, mode=c( "difference", "overlay"),
					label="", sublabel=NULL) {

	type <- match.arg( type)
	mode <- match.arg( mode)

	if ( is.null( results.path)) results.path <- getOptionValue( optionsFile, "results.path", 
				notfound=".", verbose=F)
	velvet.path <- file.path( results.path, "VelvetContigs", sampleID1, "Lymphocytes")
	profileFile1 <- file.path( velvet.path, paste( "IgBlast", type, "Profile.txt", sep="."))
	velvet.path <- file.path( results.path, "VelvetContigs", sampleID2, "Lymphocytes")
	profileFile2 <- file.path( velvet.path, paste( "IgBlast", type, "Profile.txt", sep="."))

	# the data is stored as a cbind of the profile and then the germline percentages 
	ans1 <- getIgBlastProfiles( profileFile1)
	ans2 <- getIgBlastProfiles( profileFile2)

	if ( mode == "difference") {
		ans <- diffIgBlastProfiles( ans1, ans2)
		plotIgBlastDiffProfile( ans, sampleID1=sampleID1, sampleID2=sampleID2, label=label, 
					sublabel=sublabel, min.read.pct=min.read.pct)
	} else {
		plotIgBlastOverlayProfile( ans1, ans2, sampleID1=sampleID1, sampleID2=sampleID2, label=label, 
					sublabel=sublabel, min.read.pct=min.read.pct)
	}
}


`pipe.GatherLymphocyteConstructs` <- function( sampleID=NULL, optionsFile="Options.txt", results.path=NULL, 
					type=c("IG", "TCR"), min.score=200) {

	type <- match.arg( type)

	# path for all results
	if ( is.null( results.path)) results.path <- getOptionValue( optionsFile, "results.path", 
				notfound=".", verbose=F)
	velvet.path <- file.path( results.path, "VelvetContigs", sampleID, "Lymphocytes")
	contigFile <- file.path( velvet.path, "contigs.fa")
	contigFA <- loadFasta( contigFile)

	# see how all the contigs scored
	summaryFile <- file.path( velvet.path, paste( "IgBlast", type, "Summary.txt", sep="."))
	tbl <- read.delim( summaryFile, as.is=TRUE)
	vScore <- tbl$V_SCORE
	djScore <- pmax( tbl$D_SCORE, tbl$J_SCORE)

	# only keep those with good V score and some iD/J constuct
	keep <- which( vScore >= min.score & djScore > 0)
	if (length(keep) < 1) stop( paste( "No contigs with V_SCORE >= ", min.score))
	tbl <- tbl[ keep, ]
	ids <- tbl$CONTIG_ID
	grpID <- gsub( "_", "", paste( tbl$V_GROUP, tbl$J_GROUP, sep="-"))
	idsOut <- paste( sub( "_length.+","",ids), "_", tbl$V_NAME, "_", tbl$D_NAME, "_", tbl$J_NAME, 
			"_Score:", tbl$V_SCORE, sep="")
	where <- match( ids, contigFA$desc, nomatch=0)
	if ( any( where == 0)) stop( "Missing ContigIDs!")
	seqDNA <- contigFA$seq[ where]
	lenDNA <- nchar( seqDNA)

	# turn these to AA, and only keep those that make decent length proteins
	seqAA <- DNAtoBestPeptide( seqDNA)
	lenAA <- nchar( seqAA)
	keep <- which( lenAA >= (0.85*lenDNA/3))
	if (length(keep) < 1) stop( paste( "No contigs with long enough AA translation."))

	out <- data.frame( "ID"=idsOut[keep], "GROUP"=grpID[keep], "SCORE"=tbl$V_SCORE[keep], 
				"AA_SEQ"=seqAA[keep], "DNA_SEQ"=seqDNA[keep], stringsAsFactors=F)
	rownames(out) <- 1:nrow(out)
	return( out)
}
