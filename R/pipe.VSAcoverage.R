# pipe.VSAcoverage.R -- align a sample against the VSA var gene 'genomes'


`pipe.VSAcoverage` <- function( sampleID=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt",
				results.path=NULL, doBowtie=FALSE, doTabulate=doBowtie,
				#genomes=c("3D7","DD2","HB3","IGH","IT4","PFCLIN","RAJ116","REICH"),
				genomes=c("3D7","DD2","HB3","IGH","IT4","PFCLIN","RAJ116","JOS"),
				verbose=TRUE) {

	# file to process comes from annotation and options...
	annT <- readAnnotationTable( annotationFile)
	inputFastqFile <- getAnnotationValue( annT, sampleID, "Filename", verbose=F)
	fastq.path <- getOptionValue( optionsFile, "fastqData.path", verbose=F)
	inputFastqFileSet <- strsplit( inputFastqFile, split=", *")[[1]]
	nFastq <- length( inputFastqFileSet)
	inputFastqFileSet <- file.path( fastq.path, inputFastqFileSet)
	for ( i in 1:nFastq) inputFastqFileSet[i] <- allowCompressedFileName( inputFastqFileSet[i])
	

	# path for all results
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	}
	vsa.path <- file.path( results.path, "VSA", sampleID)
	if ( ! file.exists(vsa.path)) dir.create( vsa.path, recursive=TRUE)

	vsagmap <- getVSAgeneMap()
	vsadmap <- getVSAdomainMap()

	NGenomes <- length( genomes)
	for ( ig in 1:NGenomes) {
		thisGenome <- genomes[ig]
		thisSpecies <- paste( "VSA", thisGenome, "vargenes", sep="_")
		thisPrefix <- paste( "VSA", thisGenome, sep="_")
		if ( thisGenome == "JOS") thisSpecies <- thisPrefix <- "JOS"
		cat( "\nSetting target genome: ", thisGenome)
		addTarget( targetID=thisSpecies, prefix=thisPrefix)
		if ( getCurrentSpecies() != thisSpecies) {
			cat( "\nFailed to find VSA genome MapSet for genome:  ", thisGenome)
			cat( "\nSkipping...")
			next
		}
		smap <- getCurrentSeqMap()
		ord <- order( smap$SEQ_ID)
		smap <- smap[ ord, ]

		# step 1:  align all raw reads against this tiny var genome
		bamFile <- file.path( vsa.path, paste( sampleID, thisGenome, "bam", sep="."))
		if ( doBowtie || !file.exists(bamFile)) {
			cat( "\nAligning ", sampleID, " against: ", thisGenome)
			bowtieIndexFile <- paste( "VSA", thisGenome, "vargenes", sep="_")
			ans <- fastqToBAM( inputFastqFileSet, bamFile, k=1, sampleID=sampleID, optionsFile=optionsFile,
					annotationFile=annotationFile, noHitsFile=NULL, alignIndex=bowtieIndexFile, 
					keepUnaligned=FALSE, verbose=F)
			cat( "\nBowtie done for genome: ", thisGenome, "\n")
			print( ans)
			if ( !file.exists(bamFile)) {
				cat( "\nBowtie alignment to VSA genome made no alignments...")
				next
			}
		}
	
		# step 2:  extract the hits/locations/lengths
		hitFile <- file.path( vsa.path, paste( sampleID, thisGenome, "ReadHits.txt", sep="."))
		if ( doBowtie || !file.exists( hitFile)) {
			cat( "\nExtracting BAM hits: ", thisGenome)
			reader <- bamReader( bamFile)
			refData <- getRefData( reader)
			chunkSize <- 100000
			hasMore <- TRUE
			nalign <- 0
			sidOut <- posOut <- lenOut <- vector()
			repeat {
				if ( ! hasMore) break
				chunk <- getNextChunk( reader, n=chunkSize)
				if ( (nNow <- size(chunk)) < 1) break
				if ( nNow < chunkSize) hasMore <- FALSE
				sid <- refID2seqID( refID( chunk), refData=refData)
				pos <- position(chunk)
				len <- alignLength( chunk)
				now <- (nalign+1) : (nalign + nNow)
				sidOut[now] <- sid
				posOut[now] <- pos
				lenOut[now] <- len
				nalign <- nalign + nNow
			}
			tbl <- data.frame( "GENE_ID"=sidOut, "POSITION"=posOut, "END"=(posOut+lenOut-1), stringsAsFactors=FALSE)
			write.table( tbl, hitFile, sep="\t", quote=F, row.names=F)
		} else {
			tbl <- read.delim( hitFile, as.is=T)
		}

		# step 3:  tablulate for this genome
		depthFile <- file.path( vsa.path, paste( sampleID, thisGenome, "ReadDepth.txt", sep="."))
		countFile <- file.path( vsa.path, paste( sampleID, thisGenome, "ReadCounts.txt", sep="."))
		if ( doBowtie || doTabulate || !file.exists( countFile)) {
			cat( "\nTabulating Read Depths: ", thisGenome)
			nc <- nrow( smap)
			nr <- max( smap$LENGTH)
			m <- matrix( 0, nrow=nr, ncol=nc)
			# let's use the original 3D7 naming, to be like the LAX tool
			if ( thisGenome == "3D7") {
				setCurrentSpecies( "Pf3D7")
				smap$SEQ_ID <- gene2OrigID( smap$SEQ_ID)
				tbl$GENE_ID <- gene2OrigID( tbl$GENE_ID)
			}
			colnames(m) <- smap$SEQ_ID
			who <- match( tbl$GENE_ID, colnames(m), nomatch=0)
			if ( any( who == 0)) {
				tbl <- tbl[ who > 0, ]
				who <- who[ who > 0]
			}
			for ( i in 1:nrow(tbl)) {
				j <- who[i]
				from <- tbl$POSITION[i]
				to <- min( tbl$END[i], smap$LENGTH[j])
				m[ from:to, j] <- m[ from:to, j] + 1
			}
			write.table( m, depthFile, sep="\t", quote=F, row.names=F)
			sidFac <- factor( tbl$GENE_ID, levels=smap$SEQ_ID)
			cnts <- tapply( 1:nrow(tbl), sidFac, length)
			cntValues <- as.vector(cnts)
			cntNames <- names(cnts)
			td7Names <- vsagmap$BEST_3D7_ID[ match( cntNames, vsagmap$GENE_ID)]
			isNA <- which( is.na( td7Names))
			if ( length(isNA)) td7Names[isNA] <- cntNames[isNA]

			cntTbl <- data.frame( "GENE_ID"=cntNames, "BEST_3D7_ID"=td7Names, "READ_COUNT"=cntValues, 
						"RPK"=round(cntValues*1000/smap$LENGTH, digits=3), stringsAsFactors=F)
			write.table( cntTbl, countFile, sep="\t", quote=F, row.names=F)
			cat( "\nDone: ", thisGenome, "\n")
		}
	}

	# OK, we have all the details for each VSA genome, now combine them all together
	setCurrentSpecies( "Pf3D7")
	allCounts <- data.frame()
	allMatrix <- vector( mode="list", length=NGenomes)
	bigNgenes <- 0

	# Step 1:  load them all into memory
	for ( ig in 1:NGenomes) {
		thisGenome <- genomes[ig]
		depthFile <- file.path( vsa.path, paste( sampleID, thisGenome, "ReadDepth.txt", sep="."))
		m <- as.matrix( read.delim( depthFile, as.is=T))
		allMatrix[[ig]] <- m
		countFile <- file.path( vsa.path, paste( sampleID, thisGenome, "ReadCounts.txt", sep="."))
		tbl <- read.delim( countFile, as.is=T)
		bigNgenes <- max( bigNgenes, nrow(tbl))
		allCounts <- rbind( allCounts, data.frame( "VSA"=thisGenome, tbl, stringsAsFactors=F))
	}
	names(allMatrix) <- genomes

	# step 2: make one matrix that holds counts and gene names for all VSA
	allocSize <- bigNgenes * 5
	rpkM <- matrix( NA, nrow=allocSize, ncol=NGenomes)
	gidM <- matrix( "", nrow=allocSize, ncol=NGenomes)
	colnames(rpkM) <- colnames(gidM) <- genomes
	all3D7 <- sort( unique( allCounts$BEST_3D7_ID))
	curRow <- 0
	for ( gid3D7 in all3D7) {
		sml <- subset( allCounts, BEST_3D7_ID == gid3D7)
		sml <- sml[ order( sml$VSA, sml$GENE_ID), ]
		# visit the first of each gene. repeating until all have been seen
		while ( nrow(sml) > 0) {
			myVSA <- unique( sml$VSA)
			where <- match( myVSA, sml$VSA)
			curRow <- curRow + 1
			for ( wh in where) {
				j <- match( sml$VSA[wh], genomes)
				gidM[ curRow, j] <- sml$GENE_ID[wh]
				rpkM[ curRow, j] <- sml$RPK[wh]
			}
			sml <- sml[ -where, ]
			# if more genes left to visit, replicate this row
			if ( nrow(sml) > 0) {
				gidM[ curRow+1, ] <- gidM[ curRow, ]
				rpkM[ curRow+1, ] <- rpkM[ curRow, ]
			}
		}
	}
	# trim to true size
	gidM <- gidM[ 1:curRow, ]
	rpkM <- rpkM[ 1:curRow, ]
	# repackage into final data frame layout
	out <- data.frame( "GENE_ID"=alias2Gene(gidM[ ,1]), "GENE_3D7"=gidM[ ,1], "RPK_3D7"=rpkM[ ,1], stringsAsFactors=F)
	for ( j in 2:NGenomes) {
		nc <- ncol(out)
		out <- cbind( out, "GENE"=gidM[ ,j], "RPK"=rpkM[ ,j], stringsAsFactors=F)
		colnames(out)[nc+1] <- paste( "GENE", genomes[j], sep="_")
		colnames(out)[nc+2] <- paste( "RPK", genomes[j], sep="_")
	}
	outFile <- file.path( vsa.path, paste( sampleID, "AllVSA", "ReadCounts.txt", sep="."))
	write.table( out, outFile, sep="\t", quote=F, row.names=F, na="")
	cat( "\nWrote 'AllVSA' counts summary table.")

	# step 3:  make visuals
	png.path <- file.path( vsa.path, "pngPlots")
	if ( ! file.exists(png.path)) dir.create( png.path, recursive=T)
	# start with all the 3D7 gene names
	allGenes <- sort( unique( gidM[ ,1]))
	# for each gene, find all the VSA genes tied to this one
	cat( "\nPlotting..\n")
	for ( g in allGenes) {
		who <- which( gidM[ ,1] == g)
		myVSA <- as.vector( gidM[ who, ])
		myVSA <- setdiff( unique( myVSA), "")
		nVSA <- length( myVSA)

		# size the image based on the number of genes to draw
		if ( dev.cur() > 1) dev.off()
		checkX11( width=11, height=max( 7, (nVSA*3)))
		par( mai=c( 0.4, 0.7, 0.3, 0.2))
		par( mfrow=c( nVSA, 1))

		# find the gene's genome and its pileup trace for all these VSA
		Xmax <- 1000
		Ymax <- 0
		allTraces <- vector( mode="list")
		allGenomes <- vector()
		for (igene in 1:nVSA) {
			thisGenome <- ""
			thisGene <- myVSA[igene]
			for (ig in 1:NGenomes) {
				m <- allMatrix[[ig]]
				where <- match( thisGene, colnames(m), nomatch=0)
				if ( where > 0) {
					thisGenome <- allGenomes[igene] <- genomes[ig]
					thisTrace <- m[ ,where]
					allTraces[[igene]] <- thisTrace
					Ymax <- max( Ymax, thisTrace, na.rm=T)
					break
				}
				where <- match( alias2Gene(thisGene), alias2Gene(colnames(m)), nomatch=0)
				if ( where > 0) {
					thisGenome <- allGenomes[igene] <- genomes[ig]
					thisTrace <- m[ ,where]
					allTraces[[igene]] <- thisTrace
					Ymax <- max( Ymax, thisTrace, na.rm=T)
					break
				}
			}
			if (thisGenome == "") cat( "\nError:  could not find gene for pileup plotting: ", g)
			dmap <- subset( vsadmap, GENE_ID == myVSA[igene])
			Xmax <- max( Xmax, dmap$DNA_STOP, na.rm=T)
		}
		
		# now with all these facts, we can plot
		# lead space for the domains, and clip the highest (usually ATS)
		if (nVSA > 5) {
			Ymax <- Ymax * 0.50
		} else if ( nVSA > 2) {
			Ymax <- Ymax * 0.70
		}
		Ysmall <- -(Ymax*0.13)
		for ( igene in 1:nVSA) {
			thisGene <- myVSA[igene]
			thisGenome <- allGenomes[igene]
			thisTrace <- allTraces[[igene]]
			thisTrace[ thisTrace == 0] <- NA
			plot( thisTrace, type='h', col=2, xlim=c(1,Xmax), ylim=c(Ysmall, Ymax),
				main=paste( sampleID, thisGene, sep=":  "), xlab=NA, ylab=thisGenome,
				font.lab=2, cex.main=1.2, cex.lab=1.2)
			dmap <- subset( vsadmap, GENE_ID == thisGene)
			if ( nrow(dmap)) {
				rect( dmap$DNA_START, Ysmall, dmap$DNA_STOP, Ysmall*0.2, col="goldenrod", border=1)
				text( (dmap$DNA_START + dmap$DNA_STOP)/2, Ysmall*0.65, dmap$DOMAIN_ID, col=1, font=2)
			}
		}
		# write it as all those gene names
		cat( "\r", g, nVSA)
		for ( gname in c(myVSA, alias2Gene(myVSA[1]))) {
			plotFile <- file.path( png.path, paste( gname, "png", sep="."))
			dev.print( png, plotFile, width=1000, height=max(600, nVSA*200))
		}
	}
	cat( "\nMaking HTML..")
	htmlFile <- file.path( vsa.path, paste( sampleID, "AllVSA", "ReadCounts.html", sep="."))
	table2html( out, htmlFile, title=paste( sampleID, "RNA-seq Pileups onto all VSA var genomes", sep=":   "),
			linkColumnNames=colnames(out)[c(1,seq(2,ncol(out),by=2))], linkPaths="pngPlots")
	cat( "\nDone.\n")
}
