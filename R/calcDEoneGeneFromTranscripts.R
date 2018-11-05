# calcDEoneGeneFromTranscripts.R

# evaluate the differential expression between to samples for one gene...
`calcDEoneGeneFromTranscripts` <- function( intA, intB,
					totalReads1, totalReads2, minRPKM=2, clipFold=10.0) {

	# extract the values from the 2 transcript subsets...
	rawA <- intA$READS_U[1]
	sigA <- intA$SIGMA_U[1]
	rpkmA <- intA$RPKM_U[1]
	rawAm <- intA$READS_M[1]
	sigAm <- intA$SIGMA_M[1]
	rpkmAm <- intA$RPKM_M[1]
	rawB <- intB$READS_U[1]
	sigB <- intB$SIGMA_U[1]
	rpkmB <- intB$RPKM_U[1]
	rawBm <- intB$READS_M[1]
	sigBm <- intB$SIGMA_M[1]
	rpkmBm <- intB$RPKM_M[1]
	nBases <- intA$N_EXON_BASES[1]

	# make a Pvalue... needs normalized counts...
	facA <- 1.0
	facB <- totalReads1$Unique / totalReads2$Unique
	Pvalue <- rosettaPvalue( rawA, sigA, facA, rawB, sigB, facB, nBases)	

	# RPKM stuff
	rpkmFold <- calculateRPKMfold( rpkmA, rpkmB, minRPKM, clipFold) 

	facA <- 1.0
	facB <- (totalReads1$Unique + totalReads1$Multi) / (totalReads2$Unique + totalReads2$Multi)
	combo_Pvalue <- rosettaPvalue( rawAm, sigAm, facA, rawBm, sigBm, facB, nBases)	
	combo_rpkmFold <- calculateRPKMfold( rpkmAm, rpkmBm, minRPKM, clipFold) 


	out <- list( "Pvalue"=Pvalue, "rpkmFold"=rpkmFold, "rpkm1"=rpkmA, "rpkm2"=rpkmB, 
		"rawReads1"=rawA, "rawReads2"=rawB, 
		"Pvalue.Multi"=combo_Pvalue, "rpkmFold.Multi"=combo_rpkmFold, 
		"rpkm1.Multi"=rpkmAm, "rpkm2.Multi"=rpkmBm, 
		"rawReads1.Multi"=rawAm, "rawReads2.Multi"=rawBm, 
		"detectability"= NA, "nBases"=nBases)

	return( out)
}

