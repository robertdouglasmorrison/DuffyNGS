/* fastShortReadOverlap.c

	do the Smith-Waterman local 'overlap' alignment quick and dirty
*/

#include <stdlib.h>
#include <stdio.h>


void fastShortReadOverlap( int *id1, int *n1, int *id2, int *n2, double *misMatchPenalty, double *gapPenalty, double *score, int *offset, int *gaps) {

	//  Inputs --
      	//  id1, id2: 	vectors of integer base IDs,  1=A,2=C,etc. These represent the 2 short read sequences
	//  n1,n2:	length of each read

	//  Outputs --
	//  score:	score of best alignment ( approximately 2 * length of perfect overlap )
	//  offset:	the number of bases from the start of read 1 to  the start of read 2 when the best overlap is applied
	//		( negative means read 2 starts first...)
	//  gaps:	positions and lengths of gaps induced by overlap

	int i, j, l1, l2, size1, size2, best1, best2;
	double *F;
	int *T;
	int im1, im1offset, i1offset, ij, bigI;
	double MISMATCH, GAP_COST, thisGap1, thisGap2, thisF, bestF;
	int start1, start2, thisT, gapLoc1, gapLen1;
	int inGap1, gapLoc2, gapLen2, inGap2;

	l1 = *n1;
	l2 = *n2;
	size1 = l1 + 1;
	size2 = l2 + 1;

	F = (double *) malloc( (size_t) size1*size2*sizeof(double));
	for( i=0; i < (size1*size2); i++) F[i] = 0;

	T = (int *) malloc( (size_t) size1*size2*sizeof(int));

	// lets use string 1 going down the rows, string 2 across the columns
	#define true (1)
	#define false (0)
	#define MATCH (2)

	// start the scores and trace at zero
	for(i = 0; i < size2; i++) {
		T[i] = 0;			// across the top all Js of string 2
		F[i] = (double) 0;		// across the top all Js of string 2
	}
	for(i = 0; i < size1; i++) {
		T[i*size2] = 0; 		// down the left all Is of string 1
		F[i*size2] = (double) 0; 	// down the left all Is of string 1
	}
	
	// fill the scores in
	MISMATCH = *misMatchPenalty;
	GAP_COST = *gapPenalty;

	for(i = 1; i < size1; i++) {
		im1 = i - 1;
		im1offset = im1 * size2;
		i1offset = i * size2;
		for(j = 1; j < size2; j++) {

			// get the 3 choices
			thisGap2 = F[ i1offset + j-1] + GAP_COST;
			thisGap1 = F[ im1offset + j] + GAP_COST;
			thisF = F[im1offset + (j-1)];
			if( id1[i-1] == id2[j-1]) {
				thisF += MATCH;
			} else {
				thisF += MISMATCH;
			}

			// whos best?
			ij = i1offset + j;
			if ( thisF > thisGap1 && thisF > thisGap2) {
				F[ij] = thisF;
				T[ij] = 0;
			} else if ( thisGap1 > thisGap2) {
				F[ij] = thisGap1;
				T[ij] = 1;
			} else {
				F[ij] = thisGap2;
				T[ij] = (-1);
			}
		}
	}

	// scores are all known, now find the highest 'far edge' score
	bestF = MISMATCH;
	best1= 0;
	best2= 0;
	for( i=1; i < size1; i++) {
		ij = (i+1)*size2 - 1;
		thisF = F[ ij];
		if ( thisF > bestF) {
			bestF = thisF;
			best1 = i;
			best2 = size2-1;
		}
	}
	bigI = (size1-1) * size2;
	for( j=1; j < size2; j++) {
		ij = bigI + j;
		thisF = F[ ij]; 
		if ( thisF > bestF) {
			bestF = thisF;
			best2 = j;
			best1 = size1-1;
		}
	}

	// got it, now backtrack to the start point
	start1 = best1;
	start2 = best2;
	gapLoc1 = 0;
	gapLen1 = 0;
	inGap1 = false;
	gapLoc2 = 0;
	gapLen2 = 0;
	inGap2 = false;
	while ( start1 > 0 && start2 > 0) {
		ij = start1*size2 + start2;
		thisT = T[ ij];
		if ( thisT == 0) {		// no gaps, both step back
			start1--;
			start2--;
			inGap1 = inGap2 = false;
		} else if ( thisT == 1) {	// adding a gap along string 1 was best
			start1--;
			if ( inGap1) {
				gapLen1++;
			} else {
				inGap1 = true;
				gapLoc1 = start2 + 1;
				gapLen1 = 1;
				inGap2 = false;
			}
		} else {
			start2--;
			if ( inGap2) {
				gapLen2++;
			} else {
				inGap2 = true;
				gapLoc2 = start1 + 1;
				gapLen2 = 1;
				inGap1 = false;
			}
		}
	}

	// ok, we now know where it started... so clean up...
	if (T) free( T);
	if (F) free( F);

	offset[0] = (start1 - start2);
	score[0] = bestF;
	gaps[0] = gapLoc2;
	gaps[1] = gapLen2;
	gaps[2] = gapLoc1;
	gaps[3] = gapLen1;
}
