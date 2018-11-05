/* fastColorSpace2DNA.c

	do the conversion from colorspace to DNA
*/

#include <stdlib.h>
#include <stdio.h>


void fastColorSpace2DNA( int *nReads, int *lengths, char **txtIn, char **txtOut) {

	//  Inputs --
      	//  nReads:     the number of character strings in txtIn
	//  lengths:    the number of characters in each string of txtIn
	//  txtIn:	the array of colorspace strings

	//  Outputs --
	//  txtOut:	the array of DNA strings to return

	int i, j, nStrings, nchar, offset;
	char curBase, curColor, newChar;
	char  *curStrIn, *curStrOut;

	char A_Jump[5]={'A','C','G','T','N'};
	char C_Jump[5]={'C','A','T','G','N'};
	char G_Jump[5]={'G','T','A','C','N'};
	char T_Jump[5]={'T','G','C','A','N'};
	char *curJumpTable;

	nStrings = *nReads;
	for( i=0; i < nStrings; i++) {
		curStrIn = txtIn[i];
		curStrOut = txtOut[i];
		nchar = lengths[i];

		// in a colorspace read, the first character is one of ACGT, then the rest are integers
		// the output will be 1 character shorter than the input
		curBase = curStrIn[0];
		for ( j=1; j < nchar; j++) {
			curColor = curStrIn[j];

			// the current base determines which jump table
			switch( curBase) {
				case 'A':  curJumpTable = A_Jump; break;
				case 'C':  curJumpTable = C_Jump; break;
				case 'G':  curJumpTable = G_Jump; break;
				case 'T':  curJumpTable = T_Jump; break;
			}

			// the current color determines how we step
			switch( curColor) {
				case '0':  offset = 0; break;
				case '1':  offset = 1; break;
				case '2':  offset = 2; break;
				case '3':  offset = 3; break;
				default:   offset = 4; break;
			}
			newChar = curJumpTable[offset];
			curStrOut[j-1] = newChar;
			if ( newChar != 'N') curBase = newChar;
		}
		// null terminate and clobber the last unused character too
		curStrOut[nchar-1] = '\0';
		//curStrOut[nchar] = '\0';
	}
}
