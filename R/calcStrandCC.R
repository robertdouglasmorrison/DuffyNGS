`calcStrandCC` <-
function( goodCnt, badCnt) {

	# measure like Correlation coefficient of how many reads are on the right vs wrong strand
	tot <- goodCnt + badCnt
	cc <- ifelse ( tot == 0,  0, (goodCnt-badCnt)/tot)
	return( cc)
}

